"""Separate routing-test metadata and non-inferential cage comparisons.

Do not count repeated endpoints as independent studies or one-cage contrasts
as replicated effects. Check every value against the pinned derived JSON.
"""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

META_FILE="routing_test_metadata.csv"
CAGE_FILE="cage_descriptive_contrasts.csv"
META=("effect_id","family","site_scope","period","endpoint","response_channel",
      "p_holm","family_size","n_assignments","randomization","source_json","source_blob")
CAGE=("contrast_id","panel_id","qc_lane","focal_species","control_treatment",
      "dropout_treatment","control_cage","dropout_cage","n_control_bees","n_dropout_bees",
      "n_positive_control","n_positive_dropout","zero_control_ids_json","zero_dropout_ids_json",
      "strain_named_tv","strain_rank_tv","species_named_tv","species_rank_tv",
      "median_nonfocal_copy_ratio","delta_detected_strains","evidence_class",
      "source_json","source_blob")
FIELD_JSON="validation/ahif_v1/vansteijn2026_nontarget_temporal_results_v1.json"
CAGE_JSON="validation/ahif_v1/prasad2026_nontarget_dropout_cage_descriptive_v1.json"
FAMILY_SIZES={"primary":2,"secondary":16,"order_exchange":3}
INTS={"n_control_bees","n_dropout_bees","n_positive_control","n_positive_dropout"}
METRICS={"strain_named_tv":"named_TV_conditional_positive",
         "strain_rank_tv":"rank_TV_conditional_positive",
         "species_named_tv":"species_TV_conditional_positive",
         "species_rank_tv":"species_rank_TV_conditional_positive",
         "median_nonfocal_copy_ratio":"median_nonfocal_copy_ratio",
         "delta_detected_strains":"delta_detected_strains"}

def read_rows(path,cols,n):
    with path.open(encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f)
        if tuple(r.fieldnames or ())!=cols: raise ValueError("routing schema mismatch")
        rows=list(r)
    if len(rows)!=n or len({x[cols[0]] for x in rows})!=n:
        raise ValueError("routing row count or duplicate ID")
    if any(None in x or any(v is None for v in x.values()) for x in rows):
        raise ValueError("routing truncated row")
    return rows

def load_pinned(source,rows,rel,sha):
    if any(r["source_json"]!=rel or r["source_blob"]!=sha for r in rows):
        raise ValueError("routing source provenance mismatch")
    data=(source.resolve().parents[1]/rel).read_bytes()
    actual=hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()
    if actual!=sha: raise ValueError("routing source Git blob mismatch")
    return json.loads(data)

def close(a,b):
    return math.isfinite(float(a)) and math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-12)

def holm(ps):
    result=[0.]*len(ps)
    ordered=sorted(range(len(ps)),key=lambda i:ps[i])
    previous=0.
    for j,i in enumerate(ordered):
        previous=max(previous,min(1.,(len(ps)-j)*ps[i]))
        result[i]=previous
    return result

def insert_routing(conn,source:Path):
    rows=read_rows(source/META_FILE,META,21)
    original=load_pinned(source,rows,FIELD_JSON,"3a94069f67e34284e5598ffc9c3618e74d77b733")
    expected={(family,x.get("site","MAL_UME"),x["period"],x["endpoint"]):x
              for family in FAMILY_SIZES for x in original[family]}
    if len(expected)!=21: raise ValueError("field source expected 21 unique tests")
    effects={x[0]:x[1:] for x in conn.execute(
        "SELECT effect_id,panel_id,p_value,estimand,estimate FROM effects WHERE effect_id LIKE 'E_VSF26_ROUTING_%'")}
    if set(effects)!={r["effect_id"] for r in rows}: raise ValueError("routing effect join")
    groups=defaultdict(list); seen=set(); vals=[]
    for r in rows:
        family=r["family"]; exact=family=="order_exchange"
        key=(family,r["site_scope"],r["period"],r["endpoint"])
        if key not in expected or key in seen: raise ValueError("routing endpoint mismatch")
        seen.add(key); x=expected[key]
        panel,p,estimand,estimate=effects[r["effect_id"]]
        if panel!="P_VANSTEIJNFIELD2026" or int(r["family_size"])!=FAMILY_SIZES[family]:
            raise ValueError("routing family or source panel mismatch")
        if int(r["n_assignments"])!=(1024 if exact else 9999):
            raise ValueError("routing permutation count")
        if estimand!=("joint_site_order_exchange_statistic" if exact else "E_null_adj_R2"):
            raise ValueError("wrong routing estimand")
        if not close(estimate,x["statistic"] if exact else x["E"]) or not close(p,x["p_exact"] if exact else x["p_upper"]) or not close(r["p_holm"],x["p_holm_family"]):
            raise ValueError("routing numeric source mismatch")
        if p is None or not 0<=float(r["p_holm"])<=1:
            raise ValueError("invalid adjusted p")
        groups[family].append((p,float(r["p_holm"])))
        vals.append([float(r[k]) if k=="p_holm" else int(r[k]) if k in ("family_size","n_assignments") else r[k] for k in META])
    if {k:len(v) for k,v in groups.items()}!=FAMILY_SIZES: raise ValueError("routing incomplete Holm families")
    for g in groups.values():
        if any(not close(a,b) for a,b in zip(holm([x[0] for x in g]),[x[1] for x in g])):
            raise ValueError("Holm mismatch")
    conn.execute("""CREATE TABLE routing_test_metadata(
        effect_id TEXT PRIMARY KEY REFERENCES effects(effect_id),family TEXT NOT NULL,
        site_scope TEXT NOT NULL,period TEXT NOT NULL,endpoint TEXT NOT NULL,response_channel TEXT NOT NULL,
        p_holm REAL NOT NULL CHECK(p_holm BETWEEN 0 AND 1),family_size INTEGER NOT NULL,
        n_assignments INTEGER NOT NULL,randomization TEXT NOT NULL,source_json TEXT NOT NULL,source_blob TEXT NOT NULL)""")
    conn.executemany("INSERT INTO routing_test_metadata VALUES("+",".join("?" for _ in META)+")",vals)
    conn.execute("""CREATE VIEW v_routing_tests AS
        SELECT m.*,e.estimate,e.estimand,e.p_value AS p_unadjusted,s.publication_cluster
        FROM routing_test_metadata m JOIN effects e USING(effect_id)
        JOIN panels p USING(panel_id) JOIN studies s USING(study_id)""")

    cages=read_rows(source/CAGE_FILE,CAGE,12)
    original=load_pinned(source,cages,CAGE_JSON,"da0670ac48fffffeb5b79df73b52c09a24f75695")
    expected={(lane,x["drop"]):x for lane,k in (("source_qc","rows"),
        ("all_exported_sensitivity","all_exported_bees_sensitivity")) for x in original[k]}
    if len(expected)!=12: raise ValueError("cage source expected 12 records")
    if conn.execute("SELECT COUNT(*) FROM effects WHERE panel_id='P_PRASAD_DROPOUT'").fetchone()[0]:
        raise ValueError("one-cage outcomes promoted to effects")
    if conn.execute("SELECT raw_source_audited,screen_status FROM studies WHERE study_id='ST_PRASAD2026'").fetchone()!=("no","derived_matrix_descriptive"):
        raise ValueError("derived matrices not raw reads")
    vals=[]; lanes=defaultdict(list); seen=set()
    for r in cages:
        key=(r["qc_lane"],r["dropout_treatment"])
        if key not in expected or key in seen: raise ValueError("cage key mismatch")
        seen.add(key); x=expected[key]
        if r["panel_id"]!="P_PRASAD_DROPOUT" or r["evidence_class"]!="descriptive_one_cage_per_arm_no_inference":
            raise ValueError("cage evidence grade")
        for k,kx in (("focal_species","focal"),("control_treatment","ctrl"),
                     ("control_cage","control_cage"),("dropout_cage","dropout_cage")):
            if r[k]!=x[kx]: raise ValueError("cage source identity mismatch")
        for arm in ("control","dropout"):
            n=int(r["n_"+arm+"_bees"]); pos=int(r["n_positive_"+arm])
            zero=json.loads(r["zero_"+arm+"_ids_json"])
            if n!=x["n_"+arm] or pos!=x["n_positive_"+arm] or zero!=x["zero_"+arm]:
                raise ValueError("cage source count mismatch")
            if not 0<pos<=n or len(zero)!=n-pos or len(set(zero))!=len(zero):
                raise ValueError("zero-roster gate")
        for k,kx in METRICS.items():
            if not close(r[k],x[kx]): raise ValueError("cage numeric source mismatch")
            if k.endswith("_tv") and not 0<=float(r[k])<=1: raise ValueError("invalid TV")
        if float(r["median_nonfocal_copy_ratio"])<0: raise ValueError("negative copy ratio")
        lanes[r["qc_lane"]].append(r)
        vals.append([int(r[k]) if k in INTS else float(r[k]) if k in METRICS else r[k] for k in CAGE])
    for lane,g in lanes.items():
        if len(g)!=6 or len({z["dropout_cage"] for z in g})!=6 or len({z["control_cage"] for z in g})!=2:
            raise ValueError("shared cage design")
        counts={}
        for z in g:
            for arm in ("control","dropout"):
                key=z[arm+"_cage"]; n=int(z["n_"+arm+"_bees"])
                if key in counts and counts[key]!=n: raise ValueError("shared control count mismatch")
                counts[key]=n
        if len(counts)!=8 or (lane=="source_qc" and sum(counts.values())!=69):
            raise ValueError("cage unit count")
    conn.execute("""CREATE TABLE cage_descriptive_contrasts(
        contrast_id TEXT PRIMARY KEY,panel_id TEXT NOT NULL REFERENCES panels(panel_id),
        qc_lane TEXT NOT NULL,focal_species TEXT NOT NULL,control_treatment TEXT NOT NULL,
        dropout_treatment TEXT NOT NULL,control_cage TEXT NOT NULL,dropout_cage TEXT NOT NULL,
        n_control_bees INTEGER NOT NULL,n_dropout_bees INTEGER NOT NULL,
        n_positive_control INTEGER NOT NULL,n_positive_dropout INTEGER NOT NULL,
        zero_control_ids_json TEXT NOT NULL,zero_dropout_ids_json TEXT NOT NULL,
        strain_named_tv REAL NOT NULL,strain_rank_tv REAL NOT NULL,species_named_tv REAL NOT NULL,species_rank_tv REAL NOT NULL,
        median_nonfocal_copy_ratio REAL NOT NULL,delta_detected_strains REAL NOT NULL,
        evidence_class TEXT NOT NULL,source_json TEXT NOT NULL,source_blob TEXT NOT NULL)""")
    conn.executemany("INSERT INTO cage_descriptive_contrasts VALUES("+",".join("?" for _ in CAGE)+")",vals)
    return {"routing_field_test_rows":21,"routing_field_test_publication_clusters":1,
            "routing_cage_descriptive_contrast_rows":12,"routing_cage_qc_contrasts":6,
            "routing_cage_sensitivity_contrasts":6,
            "routing_cage_independent_treatment_replications_per_arm":1,"routing_cage_meta_eligible":False}
