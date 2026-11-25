"""Finite combinatorial prediction benchmarks; no independent ecological effect promotion."""
import csv,json,math,statistics
from collections import Counter
MODEL_FILE="low_order_model_errors.csv"
CONTRAST_FILE="low_order_model_contrasts.csv"
SRC="validation/ahif_v1/ishizawa2024_equal_information_v1.json"
RUN="https://github.com/zuizui0223/knk/actions/runs/45526767902"
PANEL="P_ISHIZAWA2024"
MODELS=("mono","pair_mean","pair","trio_mean","pair_on_trio_rhs",
        "global_trio_coeff","trio","additive_log","inclusion_exclusion")
ENDPOINTS=("named_TV","rank_TV","identity_assignment_TV","simpson_abs_error",
           "top1_abs_error","log_sum_geometric_cfu_abs_error")
CONTRASTS=(("mono","pair"),("pair","trio"),("pair_mean","trio_mean"),
           ("trio_mean","trio"),("pair_on_trio_rhs","trio"),
           ("global_trio_coeff","trio"),("additive_log","trio"),
           ("inclusion_exclusion","trio"))
M=("panel_id","model","endpoint","mean_error","n_target_types","source_json","source_run")
C=("panel_id","baseline","enhanced","endpoint","absolute_improvement","relative_gain",
   "equal_size_weighted_improvement","improved","worsened","tied","n_target_types","source_json","source_run")

def close(a,b):
    return math.isfinite(float(a)) and math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-11)

def read(path,fields,n):
    with path.open(encoding="utf-8-sig",newline="") as f:
        d=csv.DictReader(f)
        if tuple(d.fieldnames or ())!=fields:raise ValueError("low-order schema")
        rows=list(d)
    if len(rows)!=n or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("low-order count")
    return rows

def insert_low_order(conn,source):
    data=json.loads((source.resolve().parents[1]/SRC).read_text())
    if (data["source_run"],data["source_commit"],data["source_doi"],data["status"])!=(RUN,"b9a81fbcf31bacce8b2db518bec41f09c2f7a708","10.5281/zenodo.8433832","SOURCE_LOCKED_EQUAL_INFORMATION_DESCRIPTIVE_NOT_HISTORY_MEDIATION"):
        raise ValueError("low-order lineage")
    if data["audit"]["raw_rows"]!=2377 or data["audit"]["strain_subset_cells"]!=448 or max(data["max_author_point_difference"].values())>1e-6:
        raise ValueError("low-order source grain")
    rep=data["replicate_label_audit"]
    if (rep["n_system_rep_label_groups"],rep["duplicate_system_rep_strain_rows"],rep["incomplete_or_duplicated_roster_groups"],rep["physical_flask_pairing_independently_verified"])!=(692,0,0,False):
        raise ValueError("replicate grain")
    targets=data["targets"]
    if len(targets)!=64 or len({x["system"] for x in targets})!=64 or Counter(x["size"] for x in targets)!={4:35,5:21,6:7,7:1}:
        raise ValueError("low-order targets")
    for x in targets:
        if set(x["errors"])!=set(MODELS):raise ValueError("low-order models")
        for m,e in x["errors"].items():
            if set(e)!=set(ENDPOINTS) or any(not math.isfinite(v) or v<0 for v in e.values()):
                raise ValueError("low-order endpoint")
            if not close(e["named_TV"],e["rank_TV"]+e["identity_assignment_TV"]):
                raise ValueError("TV partition")
    mean={(m,e):statistics.mean(x["errors"][m][e] for x in targets) for m in MODELS for e in ENDPOINTS}
    sm=data["summary"];lookup={(x["baseline"],x["enhanced"]):x for x in sm["comparisons"]}
    eql={(x["baseline"],x["enhanced"]):x for x in data["equal_size_weighted_comparisons"]}
    if set(lookup)!=set(CONTRASTS) or set(eql)!=set(CONTRASTS):raise ValueError("low-order contrast coverage")
    models=read(source/MODEL_FILE,M,54);contrasts=read(source/CONTRAST_FILE,C,48)
    if {(x["model"],x["endpoint"]) for x in models}!=set(mean):
        raise ValueError("model coverage")
    if {(x["baseline"],x["enhanced"],x["endpoint"]) for x in contrasts}!={(b,c,e) for b,c in CONTRASTS for e in ENDPOINTS}:
        raise ValueError("contrast coverage")
    if conn.execute("SELECT COUNT(*) FROM effects WHERE panel_id=?",(PANEL,)).fetchone()[0]!=0:
        raise ValueError("finite predictions promoted to effects")
    for x in models+contrasts:
        if x["panel_id"]!=PANEL or x["source_json"]!=SRC or x["source_run"]!=RUN or int(x["n_target_types"])!=64:
            raise ValueError("low-order row lineage")
    mv=[];cv=[]
    for x in models:
        m,e=x["model"],x["endpoint"];v=mean[m,e]
        if not close(x["mean_error"],v) or not close(sm["mean_errors"][m][e],v):
            raise ValueError("model mean mismatch")
        mv.append((PANEL,m,e,v,64,SRC,RUN))
    for x in contrasts:
        b,c,e=x["baseline"],x["enhanced"],x["endpoint"]
        ds=[t["errors"][b][e]-t["errors"][c][e] for t in targets]
        d=statistics.mean(ds);eq=statistics.mean(statistics.mean(v for v,t in zip(ds,targets) if t["size"]==k) for k in (4,5,6,7))
        rel=d/mean[b,e] if mean[b,e]>1e-12 else None
        counts=(sum(v>1e-12 for v in ds),sum(v<-1e-12 for v in ds),sum(abs(v)<=1e-12 for v in ds))
        q=lookup[b,c]
        if not close(x["absolute_improvement"],d) or not close(q["absolute_improvement"][e],d) or not close(x["equal_size_weighted_improvement"],eq) or not close(eql[b,c]["absolute_improvement"][e],eq):
            raise ValueError("contrast mismatch")
        if rel is None:
            if x["relative_gain"]!="" or q["relative_gain"][e] is not None:raise ValueError("zero baseline")
        elif not close(x["relative_gain"],rel) or not close(q["relative_gain"][e],rel):
            raise ValueError("gain mismatch")
        if tuple(int(x[k]) for k in ("improved","worsened","tied"))!=counts or tuple(q["counts"][e][k] for k in ("improved","worsened","tied"))!=counts:
            raise ValueError("target counts")
        cv.append((PANEL,b,c,e,d,rel,eq,*counts,64,SRC,RUN))
    conn.execute("""CREATE TABLE low_order_model_errors(
      panel_id TEXT NOT NULL REFERENCES panels(panel_id),model TEXT NOT NULL,endpoint TEXT NOT NULL,
      mean_error REAL NOT NULL CHECK(mean_error>=0),n_target_types INTEGER NOT NULL,
      source_json TEXT NOT NULL,source_run TEXT NOT NULL,PRIMARY KEY(panel_id,model,endpoint))""")
    conn.executemany("INSERT INTO low_order_model_errors VALUES("+",".join("?" for _ in M)+")",mv)
    conn.execute("""CREATE TABLE low_order_model_contrasts(
      panel_id TEXT NOT NULL REFERENCES panels(panel_id),baseline TEXT NOT NULL,enhanced TEXT NOT NULL,
      endpoint TEXT NOT NULL,absolute_improvement REAL NOT NULL,relative_gain REAL,
      equal_size_weighted_improvement REAL NOT NULL,improved INTEGER NOT NULL,worsened INTEGER NOT NULL,
      tied INTEGER NOT NULL,n_target_types INTEGER NOT NULL,source_json TEXT NOT NULL,source_run TEXT NOT NULL,
      PRIMARY KEY(panel_id,baseline,enhanced,endpoint),CHECK(improved+worsened+tied=n_target_types))""")
    conn.executemany("INSERT INTO low_order_model_contrasts VALUES("+",".join("?" for _ in C)+")",cv)
    return {"low_order_model_error_rows":54,"low_order_contrast_rows":48,
      "low_order_experiment_families":1,"low_order_fixed_target_types":64,
      "low_order_meta_eligible":False}
