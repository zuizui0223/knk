"""Validate Saine's complete source-linked outcome families, not eligibility tests."""
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
SAINE_FILE="saine_test_metadata.csv"
SRC="validation/ahif_v1/saine2025_nontarget_succession_v1.json"
COLS=("effect_id","lane","year","endpoint","p_holm","n_logs","computed_features",
      "logical_features","n_permutations","source_json")
SIZES={"ten_group":4,"any_introduction":4,"identity_among_introduced":4,
       "ten_group_scalar_secondary":6,"baseline_placebo":2}

def close(a,b):
    return math.isfinite(float(a)) and math.isclose(float(a),float(b),rel_tol=1e-10,abs_tol=1e-12)

def insert_saine(conn,source:Path):
    with (source/SAINE_FILE).open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=COLS:raise ValueError("Saine schema mismatch")
        rows=list(reader)
    if len(rows)!=20 or len({r["effect_id"] for r in rows})!=20:raise ValueError("Saine row count")
    src=json.loads((source.resolve().parents[1]/SRC).read_text())
    if (src["status"],src["source_physical_logs"],src["complete_triplets"],src["source_nontarget_otus"],src["source_qc_active_nontarget_otus"])!=("SOURCE_LOCKED_CONDITIONAL_COMPLETE_LOG_ANALYSIS_NOT_MEDIATION",460,446,4653,4646):
        raise ValueError("Saine source status or grain")
    if src["source_commit"]!="eb3f43a3f3b242705b192b2bd876940f34f90a50" or src["source_run"]!="https://github.com/zuizui0223/knk/actions/runs/45385872568":
        raise ValueError("Saine source run mismatch")
    if len(src["excluded_target_otus"])!=9 or len(src["eligibility_gate"])!=2:
        raise ValueError("Saine target or eligibility audit")
    if any(x["p_holm_family"]<=.05 for x in src["eligibility_gate"]):
        raise ValueError("Saine conditional eligibility gate failed")
    families={**src["primary_and_distinct_assignment_lanes"],
              "ten_group_scalar_secondary":src["ten_group_scalar_secondary"],
              "baseline_placebo":src["baseline_placebo"]}
    if {k:len(v) for k,v in families.items()}!=SIZES:raise ValueError("Saine incomplete test families")
    expected={(k,x["year"],x["endpoint"]):x for k,rs in families.items() for x in rs}
    effects={x[0]:x[1:] for x in conn.execute(
        "SELECT effect_id,estimate,p_value,estimand,evidence_url FROM effects WHERE panel_id='P_SAINE2025'")}
    if set(effects)!={r["effect_id"] for r in rows}:raise ValueError("Saine effect join")
    vals=[];seen=set();groups=defaultdict(list)
    for r in rows:
        k=(r["lane"],r["year"],r["endpoint"])
        if k not in expected or k in seen or r["source_json"]!=SRC:raise ValueError("Saine endpoint provenance")
        seen.add(k);x=expected[k];e,p,est,url=effects[r["effect_id"]]
        if est!="E_null_adj_R2" or url!=src["source_run"] or not close(e,x["E"]) or not close(p,x["p"]) or not close(r["p_holm"],x["p_holm_family"]):
            raise ValueError("Saine numeric source mismatch")
        for a,b in (("n_logs","n"),("computed_features","features"),("n_permutations","permutations")):
            if int(r[a])!=x[b]:raise ValueError("Saine endpoint grain")
        if int(r["logical_features"])!=x.get("logical_features",x["features"]):raise ValueError("Saine logical features")
        groups[r["lane"]].append((p,float(r["p_holm"])))
        vals.append([float(r[c]) if c=="p_holm" else int(r[c]) if c in ("n_logs","computed_features","logical_features","n_permutations") else r[c] for c in COLS])
    for k,g in groups.items():
        last=0.
        for i,(p,a) in enumerate(sorted(g,key=lambda z:z[0])):
            last=max(last,min(1.,p*(len(g)-i)))
            if not close(last,a):raise ValueError("Saine Holm mismatch")
    conn.execute("""CREATE TABLE saine_test_metadata(
      effect_id TEXT PRIMARY KEY REFERENCES effects(effect_id),lane TEXT NOT NULL,
      year TEXT NOT NULL,endpoint TEXT NOT NULL,p_holm REAL NOT NULL CHECK(p_holm BETWEEN 0 AND 1),
      n_logs INTEGER NOT NULL,computed_features INTEGER NOT NULL,logical_features INTEGER NOT NULL,
      n_permutations INTEGER NOT NULL,source_json TEXT NOT NULL)""")
    conn.executemany("INSERT INTO saine_test_metadata VALUES("+",".join("?" for _ in COLS)+")",vals)
    conn.execute("""CREATE VIEW v_saine_tests AS
       SELECT m.*,e.estimate,e.p_value,e.estimand,s.publication_cluster
       FROM saine_test_metadata m JOIN effects e USING(effect_id)
       JOIN panels p USING(panel_id) JOIN studies s USING(study_id)""")
    return {"saine_test_rows":20,"saine_ecological_postintroduction_rows":18,
       "saine_baseline_placebo_rows":2,"saine_eligibility_tests_not_effects":2,
       "saine_source_logs":460,"saine_complete_logs":446,
       "saine_publication_families":1,"saine_meta_eligible":False}
