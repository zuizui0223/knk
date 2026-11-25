"""Separate source-linked held-origin prediction records from ecological effects."""
import csv,json,math
from pathlib import Path
from collections import defaultdict
PRED_FILE="prediction_outcomes.csv"
FOLD_FILE="prediction_folds.csv"
PRED=("prediction_id","panel_id","estimand_lane","training_window","endpoint",
      "gain_fraction","baseline_loss","enhanced_loss","p_approx","p_holm","n_units",
      "n_held_locations","n_resamples","source_json","source_run","inference_scope")
FOLD=("prediction_id","held_location","n_units_or_unknown","baseline_loss","enhanced_loss",
      "gain_fraction","shrinkage_or_unknown")
SAINE="validation/ahif_v1/saine2025_depth_biology_v1.json"
REPLAY="validation/ahif_v1/replay2025_state_v1.json"

def close(a,b):
    return math.isfinite(float(a)) and math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-11)

def read(path,cols,n):
    with path.open(encoding="utf-8-sig",newline="") as f:
        reader=csv.DictReader(f)
        if tuple(reader.fieldnames or ())!=cols:raise ValueError("prediction schema")
        rows=list(reader)
    if len(rows)!=n or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("prediction row count or truncated row")
    return rows

def insert_predictions(conn,source:Path):
    root=source.resolve().parents[1]
    s=json.loads((root/SAINE).read_text());r=json.loads((root/REPLAY).read_text())
    if s["source_run"]!="https://github.com/zuizui0223/knk/actions/runs/45389204461" or r["source_run"]!="https://github.com/zuizui0223/knk/actions/runs/45427172207":
        raise ValueError("prediction provenance")
    if s["complete_logs"]!=446 or r["n_parents"]!=275 or r["n_final_replicates"]!=1100 or r["n_locations"]!=14 or r["source_asv_matrix"]["n_asvs"]!=1454:
        raise ValueError("prediction unit counts")
    expected={};ef={}
    for x in s["nested_prediction"]:
        id="PRED_SAINE_"+str(x["train_year"])+"_"+x["endpoint"].upper()
        expected[id]=("P_SAINE2025","intervention_R_loss",str(x["train_year"])+"_to_2021",x["endpoint"],x["biology_gain_over_depth"],x["losses"]["depth"],x["losses"]["nested_biology"],x["p_approx"],x["p_holm4_approx"],446,5,999,SAINE,s["source_run"],"restricted_wild_nonzero_depth_context_not_causal_moderator")
        for a in x["folds"]:
            ef[id,a["site"]]=(None,a["depth"],a["nested_biology"],a["alpha"])
    for x in r["results"]:
        id="PRED_REPLAY_"+x["endpoint"].upper()
        expected[id]=("P_REPLAY2025","spontaneous_terminal_SSE","day0_to_day7_mean_of_four",x["endpoint"],x["held_location_gain"],x["depth_sse"],x["biology_sse"],x["p_approx"],x["p_holm2_approx"],275,14,999,REPLAY,r["source_run"],"restricted_wild_conditional_prediction_not_intervention")
        for a in x["folds"]:
            ef[id,a["location"]]=(a["n_parents"],a["depth_sse"],a["biology_sse"],None)
    rows=read(source/PRED_FILE,PRED,6);folds=read(source/FOLD_FILE,FOLD,48)
    if len(expected)!=6 or {x["prediction_id"] for x in rows}!=set(expected):
        raise ValueError("prediction endpoint coverage")
    if len(ef)!=48 or {(x["prediction_id"],x["held_location"]) for x in folds}!=set(ef):
        raise ValueError("prediction fold coverage")
    if conn.execute("SELECT COUNT(*) FROM effects WHERE panel_id='P_REPLAY2025'").fetchone()[0]:
        raise ValueError("spontaneous prediction promoted to intervention effects")
    ps=[];fs=[];group=defaultdict(list)
    for row in rows:
        id=row["prediction_id"]
        v=expected[id]
        for i,key in enumerate(PRED[1:]):
            if 4<=i<=8:
                if not close(row[key],v[i]):raise ValueError("prediction numeric mismatch")
            elif 9<=i<=11:
                if int(row[key])!=v[i]:raise ValueError("prediction grain mismatch")
            elif row[key]!=v[i]:raise ValueError("prediction provenance mismatch")
        if float(row["baseline_loss"])<=0 or float(row["enhanced_loss"])<0 or not close(row["gain_fraction"],1-float(row["enhanced_loss"])/float(row["baseline_loss"])):
            raise ValueError("prediction loss or gain")
        ps.append([id,*v])
    for a in folds:
        id=a["prediction_id"];site=a["held_location"];n,b,e,alpha=ef[id,site]
        if (a["n_units_or_unknown"]=="" if n is None else int(a["n_units_or_unknown"])==n) is False:
            raise ValueError("prediction fold count")
        if (a["shrinkage_or_unknown"]=="" if alpha is None else close(a["shrinkage_or_unknown"],alpha)) is False:
            raise ValueError("prediction fold shrinkage")
        if not close(a["baseline_loss"],b) or not close(a["enhanced_loss"],e) or not close(a["gain_fraction"],1-e/b):
            raise ValueError("prediction fold loss")
        fs.append([id,site,n,b,e,1-e/b,alpha]);group[id].append((b,e,n))
    for id,g in group.items():
        v=expected[id]
        if len(g)!=v[10] or not close(sum(x[0] for x in g),v[5]) or not close(sum(x[1] for x in g),v[6]):
            raise ValueError("prediction fold sums")
        if v[0]=="P_REPLAY2025" and sum(x[2] for x in g)!=275:raise ValueError("parent repeat inflation")
    conn.execute("""CREATE TABLE prediction_outcomes(
      prediction_id TEXT PRIMARY KEY,panel_id TEXT NOT NULL REFERENCES panels(panel_id),
      estimand_lane TEXT NOT NULL,training_window TEXT NOT NULL,endpoint TEXT NOT NULL,
      gain_fraction REAL NOT NULL,baseline_loss REAL NOT NULL,enhanced_loss REAL NOT NULL,
      p_approx REAL NOT NULL CHECK(p_approx BETWEEN 0 AND 1),
      p_holm REAL NOT NULL CHECK(p_holm BETWEEN 0 AND 1),n_units INTEGER NOT NULL,
      n_held_locations INTEGER NOT NULL,n_resamples INTEGER NOT NULL,
      source_json TEXT NOT NULL,source_run TEXT NOT NULL,inference_scope TEXT NOT NULL)""")
    conn.executemany("INSERT INTO prediction_outcomes VALUES("+",".join("?" for _ in PRED)+")",ps)
    conn.execute("""CREATE TABLE prediction_folds(
      prediction_id TEXT NOT NULL REFERENCES prediction_outcomes(prediction_id),held_location TEXT NOT NULL,
      n_units_or_unknown INTEGER,baseline_loss REAL NOT NULL,enhanced_loss REAL NOT NULL,
      gain_fraction REAL NOT NULL,shrinkage_or_unknown REAL,
      PRIMARY KEY(prediction_id,held_location))""")
    conn.executemany("INSERT INTO prediction_folds VALUES("+",".join("?" for _ in FOLD)+")",fs)
    return {"prediction_outcome_rows":6,"prediction_fold_rows":48,
            "prediction_source_experiment_families":2,"prediction_distinct_estimand_lanes":2,
            "prediction_only_source_panels":1,"prediction_meta_eligible":False}
