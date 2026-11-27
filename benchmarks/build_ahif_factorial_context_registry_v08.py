"""Extend exact v0.7 with Shen factorial context; preserve separate response estimands."""
import argparse,csv,hashlib,json,sqlite3,zipfile
from pathlib import Path,PurePosixPath
BASE="16baafefe72853bd6880921dfb043b2e198b92d5b5b4d6815b0b7f80114f305c"
VAL=Path("validation/ahif_v1")
CSV=Path("database/ahif_mechanism_v0/shen_context_factorial_effects_20261127.csv")
SOURCE="shen2023";RUN=45940785495
def check(a,b):assert abs(float(a)-float(b))<1e-10,(a,b)
def main():
    p=argparse.ArgumentParser();p.add_argument("--base-zip",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    assert hashlib.sha256(a.base_zip.read_bytes()).hexdigest()==BASE
    a.out.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.base_zip) as z:
        assert len(z.infolist())<100 and sum(x.file_size for x in z.infolist())<1_000_000_000
        for x in z.infolist():
            q=PurePosixPath(x.filename);assert not q.is_absolute() and ".." not in q.parts
        z.extractall(a.out)
    old=json.loads((a.out/"finite_context_audit.json").read_text());ref=json.loads((VAL/"finite_context_registry_v0_7_receipt_20261127.json").read_text())
    assert all(old[k]==ref[k] for k in ("version","counts","sign_by_background","primary_identity_MAE","integrity","foreign_key_check"))
    z=json.loads((VAL/"shen2023_context_summary_20261127.json").read_text())
    assert z["original_scores_replayed"] and z["source_run"].endswith("/"+str(RUN))
    audit=z["audit"];assert audit["culture_n"]==24 and audit["source_rows"]==1322 and audit["focal_observed_rows"]==807
    assert len(audit["missing_entire_culture_days"])==172 and audit["zero_fills_common"]==52
    con=sqlite3.connect(a.out/"ahif_mechanism_v0.sqlite");con.execute("PRAGMA foreign_keys=ON")
    before={k:con.execute('SELECT COUNT(*) FROM "'+k+'"').fetchone()[0] for (k,) in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    con.executescript("""
CREATE TABLE factorial_context_sources(
 source TEXT PRIMARY KEY,publication_doi TEXT UNIQUE NOT NULL,data_doi TEXT,source_commit TEXT NOT NULL,
 measurement TEXT NOT NULL,design TEXT NOT NULL,n_cultures INTEGER NOT NULL,
 source_run INTEGER NOT NULL,inference_limit TEXT NOT NULL,meta_eligible TEXT CHECK(meta_eligible='no'));
CREATE TABLE factorial_context_lanes(
 source TEXT NOT NULL REFERENCES factorial_context_sources(source),lane TEXT NOT NULL,status TEXT NOT NULL,
 dates_json TEXT NOT NULL,n_dates INTEGER NOT NULL,role TEXT NOT NULL,duplicate_of TEXT,
 PRIMARY KEY(source,lane),FOREIGN KEY(source,duplicate_of) REFERENCES factorial_context_lanes(source,lane));
CREATE TABLE factorial_context_marginals(
 source TEXT NOT NULL,lane TEXT NOT NULL,focal TEXT NOT NULL,competitor TEXT NOT NULL,added_agent TEXT NOT NULL,
 baseline_CD_density REAL NOT NULL CHECK(baseline_CD_density>0),base_marginal REAL NOT NULL,background_marginal REAL NOT NULL,difference REAL NOT NULL,
 direction TEXT NOT NULL,direction_threshold_005 TEXT NOT NULL,
 PRIMARY KEY(source,lane,focal,competitor,added_agent),FOREIGN KEY(source,lane) REFERENCES factorial_context_lanes(source,lane),
 CHECK(abs(background_marginal-base_marginal-difference)<1e-10),
 CHECK(direction IN ('same','opposite','unresolved')),CHECK(direction_threshold_005 IN ('same','opposite','unresolved')));
CREATE TABLE factorial_context_forecasts(
 source TEXT NOT NULL,lane TEXT NOT NULL,focal TEXT NOT NULL,competitor TEXT NOT NULL,added_agent TEXT NOT NULL,
 n_common_folds INTEGER NOT NULL,n_total_folds INTEGER NOT NULL,no_change REAL NOT NULL,additive REAL NOT NULL,multiplicative REAL NOT NULL,
 PRIMARY KEY(source,lane,focal,competitor,added_agent),
 FOREIGN KEY(source,lane,focal,competitor,added_agent) REFERENCES factorial_context_marginals(source,lane,focal,competitor,added_agent));
CREATE VIEW factorial_context_summary AS
 SELECT source,lane,COUNT(*) AS dependent_marginals,SUM(direction='opposite') AS opposite,
 SUM(direction_threshold_005='opposite') AS opposite_005,SUM(direction_threshold_005='unresolved') AS unresolved_005
 FROM factorial_context_marginals GROUP BY source,lane;
""")
    con.execute("INSERT INTO factorial_context_sources VALUES (?,?,?,?,?,?,?,?,?,?)",(SOURCE,z["publication"],"10.5281/zenodo.7896204",z["source_commit"],"replicate_mean_individuals_per_mL","CD fixed focal pair × competitor none/P/S × predator assigned absent/present",24,RUN,"Finite marginal arm contrasts; day1 imbalance, common dates only, shared controls, no monoculture or direct interaction coefficient","no"))
    for lane,r in z["results"].items():
        role="primary" if lane=="primary_common" else "sensitivity"
        con.execute("INSERT INTO factorial_context_lanes VALUES (?,?,?,?,?,?,?)",(SOURCE,lane,r["status"],json.dumps(r["dates"]),len(r["dates"]),role,"primary_common" if lane=="early" else None))
    rows=list(csv.DictReader(CSV.open()));assert len(rows)==32
    for row in rows:
        lane=row["lane"];q=z["results"][lane];key=tuple(row[k] for k in ("focal","competitor","added_agent"))
        m=[r for r in q["marginals"] if tuple(r[k] for k in ("focal","competitor","added_agent"))==key]
        assert len(m)==1;m=m[0]
        assert row["source"]==SOURCE and row["meta_eligible"]=="no" and int(row["source_run"])==RUN and int(row["n_dates"])==len(q["dates"])
        for k in ("baseline_CD_density","base_marginal","background_marginal","difference"):check(row[k],m[k])
        for k in ("direction","direction_threshold_005"):assert row[k]==m[k]
        con.execute("INSERT INTO factorial_context_marginals VALUES (?,?,?,?,?,?,?,?,?,?,?)",(SOURCE,lane,*key,*[float(row[k]) for k in ("baseline_CD_density","base_marginal","background_marginal","difference")],row["direction"],row["direction_threshold_005"]))
    for lane,q in z["results"].items():
        f=q.get("forecast")
        if not f:continue
        assert f["n_common_folds"]==f["n_total_folds"]==512 and len(f["by_panel"])==8
        for r in f["by_panel"]:
            con.execute("INSERT INTO factorial_context_forecasts VALUES (?,?,?,?,?,?,?,?,?,?)",(SOURCE,lane,r["focal"],r["competitor"],r["added_agent"],512,512,*[r["MAE"][k] for k in ("no_change","additive","multiplicative")]))
        rec=con.execute("SELECT AVG(no_change),AVG(additive),AVG(multiplicative) FROM factorial_context_forecasts WHERE source=? AND lane=?",(SOURCE,lane)).fetchone()
        for v,k in zip(rec,("no_change","additive","multiplicative")):check(v,f["mean_MAE"][k])
    for _,lane,n,op,op05,un in con.execute("SELECT * FROM factorial_context_summary"):
        r=z["results"][lane];assert n==8 and op==r["counts"].get("opposite",0) and op05==r["counts_threshold_005"].get("opposite",0) and un==r["counts_threshold_005"].get("unresolved",0)
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0]=="ok" and con.execute("PRAGMA foreign_key_check").fetchall()==[]
    assert all(con.execute('SELECT COUNT(*) FROM "'+k+'"').fetchone()[0]==v for k,v in before.items())
    assert con.execute("SELECT COUNT(*) FROM opportunity_signed_interaction_effects").fetchone()[0]==0
    tables=["factorial_context_sources","factorial_context_lanes","factorial_context_marginals","factorial_context_forecasts"]
    counts={t:con.execute("SELECT COUNT(*) FROM "+t).fetchone()[0] for t in tables};assert list(counts.values())==[1,5,32,24],counts
    for t in tables+["factorial_context_summary"]:
        cur=con.execute("SELECT * FROM "+t)
        with (a.out/(t+".csv")).open("w",newline="") as f:
            wr=csv.writer(f);wr.writerow([r[0] for r in cur.description]);wr.writerows(cur)
    con.close()
    result={"status":"VALIDATED_FACTORIAL_CONTEXT_REGISTRY_NOT_META","version":"0.8.0","base_v07_sha256":BASE,"counts":counts,"source_run":RUN,"primary_marginals":8,"early_duplicate_marginals":8,"late_stop":True,"forecast_lanes":3,"prior_table_counts_unchanged":True,"integrity":"ok","foreign_key_check":"ok","signed_coefficient_rows":0,"meta_eligible":0}
    (a.out/"factorial_context_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    a0=json.loads((a.out/"registry_audit.json").read_text());a0["version"]="0.8.0";a0["counts"].update(counts);a0["factorial_context"]=result
    (a.out/"registry_audit.json").write_text(json.dumps(a0,indent=2)+"\n")
    print("AHIF_FACTORIAL_CONTEXT_V08_RESULT",json.dumps(result),flush=True)
if __name__=="__main__":main()
