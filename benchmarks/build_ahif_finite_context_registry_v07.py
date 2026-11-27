"""Append observed finite relative-density context records to exact v0.6 registry."""
import argparse,csv,hashlib,json,sqlite3,subprocess,sys,zipfile
from pathlib import Path,PurePosixPath
VAL=Path("validation/ahif_v1")
SOURCE="castledine2024";RUN=45939572589
BASE_DIGEST="d50d6d93d28a7c25e6b7671150266f6fef8e00b11c7f9b35c845c08cd0e7def6"
def readcsv(p):return list(csv.DictReader(p.open()))
def close(a,b):
    assert abs(float(a)-float(b))<1e-10,(a,b)
def main():
    p=argparse.ArgumentParser();p.add_argument("--base-zip",type=Path,required=True);p.add_argument("--out",type=Path,required=True);p.add_argument("--scratch",type=Path,required=True);a=p.parse_args()
    assert hashlib.sha256(a.base_zip.read_bytes()).hexdigest()==BASE_DIGEST
    a.out.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(a.base_zip) as z:
        assert len(z.infolist())<100 and sum(f.file_size for f in z.infolist())<1_000_000_000
        for f in z.infolist():
            p=PurePosixPath(f.filename);assert not p.is_absolute() and ".." not in p.parts
        z.extractall(a.out)
    base=json.loads((a.out/"opportunity_audit.json").read_text())
    old=json.loads((VAL/"opportunity_registry_v0_6_receipt_20261126.json").read_text())
    assert all(base[k]==old[k] for k in ("status","version","counts","coverage","assay_summary","source_receipts","integrity","foreign_key_check"))
    a.scratch.mkdir(parents=True,exist_ok=True);res=a.scratch/"castledine"
    subprocess.run([sys.executable,"-m","benchmarks.score_ahif_castledine2024_context","--out",str(res)],check=True)
    z=json.loads((res/"summary.json").read_text());ref=json.loads((VAL/"castledine2024_context_summary_20261127.json").read_text())
    assert z=={k:ref[k] for k in z},"Frozen complete score replay"
    ends=readcsv(res/"finite_density_endpoints.csv");ctx=readcsv(res/"context_marginal_responses.csv");full=readcsv(res/"full_endpoint_forecasts.csv")
    con=sqlite3.connect(a.out/"ahif_mechanism_v0.sqlite");con.execute("PRAGMA foreign_keys=ON")
    before={k:con.execute('SELECT COUNT(*) FROM "'+k+'"').fetchone()[0] for (k,) in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    con.executescript("""
CREATE TABLE finite_response_sources (
 source TEXT PRIMARY KEY, publication_doi TEXT UNIQUE NOT NULL, source_commit TEXT NOT NULL,
 measurement TEXT NOT NULL, design TEXT NOT NULL, inference_limit TEXT NOT NULL, source_run INTEGER NOT NULL,
 meta_eligible TEXT NOT NULL CHECK(meta_eligible='no'));
CREATE TABLE finite_density_endpoints (
 source TEXT NOT NULL REFERENCES finite_response_sources(source), focal TEXT NOT NULL, partners TEXT NOT NULL,
 n_partners INTEGER NOT NULL CHECK(n_partners>=0), relative_density REAL NOT NULL CHECK(relative_density>=0),
 finite_response_vs_mono REAL NOT NULL, evidence_class TEXT NOT NULL CHECK(evidence_class IN ('author_processed_CFU_relative_density','definition_not_observation')),
 PRIMARY KEY(source,focal,partners), CHECK(abs(relative_density-1-finite_response_vs_mono)<1e-10),
 CHECK(n_partners=length(partners)), CHECK((n_partners=0)=(evidence_class='definition_not_observation')));
CREATE TABLE finite_context_marginals (
 source TEXT NOT NULL, focal TEXT NOT NULL, added TEXT NOT NULL, background TEXT NOT NULL, target_partners TEXT NOT NULL,
 background_size INTEGER NOT NULL CHECK(background_size>=0), pair_marginal REAL NOT NULL,
 background_density REAL NOT NULL, target_density REAL NOT NULL, context_marginal REAL NOT NULL,
 zero_forecast REAL NOT NULL, additive_forecast REAL NOT NULL, multiplicative_forecast REAL NOT NULL, anonymous_forecast REAL NOT NULL,
 PRIMARY KEY(source,focal,added,background),
 FOREIGN KEY(source,focal,added) REFERENCES finite_density_endpoints(source,focal,partners),
 FOREIGN KEY(source,focal,background) REFERENCES finite_density_endpoints(source,focal,partners),
 FOREIGN KEY(source,focal,target_partners) REFERENCES finite_density_endpoints(source,focal,partners),
 CHECK(abs(target_density-background_density-context_marginal)<1e-10),
 CHECK(length(background)=background_size AND length(target_partners)=background_size+1),
 CHECK(instr(background,added)=0 AND instr(target_partners,added)>0 AND focal<>added));
CREATE TABLE finite_full_endpoint_forecasts (
 source TEXT NOT NULL, focal TEXT NOT NULL, partners TEXT NOT NULL,n_partners INTEGER NOT NULL,
 target REAL NOT NULL,no_response REAL NOT NULL, actual_product REAL NOT NULL,anonymous_product REAL NOT NULL,
 PRIMARY KEY(source,focal,partners),
 FOREIGN KEY(source,focal,partners) REFERENCES finite_density_endpoints(source,focal,partners));
CREATE VIEW finite_context_sign_summary AS
 SELECT source,background_size,COUNT(*) AS dependent_contrasts,
 SUM(CASE WHEN pair_marginal*context_marginal>0 THEN 1 ELSE 0 END) AS same_sign,
 SUM(CASE WHEN pair_marginal*context_marginal<0 THEN 1 ELSE 0 END) AS opposite_sign,
 SUM(CASE WHEN abs(pair_marginal)<=.05 OR abs(context_marginal)<=.05 THEN 1 ELSE 0 END) AS unresolved_threshold_005,
 SUM(CASE WHEN abs(pair_marginal)>.05 AND abs(context_marginal)>.05 AND pair_marginal*context_marginal<0 THEN 1 ELSE 0 END) AS opposite_threshold_005
 FROM finite_context_marginals WHERE background_size>0 GROUP BY source,background_size;
CREATE VIEW finite_context_identity_primary AS
 SELECT source,focal,added,COUNT(*) AS dependent_rows,AVG(abs(target_density-zero_forecast)) AS zero_mae,
 AVG(abs(target_density-additive_forecast)) AS additive_mae,AVG(abs(target_density-multiplicative_forecast)) AS multiplicative_mae,
 AVG(abs(target_density-anonymous_forecast)) AS anonymous_mae
 FROM finite_context_marginals WHERE background_size IN (1,2) GROUP BY source,focal,added;
""")
    con.execute("INSERT INTO finite_response_sources VALUES (?,?,?,?,?,?,?,?)",(SOURCE,z["publication"],z["source_commit"],"author_processed_CFU_relative_to_monoculture","complete additive five-species subset series; nominal focal inoculum fixed, total inoculum varies","finite net density contrasts; no raw replicate uncertainty, no direct per-capita coefficients",RUN,"no"))
    assert len(ends)==75
    for r in ends:
        assert r["meta_eligible"]=="no"
        con.execute("INSERT INTO finite_density_endpoints VALUES (?,?,?,?,?,?,?)",(SOURCE,r["focal"],r["partners"],int(r["n_partners"]),float(r["relative_density"]),float(r["finite_response_vs_mono"]),r["measurement"]))
    for i in "aopsv":con.execute("INSERT INTO finite_density_endpoints VALUES (?,?,?,?,?,?,?)",(SOURCE,i,"",0,1.,0.,"definition_not_observation"))
    def y(i,s):return con.execute("SELECT relative_density FROM finite_density_endpoints WHERE source=? AND focal=? AND partners=?",(SOURCE,i,s)).fetchone()[0]
    for r in ctx:
        i,j,s=r["focal"],r["added"],r["background"];t="".join(sorted(s+j))
        close(y(i,s),r["background_density"]);close(y(i,t),r["target"]);close(y(i,j)-1,r["pair_marginal"])
        con.execute("INSERT INTO finite_context_marginals VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",(SOURCE,i,j,s,t,int(r["background_size"]),*[float(r[k]) for k in ("pair_marginal","background_density","target","context_marginal","zero","additive","multiplicative","anonymous")]))
    for r in full:
        close(y(r["focal"],r["partners"]),r["target"])
        con.execute("INSERT INTO finite_full_endpoint_forecasts VALUES (?,?,?,?,?,?,?,?)",(SOURCE,r["focal"],r["partners"],int(r["n_partners"]),*[float(r[k]) for k in ("target","no_response","actual_product","anonymous_product")]))
    assert len(ctx)==160 and len(full)==55
    sign=con.execute("SELECT background_size,dependent_contrasts,same_sign,opposite_sign,unresolved_threshold_005,opposite_threshold_005 FROM finite_context_sign_summary ORDER BY background_size").fetchall()
    for n,total,same,oppo,unresolved,opp05 in sign:
        r=z["sign"]["0"]["by_background_size"][str(n)];q=z["sign"]["0.05"]["by_background_size"][str(n)]
        assert (total,same,oppo,unresolved,opp05)==(r["n_rows"],r["same"],r["opposite"],q["unresolved"],q["opposite"])
    primary=con.execute("SELECT COUNT(*),SUM(dependent_rows),AVG(zero_mae),AVG(additive_mae),AVG(multiplicative_mae),AVG(anonymous_mae) FROM finite_context_identity_primary").fetchone()
    assert primary[:2]==(20,120)
    for v,k in zip(primary[2:],("zero","additive","multiplicative","anonymous")):close(v,z["marginal_forecasts"]["primary_identity"]["MAE"][k])
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0]=="ok" and con.execute("PRAGMA foreign_key_check").fetchall()==[]
    assert all(con.execute('SELECT COUNT(*) FROM "'+k+'"').fetchone()[0]==v for k,v in before.items())
    assert con.execute("SELECT COUNT(*) FROM opportunity_signed_interaction_effects").fetchone()[0]==0
    newtables=["finite_response_sources","finite_density_endpoints","finite_context_marginals","finite_full_endpoint_forecasts"]
    counts={t:con.execute("SELECT COUNT(*) FROM "+t).fetchone()[0] for t in newtables}
    assert list(counts.values())==[1,80,160,55]
    for t in newtables+["finite_context_sign_summary","finite_context_identity_primary"]:
        cur=con.execute("SELECT * FROM "+t)
        with (a.out/(t+".csv")).open("w",newline="") as f:
            wr=csv.writer(f);wr.writerow([r[0] for r in cur.description]);wr.writerows(cur)
    con.close()
    audit=json.loads((a.out/"registry_audit.json").read_text());audit["version"]="0.7.0";audit["counts"].update(counts)
    audit["finite_context"]={"source":SOURCE,"source_run":RUN,"observed_endpoint_rows":75,"definition_only_rows":5,"directed_pair_alone":20,"background_rows":140,"opposite_sign":39,"primary_identity_rows":120,"meta_eligible":0,"signed_coefficient_rows":0}
    (a.out/"registry_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    result={"status":"VALIDATED_FINITE_CONTEXT_REGISTRY_NOT_META","version":"0.7.0","base_v06_sha256":BASE_DIGEST,"prior_table_counts_unchanged":True,"counts":counts,"observed_endpoint_rows":75,"definition_only_rows":5,"sign_by_background":sign,"primary_identity_MAE":dict(zip(("zero","additive","multiplicative","anonymous"),primary[2:])),"source_run":RUN,"integrity":"ok","foreign_key_check":"ok","meta_eligible":0,"signed_coefficient_rows":0}
    (a.out/"finite_context_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print("AHIF_FINITE_CONTEXT_V07_RESULT",json.dumps(result),flush=True)
if __name__=="__main__":main()
