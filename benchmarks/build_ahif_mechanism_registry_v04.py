"""Materialize separate factorial comparisons on validated mechanism registry v0.3."""
import argparse,csv,json,sqlite3,subprocess,sys
from pathlib import Path
BASE=Path("database/ahif_mechanism_v0")
VAL=Path("validation/ahif_v1/mass_direction_transfer_summary_20261126.json")
RUN=45659798054
KEYS=("OO","RO","OR","RR","old_mass_TV","recent_mass_TV","mass_first","mass_second","mass_symmetric","sign_first","sign_second","sign_symmetric","total_gain","interaction","hybrid_over_lag","hybrid_over_old","mass_only_gain")
COLS=("source","comparison","panels","n_common","n_total",*KEYS,"meta_eligible","source_run")
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    subprocess.run([sys.executable,"-m","benchmarks.build_ahif_mechanism_registry_v03","--out",str(a.out)],check=True)
    z=json.loads(VAL.read_text())
    assert z["source_run"].endswith("/"+str(RUN)) and z["source_families"]==6 and z["panels"]==33
    assert z["saine"]["status"]=="STOP_ONLY_TWO_POST_INTRODUCTION_DATES"
    src=BASE/"mass_direction_transfer_effects_20261126.csv"
    with src.open() as f:
        rd=csv.DictReader(f);assert tuple(rd.fieldnames)==COLS;rows=list(rd)
    assert len(rows)==25
    expected={(s,c) for s,x in z["results"].items() for c,r in x["comparisons"].items() if r["status"]=="SCORED_DESCRIPTIVE"}
    seen=set();records=[]
    for row in rows:
        s,c=row["source"],row["comparison"]
        assert (s,c) not in seen;seen.add((s,c))
        r=z["results"][s]["comparisons"][c];m=r["mean"]
        assert int(row["panels"])==r["panels"]
        assert int(row["n_common"])==r["n_common"]==r["n_total"]==int(row["n_total"])
        assert row["meta_eligible"]=="no" and int(row["source_run"])==RUN
        assert all(abs(float(row[k])-m[k])<1e-12 for k in KEYS)
        assert abs(m["mass_symmetric"]+m["sign_symmetric"]-m["total_gain"])<1e-10
        assert abs(m["OO"]-m["RR"]-m["total_gain"])<1e-10
        records.append(tuple(row[k] for k in COLS))
    assert seen==expected
    con=sqlite3.connect(a.out/"ahif_mechanism_v0.sqlite")
    con.execute("PRAGMA foreign_keys=ON")
    con.execute("""CREATE TABLE mass_direction_transfer_effects(
      source TEXT REFERENCES studies(source), comparison TEXT, panels INTEGER, n_common INTEGER,n_total INTEGER,
      """+",".join('"'+k+'" REAL NOT NULL' for k in KEYS)+""",
      meta_eligible TEXT CHECK(meta_eligible='no'),source_run INTEGER,
      PRIMARY KEY(source,comparison))""")
    con.executemany("INSERT INTO mass_direction_transfer_effects VALUES ("+",".join("?" for _ in COLS)+")",records)
    con.execute("CREATE VIEW equal_one_mass_direction AS SELECT s.domain,s.measurement_lane,m.* FROM mass_direction_transfer_effects m JOIN studies s USING(source) WHERE comparison='one_observation'")
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0]=="ok"
    assert con.execute("PRAGMA foreign_key_check").fetchall()==[]
    n=con.execute("SELECT COUNT(*) FROM mass_direction_transfer_effects").fetchone()[0]
    p=con.execute("""SELECT COUNT(*),SUM(panels),SUM(n_common),SUM(mass_symmetric>0),SUM(mass_symmetric>sign_symmetric),SUM(total_gain>0),SUM(hybrid_over_lag>0 AND hybrid_over_old>0) FROM equal_one_mass_direction""").fetchone()
    assert n==25 and p==(6,33,67072,6,6,6,0),p
    con.close()
    (a.out/"mass_direction_transfer_effects.csv").write_bytes(src.read_bytes())
    audit=json.loads((a.out/"registry_audit.json").read_text())
    audit["version"]="0.4.0";audit["counts"]["mass_direction_transfer_effects"]=25
    audit["equal_one_mass_direction"]={"sources":6,"panels":33,"dependent_cells":67072,"mass_positive":6,"mass_exceeds_sign":6,"total_positive":6,"hybrid_beats_both":0,"meta_eligible":0,"source_run":RUN,"saine":"STOP_ONLY_TWO_POST_INTRODUCTION_DATES"}
    (a.out/"registry_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print("AHIF_MECHANISM_V04_RESULT",json.dumps(audit),flush=True)
if __name__=="__main__":main()
