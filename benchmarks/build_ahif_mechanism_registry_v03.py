"""Build mechanism registry v0.3 from frozen v0.2 and verified replication-scale receipts."""
import argparse,csv,json,sqlite3,subprocess,sys
from pathlib import Path
BASE=Path("database/ahif_mechanism_v0")
VAL=Path("validation/ahif_v1/replication_scale_summary_20261125.json")
EXPECTED_RUN=45589332743
COLS=("source","n_train","n_held","panels","dependent_cells","G","H","hedging_advantage","allocation_overlap","overlap_weighted_sign_agreement","concurrent_G","concurrent_H","concurrent_minus_lag_G","fixed_mass_concurrent_advantage","held_gross_same","meta_eligible","source_run")
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    subprocess.run([sys.executable,"-m","benchmarks.build_ahif_mechanism_registry","--out",str(a.out)],check=True)
    z=json.loads(VAL.read_text())
    assert z["status"]=="MATCHED_REPLICATION_DESCRIPTIVE_NOT_META" and z["source_run"].endswith("/"+str(EXPECTED_RUN))
    assert len(z["results"])==7
    with (BASE/"replication_scale_effects_20261125.csv").open() as f:
        reader=csv.DictReader(f);assert tuple(reader.fieldnames)==COLS;rows=list(reader)
    assert len(rows)==53
    keys=set();records=[]
    for row in rows:
        s=row["source"];n=int(row["n_train"]);h=int(row["n_held"])
        assert (s,n,h) not in keys;keys.add((s,n,h))
        r=z["results"][s];q=r["configs"][f"{n}+{h}"];m=q["common"]
        assert int(row["panels"])==r["n_panels"] and int(row["dependent_cells"])==q["n_common"]==q["n_total"]==q["n_defined"]
        assert row["meta_eligible"]=="no" and int(row["source_run"])==EXPECTED_RUN
        o=1-m["allocation_tv"];assert o>0
        expected={"G":m["G"],"H":m["H"],"hedging_advantage":m["H"]-m["G"],"allocation_overlap":o,
         "overlap_weighted_sign_agreement":.5+m["H"]/o,"concurrent_G":m["concurrent_G"],"concurrent_H":m["concurrent_H"],
         "concurrent_minus_lag_G":m["concurrent_minus_lag_G"],"fixed_mass_concurrent_advantage":m["fixed_mass_concurrent_advantage"],"held_gross_same":m["same"]}
        assert all(abs(float(row[k])-v)<1e-12 for k,v in expected.items())
        assert 0<=expected["overlap_weighted_sign_agreement"]<=1 and expected["hedging_advantage"]>=-1e-12
        records.append(tuple(row[k] for k in COLS))
    assert keys=={(s,*map(int,k.split("+"))) for s,r in z["results"].items() for k in r["configs"]}
    assert len({s for s,n,h in keys if n==h==2})==7
    con=sqlite3.connect(a.out/"ahif_mechanism_v0.sqlite");con.execute("PRAGMA foreign_keys=ON")
    con.execute("""CREATE TABLE replication_scale_effects (
     source TEXT REFERENCES studies(source), n_train INTEGER CHECK(n_train>0), n_held INTEGER CHECK(n_held>0),
     panels INTEGER CHECK(panels>0), dependent_cells INTEGER CHECK(dependent_cells>0),
     G REAL,H REAL,hedging_advantage REAL CHECK(hedging_advantage>=-1e-12),
     allocation_overlap REAL CHECK(allocation_overlap>0 AND allocation_overlap<=1),
     overlap_weighted_sign_agreement REAL CHECK(overlap_weighted_sign_agreement>=0 AND overlap_weighted_sign_agreement<=1),
     concurrent_G REAL,concurrent_H REAL,concurrent_minus_lag_G REAL,fixed_mass_concurrent_advantage REAL,held_gross_same REAL,
     meta_eligible TEXT CHECK(meta_eligible='no'),source_run INTEGER,PRIMARY KEY(source,n_train,n_held))""")
    con.executemany("INSERT INTO replication_scale_effects VALUES ("+",".join("?" for _ in COLS)+")",records)
    con.execute("""CREATE VIEW matched_two_plus_two AS
      SELECT s.domain,s.measurement_lane,s.cluster_unit,r.*
      FROM replication_scale_effects r JOIN studies s USING(source)
      WHERE n_train=2 AND n_held=2""")
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0]=="ok" and con.execute("PRAGMA foreign_key_check").fetchall()==[]
    n=con.execute("SELECT COUNT(*) FROM replication_scale_effects").fetchone()[0]
    primary=con.execute("SELECT COUNT(*),SUM(G>0),SUM(H>0),SUM(dependent_cells) FROM matched_two_plus_two").fetchone()
    assert n==53 and primary==(7,4,7,72192)
    assert con.execute("SELECT COUNT(*) FROM replication_scale_effects WHERE meta_eligible='yes'").fetchone()[0]==0
    con.close()
    (a.out/"replication_scale_effects.csv").write_bytes((BASE/"replication_scale_effects_20261125.csv").read_bytes())
    base=json.loads((a.out/"registry_audit.json").read_text())
    base["version"]="0.3.0";base["counts"]["replication_scale_effects"]=n
    base["matched_2plus2"]={"sources":primary[0],"G_positive":primary[1],"H_positive":primary[2],"dependent_cells":primary[3],"source_run":EXPECTED_RUN,"meta_eligible":0,"note":"Descriptive dependent source-size effects; Saine strata not prior site panels"}
    (a.out/"registry_audit.json").write_text(json.dumps(base,indent=2)+"\n")
    print("AHIF_MECHANISM_V03_RESULT",json.dumps(base),flush=True)
if __name__=="__main__":main()
