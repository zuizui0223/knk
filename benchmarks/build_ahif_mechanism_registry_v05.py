"""Append an external validation lane to v0.4 without altering development rows."""
import argparse,csv,json,sqlite3,subprocess,sys
from pathlib import Path
BASE=Path("database/ahif_mechanism_v0")
VAL=Path("validation/ahif_v1/villarreal2025_external_summary_20261126.json")
RUN=45714536073
KEYS=("OO","RO","OR","RR","old_mass_TV","recent_mass_TV","mass_first","mass_second","mass_symmetric","sign_first","sign_second","sign_symmetric","total_gain","interaction","hybrid_over_lag","hybrid_over_old","mass_only_gain")
LAG=("G","H","O","Astar","J","concurrent_G","concurrent_H","fixed_mass_concurrent_advantage")
COLS=("source","panel","role","n_train","n_held","n_common","n_total",*KEYS,*LAG,"meta_eligible","source_run")
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    subprocess.run([sys.executable,"-m","benchmarks.build_ahif_mechanism_registry_v04","--out",str(a.out)],check=True)
    z=json.loads(VAL.read_text());audit=z["source_audit"]
    assert z["status"]=="INDEPENDENT_ZOOPLANKTON_EXTERNAL_VALIDATION_NOT_META" and z["source_run"].endswith("/"+str(RUN))
    assert (audit["n_tanks"],audit["n_linked_pairs"],audit["n_profiles"],audit["n_taxa"],audit["individuals"])==(40,20,200,63,88774)
    source=BASE/"external_villarreal_effects_20261126.csv"
    with source.open() as f:
        rd=csv.DictReader(f);assert tuple(rd.fieldnames)==COLS;rows=list(rd)
    assert len(rows)==6;seen=set();records=[]
    for r in rows:
        n=r["n_train"];panel=r["panel"];key=(panel,n);assert key not in seen;seen.add(key)
        q=z["by_panel"][panel][n];m=q["common"]["factorial"];y=q["common"]["one_lag"]
        assert n==r["n_held"] and r["source"]=="villarreal2025" and r["meta_eligible"]=="no"
        assert int(r["source_run"])==RUN and int(r["n_common"])==q["common"]["n"]==512 and int(r["n_total"])==q["n_total"]==512
        assert r["role"]==("secondary" if panel=="Ab-Aa" else "primary")
        assert all(abs(float(r[k])-m[k])<1e-12 for k in KEYS)
        assert all(abs(float(r[k])-y[k])<1e-12 for k in LAG)
        assert abs(m["mass_symmetric"]+m["sign_symmetric"]-m["total_gain"])<1e-10
        records.append(tuple(r[k] for k in COLS))
    assert seen=={(p,n) for p in ("Aa-A0","Ab-A0","Ab-Aa") for n in ("2","4")}
    con=sqlite3.connect(a.out/"ahif_mechanism_v0.sqlite");con.execute("PRAGMA foreign_keys=ON")
    con.execute("""CREATE TABLE external_validation_sources(
      source TEXT PRIMARY KEY, publication_doi TEXT UNIQUE NOT NULL, data_doi TEXT,
      domain TEXT, measurement_lane TEXT, cluster_unit TEXT, n_tanks INTEGER,n_clusters INTEGER,n_taxa INTEGER,
      source_receipt TEXT, meta_eligible TEXT CHECK(meta_eligible='no'))""")
    con.execute("INSERT INTO external_validation_sources VALUES (?,?,?,?,?,?,?,?,?,?,?)",("villarreal2025","10.1002/ecy.70256","10.6084/m9.figshare.28651985.v1","zooplankton","fixed_volume_count","linked_donor_recipient_pair",40,20,63,str(VAL),"no"))
    con.execute("""CREATE TABLE external_factorial_effects(
      source TEXT REFERENCES external_validation_sources(source), panel TEXT, role TEXT,
      n_train INTEGER,n_held INTEGER,n_common INTEGER,n_total INTEGER,
      """+",".join('"'+k+'" REAL NOT NULL' for k in KEYS+LAG)+""",
      meta_eligible TEXT CHECK(meta_eligible='no'),source_run INTEGER,PRIMARY KEY(source,panel,n_train,n_held))""")
    con.executemany("INSERT INTO external_factorial_effects VALUES ("+",".join("?" for _ in COLS)+")",records)
    con.execute("""CREATE VIEW external_primary_source_effects AS SELECT source,n_train,n_held,
      COUNT(*) AS panels,SUM(n_common) AS dependent_cells,AVG(mass_symmetric) AS mass_symmetric,
      AVG(sign_symmetric) AS sign_symmetric,AVG(total_gain) AS total_gain,AVG(hybrid_over_lag) AS hybrid_over_lag
      FROM external_factorial_effects WHERE role='primary' GROUP BY source,n_train,n_held""")
    for n in (2,4):
        v=con.execute("SELECT panels,dependent_cells,mass_symmetric,sign_symmetric,total_gain,hybrid_over_lag FROM external_primary_source_effects WHERE n_train=? AND n_held=?",(n,n)).fetchone()
        ref=z["primary"][str(n)]
        assert v[:2]==(2,1024)
        assert all(abs(x-ref["factorial"][k])<1e-12 for x,k in zip(v[2:],("mass_symmetric","sign_symmetric","total_gain","hybrid_over_lag")))
        assert v[2]>0 and v[2]>v[3] and v[5]<0
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0]=="ok" and con.execute("PRAGMA foreign_key_check").fetchall()==[]
    assert con.execute("SELECT COUNT(*) FROM external_validation_sources").fetchone()[0]==1
    assert con.execute("SELECT COUNT(*) FROM external_factorial_effects").fetchone()[0]==6
    assert con.execute("SELECT COUNT(*) FROM external_factorial_effects WHERE role='primary' AND mass_symmetric<0").fetchone()[0]==2
    con.close()
    (a.out/"external_villarreal_effects.csv").write_bytes(source.read_bytes())
    (a.out/"external_recruitment_screen.csv").write_bytes((BASE/"external_recruitment_screen_20261126.csv").read_bytes())
    audit=json.loads((a.out/"registry_audit.json").read_text());audit["version"]="0.5.0"
    audit["counts"]["external_validation_sources"]=1;audit["counts"]["external_factorial_effects"]=6
    audit["external_villarreal"]={"source_run":RUN,"primary_panels":2,"primary_2plus2_cells":1024,"mass_positive":True,"mass_exceeds_sign":True,"hybrid_beats_both":False,"negative_same_source_mass_both_sizes":True,"meta_eligible":0}
    (a.out/"registry_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    print("AHIF_MECHANISM_V05_RESULT",json.dumps(audit),flush=True)
if __name__=="__main__":main()
