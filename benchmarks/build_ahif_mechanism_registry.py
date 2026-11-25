"""Materialize separate AHIF mechanism registry, no effect pooling."""
import argparse,csv,json,sqlite3
from pathlib import Path
BASE=Path("database/ahif_mechanism_v0")
VAL=Path("validation/ahif_v1")
LEDGER=Path("database/ahif_v1/one_lag_study_effects_exploratory_20261125.csv")
def readjson(p):return json.loads(p.read_text())
def csvwrite(path,rows):
    assert rows
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    c=readjson(BASE/"catalog.json")
    studies=[dict(zip(c["study_columns"],v)) for v in c["studies"]]
    S={r["source"]:r for r in studies};assert len(S)==len(studies)==7
    previous=readjson(VAL/"one_lag_fungal_summary_20261125.json")
    transfer=readjson(VAL/"signal_transfer_summary_20261125.json")
    assert len(transfer["panels"])==27
    panels=[]
    for s,key,n in transfer["panels"]:
        assert s in S and len(S[s]["observed_times"])==n
        exposed=key;reference="source-defined"
        context=""
        if s=="pond":
            exposed,reference=key.split("-")
        elif s=="kansas":exposed,reference="Seeds","Control"
        elif s=="poem":context,exposed,reference=key.split("_")
        elif s=="field":context,exposed,reference=key.split("_")
        elif s=="mancuso":
            rate,freq=key.split(":");context="mean_dilution="+rate
            exposed=freq+"_pulses/day";reference="600_continuous_code"
        panels.append(dict(source=s,panel=key,context=context,treated=exposed,reference=reference,n_observed_times=n,primary_targets=n-2,evidence="validation/ahif_v1/signal_transfer_summary_20261125.json"))
    for site in sorted(previous["saine"]["binary"]["coverage"]):
        panels.append(dict(source="saine",panel=site,context="natural/felled equal-stratum contrasts",treated="any introduction",reference="sterile control",n_observed_times=3,primary_targets=1,evidence="validation/ahif_v1/one_lag_fungal_summary_20261125.json"))
    hy=readjson(VAL/"hayashi2026_direction_summary_20261125.json")
    assert hy["status"]=="DESCRIPTIVE_HELD_WELL_FORECAST_POSITION_ALIASED_NOT_META" and hy["n_complete_wells"]==760
    hp=hy["otu99_full_half"]["primary"]
    for key in sorted(hp["by_panel"]):
        source,density=key.split(":")
        panels.append(dict(source="hayashi",panel=key,context=source+"; density-position aliased",treated="founder fraction "+density,reference="undiluted founder",n_observed_times=4,primary_targets=2,evidence="validation/ahif_v1/hayashi2026_direction_summary_20261125.json"))
    assert len(panels)==38 and len({(r["source"],r["panel"]) for r in panels})==38
    for s,row in S.items():
        assert sum(r["source"]==s for r in panels)==row["n_panels"]
        row["observed_times"]=json.dumps(row["observed_times"])
    with LEDGER.open() as f: effects=list(csv.DictReader(f))
    assert len(effects)==6 and {r["source"] for r in effects}==set(S)-{"hayashi"}
    hymean=hp["equal_panel"]
    effects.append(dict(source="hayashi",study="Hayashi et al. 2026",measurement_lane="relative_read_composition",independent_unit_candidates="760",panels="6",dependent_cells=str(hp["n_defined"]),gain=str(hymean["gain"]),same=str(hymean["same"]),opposite=str(hymean["opposite"]),abstain=str(hymean["abstain"]),allocation_tv=str(hymean["allocation_tv"]),meta_eligible="no",note="760 repeated wells; treatment aliases plate position; large half-splits; no calibrated SE"))
    assert len(effects)==7 and {r["source"] for r in effects}==set(S)
    for r in effects:
        s=r["source"];assert r["measurement_lane"]==S[s]["measurement_lane"] and r["meta_eligible"]=="no"
        for key in ("gain","same","opposite","abstain","allocation_tv"):
            r[key]=float(r[key])
        assert -1<=r["gain"]<=1 and abs(r["same"]+r["opposite"]+r["abstain"]-1)<1e-9
        r["calibrated_se"]=None
        r["effect_definition"]="one_lag_signed_vs_equal_mass_signblind_TV"
        r["source_receipt"]="validation/ahif_v1/hayashi2026_direction_summary_20261125.json" if s=="hayashi" else "validation/ahif_v1/one_lag_fungal_summary_20261125.json"
    vars=[dict(zip(c["variable_columns"],v)) for v in c["variables"]]
    assert len({r["variable"] for r in vars})==len(vars)
    candidates=[{"candidate_id":r["candidate_id"],"publication_doi":r["publication_doi"],"data_doi":r["data_doi"],"status":r["status"],"note":r["note"]} for r in c["candidate_sources"]]
    db=a.out/"ahif_mechanism_v0.sqlite"
    if db.exists():db.unlink()
    con=sqlite3.connect(db);con.execute("PRAGMA foreign_keys=ON")
    con.executescript("""
    CREATE TABLE studies(source TEXT PRIMARY KEY,study_label TEXT,domain TEXT,ecosystem TEXT,measurement_lane TEXT,taxonomic_resolution TEXT,cluster_unit TEXT,n_clusters INTEGER CHECK(n_clusters>0),n_panels INTEGER CHECK(n_panels>0),exposure TEXT,exposure_class TEXT,exposure_schedule TEXT,immigration_context TEXT,time_unit TEXT,observed_times TEXT,publication_overlap_group TEXT,source_status TEXT,source_evidence TEXT);
    CREATE TABLE panels(source TEXT REFERENCES studies(source),panel TEXT,context TEXT,treated TEXT,reference TEXT,n_observed_times INTEGER,primary_targets INTEGER,evidence TEXT,PRIMARY KEY(source,panel));
    CREATE TABLE study_effects(source TEXT PRIMARY KEY REFERENCES studies(source),study TEXT,measurement_lane TEXT,independent_unit_candidates INTEGER,panels INTEGER,dependent_cells INTEGER,gain REAL,same REAL,opposite REAL,abstain REAL,allocation_tv REAL,meta_eligible TEXT CHECK(meta_eligible='no'),note TEXT,calibrated_se REAL,effect_definition TEXT,source_receipt TEXT);
    CREATE TABLE variable_dictionary(variable TEXT PRIMARY KEY,role TEXT,definition TEXT,status TEXT,caveat TEXT);
    CREATE TABLE candidates(candidate_id TEXT PRIMARY KEY,publication_doi TEXT,data_doi TEXT,status TEXT,note TEXT);
    CREATE VIEW one_lag_mechanism AS SELECT s.source,s.domain,s.ecosystem,s.exposure_class,s.exposure_schedule,s.immigration_context,s.time_unit,s.measurement_lane,s.taxonomic_resolution,s.n_clusters,s.publication_overlap_group,e.gain,e.same,e.allocation_tv,e.calibrated_se,e.meta_eligible FROM studies s JOIN study_effects e USING(source);
    """)
    for name,rows in [("studies",studies),("panels",panels),("study_effects",effects),("variable_dictionary",vars),("candidates",candidates)]:
        keys=list(rows[0]);assert all(list(r)==keys for r in rows)
        con.executemany(f"INSERT INTO {name} ({','.join(keys)}) VALUES ({','.join('?' for k in keys)})",[tuple(r.values()) for r in rows])
        csvwrite(a.out/(name+".csv"),rows)
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0]=="ok"
    assert con.execute("PRAGMA foreign_key_check").fetchall()==[]
    q={"status":"VALIDATED_DESCRIPTIVE_REGISTRY_NOT_META","version":c["version"],"counts":{t:con.execute(f"SELECT COUNT(*) FROM {t}").fetchone()[0] for t in ("studies","panels","study_effects","variable_dictionary","candidates")},"effect_directions":{"positive":con.execute("SELECT COUNT(*) FROM study_effects WHERE gain>0").fetchone()[0],"negative":con.execute("SELECT COUNT(*) FROM study_effects WHERE gain<0").fetchone()[0]},"meta_eligible":con.execute("SELECT COUNT(*) FROM study_effects WHERE meta_eligible='yes'").fetchone()[0],"by_measurement_lane":con.execute("SELECT measurement_lane,COUNT(*) FROM studies GROUP BY measurement_lane").fetchall(),"mechanism_unmeasured":con.execute("SELECT variable FROM variable_dictionary WHERE role='unmeasured' ORDER BY variable").fetchall()}
    assert q["counts"]["studies"]==7 and q["counts"]["panels"]==38 and q["meta_eligible"]==0
    (a.out/"registry_audit.json").write_text(json.dumps(q,indent=2)+"\n")
    print("AHIF_MECHANISM_DB_RESULT",json.dumps(q),flush=True);con.close()
if __name__=="__main__":main()
