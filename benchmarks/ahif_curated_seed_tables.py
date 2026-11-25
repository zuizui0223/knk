"""Curated seed-addition synthesis diagnostics, isolated from original-source effect ledger."""
import csv,json,math,statistics
from collections import Counter
from pathlib import Path

SOURCE="CURATED_LADOUCEUR2020"
AUDIT="validation/ahif_v1/ladouceur2020_synthesis_source_audit_v1.json"
DIRECTION="validation/ahif_v1/ladouceur2020_two_direction_validation_v1.json"
KANSAS="validation/ahif_v1/kansas2020_time_block_forecast_v1.json"
KRUN="https://github.com/zuizui0223/knk/actions/runs/45541282653"
DRUN="https://github.com/zuizui0223/knk/actions/runs/45539576248"
MODES=("direct","indirect")
MODELS=("zero","previous","cumulative","history_ridge","state_ridge","concurrent_sign_reference")
E=("source_id","experiment","n_plot_states","n_taxa","n_direct_states","n_indirect_states","n_mixed_states",
   "duplicate_physical_key_age_rows","drifting_keys","status","source_json")
C=("source_id","experiment","source_age_code","mode","status","n_matched_strata","n_control_states","n_seeded_states",
   "excluded_unmatched_states","final_matched_age","gross","net","cancellation","fraction","signed_total",
   "reported_total_delta","reported_richness_delta","direction_status","held_increase","held_decrease","held_min",
   "fraction_both_positive","strong_descriptive_gate","source_json","source_run")
F=("source_id","experiment","time_scope","source_age_code","model","n_orientations","sse","held_increase","held_decrease",
   "held_min","fraction_both_positive","strong_descriptive_gate","source_json","source_run")

def close(a,b):
    return math.isfinite(float(a)) and math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-10)
def finite(v):
    if v is not None and not math.isfinite(float(v)):raise ValueError("nonfinite curated score")
def export(out,name,headers,rows):
    with (out/name).open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(rows)
def insert_curated(conn,source,out):
    root=source.resolve().parents[1]
    a=json.loads((root/AUDIT).read_text());d=json.loads((root/DIRECTION).read_text());k=json.loads((root/KANSAS).read_text())
    if (a["status"],a["plot_rows"],a["species_rows"],a["n_experiments"])!=("SOURCE_BYTES_VERIFIED_SCHEMA_AUDITED_NOT_ECOLOGICAL_SCORING",3393,58426,12):
        raise ValueError("curated source audit")
    if d["status"]!="CURATED_SYNTHESIS_TWO_DIRECTION_DIAGNOSTIC_NOT_CAUSAL" or d["provenance"]["run_id"]!=45539576248 or d["n_experiment_labels"]!=12 or d["n_plot_states"]!=3393:
        raise ValueError("curated directional lineage")
    if k["status"]!="CURATED_DIRECT_KANSAS_TIME_BLOCK_FORECAST_DESCRIPTIVE" or k["source_run"]!=KRUN or k["source_commit"]!="a7cdd253cf0441d7d1e2bbd3eff4f2e19108c97d":
        raise ValueError("Kansas lineage")
    for key,sha in (("plot","0455c64afd914af5c7f85273ae950104ece31bf94fe9737f1a172aa6d9f437a3"),("species","d9f09ea90b3bac6e9913000d93acf24e1c8b0293bebac563477c9e90a9fd590e")):
        if any(z[field]["files"][key]["sha256"]!=sha for z,field in ((a,"manifest"),(d,"source_manifest"),(k,"source_manifest"))):
            raise ValueError("curated source byte drift")
    ar={v["experiment"]:v for v in a["experiments"]}
    if len(ar)!=12 or {e["experiment"] for e in d["experiments"]}!=set(ar):raise ValueError("curated roster")
    er=[];cr=[];summary=Counter();stops=set()
    for e in d["experiments"]:
        name=e["experiment"];n=e["n_states"];m=e["measurement_mode_counts"]
        if n!=ar[name]["plot_rows"] or sum(m.values())!=n:raise ValueError("state grain")
        er.append((SOURCE,name,n,e["n_taxa"],m.get("direct",0),m.get("indirect",0),m.get("mixed",0),
                   e["duplicate_physical_key_age_rows"],e["physical_keys_with_treatment_or_dose_drift"],e["status"],DIRECTION))
        if e["status"].startswith("STOP"):
            if e["comparisons"]:raise ValueError("ambiguous experiment scored")
            stops.add(name)
        for c in e["comparisons"]:
            mode=c["mode"];status=c["status"];age=c["age"]
            if mode not in MODES:raise ValueError("curated mode")
            if status!="SCORED_DESCRIPTIVE":
                cr.append((SOURCE,name,age,mode,status,None,None,None,c["excluded_unmatched_states"],None,
                           None,None,None,None,None,None,None,status,None,None,None,None,None,DIRECTION,DRUN))
                continue
            summary[mode,"scored"]+=1
            q=c["pooled_species_metrics"];g,n1,v,f,sg=(q[x] for x in ("gross","net","cancellation","fraction","signed_total"))
            if g<0 or n1<0 or v<0 or not close(g-n1,v) or not close(abs(sg),n1) or (g>1e-12 and not close(f,v/g)):
                raise ValueError("curated opposition arithmetic")
            j=c["directional_validation"];eligible=j["status"]=="DESCRIPTIVE_TWO_DIRECTION_HELD_BLOCKS"
            u=down=mn=fr=gate=None
            if eligible:
                summary[mode,"eligible"]+=1
                if c["final_matched_age"]:summary[mode,"final_eligible"]+=1
                u,down,mn,fr=(j[x] for x in ("mean_held_increase","mean_held_decrease","mean_held_min_signed_component","fraction_orientations_both_positive"))
                gate=int(u>0 and down>0 and mn>0)
                if gate!=int(j["strong_direction_descriptive_gate"]) or not close(u+down,c["crossfit"]["signed_gross"]) or not 0<=fr<=1 or j["n_blocks"]!=c["n_matched_strata"]:
                    raise ValueError("directional score mismatch")
                summary[mode,"strong"]+=gate
                if c["final_matched_age"]:summary[mode,"final_strong"]+=gate
            for x in (g,n1,v,f,sg,u,down,mn,fr):finite(x)
            cr.append((SOURCE,name,age,mode,status,c["n_matched_strata"],c["n_control_states"],c["n_seeded_states"],
                       c["excluded_unmatched_states"],int(c["final_matched_age"]),g,n1,v,f,sg,c["reported_total_delta"],c["reported_richness_delta"],j["status"],u,down,mn,fr,gate,DIRECTION,DRUN))
    if len(er)!=12 or len(cr)!=59 or len({x[:4] for x in cr})!=59 or stops!={"Germany_Montane","Texas_Temple_Prarie"}:
        raise ValueError("curated cell coverage")
    for mode,expected in (("direct",(24,23,5,3,0)),("indirect",(34,30,20,5,4))):
        if tuple(summary[mode,x] for x in ("scored","eligible","strong","final_eligible","final_strong"))!=expected:
            raise ValueError("curated gate count")
    if (k["n_states"],k["n_taxa"],k["n_blocks"],k["n_oriented_partitions"],k["n_forecast_times"])!=(120,100,4,6,10) or len(k["all_forecast_cells"])!=60:
        raise ValueError("Kansas source grid")
    cells=k["all_forecast_cells"];rows=[];exp=k["experiment"]
    for c in cells:
        if len(c["train_blocks"])!=2 or len(c["test_blocks"])!=2 or any(x in c["test_blocks"] for x in c["train_blocks"]) or set(c["scores"])!=set(MODELS):
            raise ValueError("Kansas partition or model mismatch")
    ages=k["by_age"]
    if [q["age"] for q in ages]!=list(range(6,16)) or [v["age"] for v in k["exact_concurrent"]]!=list(range(1,16)):
        raise ValueError("Kansas time coverage")
    def add(m,scope,age,n,v):
        u,dd,mn,fr=(v[x] for x in ("mean_held_increase","mean_held_decrease","mean_held_min","fraction_both_positive"))
        se=v.get("SSE")
        gate=int(u>0 and dd>0 and mn>0)
        if gate!=int(v["strong_descriptive_gate"]) or not 0<=fr<=1 or (se is not None and se<0):
            raise ValueError("Kansas score")
        for x in (u,dd,mn,fr,se):finite(x)
        rows.append((SOURCE,exp,scope,age,m,n,se,u,dd,mn,fr,gate,KANSAS,KRUN))
    for q,sel,scope,age in [(k["aggregate_forecast"],cells,"all_forecast_times",0)]+[(q,[v for v in cells if v["age"]==q["age"]],"forecast_time",q["age"]) for q in ages]:
        if len(sel)!=q["n_orientations"] or set(q["models"])!=set(MODELS):raise ValueError("Kansas aggregate grain")
        for m,v in q["models"].items():
            ds=[c["scores"][m] for c in sel]
            for f,g in (("mean_held_increase","held_increase"),("mean_held_decrease","held_decrease"),("mean_held_min","held_min"),("fraction_both_positive","both_positive")):
                if not close(v[f],statistics.mean(x[g] for x in ds)):raise ValueError("Kansas direction aggregation")
            if m!="concurrent_sign_reference" and not close(v["SSE"],sum(x["SSE"] for x in ds)):raise ValueError("Kansas loss aggregation")
            add(m,scope,age,len(sel),v)
    for v in k["exact_concurrent"]:
        add("concurrent_sign_reference","same_time_nonforecast",v["age"],6,v)
    if len(rows)!=81:raise ValueError("Kansas summary coverage")
    conn.executescript("""CREATE TABLE curated_synthesis_sources(source_id TEXT PRIMARY KEY,article_doi TEXT NOT NULL,data_doi TEXT NOT NULL,n_experiment_labels INTEGER,n_plot_states INTEGER,n_species_rows INTEGER,scope TEXT NOT NULL);
    CREATE TABLE curated_seed_experiments(source_id TEXT REFERENCES curated_synthesis_sources(source_id),experiment TEXT,n_plot_states INTEGER,n_taxa INTEGER,n_direct_states INTEGER,n_indirect_states INTEGER,n_mixed_states INTEGER,duplicate_physical_key_age_rows INTEGER,drifting_keys INTEGER,status TEXT,source_json TEXT,PRIMARY KEY(source_id,experiment));
    CREATE TABLE curated_seed_comparisons(source_id TEXT,experiment TEXT,source_age_code INTEGER,mode TEXT,status TEXT,n_matched_strata INTEGER,n_control_states INTEGER,n_seeded_states INTEGER,excluded_unmatched_states INTEGER,final_matched_age INTEGER,gross REAL,net REAL,cancellation REAL,fraction REAL,signed_total REAL,reported_total_delta REAL,reported_richness_delta REAL,direction_status TEXT,held_increase REAL,held_decrease REAL,held_min REAL,fraction_both_positive REAL,strong_descriptive_gate INTEGER,source_json TEXT,source_run TEXT,PRIMARY KEY(source_id,experiment,source_age_code,mode),FOREIGN KEY(source_id,experiment) REFERENCES curated_seed_experiments(source_id,experiment));
    CREATE TABLE curated_kansas_forecast(source_id TEXT,experiment TEXT,time_scope TEXT,source_age_code INTEGER,model TEXT,n_orientations INTEGER,sse REAL,held_increase REAL,held_decrease REAL,held_min REAL,fraction_both_positive REAL,strong_descriptive_gate INTEGER,source_json TEXT,source_run TEXT,PRIMARY KEY(source_id,experiment,time_scope,source_age_code,model),FOREIGN KEY(source_id,experiment) REFERENCES curated_seed_experiments(source_id,experiment));""")
    conn.execute("INSERT INTO curated_synthesis_sources VALUES(?,?,?,?,?,?,?)",(SOURCE,"10.1111/ele.13566","10.6084/m9.figshare.12319682.v1",12,3393,58426,"secondary_synthesis_not_raw_independent_effects"))
    for table,fields,values in (("curated_seed_experiments",E,er),("curated_seed_comparisons",C,cr),("curated_kansas_forecast",F,rows)):
        conn.executemany("INSERT INTO "+table+" VALUES("+",".join("?" for _ in fields)+")",values)
        export(out,table+".csv",fields,values)
    return {"curated_seed_synthesis_sources":1,"curated_seed_experiment_labels":12,
       "curated_seed_comparison_rows":59,"curated_seed_scored_comparisons":58,
       "curated_kansas_forecast_summary_rows":81,"curated_seed_direct_strong_cells":5,
       "curated_seed_indirect_strong_cells":20,"curated_synthesis_meta_eligible":False}
