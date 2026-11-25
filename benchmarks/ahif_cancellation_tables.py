"""Source-validated absolute response cancellation ledger; separate from effects."""
import csv,json,math,statistics
from pathlib import Path
POEM="validation/ahif_v1/poem_field_absolute_compensation_v1.json"
LEAR="validation/ahif_v1/lear2024_resident_cancellation_v1.json"
PRUN="https://github.com/zuizui0223/knk/actions/runs/45536266086"
LRUN="https://github.com/zuizui0223/knk/actions/runs/45527948842"
C=("study_id","comparison_id","scope","unit","gross","net","cancellation","fraction",
   "signed_total","crossfit_gap","interpretation_scope","source_json","source_run")
I=("study_id","comparison_id","endpoint","unit","estimate","ci95_low","ci95_high",
   "interval_family_size","interpretation_scope","source_json","source_run")
def close(a,b):
    return math.isclose(float(a),float(b),rel_tol=1e-9,abs_tol=1e-10)
def metrics(d):
    g=sum(abs(x) for x in d);n=abs(sum(d));c=max(0.,g-n)
    return g,n,c,c/g if g>1e-12 else None
def check(m,d=None):
    g,n,c,f=(m[k] for k in ("gross","net","cancellation","fraction"))
    if not all(math.isfinite(v) and v>=0 for v in (g,n,c)) or not close(g-n,c):
        raise ValueError("cancellation arithmetic")
    if g>1e-12:
        if f is None or not close(f,c/g):raise ValueError("cancellation fraction")
    elif f is not None:raise ValueError("zero denominator")
    if d is not None:
        vals=metrics(d)
        if not all(close(a,b) for a,b in zip((g,n,c),vals[:3])):raise ValueError("vector mismatch")
    return (g,n,c,f)
def export(path,headers,rows):
    with path.open("w",encoding="utf-8",newline="") as f:
        w=csv.writer(f);w.writerow(headers);w.writerows(rows)
def insert_cancellation(conn,source,out):
    root=source.resolve().parents[1];p=json.loads((root/POEM).read_text());l=json.loads((root/LEAR).read_text())
    if (p["source_run"],p["source_commit"],p["source_doi"],p["status"])!=(PRUN,"551b2034dbed049e8120764c672654f1bf7a91e5","10.5281/zenodo.10119982","SOURCE_LOCKED_POEM_ABSOLUTE_BIOMASS_COMPENSATION_NOT_MEDIATION"):
        raise ValueError("POEM lineage")
    if (l["source_run"],l["source_commit"],l["source_doi"],l["status"])!=(LRUN,"2267688c5d6fd7468ad501bdac969278ba1c4b4b","10.5281/zenodo.12528023","SOURCE_LOCKED_RESIDENT_COUNT_CANCELLATION_NOT_INTERACTION_MEDIATION"):
        raise ValueError("Lear lineage")
    if p["n_primary_plots"]!=40 or p["audit"]["total_source_rows"]!=5210 or len(p["taxon_names"])!=49 or p["n_split_draws"]!=256 or l["audit"]["source_rows"]!=197 or len(l["audit"]["ambiguous_rep_cells"])!=2:
        raise ValueError("source grain")
    rows=[];intervals=[]
    def add(st,id,scope,unit,m,sign,cf,interpretation,path,run,d=None):
        vals=check(m,d)
        if sign is not None and not close(abs(sign),m["net"]):raise ValueError("signed total mismatch")
        if cf is not None and not math.isfinite(cf):raise ValueError("crossfit nonfinite")
        rows.append((st,id,scope,unit,*vals,sign,cf,interpretation,path,run))
    def iv(st,id,ep,unit,est,ci,n,scope,path,run):
        if not math.isfinite(est):raise ValueError("nonfinite endpoint")
        if ci is not None and (not all(math.isfinite(v) for v in ci) or not ci[0]<=est<=ci[1]):raise ValueError("interval mismatch")
        intervals.append((st,id,ep,unit,est,*(ci or (None,None)),n,scope,path,run))
    ps="ST_POEM2025";ls="ST_LEAR2024";pu="g_identified_shoot_per_quadrat";lu="resident_colony_count"
    keys={(x["cohort"],x["age"],x["early"]) for x in p["contrasts"]}
    if len(p["contrasts"])!=18 or keys!={(c,a,g) for c in ("2020","2021") for a in (1,2,3) for g in ("F","G","L")}:
        raise ValueError("POEM contrasts")
    for x in p["contrasts"]:
        d=x["species_delta"];m=x["taxon_metrics"]
        if len(d)!=49 or not close(sum(d),m["signed_total"]):raise ValueError("taxon grain")
        add(ps,f'{x["cohort"]}_age{x["age"]}_{x["early"]}_vs_S',"within_cohort_taxon",pu,m,m["signed_total"],x["crossfit_taxon_signed_gap"],"descriptive_reused_control",POEM,PRUN,d)
    if len(p["ages"])!=3:raise ValueError("POEM ages")
    for a in p["ages"]:
        age=a["age"];id=f"pooled_age{age}"
        roles=a["role_estimates"]
        if [x["role"] for x in roles]!=p["role_names"] or not close(sum(x["estimate"] for x in roles[:3]),roles[3]["estimate"]):raise ValueError("role sum")
        role=a["role_metrics"];d=[x["estimate"] for x in roles[:3]]
        add(ps,id,"pooled_roles",pu,role,role["signed_total"],a["pooled_role_crossfit_gap"],"descriptive_role_aligned",POEM,PRUN,d)
        mm=a["mean_within_contrast_taxon_metrics"];sel=[x for x in p["contrasts"] if x["age"]==age]
        if len(sel)!=6 or any(not close(mm[k],statistics.mean(x["taxon_metrics"][k] for x in sel)) for k in ("gross","net","cancellation")) or not close(mm["mean_crossfit_taxon_signed_gap"],statistics.mean(x["crossfit_taxon_signed_gap"] for x in sel)):
            raise ValueError("age aggregate")
        add(ps,id,"mean_within_contrast_taxon",pu,mm,None,mm["mean_crossfit_taxon_signed_gap"],"descriptive_ratio_of_means",POEM,PRUN)
        n=a["interval_family_size"]
        if n!=(4 if age==2 else 1):raise ValueError("interval family")
        for x in roles:iv(ps,id,x["role"],pu,x["estimate"],x["ci95"],n,"approx_simultaneous_exploratory" if age==2 else "unadjusted_descriptive",POEM,PRUN)
    if len(l["contrasts"])!=4:raise ValueError("Lear contrasts")
    for x in l["contrasts"]:
        id=f'{x["after"]}_minus_{x["before"]}';m=x["pooled"]
        scope="descriptive_ambiguous_rep" if x["before"]=="6" else "descriptive_pooled"
        cf=x["crossfit"];ints=x["intervals"]
        if (cf is None)!=(x["before"]=="6") or (ints is None)!=(x["before"]=="6"):
            raise ValueError("ambiguous label inference")
        add(ls,id,"pooled_resident",lu,m,m["signed_total_delta"],None if cf is None else cf["pooled_validated_gap"],scope,LEAR,LRUN,m["delta"])
        if len(x["by_level"])!=4:raise ValueError("source levels")
        for v in x["by_level"]:
            add(ls,id+"_"+v["level"],"within_source_level",lu,v,v["signed_total_delta"],None,"descriptive_only",LEAR,LRUN,v["delta"])
        if ints:
            if len(ints["rows"])!=5 or ints["n_family"]!=(5 if x["before"]=="7.9" else 10):raise ValueError("resident interval family")
            for v,d in zip(ints["rows"],m["delta"]):
                if not close(v["delta"],d):raise ValueError("resident delta")
                iv(ls,id,v["resident"],lu,v["delta"],v["simultaneous_ci95"],ints["n_family"],"approx_simultaneous_exploratory",LEAR,LRUN)
    # Distinct already-registered field source; never pool hits with biomass/counts.
    fp="validation/ahif_v1/vansteijn2026_field_opposition_v1.json"
    fr="https://github.com/zuizui0223/knk/actions/runs/45537486831"
    d=json.loads((root/fp).read_text())
    if (d["source_run"],d["source_commit"],d["source_doi"],d["status"])!=(fr,"5143c3e29c60caa2fe738acfcaf8da670246a23e","10.5281/zenodo.15363824","SOURCE_LOCKED_FIELD_HIT_OPPOSITION_NOT_BIOMASS_OR_MEDIATION"):
        raise ValueError("field opposition lineage")
    if d["n_physical_plots"]!=40 or d["n_blocks"]!=10 or d["n_sown"]!=20 or d["n_nonsown"]!=58 or len(d["taxon_names"])!=78 or d["n_splits"]!=256 or d["audit"]["filter_rows"]!=2832:
        raise ValueError("field opposition grain")
    if len(d["taxon_contrasts"])!=32 or {(x["site"],x["year"],x["reference"],x["history"]) for x in d["taxon_contrasts"]}!={(s,y,r,h) for s in ("MAL","UME") for y in ("2021","2022","2023","2024") for r in ("SFearly","SFlate") for h in ("FtS","StF")}:
        raise ValueError("field opposition contrast coverage")
    if len(d["comparisons"])!=8 or {(x["year"],x["reference"]) for x in d["comparisons"]}!={(y,r) for y in ("2021","2022","2023","2024") for r in ("SFearly","SFlate")}:
        raise ValueError("field opposition comparison coverage")
    st="ST_VANSTEIJNFIELD2026";unit="vegetation_pin_hits"
    for x in d["taxon_contrasts"]:
        v=x["species_delta"];m=x["taxon_metrics"]
        if len(v)!=78 or not close(sum(v),m["signed_total"]) or len(x["role_delta"])!=4 or not close(sum(x["role_delta"][:3]),x["role_delta"][3]) or not close(x["role_delta"][3],sum(v)):
            raise ValueError("field taxon/role grain")
        add(st,f'{x["site"]}_{x["year"]}_{x["history"]}_vs_{x["reference"]}',"within_site_taxon",unit,m,m["signed_total"],x["crossfit_taxon_signed_gap"],"descriptive_reused_block_control",fp,fr,v)
    for x in d["comparisons"]:
        y=x["year"];ref=x["reference"];key=f"pooled_{y}_vs_{ref}"
        primary=y=="2022" and ref=="SFlate"
        if x["primary"]!=primary or x["interval_family_size"]!=(4 if primary else 1):
            raise ValueError("field primary gate")
        rr=x["role_estimates"]
        if [a["role"] for a in rr]!=d["role_names"] or len(rr)!=4:
            raise ValueError("field role roster")
        v=[a["estimate"] for a in rr]
        sel=[a for a in d["taxon_contrasts"] if a["year"]==y and a["reference"]==ref]
        if len(sel)!=4 or not close(sum(v[:3]),v[3]) or any(not close(v[j],statistics.mean(a["role_delta"][j] for a in sel)) for j in range(4)):
            raise ValueError("field role aggregation")
        m=x["role_metrics"]
        add(st,key,"pooled_roles",unit,m,m["signed_total"],x["pooled_role_crossfit_gap"],"descriptive_role_aligned",fp,fr,v[:3])
        mm=x["mean_within_contrast_taxon_metrics"]
        if any(not close(mm[k],statistics.mean(a["taxon_metrics"][k] for a in sel)) for k in ("gross","net","cancellation")) or not close(mm["mean_crossfit_taxon_signed_gap"],statistics.mean(a["crossfit_taxon_signed_gap"] for a in sel)):
            raise ValueError("field aggregate")
        add(st,key,"mean_within_contrast_taxon",unit,mm,None,mm["mean_crossfit_taxon_signed_gap"],"descriptive_ratio_of_means",fp,fr)
        for a in rr:
            iv(st,key,a["role"],unit,a["estimate"],a["ci95"],x["interval_family_size"],"approx_simultaneous_exploratory" if primary else "unadjusted_descriptive",fp,fr)
        if primary:
            if x["primary_directional_gate"]!={"other_sown_negative_nonsown_positive":False,"strict_early_positive_other_or_nonsown_negative":True}:
                raise ValueError("field directional gate changed")
    if len(rows)!=92 or len(intervals)!=59:raise ValueError("ledger coverage")
    if conn.execute("SELECT COUNT(*) FROM effects WHERE panel_id='P_LEAR2024'").fetchone()[0]!=0:
        raise ValueError("cancellation-only rows promoted to effects")
    conn.execute("""CREATE TABLE cancellation_contrasts(
      study_id TEXT NOT NULL REFERENCES studies(study_id),comparison_id TEXT NOT NULL,scope TEXT NOT NULL,
      unit TEXT NOT NULL,gross REAL NOT NULL,net REAL NOT NULL,cancellation REAL NOT NULL,fraction REAL,
      signed_total REAL,crossfit_gap REAL,interpretation_scope TEXT NOT NULL,source_json TEXT NOT NULL,
      source_run TEXT NOT NULL,PRIMARY KEY(study_id,comparison_id,scope))""")
    conn.executemany("INSERT INTO cancellation_contrasts VALUES("+",".join("?" for _ in C)+")",rows)
    conn.execute("""CREATE TABLE cancellation_intervals(
      study_id TEXT NOT NULL REFERENCES studies(study_id),comparison_id TEXT NOT NULL,endpoint TEXT NOT NULL,
      unit TEXT NOT NULL,estimate REAL NOT NULL,ci95_low REAL,ci95_high REAL,interval_family_size INTEGER NOT NULL,
      interpretation_scope TEXT NOT NULL,source_json TEXT NOT NULL,source_run TEXT NOT NULL,
      PRIMARY KEY(study_id,comparison_id,endpoint))""")
    conn.executemany("INSERT INTO cancellation_intervals VALUES("+",".join("?" for _ in I)+")",intervals)
    export(out/"cancellation_contrasts.csv",C,rows);export(out/"cancellation_intervals.csv",I,intervals)
    return {"cancellation_contrast_rows":92,"cancellation_interval_rows":59,
            "cancellation_source_families":3,"cancellation_meta_eligible":False,
            "poem_strict_early_role_primary_supported":False}
