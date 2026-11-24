#!/usr/bin/env python3
"""2026 Swedish 20-species randomized field priority experiment: named vs rank.

Original source: van Steijn et al. 2026 Ecology DOI 10.1002/ecy.70291
Data DOI 10.5281/zenodo.15363824, original 7 source MD5 checks pinned.
Pre-score contract committed before this script:
docs/AHIF_VANSTEIJN2026_FIELD_PRE_SCORE_CONTRACT_20261124.md

No source datasets redistributed. Original authors already identified germination-
speed and site-related priority effects; this is an explicitly exploratory
independent reanalysis of *representation-dependent historical retention*.
"""
from __future__ import annotations
import argparse
from collections import Counter,defaultdict
import json
from pathlib import Path
import numpy as np

from benchmarks.audit_ahif_vansteijn2026_grain import audit as source_grain, read as source_csv

SITES=("MAL","UME")
YEARS=("2021","2022","2023","2024")
TREATMENTS=("FtS","SFearly","SFlate","StF")
SEED=20261124
N_PERM=9999
N_BOOT=2000


class SourceStop(ValueError):
    pass


def source_panels(directory:Path):
    a=source_grain(directory)
    if a["filter_rows"]!=2832 or a["unique_plot_years"]!=160 or \
        a["source_trmt_disagree_with_assignment"] or \
        a["site_block_plot_keys_in_survey_not_in_assignment"] or \
        a["assigned_plot_keys_not_in_any_survey"] or \
        a["n_duplicate_same_species_per_plot_year"]!=0:
        raise SourceStop("Original experimental grain audit failed")
    plot=source_csv(directory/"plotdata.csv")
    taxa=source_csv(directory/"speciesinfo.csv")
    raw=source_csv(directory/"Pinpoint_all_years.csv")
    target=sorted({z["Species"] for z in taxa})
    if len(target)!=20:raise SourceStop("Expected 20 originally planted taxa")
    early=sorted(z["Species"] for z in taxa if z["Speed"]=="Fast")
    late=sorted(z["Species"] for z in taxa if z["Speed"]=="Slow")
    if len(early)!=10 or len(late)!=10 or set(early)&set(late):
        raise SourceStop("Original fast/slow germination guilds not balanced")
    assigned={(r["Site"],r["Block"],r["Plot"]):r["Trmt"] for r in plot}
    if len(assigned)!=40:raise SourceStop("40 independently assigned source plots expected")
    observed=[r for r in raw if r["Site"]!="GAL" and r["Type"]=="Transect"]
    names=sorted(set(r["Species"] for r in observed))
    # Two sown taxa have no records in any of the original observed transects.
    recorded_target=set(names)&set(target)
    if len(names)!=76 or len(recorded_target)!=18 or len(set(names)-set(target))!=58:
        raise SourceStop("Source species roster changed or target join failed")
    group=defaultdict(dict)
    for r in observed:
        key=(r["Site"],r["Year"],r["Block"],r["Plot"])
        plotkey=(r["Site"],r["Block"],r["Plot"])
        if plotkey not in assigned or assigned[plotkey]!=r["Trmt"]:
            raise SourceStop("Original assigned history mismatch")
        try:raw_count=float(r["Abundance"])
        except ValueError as e:raise SourceStop("Invalid positive original abundance") from e
        if not np.isfinite(raw_count) or raw_count<=0:
            raise SourceStop("Original positive count invalid; do not impute")
        if r["Species"] in group[key]:
            raise SourceStop("Duplicate taxon in same plot year")
        group[key][r["Species"]]=raw_count
    if len(group)!=160:raise SourceStop("Not 40 plots sampled four years")
    site_plotids={}
    for site in SITES:
        site_plotids[site]=sorted([k for k in assigned if k[0]==site],
                                  key=lambda x:(int(x[1]),x[2]))
        if len(site_plotids[site])!=20:
            raise SourceStop("Not exactly 20 plots per site")
        for block in range(1,6):
            x=[k for k in site_plotids[site] if int(k[1])==block]
            if len(x)!=4 or {assigned[k] for k in x}!=set(TREATMENTS):
                raise SourceStop("Randomized block structure violated")
    panel={}
    for site in SITES:
        plots=site_plotids[site]
        for year in YEARS:
            x=np.zeros((20,len(names)),float)
            for i,p in enumerate(plots):
                d=group.get((site,year,p[1],p[2]))
                if not d:raise SourceStop("Missing observational plot year")
                for j,name in enumerate(names):
                    x[i,j]=d.get(name,0.)
            if np.any(x.sum(axis=1)<=0):
                raise SourceStop("All-zero plot: do not silently drop")
            panel[site,year]={
                "X":x,
                "plots":plots,
                "history":np.asarray([assigned[p] for p in plots]),
                "names":names,
                "sown_target":target,"fast":early,"slow":late,
                "site":site,"year":year
            }
    a["source_species_union_observed"]=len(names)
    a["n_sown_species_with_actual_transect_records"]=18
    a["n_unsown_species_in_recorded_union"]=58
    a["n_sown_species_NEVER_recorded"]=2
    a["source_names_all_zero_sown"]=sorted(set(target)-set(names))
    a["correct_scoring_unit"]="independent whole plot in a randomized block; same plot repeated 4 years"
    return panel,a


def bray(x:np.ndarray)->np.ndarray:
    diff=np.sum(np.abs(x[:,None,:]-x[None,:,:]),axis=2)
    denom=np.sum(x[:,None,:]+x[None,:,:],axis=2)
    return np.divide(diff,denom,out=np.zeros_like(diff),where=denom>0)


def gower(d:np.ndarray)->np.ndarray:
    n=len(d)
    h=np.eye(n)-np.ones((n,n))/n
    return -.5*h@(d*d)@h


def treatment_r2(g:np.ndarray,labels:np.ndarray)->float:
    trace=float(np.trace(g))
    if trace<=1e-12:return 0.
    groups=np.unique(labels)
    if len(groups)!=4:
        raise SourceStop("Expected four history groups")
    explained=sum(float(g[np.ix_(i,i)].sum())/len(i)
                  for v in groups for i in [np.flatnonzero(labels==v)])
    return float(explained/trace)


def blocked_assignment_stream(labels:np.ndarray,n_perm:int=9999,seed:int=SEED)->np.ndarray:
    if len(labels)!=20 or set(labels)!=set(TREATMENTS) or \
       any(set(labels[i:i+4])!=set(TREATMENTS) for i in range(0,20,4)):
        raise SourceStop("Source four-per-block assignment ordering invalid")
    rng=np.random.default_rng(seed)
    per=np.tile(labels,(n_perm,1))
    for trial in per:
        for i in range(0,20,4):
            trial[i:i+4]=rng.permutation(trial[i:i+4])
    return per


def rank_vectors(x):
    return np.sort(x,axis=1)[:,::-1]


def score_matrix(x:np.ndarray,labels:np.ndarray,permutations:np.ndarray,
                 metric:str="bray")->dict:
    if metric=="bray":d=bray(x)
    elif metric=="euclidean":
        d=np.sqrt(np.sum((x[:,None,:]-x[None,:,:])**2,axis=2))
    else:raise ValueError("Wrong metric")
    g=gower(d)
    obs=treatment_r2(g,labels)
    nul=np.asarray([treatment_r2(g,y) for y in permutations])
    if not np.isfinite(nul).all() or not np.isfinite(obs):
        raise SourceStop("Nonfinite treatment R2")
    return {
        "R2":obs,
        "null_median":float(np.median(nul)),
        "null95":np.quantile(nul,[.025,.975]).tolist(),
        "E":obs-float(np.median(nul)),
        "p_upper":float((1+np.count_nonzero(nul>=obs-1e-12))/(len(nul)+1)),
        "n_units":len(labels),
        "n_columns":x.shape[1],
        "n_permutations":len(nul),
        "metric":metric
    }


def within_roster(panel:dict,selected:list[str],permutation:np.ndarray):
    columns={s:i for i,s in enumerate(panel["names"])}
    # The complete 20-species sown target roster contains 2 species with
    # no records anywhere: represented by true zeros, never invented biomass.
    sub=np.asarray([
        [panel["X"][row,columns[n]] if n in columns else 0.
         for n in selected] for row in range(len(panel["X"]))
    ],float)
    zeros=np.flatnonzero(sub.sum(axis=1)<=0)
    if len(zeros):
        return {"status":"STOP_ZERO_SUBCOMMUNITY_TOTAL",
                "n_unscorable_whole_plots":int(len(zeros)),
                "n_species_columns":len(selected)}
    frac=sub/sub.sum(axis=1,keepdims=True)
    return {
        "status":"SCORED_CONDITIONAL_ON_GUILD_MEMBERSHIP",
        "n_species_columns":len(selected),
        "n_whole_plot_units":len(frac),
        "named":score_matrix(frac,panel["history"],permutation),
        "sorted":score_matrix(rank_vectors(frac),panel["history"],permutation)
    }


def holm(primary:list[dict]):
    numbered=[(i,r) for i,z in enumerate(primary) for r in ("named76","rank76")]
    order=sorted(numbered,key=lambda k:primary[k[0]]["representations"][k[1]]["p_upper"])
    m=len(order)
    last=0.
    for k,(i,r) in enumerate(order):
        p=primary[i]["representations"][r]["p_upper"]
        last=max(last,min(1.,p*(m-k)))
        primary[i]["representations"][r]["p_holm_16_site_year_primary"]=float(last)


def site_plot_resampling_delta(named:np.ndarray,rank:np.ndarray,labels:np.ndarray,
                               B:int=N_BOOT,seed:int=SEED):
    # Cluster sampling entire 4-treatment blocks: one physical plot per class
    # per block; keep all copied vectors and whole-block label arrays together.
    rng=np.random.default_rng(seed)
    out=np.empty(B,float)
    for b in range(B):
        chosen=rng.choice(np.arange(5),size=5,replace=True)
        ind=np.concatenate([np.arange(4*i,4*i+4) for i in chosen])
        yn=labels[ind]
        out[b]=(treatment_r2(gower(bray(named[ind])),yn)-
                treatment_r2(gower(bray(rank[ind])),yn))
    return {
        "delta_R2_named_minus_rank":float(
            treatment_r2(gower(bray(named)),labels)-
            treatment_r2(gower(bray(rank)),labels)),
        "block_resample_bootstrap95":np.quantile(out,[.025,.975]).tolist(),
        "n_source_blocks":5,
        "n_cluster_bootstraps":B,
        "note":"Five whole source blocks per site only; percentile uncertainty descriptive, not independent weather year CI"
    }


def analyse(panel:dict,audit:dict,n_perm:int=N_PERM,n_boot:int=N_BOOT)->dict:
    primary=[]
    sens=[]
    longitudinal=[]
    for site in SITES:
        reference=panel[site,YEARS[0]]
        y=reference["history"]
        plots=reference["plots"]
        per=blocked_assignment_stream(y,n_perm=n_perm,seed=SEED+(1 if site=="MAL" else 2))
        series_named=[]
        series_rank=[]
        for year in YEARS:
            state=panel[site,year]
            if state["plots"]!=plots or not np.array_equal(state["history"],y):
                raise SourceStop("Same 20 plots and assignments not held across 4 years")
            x=state["X"];names=state["names"]
            p=x/x.sum(axis=1,keepdims=True)
            ranked=rank_vectors(p)
            named=score_matrix(p,y,per)
            rank=score_matrix(ranked,y,per)
            primary.append({
                "site":site,"year":year,"n_source_plot_units":len(x),
                "n_obs_recorded_species":len(names),
                "representations":{"named76":named,"rank76":rank},
                "delta_E_named_minus_rank":named["E"]-rank["E"],
                "delta_R2_named_minus_rank":named["R2"]-rank["R2"],
                "mean_total_pinpoint_hits":float(x.sum(axis=1).mean()),
                "mean_n_detected_species":float(np.mean((x>0).sum(axis=1))),
                "number_of_implicitly_zero_filled_taxon_by_plot_cells":int(np.sum(x==0))
            })
            present=set(names)
            selected={
                "sown20":state["sown_target"],
                "fast_sown10":state["fast"],
                "slow_sown10":state["slow"],
                "spontaneous_unsown58":sorted(present-set(state["sown_target"]))
            }
            srow={"site":site,"year":year}
            for key,taxa in selected.items():
                srow[key]=within_roster(state,taxa,per)
            srow["top1_share"]=score_matrix(ranked[:,:1],y,per,metric="euclidean")
            srow["total_pinpoint_hits"]=score_matrix(x.sum(axis=1,keepdims=True),y,per,
                                                    metric="euclidean")
            sens.append(srow)
            series_named.append(p)
            series_rank.append(ranked)
        # Exactly 20 physical plots, 4 annual repeated species compositions.
        seq_named=np.concatenate(series_named,axis=1)
        seq_rank=np.concatenate(series_rank,axis=1)
        if seq_named.shape!=(20,4*76) or seq_rank.shape!=(20,4*76):
            raise SourceStop("Year-concatenated source vectors not 304 columns")
        name_score=score_matrix(seq_named,y,per)
        rank_score=score_matrix(seq_rank,y,per)
        longitudinal.append({
            "site":site,"n_independent_whole_plots":20,
            "n_repeated_years":4,"n_features":304,
            "named_4yr":name_score,"sorted_4yr":rank_score,
            "delta_E_named_minus_rank":name_score["E"]-rank_score["E"],
            "delta_R2_named_minus_rank":name_score["R2"]-rank_score["R2"],
            "blocked_source_plot_bootstrap":site_plot_resampling_delta(
                seq_named,seq_rank,y,n_boot,seed=SEED+(51 if site=="MAL" else 52))
        })
    if len(primary)!=8 or len(sens)!=8 or len(longitudinal)!=2:
        raise SourceStop("Expected all 8 site-years and two site-longitudinal states")
    holm(primary)
    return {
        "status":"INDEPENDENT_REAL_DATA_FIELD_HISTORY_RANK_PILOT_NOT_MECHANISM_MEDIATION",
        "article_doi":"10.1002/ecy.70291",
        "data_doi":"10.5281/zenodo.15363824",
        "analysis_contract":"docs/AHIF_VANSTEIJN2026_FIELD_PRE_SCORE_CONTRACT_20261124.md",
        "source_original_method_filter":"Site != GAL and Type == Transect",
        "source_plot_grain_audit":audit,
        "analysis_design":{
            "n_independent_original_field_plots":40,
            "site":list(SITES),"year":list(YEARS),
            "n_plot_year_states":160,"original_4_history_labels":list(TREATMENTS),
            "permutation_unit":"four labels permuted WITHIN each original randomized 4-plot source block",
            "same_assignments_reused_across_years":True,
            "n_permutations_per_site":n_perm,
            "n_cluster_bootstrap":n_boot,
            "n_target_sown_taxa":20,"n_target_taxa_recorded":18,
            "n_source_union_recorded_taxa":76,
            "n_spontaneous_non_target_taxa":58,
            "original_grain":"positive vegetation survey species pin hits, not individuals or biomass"
        },
        "primary_full76_4_history_by_site_year":primary,
        "source_predeclared_secondary_roster_and_scalar_sensitivities":sens,
        "two_site_4yr_whole_plot_trajectories":longitudinal,
        "scientific_boundaries":[
            "The original source paper already published how germination speed, sowing schedule and site affect priority effects.",
            "The survey's Type=Plot and GAL source fields are excluded in source authors' code; 2832 Transect rows only.",
            "Zeros are only the unlisted species in an otherwise observed source plot-year; 2 sown never-recorded taxa are full zero vector only in sown roster.",
            "Eight site-year states reuse the SAME 40 plots; they are not eight independent studies or year-randomized climate experiments.",
            "Within-site block permutations reproduce the original 5 blocks and four labels; do not use arbitrary plot label permutations.",
            "Relative positive pin-hit profiles are not biomass, individual demographic growth rates, or direct competition strength.",
            "Distance-geometry E is NOT entropy/mutual information and effect shape changes do not prove named-vs-rank mechanism.",
            "This original study did NOT quantify separate arrival-conditioned neighbor-only coefficients alpha_ij, resource depletion or soil-conditioned mediation.",
            "The original 2025 van Steijn directed-pair study shares researcher overlap but is not a same-unit causal mediator replicate.",
            "No changes to frozen KNK v28-v38 or submitted-candidate Ecology Letters manuscript; one source cluster remains one study."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-dir",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=N_PERM)
    p.add_argument("--bootstrap",type=int,default=N_BOOT)
    args=p.parse_args()
    P,A=source_panels(args.original_dir)
    result=analyse(P,A,n_perm=args.permutations,n_boot=args.bootstrap)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("AHIF_FIELD_RANK_SCORE_START")
    print("SOURCE_EXPERIMENT",json.dumps({
        "n_independent_plots":40,"n_plot_years":160,
        "n_observed_plant_taxa":76,"n_sown_taxa_recorded":18}))
    for row in result["primary_full76_4_history_by_site_year"]:
        print("SITE_YEAR",json.dumps(row))
    for row in result["two_site_4yr_whole_plot_trajectories"]:
        print("SITE_LONGITUDINAL",json.dumps(row))
    for row in result["source_predeclared_secondary_roster_and_scalar_sensitivities"]:
        print("SITE_SECONDARY",json.dumps({
            "site":row["site"],"year":row["year"],
            **{k:(
                {"status":z["status"],"named_E":z["named"]["E"],
                 "named_p":z["named"]["p_upper"],
                 "rank_E":z["sorted"]["E"],"rank_p":z["sorted"]["p_upper"]}
                 if z["status"].startswith("SCORED") else z)
              for k,z in row.items() if k not in ("site","year") and "status" in z}
        }))
    print("AHIF_FIELD_RANK_SCORE_DONE",args.out)

if __name__=="__main__":main()
