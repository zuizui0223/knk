#!/usr/bin/env python3
"""POEM follow-up: study-level temporal composite and year-2 storage-channel audit.

Exploratory after six primary scores were revealed (not a fresh preregistered
confirmatory endpoint). Repeated observations are concatenated per persistent
whole plot, not treated as independent years. Raw source is never redistributed.
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path
import numpy as np
from benchmarks.audit_poem_design import SOURCE_SOWN
from benchmarks.score_poem_field_storage import (
    source_panels, matrix_for, FOUR, SEED, N_PERM,
    bray_curtis, gower, history_r2, evaluate_one
)


def score_g(g, y, perms):
    obs=history_r2(g,y)
    null=np.array([history_r2(g,p) for p in perms])
    return {"R2":obs,"null_median":float(np.median(null)),
            "E":obs-float(np.median(null)),
            "p_upper":float((1+np.count_nonzero(null>=obs-1e-12))/(len(null)+1)),
            "null_ci95":np.quantile(null,[.025,.975]).tolist(),
            "n_permutations":len(perms)}


def hellinger_score(p,y,perms):
    sq=np.sqrt(p)
    d=np.linalg.norm(sq[:,None,:]-sq[None,:,:],axis=2)
    return score_g(gower(d),y,perms)


def area_source_composites(panels, n_perm, seed):
    cohort_samples=[]
    outcomes=[]
    # Each of 20 matched S/F/G/L field plots contributes one vector
    # concatenating its same, fixed 49-taxon proportions over 3 seasons.
    for cohort in ("2020","2021"):
        series=[]
        for age in (1,2,3):
            x,y,plots,names=matrix_for(panels[(cohort,age)],FOUR)
            if np.any(x.sum(axis=1)<=0):raise ValueError("Zero total mass")
            p=x/x.sum(axis=1,keepdims=True)
            series.append((p,y,plots))
        if not all(np.array_equal(series[0][1],q[1]) and series[0][2]==q[2] for q in series):
            raise ValueError("Repeated plot identities not aligned")
        stacked=np.concatenate([q[0] for q in series],axis=1)
        ranks=np.concatenate([np.sort(q[0],axis=1)[:,::-1] for q in series],axis=1)
        rng=np.random.default_rng(seed+int(cohort))
        perms=np.array([rng.permutation(y) for _ in range(n_perm)])
        named=evaluate_one(stacked,y,perms)
        rank=evaluate_one(ranks,y,perms)
        cohort_samples.append((cohort,stacked,ranks,y))
        outcomes.append({"cohort":cohort,"n_independent_plots":20,
                         "n_repeated_seasons":3,"n_features":stacked.shape[1],
                         "named":named,"rank":rank,
                         "delta_E":named["E"]-rank["E"]})
    return outcomes


def year2_audit(panels, n_perm:int=9999,seed:int=20261126,
                bootstrap:int=3000)->dict:
    cohorts={}
    source_sown=set(SOURCE_SOWN)
    # Compare the same pre-defined 4-group histories at biological age 2.
    for cohort in ("2020","2021"):
        x,y,plots,names=matrix_for(panels[(cohort,2)],FOUR)
        p=x/x.sum(axis=1,keepdims=True)
        rank=np.sort(p,axis=1)[:,::-1]
        rng=np.random.default_rng(seed+int(cohort))
        perm=np.asarray([rng.permutation(y) for _ in range(n_perm)])
        named_h=hellinger_score(p,y,perm)
        rank_h=hellinger_score(rank,y,perm)
        # Explicitly post hoc: separate naturally recruited unsown taxa.
        unsown=np.array([s not in source_sown for s in names])
        unsown_g=x[:,unsown]
        if np.all(unsown_g.sum(axis=1)>0):
            un=unsown_g/unsown_g.sum(axis=1,keepdims=True)
            un_named=evaluate_one(un,y,perm)
            un_rank=evaluate_one(np.sort(un,axis=1)[:,::-1],y,perm)
            unsown_result={"status":"SCORED","n_unsown_taxa":int(unsown.sum()),
                           "named":un_named,"rank":un_rank}
        else:
            unsown_result={"status":"STOP_ZERO_UNSOWN_MASS",
                           "zero_unsown_plots":int(np.count_nonzero(unsown_g.sum(axis=1)==0))}
        sown_mass=x[:,~unsown].sum(axis=1)
        shares=sown_mass/x.sum(axis=1)
        means={t:{"n_plots":int(np.sum(y==t)),
                  "sown_share_mean":float(shares[y==t].mean()),
                  "top1_share_mean":float(rank[y==t,0].mean()),
                  "mean_total_grams_per_quadrat":float(x[y==t].sum(axis=1).mean())}
               for t in FOUR}
        winner=np.argmax(x,axis=1)
        common_top=Counter((t,names[w]) for t,w in zip(y,winner))
        top_treatment={
            t:sorted([{"species":sp,"count":int(ct)}
                      for (tt,sp),ct in common_top.items() if tt==t],
                     key=lambda z:-z["count"])
            for t in FOUR
        }
        cohorts[cohort]={
            "hellinger_named":named_h,"hellinger_sorted":rank_h,
            "unsown_species_only":unsown_result,
            "treatment_descriptives":means,
            "rank1_winner_counts_by_treatment":top_treatment,
            "raw_data_scoring_n_units":len(y),
        }
    # Bootstrap each cohort by resampling the same five assigned plot IDs within
    # each treatment. Bootstrap is only within-study descriptive uncertainty.
    rng=np.random.default_rng(seed)
    cases={}
    for cohort in ("2020","2021"):
        x,y,_,_=matrix_for(panels[(cohort,2)],FOUR)
        frac=x/x.sum(axis=1,keepdims=True)
        idx=[np.flatnonzero(y==t) for t in FOUR]
        draws=[]
        for _ in range(bootstrap):
            bidx=np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in idx])
            named=gower(bray_curtis(frac[bidx]))
            ranked=gower(bray_curtis(np.sort(frac[bidx],axis=1)[:,::-1]))
            draws.append(history_r2(named,y[bidx])-history_r2(ranked,y[bidx]))
        cases[cohort]=np.asarray(draws)
        cohorts[cohort]["delta_R2_named_minus_sorted_bootstrap"]={
            "point_estimate":float(
                history_r2(gower(bray_curtis(frac)),y)-
                history_r2(gower(bray_curtis(np.sort(frac,axis=1)[:,::-1])),y)
            ),
            "percentile95":np.quantile(cases[cohort],[.025,.975]).tolist(),
            "n_plot_bootstrap":bootstrap
        }
    difference=cases["2021"]-cases["2020"]
    return {"status":"POST_HOC_YEAR2_CHANNEL_MECHANISM_AUDIT",
            "two_cohorts_one_published_study":True,
            "seed":seed,"n_permutations":n_perm,
            "year2_cohorts":cohorts,
            "cross_cohort_deltaR2_difference_2021_minus_2020":{
                "point_difference":float(
                    cohorts["2021"]["delta_R2_named_minus_sorted_bootstrap"]["point_estimate"]-
                    cohorts["2020"]["delta_R2_named_minus_sorted_bootstrap"]["point_estimate"]),
                "percentile95":np.quantile(difference,[.025,.975]).tolist(),
                "bootstrap_reps":bootstrap,
                "interpretation":"Within-cohort plot resampling cannot generalize to new weather years (n=1 per establishment year)"
            },
            "limitations":["Unsown vs sown selection was motivated by the revealed main result",
                           "Cohort is not randomized; difference cannot identify weather as cause",
                           "Contrast in p values is not itself a formal interaction test",
                           "Repeated plots not independent across years",
                           "No claim of niche pre-emption or source-level meta pooled effect"]}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--permutations",type=int,default=9999)
    ap.add_argument("--bootstrap",type=int,default=3000)
    args=ap.parse_args()
    panels,audit=source_panels(args.zip)
    composites=area_source_composites(panels,args.permutations,SEED)
    y2=year2_audit(panels,n_perm=args.permutations,bootstrap=args.bootstrap)
    output={
        "status":"EXPLORATORY_SOURCE_LEVEL_REPEATED_PLOT_AND_YEAR2_SENSITIVITIES",
        "source_sha256":audit["source_zip_sha256"],
        "source_contract":"docs/POEM_FIELD_PRE_SCORE_CONTRACT_20261123.md",
        "source_level_three_age_composites":composites,
        "year2_posthoc":y2,
        "claim_boundaries":["Three years of repeated plots form one trajectory per plot",
                           "Year 2 selected for diagnostics after results: not confirmatory",
                           "Four source clusters including this study; no random-effects meta-analysis"]
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2)+"\n")
    print("POEM_SECONDARY_RESULT_START")
    for v in composites:
        print("COHORT_COMPOSITE",json.dumps(v))
    for k,v in y2["year2_cohorts"].items():
        print("YEAR2_AUDIT",json.dumps({"cohort":k,**v}))
    print("COHORT_SWITCH_DIFF",json.dumps(y2["cross_cohort_deltaR2_difference_2021_minus_2020"]))
    print("POEM_SECONDARY_RESULT_END",args.output)

if __name__=="__main__":
    main()
