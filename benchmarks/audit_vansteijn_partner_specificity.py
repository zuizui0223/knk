#!/usr/bin/env python3
"""Post hoc specificity audit for early-neighbor identity in timing response matrix.

Source-locked to Zenodo v2, uses prep() and output of original documented
van-Steijn study data. Does NOT identify alpha_ij or prove direct antagonism.
"""
from __future__ import annotations

import argparse
from collections import defaultdict,Counter
import json
from pathlib import Path

import numpy as np
from benchmarks.score_vansteijn_partner_effects import prep, neighbor_consistency


def centered_by_focal(rows, key):
    by=defaultdict(list)
    for i,r in enumerate(rows):
        if r[key] is not None:by[r["focal"]].append(i)
    obs=np.asarray([r[key] for r in rows],dtype=float)
    means={f:np.mean(obs[indices]) for f,indices in by.items()}
    return np.asarray([obs[k]-means[r["focal"]] for k,r in enumerate(rows)])


def group_score(y, labels):
    d=float(np.dot(y,y))
    if d<1e-12:raise ValueError("No residual variation")
    keys=sorted(set(labels))
    score=sum(
        len([q for q in labels if q==s]) *
        (np.mean([yi for yi,lab in zip(y,labels) if lab==s])**2)
        for s in keys
    )
    return float(score/d)


def within_focal_neighbor_guild_residual(rows,key):
    # Regression against 15 focal names and three guild labels.
    # Residualization is once per analysis: permuting labels inside these cells
    # leaves this nuisance model invariant under conditional exchangeability.
    y=np.asarray([r[key] for r in rows],float)
    focal=sorted({r["focal"] for r in rows})
    guild=sorted({r["neighbor_guild"] for r in rows})
    F=np.zeros((len(rows),len(focal)+len(guild)-1),dtype=float)
    for i,r in enumerate(rows):
        F[i,focal.index(r["focal"])]=1
        if guild.index(r["neighbor_guild"])>0:
            F[i,len(focal)+guild.index(r["neighbor_guild"])-1]=1
    coef=np.linalg.lstsq(F,y,rcond=None)[0]
    return y-F@coef


def conditional_species_identity(rows:list[dict],key:str,n_perm:int=9999,seed:int=20261124):
    dat=[r for r in rows if r.get(key) is not None]
    response=within_focal_neighbor_guild_residual(dat,key)
    labels=[r["neighbor"] for r in dat]
    obs=group_score(response,labels)
    cells=defaultdict(list)
    for idx,r in enumerate(dat):
        cells[(r["focal"],r["neighbor_guild"])].append(idx)
    rng=np.random.default_rng(seed)
    null=np.empty(n_perm)
    for b in range(n_perm):
        perm=list(labels)
        for pos in cells.values():
            shuf=rng.permutation([labels[k] for k in pos])
            for k,j in zip(pos,shuf):perm[k]=j
        null[b]=group_score(response,perm)
    return {
        "focal_and_neighbor_guild_residual_identity_R2":obs,
        "null_median":float(np.median(null)),
        "null_95pct":np.quantile(null,[.025,.975]).tolist(),
        "p_exchangeability_upper":float((1+np.count_nonzero(null>=obs-1e-12))/(n_perm+1)),
        "unit":"directed pair, permute early-neighbor names within focal × guild",
        "n_pairs":len(dat),"n_permutations":n_perm,
        "caveat":"This is an identity-label consistency diagnostic, not randomized causal treatment p.",
    }


def guild_mapping_placebo(rows:list[dict],key:str,n_perm:int=9999,seed:int=20261125):
    dat=[r for r in rows if r.get(key) is not None]
    y=centered_by_focal(dat,key)
    actual=[r["neighbor_guild"] for r in dat]
    score=group_score(y,actual)
    neighbor_names=sorted({r["neighbor"] for r in dat})
    guild_names={r["neighbor"]:r["neighbor_guild"] for r in dat}
    levels=[guild_names[x] for x in neighbor_names]
    assert sorted(Counter(levels).values())==[5,5,5]
    rng=np.random.default_rng(seed)
    vals=np.empty(n_perm)
    for i in range(n_perm):
        p=rng.permutation(levels)
        by=dict(zip(neighbor_names,p))
        vals[i]=group_score(y,[by[r["neighbor"]] for r in dat])
    return {
        "guild_R2_after_focal_centering":score,
        "random_5_5_5_mapping_null_median":float(np.median(vals)),
        "random_5_5_5_mapping_null_95pct":np.quantile(vals,[.025,.975]).tolist(),
        "all_mapping_upper_fraction_MC":float((1+np.count_nonzero(vals>=score-1e-12))/(n_perm+1)),
        "n_random_5_5_5_species_partitions":n_perm,
        "caveat":"Random guild mapping ranks are placebo specificity diagnostics, not experiment-randomized p-values.",
    }


def extreme_neighbor_report(rows, key):
    data=[r for r in rows if r.get(key) is not None]
    resid=centered_by_focal(data,key)
    by=defaultdict(list)
    for i,r in enumerate(data):
        by[r["neighbor"]].append((r["focal"],resid[i]))
    return [
        {"neighbor":n,"n_focal_species":len(pairs),
         "focal_adjusted_mean":float(np.mean([v for _,v in pairs])),
         "n_focals_better_than_average":sum(v>0 for _,v in pairs),
         "n_focals_worse_than_average":sum(v<0 for _,v in pairs)}
        for n,pairs in sorted(by.items(),
                               key=lambda kv:np.mean([q[1] for q in kv[1]]))
    ]


def audit(folder:Path):
    rows,source,_=prep(folder)
    report={
        "status":"POSTHOC_SPECIES_IDENTITY_BEYOND_FUNCTIONAL_GUILD_DIAGNOSTIC",
        "source_audit":source,
        "late_vs_simultaneous":{
            "identity_beyond_guild":conditional_species_identity(rows,"RII_late_total"),
            "functional_guild_specificity":guild_mapping_placebo(rows,"RII_late_total"),
            "neighbor_rankings":extreme_neighbor_report(rows,"RII_late_total"),
            "leave_hubs_out":{}
        },
        "simultaneous_neighbor_vs_alone":{
            "identity_beyond_guild":conditional_species_identity(rows,"RII_sim_total",seed=20261126),
            "functional_guild_specificity":guild_mapping_placebo(rows,"RII_sim_total",seed=20261127),
        },
        "positive_late_pairs":[{"focal":r["focal"],"neighbor":r["neighbor"],
                                "RII_late_total":r["RII_late_total"]}
                               for r in rows if r["RII_late_total"]>0],
    }
    for exclude in (("Antennaria dioica",),("Briza media",),
                    ("Antennaria dioica","Briza media")):
        subset=[r for r in rows if r["neighbor"] not in exclude]
        x=neighbor_consistency(subset,key="RII_late_total",n_perm=1999,seed=20261128)
        report["late_vs_simultaneous"]["leave_hubs_out"]["+".join(exclude)]={
            k:v for k,v in x.items() if k!="neighbor_adjusted_means"
        }
    report["claims"]=[
        "No independent replication: 202 directed pair pots share 15 focal and 15 neighbor species",
        "Response reflects delayed planting relative to same-species simultaneous planting, not pure alpha_ij formation",
        "Source 2025 Ecology paper already reported trait-driven competition and priority effects",
        "Guild partitions biologically grounded but post hoc specificity comparisons remain exploratory",
        "No claim of directly measured pairwise competition coefficient changes without age-matched partner-free controls",
    ]
    return report


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--source-dir",type=Path,required=True)
    a.add_argument("--output",type=Path,required=True)
    args=a.parse_args()
    result=audit(args.source_dir)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("NEIGHBOR_SPECIFICITY_START")
    for stage in ("late_vs_simultaneous","simultaneous_neighbor_vs_alone"):
        x=result[stage]
        for k in ("identity_beyond_guild","functional_guild_specificity"):
            print(stage,k,json.dumps(x[k]))
    print("HUB_DELETION",json.dumps(result["late_vs_simultaneous"]["leave_hubs_out"]))
    print("NEIGHBOR_EXTREMES",json.dumps(
        result["late_vs_simultaneous"]["neighbor_rankings"][:4]+
        result["late_vs_simultaneous"]["neighbor_rankings"][-4:]))
    print("POSITIVE_LATE",json.dumps(result["positive_late_pairs"]))
    print("NEIGHBOR_SPECIFICITY_DONE")

if __name__=="__main__":main()
