#!/usr/bin/env python3
"""Post hoc, source-locked test of earlier-neighbor SOLO growth rate transfer.

DO NOT call this a discovery of the original competitive-trait hypothesis.
Contract: docs/VANSTEIJN_NEIGHBOR_GR_POSTHOC_TRANSFER_CONTRACT_20261124.md
Source: van Steijn et al. (2025) Ecology doi 10.1002/ecy.4528
"""
from __future__ import annotations
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from benchmarks.score_vansteijn_partner_effects import prep
from benchmarks.score_vansteijn_age_matched_neighbor_effect import transform
from benchmarks.predict_vansteijn_unseen_neighbor_traits import (
    traits_for_source, ridge_prediction, GUILDS, AuditStop
)

MODELS=(
    "focal_morph6",
    "focal_morph6_plus_focal_GR",
    "focal_morph6_plus_neighbor_GR",
    "focal_morph6_plus_both_GR",
    "focal_morph6_plus_neighbor_morph6",
    "focal_morph6_plus_neighbor_morph6_plus_neighbor_GR",
    "focal_morph6_plus_neighbor_morph6_plus_abs_gap6_plus_neighbor_GR",
    "guild_focal_morph6_plus_neighbor_GR",
)
LAMBDAS=(10.,1.,100.)
MAIN=MODELS[2]
BASELINE=MODELS[0]


def guild_code(species: str,traits:dict):
    return np.asarray([float(traits[species]["guild"]==g)
                       for g in GUILDS[1:]],dtype=float)


def feat(rows:list[dict],traits:dict,train_species:list[str],mode:str,
         neighbor_gr_mapping:dict[str,str]|None=None)->np.ndarray:
    zbase=np.stack([traits[sp]["traits"] for sp in train_species])
    if zbase.shape!=(14,7):
        raise AuditStop("Expected 14 training species with six morphological traits and GR")
    mu=zbase.mean(axis=0)
    sd=zbase.std(axis=0)
    sd=np.where(sd<1e-8,1.0,sd)
    if neighbor_gr_mapping is None:
        neighbor_gr_mapping={s:s for s in traits}
    out=[]
    for r in rows:
        f,n=r["focal"],r["neighbor"]
        fv=(traits[f]["traits"]-mu)/sd
        nv=(traits[n]["traits"]-mu)/sd
        ng=(traits[neighbor_gr_mapping[n]]["traits"][6]-mu[6])/sd[6]
        morph_f=fv[:6]
        morph_n=nv[:6]
        if mode=="focal_morph6":
            x=morph_f
        elif mode=="focal_morph6_plus_focal_GR":
            x=np.r_[morph_f,fv[6]]
        elif mode=="focal_morph6_plus_neighbor_GR":
            x=np.r_[morph_f,ng]
        elif mode=="focal_morph6_plus_both_GR":
            x=np.r_[morph_f,fv[6],ng]
        elif mode=="focal_morph6_plus_neighbor_morph6":
            x=np.r_[morph_f,morph_n]
        elif mode=="focal_morph6_plus_neighbor_morph6_plus_neighbor_GR":
            x=np.r_[morph_f,morph_n,ng]
        elif mode=="focal_morph6_plus_neighbor_morph6_plus_abs_gap6_plus_neighbor_GR":
            x=np.r_[morph_f,morph_n,np.abs(morph_f-morph_n),ng]
        elif mode=="guild_focal_morph6_plus_neighbor_GR":
            x=np.r_[guild_code(f,traits),guild_code(n,traits),morph_f,ng]
        else:
            raise AuditStop("Unfrozen predictor "+mode)
        out.append(x)
    return np.stack(out)


def score(rows:list[dict],traits:dict,models=MODELS,penalty=10.,
          neighbor_gr_mapping=None)->dict:
    taxon=sorted(traits)
    if len(taxon)!=15:raise AuditStop("Expected 15 taxa")
    Y=np.asarray([r["Cdelta_total"] for r in rows],dtype=float)
    preds={m:np.full(len(rows),np.nan) for m in models}
    pscore={}
    for withholding in taxon:
        tr=np.array([i for i,r in enumerate(rows) if
                     r["focal"]!=withholding and r["neighbor"]!=withholding])
        te=np.array([i for i,r in enumerate(rows) if r["neighbor"]==withholding])
        if len(te)<10 or len(tr)<150:
            raise AuditStop("Invalid directed species holdout")
        train_species=[z for z in taxon if z!=withholding]
        for m in models:
            a=feat([rows[i] for i in tr],traits,train_species,m,neighbor_gr_mapping)
            b=feat([rows[i] for i in te],traits,train_species,m,neighbor_gr_mapping)
            preds[m][te]=ridge_prediction(a,Y[tr],b,penalty)
    if any(not np.isfinite(p).all() for p in preds.values()):
        raise AuditStop("Incomplete species-out predictions")
    species_results={}
    for name in taxon:
        mask=np.asarray([r["neighbor"]==name for r in rows])
        species_results[name]={
            "n_directed_test_pots":int(mask.sum()),
            "MSE":{m:float(np.mean((preds[m][mask]-Y[mask])**2)) for m in models}
        }
    summaries={}
    for m in models:
        mse=np.mean((preds[m]-Y)**2)
        species_mse=np.mean([species_results[s]["MSE"][m] for s in taxon])
        summaries[m]={
            "directed_pair_RMSE":float(np.sqrt(mse)),
            "species_equal_MSE":float(species_mse),
            "species_equal_RMSE":float(np.sqrt(species_mse)),
        }
    return {"lambda":penalty,"n_species":len(taxon),"n_directed_pairs":len(rows),
            "model_scores":summaries,"per_heldout_neighbor":species_results}


def paired_bootstrap(score_result:dict,candidate:str,base:str=BASELINE,
                     n_boot:int=20000,seed:int=20261124):
    names=sorted(score_result["per_heldout_neighbor"])
    dif=np.array([score_result["per_heldout_neighbor"][s]["MSE"][base]-
                  score_result["per_heldout_neighbor"][s]["MSE"][candidate]
                  for s in names])
    rng=np.random.default_rng(seed)
    res=rng.choice(dif,size=(n_boot,len(dif)),replace=True).mean(axis=1)
    return {
        "candidate":candidate,"baseline":base,
        "mean_MSE_gain_positive_is_better":float(np.mean(dif)),
        "bootstrap_95pct":np.quantile(res,[.025,.975]).tolist(),
        "n_heldout_species_improved":int((dif>0).sum()),
        "n_heldout_species_worse":int((dif<0).sum()),
        "total_heldout_species":len(names),
        "excluding_Antennaria_species_mean_MSE_gain":float(
            np.mean([v for n,v in zip(names,dif) if n!="Antennaria dioica"])),
        "species_level_gains":{n:float(v) for n,v in zip(names,dif)},
        "unit":"early-neighbor species, bootstrap is descriptive; 15 source taxa not random flora sample"
    }


def permute_GR_within_guild(traits:dict,rng:np.random.Generator):
    mapping={}
    for guild in GUILDS:
        group=sorted(s for s in traits if traits[s]["guild"]==guild)
        mapping.update(zip(group,rng.permutation(group).tolist()))
    return mapping


def test_GR_mapping(rows:list[dict],traits:dict,n_perm:int=499,seed:int=20261124):
    base=score(rows,traits,models=(MAIN,))["model_scores"][MAIN]["species_equal_MSE"]
    rng=np.random.default_rng(seed)
    null=np.empty(n_perm)
    for i in range(n_perm):
        mapping=permute_GR_within_guild(traits,rng)
        x=score(rows,traits,models=(MAIN,),neighbor_gr_mapping=mapping)
        null[i]=x["model_scores"][MAIN]["species_equal_MSE"]
    return {
        "observed_MSE":base,
        "within_guild_shuffled_neighbor_GR_MSE_median":float(np.median(null)),
        "within_guild_shuffled_neighbor_GR_MSE_95pct":np.quantile(null,[.025,.975]).tolist(),
        "fraction_shuffles_as_good_or_better":float((1+np.sum(null<=base+1e-12))/(n_perm+1)),
        "n_permutations":n_perm,
        "seed":seed,
        "caution":"Posthoc mapping specificity diagnostic, not source randomized trait intervention"
    }


def audit(source:Path,n_perm:int=499):
    pair,source_meta,_=prep(source)
    rows=transform(pair)
    t,trait_names=traits_for_source(source,with_growth_rate=True)
    s=score(rows,t,penalty=10.)
    controls=[c for c in MODELS if c!=BASELINE]
    contrast={c:paired_bootstrap(s,c,n_boot=20000,seed=20261124) for c in controls}
    permutations=test_GR_mapping(rows,t,n_perm=n_perm)
    sens={}
    for lam in (1.,100.):
        fit=score(rows,t,models=(BASELINE,MAIN),penalty=lam)
        sens[str(lam)]={
            "scores":fit["model_scores"],
            "comparison":paired_bootstrap(fit,MAIN,n_boot=20000,seed=20261124)}
    return {
        "status":"POSTHOC_PREDICTIVE_DECOMPOSITION_OF_SOURCE_GROWTH_RATE_NOT_CAUSAL",
        "source":"van Steijn et al. 2025 Ecology DOI 10.1002/ecy.4528, Zenodo 10.5281/zenodo.10512662",
        "protocol":"docs/VANSTEIJN_NEIGHBOR_GR_POSTHOC_TRANSFER_CONTRACT_20261124.md",
        "source_audit":source_meta,
        "tested_features":list(trait_names),
        "metric":"deltaC=pre-established neighbor net effect minus simultaneous neighbor effect, focal 10 weeks in all arms",
        "cv":"15 species folds: heldout neighbor species absent from ALL training pairs even if focal",
        "lambda_main":10.,
        "primary_results":s,
        "paired_15_species_bootstraps":contrast,
        "neighbor_GR_mapping_null":permutations,
        "alternative_penalty_sensitivities":sens,
        "limits":[
            "GR_j from 10-week solo-grown plant, NOT from an ex-ante field trait sample",
            "GR_j does not enter Cdelta_i control denominator for j!=i; focal GR_i could enter mathematically and is separated",
            "Ridge penalty fixed prior to this source scoring, penalty values 1 and 100 only sensitivity",
            "15 nonrandom species folds are not 15 independent experimental studies",
            "Trait modeling after finding a sensitivity signal is post hoc and must not be called preregistered confirmatory evidence",
            "No direct root/soil access or niche pre-emption mechanisms measured",
            "Frozen KNK v28–v38 and existing manuscripts unchanged",
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=499)
    args=p.parse_args()
    report=audit(args.source_dir,n_perm=args.permutations)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print("NEIGHBOR_GR_DECOMPOSITION_START")
    print("MODELS",json.dumps(report["primary_results"]["model_scores"]))
    print("MAIN_MSE_PAIR_CONTRAST",json.dumps(report["paired_15_species_bootstraps"][MAIN]))
    print("MAPPING_NULL",json.dumps(report["neighbor_GR_mapping_null"]))
    print("PENALTY_SENSITIVITY",json.dumps(report["alternative_penalty_sensitivities"]))
    print("NEIGHBOR_GR_DECOMPOSITION_END",args.output)


if __name__=="__main__":main()
