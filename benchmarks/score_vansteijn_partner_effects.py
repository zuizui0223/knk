#!/usr/bin/env python3
"""15-species directed interaction-response matrix in van Steijn et al. 2025.

Source doi:10.5281/zenodo.10512662, Ecology doi:10.1002/ecy.4528.
Contract committed earlier: docs/VANSTEIJN_PAIR_INTERACTION_PRE_SCORE_20261123.md

IMPORTANT: RII_late compares later-planted and simultaneous-planted focal
biomass in mixed pairs. It is NOT a mechanistic interspecific alpha coefficient,
because planting schedules and age also differ.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np

from benchmarks.probe_vansteijn_competition import audit as audit_sources

SEED=20261123
N_PERM=9999
N_BOOT=20000


class SourceStop(ValueError):
    pass


def as_positive_number(raw: str, source: str, allow_zero: bool = True) -> float:
    try:
        x=float(raw)
    except (ValueError,TypeError) as e:
        raise SourceStop(f"Non-numeric biomass {raw!r} in {source}") from e
    if not math.isfinite(x) or x<0 or (not allow_zero and x==0):
        raise SourceStop(f"Invalid biomass {raw!r} in {source}")
    return x


def read_tab(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f,delimiter=";"))
    if not rows:
        raise SourceStop(f"Empty data file: {path.name}")
    return rows


def prep(source_dir: Path):
    source=audit_sources(source_dir)
    alone=read_tab(source_dir/"alone.csv")
    sim=read_tab(source_dir/"simultaneous_arrival.csv")
    seq=read_tab(source_dir/"second_arrival.csv")
    trait=read_tab(source_dir/"traits.csv")
    guild={r["focal"]:r["func.group"] for r in trait}
    if len(guild)!=15 or Counter(guild.values())!=Counter({"Grass":5,"Forb":5,"Legume":5}):
        raise SourceStop("15-species trait guild balance mismatch")
    if len(alone)!=74 or len(sim)!=206 or len(seq)!=206:
        raise SourceStop("Original dataset row count changed")
    organ={"total":("biomass.total.CA","biomass.total.CT","biomass.total.PE"),
           "above":("biomass.above.CA","biomass.above.CT","biomass.above.PE"),
           "below":("biomass.below.CA","biomass.below.CT","biomass.below.PE")}
    controls=defaultdict(list)
    for r in alone:
        focal=r["focal"]
        if focal not in guild: raise SourceStop(f"Unknown control focal {focal}")
        controls[focal].append({k:as_positive_number(r[col[0]],"alone.csv")
                                for k,col in organ.items()})
    if len(controls)!=15 or min(map(len,controls.values()))<4:
        raise SourceStop("Controls missing focal or insufficient replication")
    controls_mean={f:{k:float(np.mean([q[k] for q in v])) for k in organ}
                   for f,v in controls.items()}
    def index_source(rows,kind):
        out={}
        for r in rows:
            f=r["focal"]; n=r["neighbour"]
            if f not in guild or n not in guild or f==n:
                raise SourceStop("Unexpected taxa or self pair in "+kind)
            key=f,n
            if key in out:
                raise SourceStop("Duplicate directed focal-neighbor combination in "+kind)
            out[key]={k:as_positive_number(r[organ[k][1 if kind=="sim" else 2]],kind)
                      for k in organ}
        return out
    S=index_source(sim,"sim")
    Q=index_source(seq,"late")
    both=sorted(set(S)&set(Q))
    if len(both)<190:
        raise SourceStop(f"Insufficient directed-pair alignment: {len(both)}")
    unmatched={"sim_only":[list(x) for x in sorted(set(S)-set(Q))],
               "late_only":[list(x) for x in sorted(set(Q)-set(S))]}
    scored=[]
    undefined=Counter()
    for f,n in both:
        s=S[(f,n)];q=Q[(f,n)]; c=controls_mean[f]
        rec={"focal":f,"neighbor":n,"focal_guild":guild[f],
             "neighbor_guild":guild[n]}
        for k in organ:
            d1=s[k]+c[k]
            d2=q[k]+s[k]
            if d1<=0 or d2<=0:
                undefined[f"RII_{k}"]+=1
                rec[f"RII_sim_{k}"]=None if d1<=0 else (s[k]-c[k])/d1
                rec[f"RII_late_{k}"]=None if d2<=0 else (q[k]-s[k])/d2
            else:
                rec[f"RII_sim_{k}"]=(s[k]-c[k])/d1
                rec[f"RII_late_{k}"]=(q[k]-s[k])/d2
        scored.append(rec)
    audit={
        "status":"SOURCE_ALIGN_PASS",
        "n_species":15,"n_directed_sim":len(S),"n_directed_late":len(Q),
        "n_joined_directed_pairs":len(scored),
        "controls_per_species":{k:len(v) for k,v in sorted(controls.items())},
        "unmatched_pairs":unmatched,
        "undefined_rii":dict(undefined),
        "source_files":{
            f:v["md5"] for f,v in source["files"].items()
        },
        "unit":"single directed-pair pot; 15 focal species as higher-order resampling clusters",
    }
    return scored, audit, controls_mean


def mean_ci_cluster(x:np.ndarray, seed:int=SEED, n_boot:int=N_BOOT):
    x=np.asarray(x,dtype=float)
    rng=np.random.default_rng(seed)
    draw=rng.choice(x,size=(n_boot,len(x)),replace=True).mean(axis=1)
    return {"mean":float(x.mean()),
            "bootstrap_percentile95":np.quantile(draw,[.025,.975]).tolist(),
            "n_clusters":len(x),"n_boot":n_boot}


def exact_signflip(x:np.ndarray)->float:
    """Exact two-sided sign randomization over independent focal clusters.

    This is a diagnostic under symmetric species-cluster null, not a claim of
    15 independent ecological taxa sampled randomly from a superpopulation.
    """
    x=np.asarray(x,dtype=float)
    n=len(x)
    if n>17: raise ValueError("Too many exact flips")
    target=abs(x.mean())
    count=0
    total=2**n
    for bits in range(total):
        signed=sum(x[i]*(1 if bits & (1<<i) else -1) for i in range(n))/n
        if abs(signed)>=target-1e-13: count+=1
    return float(count/total)


def focal_cluster_test(rows:list[dict],response:str)->dict:
    by=defaultdict(list)
    for r in rows:
        v=r.get(response)
        if v is not None and math.isfinite(v):by[r["focal"]].append(v)
    if len(by)!=15:
        raise SourceStop(f"Not all 15 focal species have {response}")
    x=np.array([np.mean(by[s]) for s in sorted(by)])
    d=mean_ci_cluster(x)
    d.update({"exact_focal_signflip_p_two_sided":exact_signflip(x),
              "n_obs_nonindependent":sum(len(a) for a in by.values()),
              "n_focal_species":len(by),
              "focal_species_means":{k:float(np.mean(v)) for k,v in sorted(by.items())},
              "sign_of_focal_means":{
                  "negative":int(np.sum(x<0)),"positive":int(np.sum(x>0)),
                  "zero":int(np.sum(x==0))}})
    return d


def neighbor_consistency(rows:list[dict],key:str="RII_late_total",n_perm:int=N_PERM,seed:int=SEED):
    """Assess whether early-neighbor name gives consistent response across focal species.

    Shuffle neighbors *within focal species*. Not an exact randomization p
    unless actual species-pair treatment assignment was randomized/exchangeable.
    """
    clean=[r for r in rows if r.get(key) is not None]
    foci=sorted({r["focal"] for r in clean})
    neigh=sorted({r["neighbor"] for r in clean})
    fi={t:i for i,t in enumerate(foci)}
    ni={t:i for i,t in enumerate(neigh)}
    yi=np.array([fi[r["focal"]] for r in clean])
    xi=np.array([ni[r["neighbor"]] for r in clean])
    response=np.array([r[key] for r in clean],float)
    residual=response-np.array([response[yi==i].mean() for i in range(len(foci))])[yi]
    denominator=float(np.sum(residual**2))
    if denominator<=1e-14:raise SourceStop("No variation after focal centering")
    def conditional_nbr_r2(labels):
        n=np.bincount(labels,minlength=len(neigh))
        sums=np.bincount(labels,weights=residual,minlength=len(neigh))
        effect=np.divide(sums,n,out=np.zeros_like(sums),where=n>0)
        return float(np.sum(n*effect**2)/denominator)
    actual=conditional_nbr_r2(xi)
    rng=np.random.default_rng(seed)
    positions=[np.flatnonzero(yi==i) for i in range(len(foci))]
    null=np.empty(n_perm)
    for b in range(n_perm):
        labels=xi.copy()
        for pos in positions:
            labels[pos]=rng.permutation(labels[pos])
        null[b]=conditional_nbr_r2(labels)
    neighbor_means={}
    for n,idx in ni.items():
        neighbor_means[n]={
            "n_dyads":int(np.sum(xi==idx)),
            "mean_focal_adjusted_late_RII":float(np.mean(residual[xi==idx])),
        }
    return {"R2_neighbor_conditioned_on_focal":actual,
            "null_median":float(np.median(null)),
            "null_95pct":np.quantile(null,[.025,.975]).tolist(),
            "p_upper_pseudorandomization":float((1+np.count_nonzero(null>=actual-1e-12))/(n_perm+1)),
            "permutations":n_perm,"seed":seed,
            "neighbor_adjusted_means":neighbor_means,
            "caution":"Not an exact treatment-randomization p; identity-label-shuffling conditional diagnostic."}


def guild_matrix(rows:list[dict],key:str):
    out=[]
    for f in ("Grass","Forb","Legume"):
        for n in ("Grass","Forb","Legume"):
            v=np.array([r[key] for r in rows if r.get(key) is not None and
                         r["focal_guild"]==f and r["neighbor_guild"]==n])
            out.append({"focal_guild":f,"neighbor_guild":n,"n_directed_pair_pots":len(v),
                        "mean":float(v.mean()) if len(v) else None,
                        "median":float(np.median(v)) if len(v) else None,
                        "proportion_negative":float(np.mean(v<0)) if len(v) else None})
    return out


def reciprocal_stats(rows:list[dict],key:str):
    score={(r["focal"],r["neighbor"]):r.get(key) for r in rows}
    pairs=[]
    for f,n in sorted(score):
        if f>=n or (n,f) not in score:continue
        a,b=score[(f,n)],score[(n,f)]
        if a is None or b is None:continue
        pairs.append((f,n,a,b))
    v=np.array([[p[2],p[3]] for p in pairs])
    return {
        "unordered_pairs_with_both_directions":len(v),
        "mean_abs_dyadic_asymmetry":float(np.mean(np.abs(v[:,0]-v[:,1]))),
        "n_opposite_signs":int(np.sum((v[:,0]>0)&(v[:,1]<0))+np.sum((v[:,0]<0)&(v[:,1]>0))),
        "fraction_opposite_sign":float(np.mean((v[:,0]*v[:,1])<0)),
        "correlation_between_directions":float(np.corrcoef(v[:,0],v[:,1])[0,1]),
        "caution":"Each ordered pair measured once; asymmetry is descriptive and can be noisy."}


def predictive_loo_unordered(rows:list[dict],key:str="RII_late_total",penalty:float=4.0):
    """Leave both directed orientations of a species pair out of training.

    Does neighbor identity help predict partner-specific schedule costs beyond
    which focal species is responding? No model tuning or external validation.
    """
    clean=[r for r in rows if r.get(key) is not None]
    taxa=sorted({r["focal"] for r in clean}|{r["neighbor"] for r in clean})
    ix={s:i for i,s in enumerate(taxa)}
    y=np.array([r[key] for r in clean],float)
    f=np.array([ix[r["focal"]] for r in clean])
    n=np.array([ix[r["neighbor"]] for r in clean])
    ids=[tuple(sorted((r["focal"],r["neighbor"]))) for r in clean]
    by=defaultdict(list)
    for i,k in enumerate(ids):by[k].append(i)
    def features(include_neighbor):
        dim=len(taxa) * (2 if include_neighbor else 1)
        X=np.zeros((len(y),dim))
        X[np.arange(len(y)),f]=1
        if include_neighbor:X[np.arange(len(y)),n+len(taxa)]=1
        return X
    def estimate(X,tr,te):
        A=np.column_stack([np.ones(len(tr)),X[tr]])
        B=np.column_stack([np.ones(len(te)),X[te]])
        R=np.diag(np.r_[0.,np.full(X.shape[1],penalty)])
        return B@np.linalg.solve(A.T@A+R,A.T@y[tr])
    out={}
    keys=sorted(by)
    for mode in ("focal_only","focal_and_neighbor"):
        X=features(mode=="focal_and_neighbor")
        pred=np.zeros(len(y))
        for key_pair in keys:
            te=np.asarray(by[key_pair])
            train=np.asarray([i for i in range(len(y)) if i not in by[key_pair]])
            pred[te]=estimate(X,train,te)
        errors=y-pred
        out[mode]={"mse":float(np.mean(errors**2)),
                   "rmse":float(np.sqrt(np.mean(errors**2))),
                   "r2_against_training_global_mean":float(1-np.sum(errors**2)/np.sum((y-y.mean())**2))}
    out["gain_in_MSE_of_focal_plus_neighbor"]=float(out["focal_only"]["mse"]-out["focal_and_neighbor"]["mse"])
    out["heldout_unordered_pairs"]=len(keys)
    out["penalty_fixed"]=penalty
    out["caution"]="Source has one pot per ordered dyad; cross validation does not validate competition causality."
    return out


def run(folder:Path,n_perm:int=N_PERM):
    rows,source,controls=prep(folder)
    primary=focal_cluster_test(rows,"RII_late_total")
    # The original publication analyses above and below responses separately;
    # this paired contrast checks whether the same focal species changes organ.
    diffrows=[]
    for r in rows:
        x,y=r["RII_late_above"],r["RII_late_below"]
        if x is not None and y is not None:
            diffrows.append({**r,"late_above_minus_below":x-y})
    diff=focal_cluster_test(diffrows,"late_above_minus_below")
    sorted_tot=np.array([r["RII_late_total"] for r in rows if r["RII_late_total"] is not None])
    refs={
        "total_late_pairwise_rows_nonindependent":len(sorted_tot),
        "proportion_negative":float(np.mean(sorted_tot<0)),
        "proportion_positive":float(np.mean(sorted_tot>0)),
        "proportion_zero":float(np.mean(sorted_tot==0)),
        "mean_unweighted_pairs":float(np.mean(sorted_tot)),
        "median_unweighted_pairs":float(np.median(sorted_tot)),
    }
    return {
        "status":"SOURCE_VERIFIED_INTERACTION_RESPONSE_MATRIX_NOT_DIRECT_CAUSAL_COEFFICIENTS",
        "source_audit":source,
        "analysis_contract":"docs/VANSTEIJN_PAIR_INTERACTION_PRE_SCORE_20261123.md",
        "focal_cluster_mean_late_response":primary,
        "focal_cluster_shoot_minus_root_response":diff,
        "late_response_distribution":refs,
        "neighbor_identity_consistency":neighbor_consistency(rows,n_perm=n_perm),
        "reciprocal_pair_asymmetry":reciprocal_stats(rows,"RII_late_total"),
        "directed_guild_3x3_late_response":guild_matrix(rows,"RII_late_total"),
        "directed_guild_3x3_simultaneous_competition":guild_matrix(rows,"RII_sim_total"),
        "heldout_pair_prediction":predictive_loo_unordered(rows),
        "caveats":[
            "These are original single-pot directed pair effects, not 206 experimental replications of one contrast.",
            "Original 2025 Ecology paper already interpreted competitive response, priority effects and traits.",
            "RII_late compares later-versus-simultaneous planting with different timing; not an unbiased partner-only interaction coefficient.",
            "Only 15 nonrandom species sampled; sign flips assume symmetric focal-species variation and no complete phylogenetic independence.",
            "Relative species-pair biomass responses do not establish stability of coexistence.",
            "No changes to frozen KNK v28-v38 or Draft PR #3/#4."
        ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source-dir",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    data=run(args.source_dir)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("VANSTEIJN_RESULT_START")
    print("AUDIT",json.dumps(data["source_audit"]))
    for key in ("focal_cluster_mean_late_response","focal_cluster_shoot_minus_root_response",
                "late_response_distribution","neighbor_identity_consistency",
                "reciprocal_pair_asymmetry","heldout_pair_prediction"):
        print("SUMMARY_"+key,json.dumps(
            data[key] if key!="neighbor_identity_consistency" else
            {k:v for k,v in data[key].items() if k!="neighbor_adjusted_means"}))
    print("GUILD_3x3",json.dumps(data["directed_guild_3x3_late_response"]))
    print("VANSTEIJN_RESULT_END",args.output)


if __name__=="__main__":
    main()
