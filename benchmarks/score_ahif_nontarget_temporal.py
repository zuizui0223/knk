#!/usr/bin/env python3
"""Posthoc field non-target temporal dynamics; source-lock, no mediation claim."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_vansteijn2026_field import (
    source_panels, blocked_assignment_stream, score_matrix, rank_vectors,
    SourceStop, SITES, SEED
)
from benchmarks.probe_ahif_vansteijn2026_field import MD5

def holm_values(p):
    p = np.asarray(p, dtype=float)
    if np.any(~np.isfinite(p)) or np.any((p < 0) | (p > 1)):
        raise ValueError("invalid p")
    order = np.argsort(p, kind="stable")
    result = np.empty(len(p), dtype=float)
    result[order] = np.minimum(1., np.maximum.accumulate(p[order] * np.arange(len(p), 0, -1)))
    return result.tolist()

def roster_data(state):
    names = state["names"]
    sown = set(state["sown_target"])
    x = np.asarray(state["X"],float)
    if x.shape != (20,76) or len(names)!=76 or len(sown)!=20 or len(set(names)&sown)!=18:
        raise SourceStop("Unexpected fixed source roster")
    other = x[:,[i for i,n in enumerate(names) if n not in sown]]
    target = np.column_stack([x[:,names.index(n)] if n in names else np.zeros(20) for n in sorted(sown)])
    if other.shape!=(20,58) or target.shape!=(20,20):
        raise SourceStop("Roster mismatch")
    n = other.sum(axis=1); t = target.sum(axis=1)
    if np.any(n<=0) or np.any(t<=0):
        raise SourceStop("Zero roster total; do not drop plots or add pseudocount")
    pn=other/n[:,None]; pt=target/t[:,None]
    rank=rank_vectors(pn)
    return {
        "nontarget_rank":rank,
        "nontarget_named":pn,
        "sown_rank":rank_vectors(pt),
        "sown_named":pt,
        "nontarget_richness":(other>0).sum(axis=1)[:,None],
        "nontarget_top1":rank[:,:1],
        "nontarget_simpson_concentration":np.sum(pn**2,axis=1)[:,None],
        "nontarget_total_pin_hits":n[:,None],
        "nontarget_share_all_hits":(n/(n+t))[:,None],
    }

def changes(panel,site,end="2024"):
    a=panel[site,"2021"]; b=panel[site,end]
    if a["plots"]!=b["plots"] or a["names"]!=b["names"] or not np.array_equal(a["history"],b["history"]):
        raise SourceStop("Repeated physical plot linkage mismatch")
    x=roster_data(a); y=roster_data(b)
    return {k:y[k]-x[k] for k in x},a["history"]

def pair_differences(x,labels):
    if len(labels)!=20 or len(x)!=20: raise SourceStop("Invalid complete blocks")
    rows=[]
    for i in range(0,20,4):
        g=list(labels[i:i+4])
        if set(g)!={"FtS","StF","SFearly","SFlate"}: raise SourceStop("Invalid four-treatment block")
        rows.append(x[i+g.index("FtS")]-x[i+g.index("StF")])
    return np.asarray(rows)

def exact_joint_swap(diffs):
    """Conditional sharp null of no effect of exchanging two asynchronous labels."""
    if len(diffs)!=2 or any(d.ndim!=2 or d.shape[0]!=5 for d in diffs):
        raise SourceStop("Exactly five pairs at each of two source sites expected")
    def site_stat(d,s):
        den=float(np.sum(d*d))
        return 0. if den<=1e-15 else float(np.sum(np.mean(d*s[:,None],axis=0)**2)/den)
    obs=sum(site_stat(d,np.ones(5)) for d in diffs)
    vals=[]
    for sign in itertools.product((-1.,1.),repeat=10):
        s=np.asarray(sign)
        vals.append(site_stat(diffs[0],s[:5])+site_stat(diffs[1],s[5:]))
    return {"statistic":obs,"p_exact":float(np.mean(np.asarray(vals)>=obs-1e-12)),
            "n_exact_assignments":1024,"n_pairs_by_site":[5,5],
            "site_mean_difference_l2":[float(np.linalg.norm(np.mean(d,axis=0))) for d in diffs]}

def analyse(panel,n_perm=9999):
    primary=[]; secondary=[]; contrasts=[]
    delta={s:changes(panel,s) for s in SITES}
    for i,s in enumerate(SITES):
        x,y=delta[s]
        perms=blocked_assignment_stream(y,n_perm=n_perm,seed=SEED+301+i)
        for k,v in x.items():
            row={"site":s,"endpoint":k,"period":"2024_minus_2021",
                 **score_matrix(v,y,perms,metric="euclidean")}
            (primary if k=="nontarget_rank" else secondary).append(row)
    for rows in (primary,secondary):
        for row,p in zip(rows,holm_values([z["p_upper"] for z in rows])):
            row["p_holm_family"]=p
    for end,k in (("2024","nontarget_rank"),("2024","nontarget_named"),("2023","nontarget_rank")):
        pairs=[]
        for s in SITES:
            x,y=changes(panel,s,end=end)
            pairs.append(pair_differences(x[k],y))
        contrasts.append({"period":end+"_minus_2021","endpoint":k,**exact_joint_swap(pairs)})
    for row,p in zip(contrasts,holm_values([z["p_exact"] for z in contrasts])):
        row["p_holm_family"]=p
    return {"status":"SOURCE_REANALYSIS_EXPLORATORY_NOT_CAUSAL_MEDIATION",
            "n_physical_plots":40,"n_publication_clusters":1,"primary":primary,
            "secondary":secondary,"order_exchange":contrasts,
            "limits":"Same source plot trajectories, not independent years or proof of direct species interactions."}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=9999)
    a=p.parse_args()
    for name,md5 in MD5.items():
        if hashlib.md5((a.original_dir/name).read_bytes()).hexdigest()!=md5:
            raise SourceStop("Original hash mismatch: "+name)
    panel,_=source_panels(a.original_dir)
    result=analyse(panel,a.permutations)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("NON_TARGET_TEMPORAL_RESULT",json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
