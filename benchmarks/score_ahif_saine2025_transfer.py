"""Cross-site and cross-year prediction of inoculum-specific non-target response vectors."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from benchmarks.score_ahif_saine2025_succession import design,SEED,holm

GRID=np.array([0.,.25,.5,.75,1.])
N_PERM=999

def group_response(weights,changes,control):
    # weights: one row per site/context/introduced identity, equal within arm
    # changes: one matrix for each time endpoint, CSR for sparse OTU response
    return np.stack([(weights@x).toarray().reshape(control[k].shape)-control[k]
                     for k,x in enumerate(changes)])

def centered(a):
    return a-a.mean(axis=2,keepdims=True)

def loss(a):
    return float(np.sum(a*a))

def cv(a):
    """a shape (year,site,context,species,feature), year 1 training -> year 2 target."""
    if a.ndim!=5 or a.shape[0]!=2 or a.shape[1]<4 or a.shape[2]!=2 or a.shape[3]<2:
        raise ValueError("CV requires matched sites, contexts and species")
    if not np.isfinite(a).all():raise ValueError("nonfinite response")
    x,y=a
    n=len(x);xc=centered(x);yc=centered(y)
    rows=[]
    for held in range(n):
        train=[i for i in range(n) if i!=held]
        xy=xx=0.
        for v in train:
            inner=[i for i in train if i!=v]
            pred=xc[inner].mean(axis=0)
            xy+=float(np.sum(pred*yc[v]))
            xx+=loss(pred)
        alpha=float(GRID[np.argmin(GRID*GRID*xx-2*GRID*xy)])
        avg=x[train].mean(axis=0)
        pooled=avg.mean(axis=1,keepdims=True)
        specific=avg-pooled
        rows.append({"held_site_index":held,"alpha":alpha,
          "zero_sse":loss(y[held]),"pooled_sse":loss(y[held]-pooled),
          "nested_identity_sse":loss(y[held]-pooled-alpha*specific),
          "full_identity_sse":loss(y[held]-avg)})
    sums={k:sum(r[k] for r in rows) for k in ("zero_sse","pooled_sse","nested_identity_sse","full_identity_sse")}
    def gain(base,model):
        return (base-model)/base if base>1e-15 else 0.
    for r in rows:
        r["identity_gain_over_pooled"]=gain(r["pooled_sse"],r["nested_identity_sse"])
    return {"folds":rows,"total_sse":sums,
      "identity_gain_over_pooled":gain(sums["pooled_sse"],sums["nested_identity_sse"]),
      "pooled_gain_over_zero":gain(sums["zero_sse"],sums["pooled_sse"]),
      "full_identity_gain_over_pooled":gain(sums["pooled_sse"],sums["full_identity_sse"])}

def prepare(meta,x):
    ix,t,strata,sites,types,good,groups,reads=design(meta,x)
    if int(good.sum())!=446:raise ValueError("Complete cohort changed")
    raw=x[ix[good]]
    p=raw/raw.sum(2,keepdims=True)
    rank=np.sort(p,axis=2)[:,:,::-1]
    k=int(np.max(np.sum(raw>0,axis=2)))
    t=t[good];sites=sites[good];types=types[good]
    snames=sorted(set(sites));spnames=sorted(set(t)-{"control"})
    if len(snames)!=5 or len(spnames)!=9 or int((t!="control").sum())!=347:
        raise ValueError("Unexpected complete cohort")
    tr=t!="control"; n=int(tr.sum())
    grp=np.asarray([snames.index(s)*2+(0 if lt=="R" else 1) for s,lt in zip(sites,types)])
    trgrp=grp[tr]
    ids=np.array([spnames.index(a) for a in t[tr]],int)
    blocks=[np.flatnonzero(trgrp==g) for g in range(10)]
    reps={}
    for label,v in (("named",p),("rank",rank[:,:,:k])):
        changes=[v[:,j]-v[:,0] for j in (1,2)]
        controls=np.stack([np.asarray([q[(grp==g)&~tr].mean(0) for g in range(10)])
                           for q in changes]).reshape(2,5,2,1,-1)
        controls=np.broadcast_to(controls,(2,5,2,9,v.shape[2])).copy()
        reps[label]=([csr_matrix(q[tr]) for q in changes],controls)
    counts=np.bincount(trgrp*9+ids,minlength=90)
    if len(counts)!=90 or np.any(counts<=0):raise ValueError("Lost introduced identity in a stratum")
    return ids,trgrp,blocks,counts,reps,snames,spnames

def run(meta,x,nperm=N_PERM):
    ids,grp,blocks,counts,reps,sites,species=prepare(meta,x)
    col=np.arange(len(ids))
    def score(labels):
        row=grp*9+labels
        if not np.array_equal(np.bincount(row,minlength=90),counts):
            raise ValueError("Permutation changed original complete arm sizes")
        w=csr_matrix((1./counts[row],(row,col)),shape=(90,len(ids)))
        return {key:cv(group_response(w,change,control)) for key,(change,control) in reps.items()}
    obs=score(ids)
    rng=np.random.default_rng(SEED+1201); null={k:[] for k in obs}
    for _ in range(nperm):
        labels=ids.copy()
        for idx in blocks:labels[idx]=rng.permutation(labels[idx])
        out=score(labels)
        for k in obs:null[k].append(out[k]["identity_gain_over_pooled"])
    for k,o in obs.items():
        n=np.asarray(null[k]);v=o["identity_gain_over_pooled"]
        o["p_upper"]=float((1+np.count_nonzero(n>=v-1e-12))/(nperm+1))
        o["null_gain_median"]=float(np.median(n))
        o["null_gain_95"]=np.quantile(n,[.025,.975]).tolist()
        for f in o["folds"]:f["held_site"]=sites[f.pop("held_site_index")]
    for k,p in zip(obs,holm([o["p_upper"] for o in obs.values()])):
        obs[k]["p_holm_family"]=p
    return {"status":"EXPLORATORY_NESTED_SITE_HOLDOUT_NOT_DIRECT_INTERACTION_COEFFICIENT",
      "source_assigned_logs":460,"complete_logs":446,"treated_logs":347,"control_logs":99,
      "site_names":sites,"context_order":["felled","natural"],"introduced_species":species,
      "n_permutations":nperm,"alpha_grid":GRID.tolist(),"results":obs,
      "limits":["Same source experiment, five fixed sites; no new independent publication.",
                "Controls shared across species; whole trajectory assignment permutations preserve dependence.",
                "Conditional on post-treatment QC, not exact full-cohort ITT.",
                "Training site year-2 outcomes select shrinkage; held-out site excluded from all fitting.",
                "Response vectors are not independent direct interaction coefficients."]}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=N_PERM)
    a=p.parse_args()
    with (a.source/"meta.csv").open() as f:meta=list(csv.DictReader(f))
    names=(a.source/"otu_names.txt").read_text().splitlines()
    if len(names)!=4653 or len(set(names))!=4653 or len((a.source/"excluded_target_otus.txt").read_text().splitlines())!=9:
        raise ValueError("Source roster mismatch")
    x=np.fromfile(a.source/"counts.bin",dtype="<f8").reshape(1380,4653)
    z=run(meta,x,a.permutations)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SAINE_TRANSFER_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
