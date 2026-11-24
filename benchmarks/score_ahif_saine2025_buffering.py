"""Saine structural-buffering falsification: paired TV gap and fixed-depth sensitivity."""
from __future__ import annotations
import argparse,csv,json
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix
from benchmarks.score_ahif_saine2025_succession import design,BlockPermutation,adjust,SEED

DEPTHS=(1000,5000)
DRAWS=8
NPERM=1999
NTV=9999
NBOOT=4000

def tv_outcomes(x):
    x=np.asarray(x,float)
    if x.ndim!=3 or x.shape[1]!=3 or np.any(x<0) or np.any(x.sum(2)<=0):
        raise ValueError("invalid log trajectory")
    p=x/x.sum(2,keepdims=True)
    rank=np.sort(p,axis=2)[:,:,::-1]
    named=.5*np.abs(p[:,1:]-p[:,:1]).sum(2)
    sorted_tv=.5*np.abs(rank[:,1:]-rank[:,:1]).sum(2)
    gap=named-sorted_tv
    if np.any(gap<-1e-10):raise ValueError("TV rearrangement inequality")
    # Columns grouped by year, preserving same log across representations.
    return np.stack((named,sorted_tv,np.maximum(gap,0.)),axis=2).reshape(len(x),6)

def weighted_contrasts(y,t,strata,weights):
    y=np.asarray(y,float);t=np.asarray(t,bool)
    out=np.zeros(y.shape[1])
    for s,w in weights.items():
        g=strata==s
        if not np.any(g&t) or not np.any(g&~t):
            raise ValueError("missing binary arm")
        out+=w*(y[g&t].mean(0)-y[g&~t].mean(0))
    return out

def tv_score(y,t,strata,weights,seed,nperm=NTV,nboot=NBOOT):
    rng=np.random.default_rng(seed); t=np.asarray(t,bool)
    y=np.asarray(y,float)
    if y.shape!=(len(t),6) or not np.isfinite(y).all():raise ValueError("TV endpoint shape")
    if set(strata)!=set(weights) or not np.isclose(sum(weights.values()),1):
        raise ValueError("original stratum weights mismatch")
    obs=weighted_contrasts(y,t,strata,weights)
    null=np.zeros((nperm,6));boot=np.zeros((nboot,6))
    for s,w in weights.items():
        g=np.flatnonzero(strata==s);a=g[t[g]];b=g[~t[g]]
        if not len(a) or not len(b):raise ValueError("missing binary arm")
        # Randomize original observed treatment counts within each retained stratum.
        p=np.asarray([rng.permutation(g) for _ in range(nperm)])
        null+=w*(y[p[:,:len(a)]].mean(1)-y[p[:,len(a):]].mean(1))
        aa=rng.choice(a,(nboot,len(a)),replace=True);bb=rng.choice(b,(nboot,len(b)),replace=True)
        boot+=w*(y[aa].mean(1)-y[bb].mean(1))
    rows=[]
    for j,(year,name) in enumerate((y,k) for y in ("2020","2021") for k in ("named_TV","rank_TV","identity_excess_TV")):
        rows.append({"year":year,"endpoint":name,"difference_treated_minus_control":float(obs[j]),
                     "ci95_unadjusted_conditional_bootstrap":np.quantile(boot[:,j],[.025,.975]).tolist(),
                     "p":float((1+np.count_nonzero(np.abs(null[:,j])>=abs(obs[j])-1e-12))/(nperm+1)),
                     "n_logs":len(t),"n_permutations":nperm,"n_bootstrap":nboot})
    adjust(rows)
    return rows

def rarefy(raw,depth,rng):
    shape=raw.shape
    flat=raw.reshape(-1,shape[-1]).astype(np.int64)
    if np.any(flat.sum(1)<depth) or not np.all(flat==raw.reshape(flat.shape)):
        raise ValueError("rarefaction requires sufficient integer counts")
    out=np.zeros_like(flat)
    for i,row in enumerate(flat):
        idx=np.flatnonzero(row)
        out[i,idx]=rng.multivariate_hypergeometric(row[idx],depth)
    if not np.all(out.sum(1)==depth):raise ValueError("rarefaction total")
    return out.reshape(shape)

def r2_stream(des,y,sparse=False):
    y=np.asarray(y,float)
    if y.ndim!=2 or y.shape[0]!=des.n or not np.isfinite(y).all():
        raise ValueError("invalid response")
    # Gram centering avoids materializing dense centered OTU vectors for sparse matrices.
    if sparse:
        a=csr_matrix(y);g=(a@a.T).toarray()
    else:g=y@y.T
    for idx in des.groups:g[idx]-=g[idx].mean(0)
    for idx in des.groups:g[:,idx]-=g[:,idx].mean(1,keepdims=True)
    total=float(np.trace(g))
    if total<=1e-15:return np.zeros(des.nperm+1)
    # Binary comparison has one contrast column.
    w=des.W[:,:,0]
    z=np.sum((w@g)*w,axis=1)/total
    if not np.isfinite(z).all() or np.any(z<-1e-8) or np.any(z>1+1e-8):
        raise ValueError("invalid R2 stream")
    return z

def r2_arrays(raw,des):
    p=raw/raw.sum(2,keepdims=True)
    rank=np.sort(p,axis=2)[:,:,::-1]
    k=int(np.max(np.sum(raw>0,axis=2)))
    return np.asarray([r2_stream(des,v[:,j]-v[:,0],sparse=label=="named")
                       for j in (1,2) for label,v in (("named",p),("rank",rank[:,:,:k]))])

def summarize_r2(a,draws,n):
    rows=[]
    for j,(year,name) in enumerate((y,k) for y in ("2020","2021") for k in ("named","rank")):
        z=a[j];e=z[0]-np.median(z[1:])
        ds=[float(d[j,0]-np.median(d[j,1:])) for d in draws]
        rows.append({"year":year,"endpoint":name,"R2":float(z[0]),"E":float(e),
                     "p":float((1+np.count_nonzero(z[1:]>=z[0]-1e-12))/len(z)),
                     "n_logs":n,"draw_E_range":[min(ds),max(ds)],"n_draws":len(draws)})
    adjust(rows);return rows

def run(meta,x,nperm=NPERM,ntv=NTV,nboot=NBOOT,draws=DRAWS):
    ix,t,strata,sites,types,good,coverage,reads=design(meta,x)
    if int(good.sum())!=446:raise ValueError("original complete-log cohort drift")
    weights={s:float(np.sum(strata==s)/460) for s in np.unique(strata)}
    binary=t!="control"
    raw=x[ix[good]]
    result={"status":"EXPLORATORY_TV_GAP_AND_DEPTH_NOT_CAUSAL_MEDIATION",
            "source_assigned_logs":460,"original_complete_logs":446,
            "raw_TV":tv_score(tv_outcomes(raw),binary[good],strata[good],weights,SEED+701,ntv,nboot),
            "depth_sensitivity":[],"limits":["TV gap mathematically nonnegative per log; treatment contrast is the test.",
            "Conditional complete-case inference after post-treatment QC is not full ITT.",
            "TV displacement not the centroid R2 estimand; no equivalence or interaction coefficients.",
            "Depth standardization equalizes non-target sequencing effort, not biomass."]}
    totals=x.sum(1)[ix]
    result["original_non_target_read_min_quantiles"]=np.quantile(totals[good].min(1),[0,.01,.05,.25,.5,.75,1]).tolist()
    for i,depth in enumerate(DEPTHS):
        keep=good & np.all(totals>=depth,axis=1)
        rec={"depth":depth,"n_logs":int(keep.sum()),"n_treated":int((keep&binary).sum()),
             "n_controls":int((keep&~binary).sum()),"retention_vs_446":float(keep.sum()/446),
             "coverage_by_stratum_and_binary":[{"stratum":s,"treated":int((keep&binary&(strata==s)).sum()),
             "control":int((keep&~binary&(strata==s)).sum())} for s in sorted(weights)]}
        qc=[]
        for j,(lane,lab) in enumerate((("ten_group",t),("binary",binary.astype(int)))):
            d=BlockPermutation(lab,strata,nperm,SEED+801+i*10+j)
            qc.append({"lane":lane,**d.score(keep.astype(float))})
        adjust(qc);rec["eligibility_gate"]=qc
        if any(z["p_holm_family"]<=.05 for z in qc):
            rec["status"]="STOP_TREATMENT_ASSOCIATED_DEPTH_ELIGIBILITY";result["depth_sensitivity"].append(rec);continue
        if any(c["treated"]==0 or c["control"]==0 for c in rec["coverage_by_stratum_and_binary"]):
            rec["status"]="STOP_MISSING_ARM";result["depth_sensitivity"].append(rec);continue
        v=x[ix[keep]]
        d=BlockPermutation(binary[keep].astype(int),strata[keep],nperm,SEED+901+i)
        raw_a=r2_arrays(v,d);r=[];tv=[]
        rng=np.random.default_rng(SEED+1001+i)
        for _ in range(draws):
            q=rarefy(v,depth,rng)
            r.append(r2_arrays(q,d));tv.append(tv_outcomes(q))
        rec.update({"status":"SCORED_CONDITIONAL_DEPTH_COHORT",
                    "matched_original_read_R2":summarize_r2(raw_a,[raw_a],int(keep.sum())),
                    "rarefied_R2":summarize_r2(np.mean(r,axis=0),r,int(keep.sum())),
                    "matched_original_read_TV":tv_score(tv_outcomes(v),binary[keep],strata[keep],weights,SEED+1101+i,ntv,nboot),
                    "rarefied_TV":tv_score(np.mean(tv,axis=0),binary[keep],strata[keep],weights,SEED+1101+i,ntv,nboot)})
        result["depth_sensitivity"].append(rec)
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    with (a.source/"meta.csv").open() as f:meta=list(csv.DictReader(f))
    if len((a.source/"otu_names.txt").read_text().splitlines())!=4653 or len((a.source/"excluded_target_otus.txt").read_text().splitlines())!=9:
        raise ValueError("source roster mismatch")
    x=np.fromfile(a.source/"counts.bin",dtype="<f8").reshape(1380,4653)
    z=run(meta,x)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SAINE_BUFFER_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
