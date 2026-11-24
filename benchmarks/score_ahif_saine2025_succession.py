"""Source locked Saine 460-log conditional succession reanalysis."""
from __future__ import annotations
import argparse,csv,json
from collections import Counter
from pathlib import Path
import numpy as np

SEED=20261124
N_PERM=1999

def holm(ps):
    result=[0.]*len(ps);last=0.
    for j,i in enumerate(sorted(range(len(ps)),key=lambda i:ps[i])):
        last=max(last,min(1.,ps[i]*(len(ps)-j)));result[i]=last
    return result

def adjust(rows):
    for r,p in zip(rows,holm([x["p"] for x in rows])):r["p_holm_family"]=p

def design(meta,x):
    if len(meta)!=1380 or x.shape!=(1380,4646) or not np.isfinite(x).all() or (x<0).any():
        raise ValueError("Unexpected original non-target matrix")
    if len(set(m["sample"] for m in meta))!=1380:
        raise ValueError("Duplicate original sample")
    ids=sorted(set(m["log"] for m in meta));years=["2019","2020","2021"]
    if len(ids)!=460:raise ValueError("Expected 460 unique logs")
    key={(m["log"],m["year"]):i for i,m in enumerate(meta)}
    if len(key)!=1380 or set(key)!={(i,y) for i in ids for y in years}:
        raise ValueError("Missing or duplicate source log-year")
    ix=np.array([[key[i,y] for y in years] for i in ids])
    baseline=[meta[j] for j in ix[:,0]]
    for i,row in enumerate(ix):
        for j in row:
            if any(meta[j][k]!=baseline[i][k] for k in ("site","logtype","treatment")):
                raise ValueError("Source log assignment drift")
    treatments=np.array([m["treatment"] for m in baseline])
    sites=np.array([m["site"] for m in baseline]);types=np.array([m["logtype"] for m in baseline])
    if len(np.unique(sites))!=5 or set(types)!={"B","R","U"}:raise ValueError("Wrong site/type")
    strata=np.array([s+"|"+("felled" if t=="R" else "natural") for s,t in zip(sites,types)])
    labels=sorted(set(treatments)-{"control"})
    if len(labels)!=9:raise ValueError("Nine target treatment labels")
    for s in np.unique(strata):
        counts=Counter(treatments[strata==s])
        expected={l:3 if s.endswith("|felled") else 5 for l in labels}
        expected["control"]=10
        if counts!=expected:raise ValueError("Allocation mismatch "+s+" "+str(counts))
    reads=np.asarray([float(m["reads"]) for m in meta])
    if not np.isfinite(reads).all() or np.count_nonzero(reads>=10000)!=1369:
        raise ValueError("Read QC mismatch")
    valid=(reads>=10000)&(x.sum(axis=1)>0)
    good=np.all(valid[ix],axis=1)
    groups=[]
    for s in sorted(set(strata)):
        for t in sorted(set(treatments)):
            m=(strata==s)&(treatments==t)
            groups.append({"stratum":s,"treatment":t,"assigned":int(m.sum()),
                           "complete":int((m&good).sum())})
    return ix,treatments,strata,sites,types,good,groups,reads

class BlockPermutation:
    """Partial Euclidean R² for treatment after fixed-stratum adjustment."""
    def __init__(self,labels,strata,nperm,seed):
        self.labels=np.asarray(labels);self.strata=np.asarray(strata)
        self.groups=[np.flatnonzero(self.strata==s) for s in np.unique(self.strata)]
        self.levels=np.unique(self.labels);self.n=len(labels)
        if len(self.levels)<2:raise ValueError("No contrast")
        for g in self.groups:
            if len(set(self.labels[g]))!=len(self.levels):
                raise ValueError("Treatment lost in stratum after QC; no inferential scoring")
        rng=np.random.default_rng(seed); p=np.tile(self.labels,(nperm+1,1))
        for i in range(1,len(p)):
            for g in self.groups:p[i,g]=rng.permutation(p[0,g])
        # One-hot treatment design residualized within original stratum.
        z=(p[:,:,None]==self.levels[None,None,:]).astype(float)
        for g in self.groups:z[:,g,:]-=z[:,g,:].mean(axis=1,keepdims=True)
        gram=z[0].T@z[0];w,v=np.linalg.eigh(gram)
        keep=w>1e-9
        if keep.sum()!=len(self.levels)-1:raise ValueError("Treatment contrasts not identified")
        self.W=z@(v[:,keep]/np.sqrt(w[keep]))
        self.nperm=nperm
    def score(self,x):
        x=np.asarray(x,float)
        if x.ndim==1:x=x[:,None]
        if x.shape[0]!=self.n or not np.isfinite(x).all():raise ValueError("Bad response")
        y=x.copy()
        for g in self.groups:y[g]-=y[g].mean(axis=0)
        total=float(np.sum(y*y))
        if total<=1e-15:
            return {"R2":0.,"E":0.,"null_median":0.,"p":1.,"n":self.n,"features":x.shape[1],
                    "status":"CONSTANT_RESPONSE","permutations":self.nperm}
        # Use a compact response or its Gram matrix; never huge permutation x OTU array.
        k=self.W.shape[2];flat=self.W.transpose(0,2,1).reshape(-1,self.n)
        if y.shape[1]>self.n:
            g=y@y.T
            ss=np.sum((flat@g)*flat,axis=1)
        else:ss=np.sum((flat@y)**2,axis=1)
        vals=ss.reshape(self.nperm+1,k).sum(axis=1)/total
        if not np.isfinite(vals).all() or np.any(vals<-1e-8) or np.any(vals>1+1e-8):
            raise ValueError("Invalid partial R2")
        obs=float(vals[0]);nul=vals[1:];med=float(np.median(nul))
        return {"R2":obs,"E":obs-med,"null_median":med,"p":float((1+np.count_nonzero(nul>=obs-1e-12))/(self.nperm+1)),
                "n":self.n,"features":x.shape[1],"permutations":self.nperm,"status":"SCORED"}

def run(meta,x,nperm=N_PERM):
    ix,t,strata,sites,types,good,groups,reads=design(meta,x)
    z={"status":"SOURCE_LOCKED_CONDITIONAL_COMPLETE_LOG_ANALYSIS_NOT_MEDIATION",
       "source_physical_logs":460,"source_sample_states":1380,"source_nontarget_otus":4646,
       "complete_triplets":int(good.sum()),"excluded_logs":int((~good).sum()),
       "low_read_states":int((reads<10000).sum()),"zero_nontarget_states":int((x.sum(1)==0).sum()),
       "coverage_by_stratum_treatment":groups,"n_source_publication_families":1,
       "limits":["Conditioning on post-treatment QC/completeness may bias randomization; missingness test does not prove exchangeability.",
                 "Read composition not absolute biomass or direct interaction coefficients.",
                 "Randomized inoculum identity, not randomized arrival order or environment."]}
    qc=[]
    for i,(lane,labels) in enumerate([("ten_group",t),("any_introduction",(t!="control").astype(int))]):
        d=BlockPermutation(labels,strata,nperm,SEED+10+i)
        qc.append({"lane":lane,"endpoint":"complete_triplet",**d.score(good.astype(float))})
    adjust(qc);z["eligibility_gate"]=qc
    if any(r["p_holm_family"]<=.05 for r in qc):
        z["status"]="STOP_TREATMENT_ASSOCIATED_ELIGIBILITY";return z
    if not good.any():z["status"]="STOP_NO_COMPLETE_TRIPLETS";return z
    raw=x[ix[good]]
    prop=raw/raw.sum(axis=2,keepdims=True)
    rank=np.sort(prop,axis=2)[:,:,::-1]
    # Trailing all-zero columns in sorted profiles carry no information; keep logical width.
    last=int(np.max(np.sum(raw>0,axis=2)))
    compact_rank=rank[:,:,:last]
    t=t[good];s=strata[good]
    groupspec=[("ten_group",np.ones(len(t),bool),t),
               ("any_introduction",np.ones(len(t),bool),(t!="control").astype(int)),
               ("identity_among_introduced",t!="control",t[t!="control"])]
    results={};secondary=[];placebo=[]
    for j,(lane,mask,label) in enumerate(groupspec):
        try:d=BlockPermutation(label,s[mask],nperm,SEED+101+j)
        except ValueError as e:results[lane]={"status":"STOP_DESIGN","reason":str(e)};continue
        rows=[]
        for year,k in [("2020",1),("2021",2)]:
            for name,v in [("named",prop),("rank",compact_rank)]:
                row={"year":year,"endpoint":name,**d.score((v[:,k]-v[:,0])[mask])}
                if name=="rank":row["logical_features"]=4646
                rows.append(row)
        adjust(rows);results[lane]=rows
        if lane=="ten_group":
            for year,k in [("2020",1),("2021",2)]:
                for name,v in [("richness",np.sum(raw>0,axis=2)),("simpson",np.sum(prop*prop,axis=2)),("top1",compact_rank[:,:,0])]:
                    secondary.append({"year":year,"endpoint":name,**d.score(v[:,k]-v[:,0])})
            for name,v in [("named",prop),("rank",compact_rank)]:
                placebo.append({"year":"2019","endpoint":name,**d.score(v[:,0])})
    adjust(secondary);adjust(placebo)
    z.update({"primary_and_distinct_assignment_lanes":results,"ten_group_scalar_secondary":secondary,
              "baseline_placebo":placebo,"n_rank_nonzero_width":last})
    return z

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=N_PERM);a=p.parse_args()
    with (a.source/"meta.csv").open() as f:meta=list(csv.DictReader(f))
    names=(a.source/"otu_names.txt").read_text().splitlines()
    excluded=(a.source/"excluded_target_otus.txt").read_text().splitlines()
    if len(names)!=4646 or len(set(names))!=4646 or len(excluded)!=16:
        raise ValueError("Source target exclusion/roster failure")
    x=np.fromfile(a.source/"counts.bin",dtype="<f8").reshape(1380,4646)
    z=run(meta,x,a.permutations);z["excluded_target_otus"]=excluded
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SAINE_SUCCESSION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
