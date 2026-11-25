"""Baseline-conditioned cross-site introduction response prediction using R-loss."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_saine2025_succession import design,SEED,holm
from benchmarks.score_ahif_saine2025_buffering import tv_outcomes

ALPHAS=np.array([0.,.25,.5,.75,1.])
NPERM=999

def fit_ridge(a,y,penalty,weight):
    a=np.asarray(a,float);y=np.asarray(y,float);weight=np.asarray(weight,float)
    g=a.T@(weight[:,None]*a)+np.diag(penalty)
    rhs=a.T@(weight[:,None]*y)
    try:return np.linalg.solve(g,rhs)
    except np.linalg.LinAlgError:raise ValueError("Unidentified fit")

def baseline_features(raw):
    z=np.asarray(raw[:,0],float)
    n=z.sum(1)
    if np.any(n<=0):raise ValueError("zero baseline")
    p=z/n[:,None]
    return np.column_stack((np.log1p((z>0).sum(1)),np.sum(p*p,axis=1),p.max(1),np.log1p(n)))

def propensity(t,strata):
    w=np.zeros(len(t))
    for g in np.unique(strata):
        ix=strata==g;e=float(np.mean(t[ix]))
        if not 0<e<1:raise ValueError("missing treatment arm")
        w[ix]=e
    return w

def make_split(features,types,sites,y,weights,train_sites,held_site):
    tr=np.isin(sites,train_sites);te=sites==held_site
    if not tr.any() or not te.any() or np.any(tr&te):raise ValueError("bad site split")
    mu=features[tr].mean(0);sd=features[tr].std(0);sd=np.where(sd>1e-12,sd,1.)
    f=(features-mu)/sd
    c=np.column_stack((types=="R",types!="R")).astype(float)
    nuisance=np.column_stack((np.ones(len(f)),types=="R",f))
    ww=weights[tr]
    # all four response columns: 2020 named/rank, 2021 named/rank.
    coef=fit_ridge(nuisance[tr],y[tr],[0.,1.,1.,1.,1.,1.],ww)
    resid=y-nuisance@coef
    return {"tr":tr,"te":te,"c":c,"f":np.column_stack((c,f)),"resid":resid,
            "weights":weights}

def prediction_split(split,w,train_year):
    tr=split["tr"];te=split["te"];a=split["c"];b=split["f"]
    weight=split["weights"];r=split["resid"]
    yy=r[:,0:2] if train_year==2020 else r[:,2:4]
    m1=fit_ridge(w[tr,None]*a[tr],yy[tr],[0.,0.],weight[tr])
    m2=fit_ridge(w[tr,None]*b[tr],yy[tr],[0.,0.,1.,1.,1.,1.],weight[tr])
    base=w[te,None]*(a[te]@m1)
    delta=w[te,None]*(b[te]@m2)-base
    target=r[te,2:4];wt=weight[te,None]
    return {"zero":np.sum(wt*target*target,axis=0),
            "base":np.sum(wt*(target-base)**2,axis=0),
            "xy":np.sum(wt*delta*(target-base),axis=0),
            "xx":np.sum(wt*delta*delta,axis=0)}

def prepare_cache(features,types,sites,y,weight):
    names=sorted(set(sites));cache=[]
    if len(names)<4:raise ValueError("nested site design")
    for held in names:
        tr=[v for v in names if v!=held]
        outer=make_split(features,types,sites,y,weight,tr,held)
        inner=[make_split(features,types,sites,y,weight,[s for s in tr if s!=v],v) for v in tr]
        cache.append((held,outer,inner))
    return cache

def cv(cache,w,train_year):
    rows=[]
    for held,out,inner in cache:
        ins=[prediction_split(z,w,train_year) for z in inner]
        xy=sum(z["xy"] for z in ins);xx=sum(z["xx"] for z in ins)
        alpha=ALPHAS[np.argmin(ALPHAS[:,None]**2*xx[None,:]-2*ALPHAS[:,None]*xy[None,:],axis=0)]
        o=prediction_split(out,w,train_year)
        nested=o["base"]-2*alpha*o["xy"]+alpha**2*o["xx"]
        full=o["base"]-2*o["xy"]+o["xx"]
        rows.append({"site":held,"alpha":alpha,"zero":o["zero"],"pooled":o["base"],
                     "nested":nested,"full":full})
    total={k:sum(z[k] for z in rows) for k in ("zero","pooled","nested","full")}
    def gain(base,pred):return np.divide(base-pred,base,out=np.zeros_like(base),where=base>1e-15)
    return {"context_gain":gain(total["pooled"],total["nested"]),
            "pooled_gain":gain(total["zero"],total["pooled"]),
            "full_context_gain":gain(total["pooled"],total["full"]),
            "losses":total,"folds":rows}

def run(meta,x,nperm=NPERM):
    ix,t,strata,sites,types,good,coverage,reads=design(meta,x)
    if int(good.sum())!=446:raise ValueError("original complete cohort changed")
    raw=x[ix[good]]
    feat=baseline_features(raw)
    tv=tv_outcomes(raw);y=tv[:,[0,1,3,4]]
    binary=(t[good]!="control").astype(float)
    ss=strata[good];sites=sites[good];types=types[good]
    # rescale original stratum weights by average retained log count; penalty
    # scale is approximately one observation; relative loss invariant to scale.
    weight=np.array([np.sum(strata==s)/np.sum(ss==s) for s in ss],float)
    e=propensity(binary,ss);cache=prepare_cache(feat,types,sites,y,weight)
    obs={str(year):cv(cache,binary-e,year) for year in (2020,2021)}
    rng=np.random.default_rng(SEED+1401)
    blocks=[np.flatnonzero(ss==s) for s in np.unique(ss)]
    null={k:[] for k in obs}
    for _ in range(nperm):
        p=binary.copy()
        for g in blocks:p[g]=rng.permutation(p[g])
        for k in obs:null[k].append(cv(cache,p-e,int(k))["context_gain"])
    out={}
    for k,z in obs.items():
        nulls=np.array(null[k]);p=(1+np.count_nonzero(nulls>=z["context_gain"]-1e-12,axis=0))/(nperm+1)
        adj=holm(p.tolist());rows=[]
        for j,name in enumerate(("named_TV","rank_TV")):
            rows.append({"endpoint":name,
                "context_gain_over_pooled":float(z["context_gain"][j]),
                "pooled_gain_over_zero":float(z["pooled_gain"][j]),
                "full_context_gain_over_pooled":float(z["full_context_gain"][j]),
                "total_loss":{a:float(b[j]) for a,b in z["losses"].items()},
                "folds":[{"site":r["site"],**{a:float(b[j]) for a,b in r.items() if a!="site"}} for r in z["folds"]],
                "sharp_null_p_upper":float(p[j]),"sharp_null_p_holm":adj[j],
                "null_gain_median":float(np.median(nulls[:,j])),
                "null_gain_95":np.quantile(nulls[:,j],[.025,.975]).tolist()})
        out[k]=rows
    return {"status":"EXPLORATORY_BASELINE_STATE_R_LOSS_NOT_CAUSAL_MODERATION",
        "assigned_logs":460,"complete_logs":446,"treated_logs":int(binary.sum()),
        "feature_names":["log1p_richness","simpson","top1","log1p_non_target_reads"],
        "test_year":2021,"alpha_grid":ALPHAS.tolist(),"n_permutations":nperm,
        "by_training_year":out,
        "limits":["Same source experiment; five fixed sites, no independent new study.",
            "R-loss evaluates predicted treatment response, not individual causal effects.",
            "Permutation reference is sharp no treatment effect, not null of no effect modification.",
            "Post-treatment QC conditional, read counts not biomass, no direct interaction coefficients."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True);p.add_argument("--permutations",type=int,default=NPERM)
    a=p.parse_args()
    with (a.source/"meta.csv").open() as f:meta=list(csv.DictReader(f))
    names=(a.source/"otu_names.txt").read_text().splitlines()
    if len(names)!=4653 or len(set(names))!=4653 or len((a.source/"excluded_target_otus.txt").read_text().splitlines())!=9:
        raise ValueError("source roster mismatch")
    x=np.fromfile(a.source/"counts.bin",dtype="<f8").reshape(1380,4653)
    z=run(meta,x,a.permutations)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SAINE_BASELINE_CONTEXT_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
