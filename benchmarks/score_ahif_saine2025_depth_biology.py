"""Depth-adjusted biological state prediction, nonzero-effect restricted wild reference."""
from __future__ import annotations
import argparse,csv,json
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_saine2025_succession import design,SEED,holm
from benchmarks.score_ahif_saine2025_buffering import tv_outcomes
from benchmarks.score_ahif_saine2025_baseline_context import (
    baseline_features,propensity,fit_ridge,ALPHAS
)
NBOOT=999

def std_on(x,mask):
    mean=x[mask].mean(0);sd=x[mask].std(0)
    return (x-mean)/np.where(sd>1e-12,sd,1.)

def treatment_features(f,types):
    c=np.column_stack((types=="R",types!="R")).astype(float)
    base=np.column_stack((c,c*f[:,3:4]))
    return base,np.column_stack((base,f[:,:3]))

def split_design(f,types,sites,weight,train_sites,held):
    tr=np.isin(sites,train_sites);te=sites==held
    if not tr.any() or not te.any() or np.any(tr&te):
        raise ValueError("split mismatch")
    z=std_on(f,tr)
    a,b=treatment_features(z,types)
    nuisance=np.column_stack((np.ones(len(f)),types=="R",z))
    # Precompute nuisance response map; no held-out outcomes in coefficients.
    g=nuisance[tr].T@(weight[tr,None]*nuisance[tr])+np.diag([0.,1.,1.,1.,1.,1.])
    m=np.linalg.solve(g,nuisance[tr].T*weight[tr])
    return tr,te,a,b,nuisance,m

def make_cache(f,types,sites,weight,w):
    names=sorted(set(sites));out=[]
    def prepare(train,held):
        tr,te,a,b,n,m=split_design(f,types,sites,weight,train,held)
        def fitmap(x,penalty):
            d=w[tr,None]*x[tr]
            g=d.T@(weight[tr,None]*d)+np.diag(penalty)
            return w[te,None]*x[te]@np.linalg.solve(g,d.T*weight[tr])
        p1=fitmap(a,[0.,0.,1.,1.])
        p2=fitmap(b,[0.,0.,1.,1.,1.,1.,1.])
        return tr,te,n,m,p1,p2-p1
    for held in names:
        train=[x for x in names if x!=held]
        out.append((held,prepare(train,held),[prepare([s for s in train if s!=v],v) for v in train]))
    return out

def split_loss(split,y,weight,train_year):
    tr,te,n,m,p1,pdelta=split
    # Fit prognostic nuisance on training data only; linear maps reused but
    # coefficients recomputed for every simulated outcome.
    rtr=y[tr]-n[tr]@(m@y[tr])
    rte=y[te]-n[te]@(m@y[tr])
    target=rte[:,2:4]
    past=rtr[:,0:2] if train_year==2020 else rtr[:,2:4]
    base=p1@past;delta=pdelta@past
    wt=weight[te,None];err=target-base
    return (np.sum(wt*err**2,axis=0),np.sum(wt*delta*err,axis=0),
            np.sum(wt*delta*delta,axis=0),np.sum(wt*target**2,axis=0))

def cv(cache,y,weight,train_year):
    folds=[]
    for site,outer,inner in cache:
        inside=[split_loss(s,y,weight,train_year) for s in inner]
        xy=sum(v[1] for v in inside);xx=sum(v[2] for v in inside)
        alpha=ALPHAS[np.argmin(ALPHAS[:,None]**2*xx-2*ALPHAS[:,None]*xy,axis=0)]
        b,xy,xx,zero=split_loss(outer,y,weight,train_year)
        folds.append({"site":site,"alpha":alpha,"zero":zero,"depth":b,
          "nested_biology":b-2*alpha*xy+alpha**2*xx,"full_biology":b-2*xy+xx})
    sums={k:sum(f[k] for f in folds) for k in ("zero","depth","nested_biology","full_biology")}
    def gain(a,b):
        return np.divide(a-b,a,out=np.zeros_like(a),where=a>1e-15)
    return {"gain":gain(sums["depth"],sums["nested_biology"]),
            "depth_gain":gain(sums["zero"],sums["depth"]),
            "full_gain":gain(sums["depth"],sums["full_biology"]),
            "folds":folds,"losses":sums}

class LinearDesign:
    """Fixed weighted-design OLS and HC3 coefficient covariance."""
    def __init__(self,x,weight):
        self.x=np.asarray(x,float);self.sw=np.sqrt(weight)
        self.a=self.x*self.sw[:,None]
        if np.linalg.matrix_rank(self.a)!=self.a.shape[1]:
            raise ValueError("collinear linear null/full design")
        self.inv=np.linalg.inv(self.a.T@self.a)
        self.map=self.inv@self.a.T
        self.h=np.sum(self.a*self.map.T,axis=1)
        if np.any(self.h>=1-1e-9):raise ValueError("unit leverage")
    def fit(self,y):
        b=self.map@(self.sw[:,None]*y)
        pred=self.x@b
        return b,pred,y-pred
    def wald3(self,y):
        b,_,res=self.fit(y)
        robust=(self.sw[:,None]*res)/(1-self.h[:,None])
        influence=self.map[-3:,:]
        values=[]
        for i in range(y.shape[1]):
            v=influence*robust[:,i]
            cov=v@v.T
            if np.linalg.matrix_rank(cov)!=3:raise ValueError("singular HC3 covariance")
            values.append(float(b[-3:,i]@np.linalg.solve(cov,b[-3:,i])))
        return np.array(values),b[-3:,:]

def null_design(features,types,strata,w,weight):
    f=std_on(features,np.ones(len(features),bool))
    a,b=treatment_features(f,types)
    g=(strata[:,None]==np.unique(strata)[None,:]).astype(float)
    prognosis=np.column_stack((g,f,f*f))
    z0=np.column_stack((prognosis,w[:,None]*a))
    z1=np.column_stack((prognosis,w[:,None]*b))
    return LinearDesign(z0,weight),LinearDesign(z1,weight)

def evaluate(features,types,sites,strata,binary,y,weight,nboot=NBOOT):
    w=binary-propensity(binary,strata)
    cache=make_cache(features,types,sites,weight,w)
    obs={str(k):cv(cache,y,weight,k) for k in (2020,2021)}
    restricted,full=null_design(features,types,strata,w,weight)
    _,fit,res=restricted.fit(y)
    resid=res/np.sqrt(1-restricted.h[:,None])
    wald,coef=full.wald3(y)
    rng=np.random.default_rng(SEED+1701)
    null_cv={k:[] for k in obs};null_w=[]
    for _ in range(nboot):
        signs=rng.choice([-1.,1.],size=len(y))
        z=fit+signs[:,None]*resid
        for k in obs:
            null_cv[k].append(cv(cache,z,weight,int(k))["gain"])
        null_w.append(full.wald3(z)[0])
    rows=[]
    for k,v in obs.items():
        null=np.array(null_cv[k])
        for j,name in enumerate(("named_TV","rank_TV")):
            gain=float(v["gain"][j])
            rows.append({"train_year":int(k),"test_year":2021,"endpoint":name,
              "biology_gain_over_depth":gain,"depth_gain_over_zero":float(v["depth_gain"][j]),
              "full_biology_gain_over_depth":float(v["full_gain"][j]),
              "p_approx":float((1+np.count_nonzero(null[:,j]>=gain-1e-12))/(nboot+1)),
              "null_gain_median":float(np.median(null[:,j])),
              "null_gain_95":np.quantile(null[:,j],[.025,.975]).tolist(),
              "losses":{a:float(b[j]) for a,b in v["losses"].items()},
              "folds":[{"site":r["site"],**{a:float(b[j]) for a,b in r.items() if a!="site"}} for r in v["folds"]]})
    for r,p in zip(rows,holm([r["p_approx"] for r in rows])):r["p_holm4_approx"]=p
    null_w=np.asarray(null_w);direct=[]
    for i,(year,name) in enumerate((yr,ep) for yr in (2020,2021) for ep in ("named_TV","rank_TV")):
        direct.append({"year":year,"endpoint":name,"wald_HC3":float(wald[i]),
          "p_approx":float((1+np.count_nonzero(null_w[:,i]>=wald[i]-1e-12))/(nboot+1)),
          "bio_standardized_slopes":coef[:,i].tolist()})
    for r,p in zip(direct,holm([r["p_approx"] for r in direct])):r["p_holm4_approx"]=p
    return {"nested_prediction":rows,"direct_linear_diagnostic":direct,
            "n_wild_bootstrap":nboot,"restricted_design_columns":restricted.x.shape[1],
            "full_design_columns":full.x.shape[1],
            "bio_slope_order":["log1p_richness","simpson","top1"],
            "max_restricted_leverage":float(restricted.h.max()),
            "max_full_leverage":float(full.h.max())}

def run(meta,x,nboot=NBOOT):
    ix,t,strata,sites,types,good,groups,reads=design(meta,x)
    if int(good.sum())!=446:raise ValueError("complete cohort changed")
    raw=x[ix[good]]
    f=baseline_features(raw);tv=tv_outcomes(raw);y=tv[:,[0,1,3,4]]
    binary=(t[good]!="control").astype(float);ss=strata[good]
    weight=np.array([np.sum(strata==s)/np.sum(ss==s) for s in ss],float)
    result=evaluate(f,types[good],sites[good],ss,binary,y,weight,nboot)
    return {"status":"EXPLORATORY_RESTRICTED_WILD_NOT_EXACT_CAUSAL_MODERATION",
            "assigned_logs":460,"complete_logs":446,"treated_logs":int(binary.sum()),
            "n_sites":5,**result,
            "limits":["Same experiment, post-treatment complete-case QC.",
             "Restricted conditional linear model, independent log errors; not exact randomization.",
             "Initial read totals proxy detection effort; richness/dominance still detection-sensitive.",
             "No directly measured competition coefficients or independent study replication."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True);p.add_argument("--bootstrap",type=int,default=NBOOT)
    a=p.parse_args()
    with (a.source/"meta.csv").open() as f:meta=list(csv.DictReader(f))
    names=(a.source/"otu_names.txt").read_text().splitlines()
    if len(names)!=4653 or len(set(names))!=4653 or len((a.source/"excluded_target_otus.txt").read_text().splitlines())!=9:
        raise ValueError("source roster mismatch")
    x=np.fromfile(a.source/"counts.bin",dtype="<f8").reshape(1380,4653)
    z=run(meta,x,a.bootstrap)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SAINE_DEPTH_BIOLOGY_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
