"""Predator × competitor finite contextual responses, no monoculture extrapolation."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
import numpy as np
ARMS=("C+D","C+D+P","C+D+S","C+D+Pd","C+D+P+Pd","C+D+S+Pd")
FOCAL=("Colp","Dexio")
DAYS=tuple(range(1,54,2))
SEED=2026112717;N=256
def load(path):
    rows=list(csv.DictReader(path.open()))
    allobs=set();v={};counts=Counter();roster=set()
    for r in rows:
        tr=r["treatment"];rep=int(float(r["replicate"]));day=int(float(r["NumDays"]));sp=r["predict_spec"]
        assert tr in ARMS and rep in (1,2,3,4) and day in DAYS,(tr,rep,day)
        if day==1:roster.add((tr,rep,sp))
        val=r["mean.dens.ml"]
        if val in ("","NA","NaN"):counts["explicit_na"]+=1;continue
        x=float(val);assert np.isfinite(x) and x>=0
        allobs.add((tr,rep,day))
        if sp not in FOCAL:continue
        key=(tr,rep,sp,day);assert key not in v,("duplicate focal row",key)
        v[key]=x
    assert all((a,r,s) in roster and (a,r,s,1) in v for a in ARMS for r in range(1,5) for s in FOCAL)
    common=[d for d in DAYS if d!=1 and all((a,r,d) in allobs for a in ARMS for r in range(1,5))]
    strict=[d for d in common if all((a,r,s,d) in v for a in ARMS for r in range(1,5) for s in FOCAL)]
    assert len(common)>=3
    def matrix(ds):
        return np.array([[[[v.get((a,r,s,d),0.) for d in ds] for r in range(1,5)] for a in ARMS] for s in FOCAL])
    audit={"source_rows":len(rows),"focal_observed_rows":len(v),"explicit_na_rows":counts["explicit_na"],"culture_n":24,"common_dates":common,"strict_dates":strict,
           "missing_entire_culture_days":[[a,r,d] for a in ARMS for r in range(1,5) for d in DAYS if (a,r,d) not in allobs],
           "zero_fills_common":sum((a,r,s,d) not in v for a in ARMS for r in range(1,5) for s in FOCAL for d in common),
           "observed_zero_rows_common":sum(v.get((a,r,s,d),-1)==0 for a in ARMS for r in range(1,5) for s in FOCAL for d in common),
           "day1_means":{s:dict(zip(ARMS,matrix([1])[i].mean((1,2)).tolist())) for i,s in enumerate(FOCAL)}}
    return matrix,common,strict,audit
def sign(x,e=0):return np.where(x>e,1,np.where(x< -e,-1,0))
def relation(x,y,e=0):
    a,b=int(sign(x,e)),int(sign(y,e))
    return "unresolved" if not a or not b else "same" if a==b else "opposite"
def rows_for(x,ds,ref):
    means=x.mean(2);ret=[]
    for i,sp in enumerate(FOCAL):
        scale=ref[i];assert scale>0
        for comp,b,d in (("P",1,4),("S",2,5)):
            a,B,c,D=means[i,[0,b,3,d]]/scale
            for typ,u,v in (("competitor",B-a,D-c),("predator",c-a,D-B)):
                m=float(u.mean());n=float(v.mean())
                rr={"focal":sp,"competitor":comp,"added_agent":typ,"base_marginal":m,"background_marginal":n,"difference":n-m,"direction":relation(m,n),"direction_threshold_005":relation(m,n,.05),
                    "daily_direction_counts":dict(Counter(relation(k,l) for k,l in zip(u,v))),
                    "n_dates":len(ds),"baseline_CD_density":float(scale),
                    "mean_arm_densities":dict(zip(("CD","CDcomp","CDpred","CDcompPred"),(means[i,[0,b,3,d]].mean(1)).tolist()))}
                ret.append(rr)
    return ret
def forecast(x):
    # Each arm independently permuted; physical units do not overlap train/held.
    rng=np.random.default_rng(SEED)
    orders=np.stack([rng.permutation(4) for _ in range(N*6)]).reshape(N,6,4)
    vals=[];ok=[]
    for order in orders:
        for orientation in (0,1):
            tr=order[:,:2] if orientation==0 else order[:,2:]
            te=order[:,2:] if orientation==0 else order[:,:2]
            assert all(not set(a)&set(b) for a,b in zip(tr,te))
            tm=np.stack([x[:,a,ids].mean((1,2)) for a,ids in enumerate(tr)],axis=1)
            hm=np.stack([x[:,a,ids].mean((1,2)) for a,ids in enumerate(te)],axis=1)
            if (tm[:,0]<=0).any():continue
            vv=[]
            for i in range(2):
                a=tm[i,0]
                for b,d in ((1,4),(2,5)):
                    for c,background in ((b,3),(3,b)):
                        train_t=tm[i,c];h=hm[i,background];target=hm[i,d]
                        pp=[h,max(0.,h+train_t-a),h*train_t/a]
                        vv.append([abs(p-target)/a for p in pp])
            vals.append(vv)
    assert vals
    a=np.array(vals);m=a.mean(0)
    return {"n_common_folds":len(a),"n_total_folds":N*2,"dependent_panel_folds":len(a)*8,"mean_MAE":dict(zip(("no_change","additive","multiplicative"),m.mean(0).tolist())),"by_panel":[{"focal":s,"competitor":c,"added_agent":j,"MAE":dict(zip(("no_change","additive","multiplicative"),v.tolist()))} for (s,c,j),v in zip([(s,c,j) for s in FOCAL for c in ("P","S") for j in ("competitor","predator")],m)]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--csv",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    mat,common,strict,audit=load(a.csv)
    print("SHEN_DENSITY_AUDIT",json.dumps(audit),flush=True)
    x=mat(common);ref=x[:,0].mean((1,2));res={}
    for lane,ds in (("primary_common",common),("early",[d for d in common if d<=27]),("late",[d for d in common if d>=29]),("strict_complete_focals",strict),("after_second_predator_addition",[d for d in common if d>=9])):
        if len(ds)<3:res[lane]={"status":"STOP_INSUFFICIENT_DATES","dates":ds};continue
        xx=mat(ds);r=rows_for(xx,ds,ref)
        res[lane]={"status":"SCORED_DESCRIPTIVE","dates":ds,"marginals":r,"counts":dict(Counter(v["direction"] for v in r)),"counts_threshold_005":dict(Counter(v["direction_threshold_005"] for v in r))}
        if lane in ("primary_common","strict_complete_focals","after_second_predator_addition"):res[lane]["forecast"]=forecast(xx)
    z={"status":"FACTORIAL_CONTEXT_NOT_PAIR_ALONE_NOT_META","source_commit":"6df03dc4f0a8eb4efca99d729e23e917a06b1812","publication":"10.1002/ece3.10502","audit":audit,"results":res}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SHEN_CONTEXT_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
