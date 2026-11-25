"""Held-plot static history vs lag prediction, fixed pre-2019 calibration."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from benchmarks.score_ahif_portal_exclosure import prepare,score_lane,C,E,SEED,N,FILES
EPS=1e-12
SIZES=(2,4)
KEYS=("lag_loss","fixed_loss","cumulative_loss","fixed_minus_lag","cumulative_minus_lag",
      "fixed_sign_only_minus_lag","cumulative_sign_only_minus_lag","concurrent_sign_only_minus_lag",
      "lag_G","fixed_G","cumulative_G")
def norm(x):
    s=x.sum(-1,keepdims=True)
    return np.divide(x,s,out=np.zeros_like(x),where=s>EPS),s[...,0]>EPS
def signalloc(x):
    return np.concatenate((np.maximum(x,0),np.maximum(-x,0)),axis=-1)
def withsign(m,x):
    return np.concatenate((m*np.where(x>0,1.,np.where(x<0,0.,.5)),
                           m*np.where(x<0,1.,np.where(x>0,0.,.5))),axis=-1)
def tv(p,q):return np.abs(p-q).sum(-1)/2
def hedge(p):
    a,b=np.split(p,2,axis=-1);return np.concatenate(((a+b)/2,(a+b)/2),axis=-1)
def calc(d,h,cal,lag,target):
    L=d[:,:,lag];D=d[:,:,target];Y=h[:,:,target]
    F=np.broadcast_to(d[:,:,cal].mean(2,keepdims=True),L.shape)
    P=np.stack([d[:,:,:t].mean(2) for t in lag],axis=2)
    l,okl=norm(signalloc(L));f,okf=norm(signalloc(F));p,okp=norm(signalloc(P));q,okq=norm(signalloc(Y))
    mass,okm=norm(np.abs(L));ok=okl&okf&okp&okq&okm
    ll=tv(l,q);fl=tv(f,q);pl=tv(p,q)
    v=np.stack((ll,fl,pl,fl-ll,pl-ll,tv(withsign(mass,F),q)-ll,tv(withsign(mass,P),q)-ll,
      tv(withsign(mass,D),q)-ll,tv(hedge(l),q)-ll,tv(hedge(f),q)-fl,tv(hedge(p),q)-pl),axis=-1)
    return v,ok
def get_contrasts(x,oc,oe,n):
    tr=[];te=[]
    for o in (0,1):
        ic=oc[:,:n] if o==0 else oc[:,-n:];jc=oc[:,-n:] if o==0 else oc[:,:n]
        ie=oe[:,:n] if o==0 else oe[:,-n:];je=oe[:,-n:] if o==0 else oe[:,:n]
        assert not np.any(ic[:,:,None]==jc[:,None,:]) and not np.any(ie[:,:,None]==je[:,None,:])
        tr.append((x[:,len(C)+ie,:].mean(2)-x[:,ic,:].mean(2)).transpose(1,0,2))
        te.append((x[:,len(C)+je,:].mean(2)-x[:,jc,:].mean(2)).transpose(1,0,2))
    return np.stack(tr,axis=1),np.stack(te,axis=1)
def summarize(v,mask):
    valid=[i for i in range(v.shape[2]) if mask[:,:,i].any()]
    if not valid:return None
    by=np.array([v[:,:,i][mask[:,:,i]].mean(0) for i in valid])
    return {"n_defined":int(mask.sum()),"n_dates":len(valid),"mean":dict(zip(KEYS,by.mean(0).tolist())),
      "per_date":[{"date_index":i,**dict(zip(KEYS,r.tolist()))} for i,r in zip(valid,by)]}
def selftest():
    a=np.array([1.,-1.]);p,_=norm(signalloc(a));q,_=norm(signalloc(a))
    m,_=norm(np.abs(a))
    assert abs(tv(p,q))<EPS and abs(tv(withsign(m,-a),q)-1)<EPS
    assert abs(tv(withsign(m,np.zeros(2)),q)-.5)<EPS
def main():
    a=argparse.ArgumentParser();a.add_argument("--root",type=Path,required=True);a.add_argument("--out",type=Path,required=True);a=a.parse_args()
    selftest();x,x2,spp,rem,pairs,audit=prepare(a.root);assert audit["n_eligible_periods"]==70 and len(pairs)==56
    rng=np.random.default_rng(SEED);oc=np.array([rng.permutation(len(C)) for _ in range(N)]);oe=np.array([rng.permutation(len(E)) for _ in range(N)])
    x=x[:,:,[i for i,s in enumerate(spp) if s not in rem]]
    # Replay original common-mask score before follow-up.
    z=score_lane(x,pairs,oc,oe)["2"]["common"]
    old=json.loads(Path("validation/ahif_v1/portal_exclosure_summary_20261125.json").read_text())["primary_nontarget"]["2"]["common"]
    assert all(abs(z[k]-old[k])<1e-10 for k in ("G","H","O","Astar"))
    moons=pd.read_csv(a.root/FILES["moons"][0]);moons.columns=moons.columns.str.lower()
    years=dict(zip(moons.newmoonnumber.astype(int),pd.to_datetime(moons.newmoondate).dt.year))
    em=audit["eligible_moons"];yr=np.array([years[t] for t in em])
    cal=np.flatnonzero((yr>=2016)&(yr<=2018))
    pp=[(i,j) for i,j in pairs if yr[i]>=2019 and yr[j]>=2019]
    assert len(cal)>=6 and pp and max(cal)<min(i for i,j in pp)
    lag=[i for i,j in pp];target=[j for i,j in pp]
    vv={};oo={}
    for n in SIZES:
        d,h=get_contrasts(x,oc,oe,n);vv[n],oo[n]=calc(d,h,cal,lag,target)
    common=np.logical_and.reduce([oo[n] for n in SIZES])
    results={}
    for n in SIZES:
        results[str(n)]={"n_total":int(common.size),"common":summarize(vv[n],common),"all_defined":summarize(vv[n],oo[n])}
    out=dict(status="PORTAL_FIXED_HISTORICAL_DIRECTION_NOT_CAUSAL",source_commit=audit["source_commit"],calibration_moons=[em[i] for i in cal],
       target_pairs=[[em[i],em[j]] for i,j in pp],n_calibration=len(cal),n_target_pairs=len(pp),results=results)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False)+"\n")
    brief={**out,"results":{n:{k:({j:v for j,v in q.items() if j!="per_date"} if isinstance(q,dict) else q) for k,q in r.items()} for n,r in results.items()}}
    print("PORTAL_STATIC_DIRECTION_RESULT",json.dumps(brief,allow_nan=False),flush=True)
if __name__=="__main__":main()
