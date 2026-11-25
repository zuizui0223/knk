"""Fixed factorial transfer of response allocation and sign histories across six studies."""
import argparse,json
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_replication_scale as rs
EPS=1e-12
SOURCES=("pond","kansas","poem","field","mancuso","hayashi")
KEYS=("OO","RO","OR","RR","old_mass_TV","recent_mass_TV","mass_first","mass_second","mass_symmetric","sign_first","sign_second","sign_symmetric","total_gain","interaction","hybrid_over_lag","hybrid_over_old","mass_only_gain")
def norm(x):
    x=np.asarray(x,float);s=x.sum(-1,keepdims=True)
    return np.divide(x,s,out=np.zeros_like(x),where=s>EPS),s[...,0]>EPS
def signmass(m,s):
    return np.concatenate((m*np.where(s>0,1.,np.where(s<0,0.,.5)),
                           m*np.where(s<0,1.,np.where(s>0,0.,.5))),axis=-1)
def tv(p,q):return .5*np.abs(p-q).sum(-1)
def factorial(old,recent,target):
    a,oa=norm(np.abs(old));b,ob=norm(np.abs(recent));y,oy=norm(np.abs(target))
    q=signmass(y,target)
    oo=tv(signmass(a,old),q);ro=tv(signmass(b,old),q)
    orr=tv(signmass(a,recent),q);rr=tv(signmass(b,recent),q)
    om=tv(a,y);rm=tv(b,y)
    mf=oo-ro;ms=orr-rr;sf=oo-orr;ss=ro-rr
    m=(mf+ms)/2;s=(sf+ss)/2
    total=oo-rr
    assert np.max(np.abs(total-m-s))<1e-10
    return np.stack((oo,ro,orr,rr,om,rm,mf,ms,m,sf,ss,s,total,mf-ms,-ss,mf,om-rm),axis=-1),oa&ob&oy
def selftest():
    x=np.array([1.,-2.,3.]);y=np.array([-2.,4.,0.])
    v,ok=factorial(x,y,y);r=dict(zip(KEYS,v.tolist()))
    assert ok and abs(r["RR"])<EPS and abs(r["total_gain"]-r["mass_symmetric"]-r["sign_symmetric"])<EPS
    assert abs(r["hybrid_over_lag"]+r["sign_second"])<EPS
    z,k=factorial(x,x,x);assert k and np.allclose(z,0)
    assert not factorial(np.zeros(3),y,y)[1]
    # Exact sign hedge on an absent taxon with positive mass.
    assert np.allclose(signmass(np.array([.3,.7]),np.array([0.,1.])),[.15,.7,.15,0.])
    print("Factorial and loss-decomposition unit checks PASS")
def reduce(v,ok):
    assert v.shape[:-1]==ok.shape
    if not ok.any(axis=(0,1)).all():raise ValueError("STOP_PANEL_DATE_NO_COMMON_DEFINED_CELLS")
    return dict(zip(KEYS,np.mean([v[:,:,i,:][ok[:,:,i]].mean(0) for i in range(v.shape[2])],axis=0).tolist()))
def summarize_panels(p):
    return {"panels":len(p),"n_common":sum(x["n_common"] for x in p.values()),
            "n_total":sum(x["n_total"] for x in p.values()),
            "mean":{k:float(np.mean([v["mean"][k] for v in p.values()])) for k in KEYS},
            "by_panel":p}
def entry(v,k,dates):
    return {"n_common":int(k.sum()),"n_total":int(k.size),"n_dates":len(dates),"target_dates":list(map(str,dates)),"mean":reduce(v,k)}
def source(st,prior):
    ps=rs.panels[st];group_ids=sorted({k for p in ps for k in p["groups"]})
    rng=np.random.default_rng(rs.SEED+rs.SOURCES.index(st))
    orders={k:np.array([rng.permutation(rs.groups[k]) for _ in range(rs.N)]) for k in group_ids}
    out={"one_observation":{},"cumulative_older":{},"one_observation_longer_dates":{},"cumulative_older_longer_dates":{},"equal_two_nonoverlap":{}}
    oldmetrics=[];oldn=0
    for p in ps:
        d,h=rs.panel_contrasts(p,orders,2,2);T=d.shape[2];dates=p["dates"]
        assert T>=3 and dates==p["dates"]
        rv,rok=rs.scores(d[:,:,1:-1],d[:,:,2:],h[:,:,2:])
        assert rok.all(),"frozen replay expected full support"
        oldmetrics.append(np.mean(rv,axis=(0,1,2)));oldn+=rok.size
        # t=2,3,...; all predictions strictly prior to target.
        recent=d[:,:,1:-1];target=h[:,:,2:]
        old=d[:,:,:-2]
        cumulative=np.stack([d[:,:,:t-1].mean(2) for t in range(2,T)],axis=2)
        v1,k1=factorial(old,recent,target);v2,k2=factorial(cumulative,recent,target)
        common=k1&k2
        out["one_observation"][p["key"]]=entry(v1,common,dates[2:])
        out["cumulative_older"][p["key"]]=entry(v2,common,dates[2:])
        if T>=4:
            out["one_observation_longer_dates"][p["key"]]=entry(v1[:,:,1:],common[:,:,1:],dates[3:])
            out["cumulative_older_longer_dates"][p["key"]]=entry(v2[:,:,1:],common[:,:,1:],dates[3:])
        if T>=5:
            a=(d[:,:,:-4]+d[:,:,1:-3])/2
            b=(d[:,:,2:-2]+d[:,:,3:-1])/2
            v,k=factorial(a,b,h[:,:,4:])
            out["equal_two_nonoverlap"][p["key"]]=entry(v,k,dates[4:])
    z=np.mean(oldmetrics,axis=0)
    ref=prior["results"][st]["configs"]["2+2"]
    assert oldn==ref["n_common"]==ref["n_total"]
    assert all(abs(z[i]-ref["common"][k])<1e-10 for i,k in enumerate(rs.METRICS)),(st,z,ref["common"])
    res={}
    for mode,p in out.items():
        if not p:res[mode]={"status":"STOP_INSUFFICIENT_POST_TREATMENT_DATES"};continue
        assert len(p)==len(ps)
        res[mode]=summarize_panels(p)
        res[mode]["status"]="SCORED_DESCRIPTIVE"
    return {"prior_replay":{"n":oldn,"metrics":dict(zip(rs.METRICS,z.tolist()))},"comparisons":res}
def main():
    p=argparse.ArgumentParser()
    for k in ("root","poem","field","mancuso","fungal","hayashi","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();selftest();a.root.mkdir(parents=True,exist_ok=True);rs.load(a)
    prior=json.loads(Path("validation/ahif_v1/replication_scale_summary_20261125.json").read_text())
    assert sum(len(rs.panels[s]) for s in SOURCES)==33
    assert all(len(p["dates"])==3 for p in rs.panels["saine"])
    results={}
    for st in SOURCES:
        z=source(st,prior);results[st]=z
        brief={**z,"comparisons":{k:{j:v for j,v in q.items() if j!="by_panel"} for k,q in z["comparisons"].items()}}
        print("MASS_DIRECTION_TRANSFER_RESULT",st,json.dumps(brief,allow_nan=False),flush=True)
    z={"status":"RETROSPECTIVE_FACTORIAL_TRANSFER_NOT_META","source_families":len(SOURCES),"panels":33,
       "saine":{"status":"STOP_ONLY_TWO_POST_INTRODUCTION_DATES","panels":10},"results":results}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
if __name__=="__main__":main()
