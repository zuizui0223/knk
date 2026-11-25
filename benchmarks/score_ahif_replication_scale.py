"""Matched train/held sample-size curves; no independent-fold inference."""
import argparse,json
from collections import defaultdict,Counter
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_cross_source_concentration as src
from benchmarks import score_ahif_study_unit_uncertainty as units
from benchmarks import score_ahif_mancuso_response as m
from benchmarks import score_ahif_hayashi2026 as hy
from benchmarks import score_ahif_one_lag_fungal as ol
SEED=2026112508;N=256
SIZES=(1,2,3,4,8,16,32)
CONFIGS=sorted({(k,k) for k in SIZES}|{(k,2) for k in SIZES}|{(2,k) for k in SIZES})
SOURCES=("pond","kansas","poem","field","mancuso","saine","hayashi")
EXPECTED={"pond":4,"kansas":1,"poem":6,"field":8,"mancuso":8,"saine":10,"hayashi":6}
METRICS=("G","H","same","allocation_tv","concurrent_G","concurrent_H","concurrent_minus_lag_G","fixed_mass_concurrent_advantage","lag_signed_tv","concurrent_signed_tv")
panels=defaultdict(list);arms={};groups={}
EPS=1e-12
def register(st,key,dates,a,b,ga,gb,pa,paired=False):
    # arm keys identify actual arrays; group keys identify shared index partitions.
    for arm,grp,x in [(ga,pa[0],a),(gb,pa[1],b)]:
        if arm in arms:assert np.array_equal(arms[arm],x),(st,arm,"shared arm drift")
        else:arms[arm]=x
        if grp in groups:assert groups[grp]==x.shape[1]
        else:groups[grp]=x.shape[1]
    assert len(dates)==a.shape[0]==b.shape[0] and a.shape[2]==b.shape[2]
    assert not paired or (pa[0]==pa[1] and a.shape[1]==b.shape[1])
    panels[st].append(dict(key=key,dates=list(dates),a=a,b=b,groups=pa,paired=paired))
def previous(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group=False):
    # Original score function includes source shape and prior numeric contrast checks.
    # Exact archive replay below via original one-lag scorer.
    ol.register(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group)
    ka,kb,isp=units.keys(st,key,a.shape[1]);assert isp==paired
    ga=ka.rsplit("/",1)[0]+"/BLOCK" if paired else ka
    gb=kb.rsplit("/",1)[0]+"/BLOCK" if paired else kb
    register(st,key,dates,a,b,ka,kb,(ga,gb),paired)
    return []
def load(args):
    src.panel=previous
    src.get_sources(args.root,args.poem,args.field,Path("database/ahif_v1"),Path("validation/ahif_v1"))
    g,names,named,audit=m.source(args.mancuso); mats,rosters,coverage,pairs=m.panels(g,names,named)
    assert set(pairs)==m.EXPECTED
    for r,f in pairs:
        a=mats[r,600];b=mats[r,f]
        previous("mancuso",f"{r}:{f}",m.DAYS,a,b,names,np.ones((len(names),1)),["all"],False,b.mean(1)-a.mean(1))
    ref=json.loads(Path("validation/ahif_v1/one_lag_fungal_summary_20261125.json").read_text())
    rr=[r for r in ol.newrows if r["primary"]]
    v=ol.alloc.nested(rr,"gain")
    assert v.keys()==ref["previous_sources"]["primary"]["by_metric"]["gain"].keys()
    assert all(abs(v[s]-ref["previous_sources"]["primary"]["by_metric"]["gain"][s])<1e-10 for s in v)
    # Fungal source's design/QC/coverage verified by source_fungal.
    prop,t,ss,sites,_=ol.source_fungal(args.fungal)
    for st in sorted(set(ss)):
        ids=np.flatnonzero(ss==st);a=ids[t[ids]=="control"];b=ids[t[ids]!="control"]
        ka,kb=f"saine/{st}/control",f"saine/{st}/introduced"
        register("saine",st,[2019,2020,2021],prop[a].transpose(1,0,2),prop[b].transpose(1,0,2),ka,kb,(ka,kb))
    # Author's fixed rarefied OTU99 complete wells; same ASV eligibility separately checked by exporter.
    a,ids,taxa=hy.read(args.hayashi,"otu99");a2,ids2,_=hy.read(args.hayashi,"asv")
    assert ids==ids2 and len(taxa)==331 and sum(map(len,ids.values()))==760
    old=json.loads(Path("validation/ahif_v1/hayashi2026_direction_summary_20261125.json").read_text())
    assert {f"{s}:{d}":len(x) for (s,d),x in ids.items()}==old["coverage"]
    for s in hy.SOURCES:
        for d in hy.LEVELS[1:]:
            ka,kb=f"hayashi/{s}/1.0",f"hayashi/{s}/{d}"
            register("hayashi",f"{s}:{d}",hy.DAYS,a[s,1.],a[s,d],ka,kb,(ka,kb))
    assert {s:len(p) for s,p in panels.items()}==EXPECTED
def scores(lag,now,target):
    # Arrays shape [iteration,orientation,date,taxon]; TV additive by taxon.
    def signed(d):return np.concatenate((np.maximum(d,0),np.maximum(-d,0)),axis=-1)
    def norm(d):
        tot=d.sum(-1,keepdims=True)
        return np.divide(d,tot,out=np.zeros_like(d),where=tot>EPS),tot[...,0]>EPS
    def tv(p,q):return .5*np.abs(p-q).sum(-1)
    def flip(p):
        a,b=np.split(p,2,axis=-1);return np.concatenate((b,a),axis=-1)
    p,ok1=norm(signed(lag));c,ok2=norm(signed(now));q,ok3=norm(signed(target))
    mass,mok=norm(np.abs(lag))
    ab=(np.abs(lag)>EPS)
    sp=np.sign(now)
    pc=np.concatenate((mass*np.where(sp>0,1,np.where(sp<0,0,.5)),mass*np.where(sp<0,1,np.where(sp>0,0,.5))),axis=-1)
    pl=tv(p,q);cl=tv(c,q)
    G=tv((p+flip(p))/2,q)-pl;H=.5*(tv(flip(p),q)-pl)
    CG=tv((c+flip(c))/2,q)-cl;CH=.5*(tv(flip(c),q)-cl)
    gross=np.abs(target);den=gross.sum(-1)
    same=np.divide((gross*(np.sign(lag)*np.sign(target)>0)).sum(-1),den,out=np.zeros_like(den),where=den>EPS)
    tq,_=norm(gross)
    out=np.stack((G,H,same,tv(mass,tq),CG,CH,CG-G,pl-tv(pc,q),pl,cl),axis=-1)
    ok=ok1&ok2&ok3&mok
    return out,ok
def mean_for(x,order,n,end):
    # One permutation per physical arm/block, reused across all configurations.
    z=order[:,:n] if end==0 else order[:,-n:]
    return x[:,z,:].mean(2).transpose(1,0,2)
def panel_contrasts(p,orders,n,h):
    ga,gb=p["groups"];oa,ob=orders[ga],orders[gb];a,b=p["a"],p["b"]
    assert n+h<=min(a.shape[1],b.shape[1])
    ds=[];hs=[]
    for o in (0,1):
        d=mean_for(b,ob,n,o)-mean_for(a,oa,n,o)
        q=mean_for(b,ob,h,1-o)-mean_for(a,oa,h,1-o)
        ds.append(d);hs.append(q)
    return np.stack(ds,axis=1),np.stack(hs,axis=1)
def evaluate_all(st):
    rng=np.random.default_rng(SEED+SOURCES.index(st))
    ps=panels[st]
    group_ids=sorted({k for p in ps for k in p["groups"]})
    capacity=min(groups[k] for k in group_ids)
    configs=[(n,h) for n,h in CONFIGS if n+h<=capacity]
    assert (2,2) in configs
    orders={k:np.array([rng.permutation(groups[k]) for _ in range(N)]) for k in group_ids}
    vals={};ok={};sitevals={};siteok={}
    for n,h in configs:
        dct={}
        for p in ps:
            d,q=panel_contrasts(p,orders,n,h)
            # t>=2 targets, previous original target dates.
            v,k=scores(d[:,:,1:-1],d[:,:,2:],q[:,:,2:])
            vals[n,h,p["key"]]=v;ok[n,h,p["key"]]=k
            if st=="saine":dct[p["key"]]=(d,q)
        if st=="saine":
            # Secondary: preserve previous equal-stratum site average before norm.
            for site in sorted({s.split("|")[0] for s in dct}):
                q=[x for k,x in dct.items() if k.split("|")[0]==site];assert len(q)==2
                d,hv=np.mean([x[0] for x in q],axis=0),np.mean([x[1] for x in q],axis=0)
                v,k=scores(d[:,:,1:-1],d[:,:,2:],hv[:,:,2:])
                sitevals[n,h,site]=v;siteok[n,h,site]=k
    def summarize(vv,oo,keys):
        common={k:np.logical_and.reduce([oo[n,h,k] for n,h in configs]) for k in keys}
        if any(not x.any(axis=(0,1)).all() for x in common.values()):
            raise ValueError(("no common support for panel-date",st))
        out={}
        for n,h in configs:
            pp={}
            for key in keys:
                v=vv[n,h,key];k=oo[n,h,key];cm=common[key]
                def reduce(mask):
                    # Equal dates, then folds. No selective loss of whole dates.
                    return np.mean([v[:,:,t,:][mask[:,:,t]].mean(0) for t in range(v.shape[2])],axis=0).tolist()
                pp[key]={"n_total":int(k.size),"n_common":int(cm.sum()),"n_defined":int(k.sum()),
                         "common":dict(zip(METRICS,reduce(cm))),"all_defined":dict(zip(METRICS,reduce(k)))}
            out[f"{n}+{h}"]={"by_panel":pp,"common":{k:float(np.mean([q["common"][k] for q in pp.values()])) for k in METRICS},
                              "all_defined":{k:float(np.mean([q["all_defined"][k] for q in pp.values()])) for k in METRICS},
                              "n_common":sum(q["n_common"] for q in pp.values()),"n_total":sum(q["n_total"] for q in pp.values()),
                              "n_defined":sum(q["n_defined"] for q in pp.values())}
        return out
    x=summarize(vals,ok,[p["key"] for p in ps])
    result=dict(capacity=capacity,n_panels=len(ps),n_group_units=sum(groups[k] for k in group_ids),configs=x)
    if st=="saine":result["site_average_secondary"]=summarize(sitevals,siteok,sorted({p["key"].split("|")[0] for p in ps}))
    return result
def selftest():
    d=np.array([1.,-1.]).reshape(1,1,1,2)
    a,k=scores(d,d,d);assert k.all() and np.allclose(a[0,0,0],[.5,.5,1.,0.,.5,.5,0.,0.,0.,0.])
    a,k=scores(d,-d,-d);assert k.all() and abs(a[0,0,0,0]+.5)<EPS and abs(a[0,0,0,7]-1)<EPS
    # Exact expected hard-sign reference equals average over all independent sign assignments.
    x=np.array([2.,-1.,.5]);y=np.array([1.,3.,-4.])
    import itertools
    z,_=scores(x.reshape(1,1,1,3),x.reshape(1,1,1,3),y.reshape(1,1,1,3))
    losses=[ol.evaluate(np.abs(x)*np.array(s),y)["signed"] for s in itertools.product((-1,1),repeat=3)]
    assert abs(z[0,0,0,1]-(np.mean(losses)-ol.evaluate(x,y)["signed"]))<EPS
    print("Matched sample-size and fair-hard-sign synthetic checks PASS")
def main():
    p=argparse.ArgumentParser()
    for k in ("root","poem","field","mancuso","fungal","hayashi","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();selftest();a.root.mkdir(parents=True,exist_ok=True);load(a)
    results={}
    for st in SOURCES:
        r=evaluate_all(st);results[st]=r
        brief={**r,"configs":{k:{j:v for j,v in q.items() if j!="by_panel"} for k,q in r["configs"].items()}}
        if "site_average_secondary" in brief:brief["site_average_secondary"]={k:{j:v for j,v in q.items() if j!="by_panel"} for k,q in r["site_average_secondary"].items()}
        print("REPLICATION_SCALE_RESULT",st,json.dumps(brief,allow_nan=False),flush=True)
    z={"status":"MATCHED_REPLICATION_DESCRIPTIVE_NOT_META","seed":SEED,"n_splits":N,"results":results}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
if __name__=="__main__":main()
