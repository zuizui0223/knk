"""One-lag directional forecasts on fixed source rosters, plus Saine fungal logs.
No formal pooling or inference from dependent folds.
"""
import argparse,csv,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from benchmarks import score_ahif_cross_source_concentration as src
from benchmarks import score_ahif_response_specificity as spec
from benchmarks import score_ahif_response_allocation as alloc
from benchmarks import score_ahif_mancuso_response as m
from benchmarks.score_ahif_saine2025_succession import design
EPS=1e-12
SEED=2026112504
N_SPLITS=256
METRICS=("gain","signed","blind","same","opposite","abstain","allocation_tv")
EXCLUDED={"Physisporinus_crocatus","Steccherinum_collabens"}
EXPECTED={"pond":4,"kansas":1,"poem":6,"field":8,"mancuso":8}
oldrows=[];newrows=[];seen=[]
def norm(v):
    s=v.sum()
    return v/s if s>EPS else None
def signed(v):
    return np.r_[np.maximum(v,0),np.maximum(-v,0)]
def blind(p):
    a,b=np.split(p,2);return np.r_[(a+b)/2,(a+b)/2]
def tv(a,b):return float(np.abs(a-b).sum()/2)
def evaluate(d,h):
    p=norm(signed(d));q=norm(signed(h))
    if p is None or q is None:return {k:None for k in METRICS}
    v=tv(p,q);b=tv(blind(p),q)
    x,y=np.split(p,2);u,w=np.split(q,2)
    gross=np.abs(h);s=gross.sum();agree=np.sign(d)*np.sign(h)
    same=float(gross[agree>0].sum()/s);opposite=float(gross[agree<0].sum()/s);abstain=float(gross[agree==0].sum()/s)
    assert abs(same+opposite+abstain-1)<1e-10
    return dict(gain=b-v,signed=v,blind=b,same=same,opposite=opposite,abstain=abstain,allocation_tv=tv(x+y,u+w))
def fixed_mass_baseline(d,pre,h):
    p=norm(np.abs(d));q=norm(signed(h))
    if p is None or q is None:return None
    v=np.r_[p*np.where(pre>0,1.,np.where(pre<0,0.,.5)),p*np.where(pre<0,1.,np.where(pre>0,0.,.5))]
    return tv(v,q)-tv(norm(signed(d)),q)
def register(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group=False):
    old=spec.evaluate(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group)
    for r in old:
        r["sign_gain"]=r["signblind"]-r["signed"] if r["signblind"] is not None and r["signed"] is not None else None
    oldrows.extend(old);seen.append((st,key,len(dates)))
    for fold,(ia,ib,ha,hb) in enumerate(src.splits(a.shape[1],paired)):
        assert not set(ia)&set(ha) and not set(ib)&set(hb)
        d=b[:,list(ib)].mean(1)-a[:,list(ia)].mean(1)
        h=b[:,list(hb)].mean(1)-a[:,list(ha)].mean(1)
        for t in range(1,len(dates)):
            newrows.append(dict(source=st,panel=key,time=str(dates[t]),fold=fold,primary=t>=2,**evaluate(d[t-1],h[t])))
    return []
def harmonized(args):
    src.panel=register
    src.get_sources(args.root,args.poem,args.field,Path("database/ahif_v1"),Path("validation/ahif_v1"))
    g,names,named,audit=m.source(args.mancuso)
    mats,rosters,coverage,pairs=m.panels(g,names,named)
    assert set(pairs)==m.EXPECTED
    for r,f in pairs:
        A=mats[r,600];B=mats[r,f]
        register("mancuso",f"{r}:{f}",m.DAYS,A,B,names,np.ones((len(names),1)),["all"],False,B.mean(1)-A.mean(1))
    assert Counter(z[0] for z in seen)==Counter(EXPECTED)
    ref=json.loads(Path("validation/ahif_v1/response_specificity_summary_20261125.json").read_text())
    mr=json.loads(Path("validation/ahif_v1/mancuso2021_response_summary_20261125.json").read_text())
    for k in ("sign_gain","same"):
        v=alloc.nested(oldrows,k)
        exp={**ref["signed"]["by_source"][k],"mancuso":mr["primary"]["direction"]["source_equal"][k]}
        assert v.keys()==exp.keys() and all(abs(v[s]-exp[s])<1e-10 for s in v),(k,v,exp)
    out={}
    for mode,rr in (("primary",[r for r in newrows if r["primary"]]),("all_adjacent",newrows)):
        ok=[r for r in rr if all(r[k] is not None for k in METRICS)]
        out[mode]=dict(total=dict(Counter(r["source"] for r in rr)),defined=dict(Counter(r["source"] for r in ok)),by_metric={k:alloc.nested(ok,k) for k in METRICS})
    return out
def source_fungal(path):
    with (path/"meta.csv").open() as f:meta=list(csv.DictReader(f))
    names=(path/"otu_names.txt").read_text().splitlines()
    excluded=(path/"excluded_target_otus.txt").read_text().splitlines()
    assert len(names)==len(set(names))==4653 and len(excluded)==9
    x=np.fromfile(path/"counts.bin",dtype="<f8").reshape(1380,4653)
    ix,t,strata,sites,types,good,coverage,reads=design(meta,x)
    assert good.sum()==446 and (t[good]=="control").sum()==99
    raw=x[ix[good]];prop=raw/raw.sum(2,keepdims=True)
    old=json.loads(Path("validation/ahif_v1/saine2025_nontarget_succession_v1.json").read_text())
    assert old["complete_triplets"]==446 and coverage==old["coverage_by_stratum_treatment"]
    t=t[good];ss=strata[good];site=sites[good]
    assert len(set(site))==5 and len(set(ss))==10
    single={r["treatment"] for r in coverage if r["complete"]<2}
    assert single==EXCLUDED,single
    return prop,t,ss,site,coverage
def cells_for(site,t,ss,sites,mode):
    groups={}
    for s in sorted(set(ss[sites==site])):
        ix=np.flatnonzero(ss==s)
        ctrl=ix[t[ix]=="control"]
        labels=sorted(set(t[ix])-{"control"}-EXCLUDED)
        assert len(labels)==7
        if mode=="binary":
            groups[s]=[(ctrl,-1.),(ix[t[ix]!="control"],1.)]
        else:
            groups[s]=[(ctrl,-1.)]+[(ix[t[ix]==label],1/7) for label in labels]
        assert all(len(ids)>=2 for ids,c in groups[s])
    assert len(groups)==2
    return groups
def weights(groups,n,rng=None):
    W=np.zeros((2*N_SPLITS if rng is not None else 1,n))
    for s,groups_s in groups.items():
        for ids,coef in groups_s:
            if rng is None:
                W[0,ids]+=coef/(len(ids)*len(groups))
            else:
                for j in range(N_SPLITS):
                    z=rng.permutation(ids);a,b=np.array_split(z,2)
                    assert len(a)>0 and len(b)>0 and not set(a)&set(b)
                    W[2*j,a]+=coef/(len(a)*len(groups))
                    W[2*j+1,b]+=coef/(len(b)*len(groups))
    assert np.max(np.abs(W.sum(1)))<1e-10
    assert np.max(np.abs(np.abs(W).sum(1)-2))<1e-10
    if rng is not None:assert not np.any((W[::2]!=0)&(W[1::2]!=0))
    return W
def contrast(W,prop):
    # One physical-log weight across all dates, no log duplication.
    return (W@prop.reshape(len(prop),-1)).reshape(len(W),3,prop.shape[2])
def diagnostic(tr,te):
    lag=evaluate(tr[1],te[2])
    placebo=evaluate(tr[0],te[2])
    concurrent=evaluate(tr[2],te[2])
    change=evaluate(tr[1]-tr[0],te[2]-te[0])
    return {**{f"lag_{k}":v for k,v in lag.items()},**{f"pre_{k}":v for k,v in placebo.items()},
            **{f"concurrent_{k}":v for k,v in concurrent.items()},**{f"change_{k}":v for k,v in change.items()},
            "lag_vs_pre_sign_fixed_mass":fixed_mass_baseline(tr[1],tr[0],te[2]),
            "lag_vs_pre_gain":lag["gain"]-placebo["gain"] if lag["gain"] is not None and placebo["gain"] is not None else None}
def summarize(rows):
    names=[k for k in rows[0] if k not in ("site","split","orientation")]
    by={}
    for s in sorted(set(r["site"] for r in rows)):
        rr=[r for r in rows if r["site"]==s]
        by[s]={k:float(np.mean([r[k] for r in rr if r[k] is not None])) if any(r[k] is not None for r in rr) else None for k in names}
        by[s]["n_cells"]=len(rr)
        by[s]["defined"]={k:sum(r[k] is not None for r in rr) for k in names}
    # do not silently remove a site
    mean={k:float(np.mean([r[k] for r in by.values()])) if all(r[k] is not None for r in by.values()) else None for k in names}
    return dict(n_sites=len(by),n_cells=len(rows),site_equal=mean,by_site=by)
def fungal(path):
    prop,t,ss,sites,coverage=source_fungal(path)
    modes={};cells={}
    for j,mode in enumerate(("binary","balanced7")):
        rng=np.random.default_rng(SEED+j)
        site_names=sorted(set(sites));rw=[];site_vectors={};site_coverage={}
        for s in site_names:
            groups=cells_for(s,t,ss,sites,mode)
            used=np.unique(np.concatenate([ids for gg in groups.values() for ids,_ in gg]))
            site_coverage[s]={"n_logs":len(used),"n_intro":int((t[used]!="control").sum()),"n_control":int((t[used]=="control").sum()),"n_excluded":int((sites==s).sum()-len(used))}
            W=weights(groups,len(t),rng)
            D=contrast(W,prop)
            for k in range(N_SPLITS):
                for o in (0,1):
                    rw.append(dict(site=s,split=k,orientation=o,**diagnostic(D[2*k+o],D[2*k+1-o])))
            site_vectors[s]=contrast(weights(groups,len(t)),prop)[0]
        held=[]
        for s in site_names:
            tr=np.mean([site_vectors[z] for z in site_names if z!=s],axis=0)
            held.append(dict(site=s,**diagnostic(tr,site_vectors[s])))
        modes[mode]=dict(coverage=site_coverage,n_logs=sum(v["n_logs"] for v in site_coverage.values()),
                         within_site=summarize(rw),held_site=summarize(held))
        cells[mode]={"within_site":rw,"held_site":held}
    assert modes["binary"]["n_logs"]==446 and modes["binary"]["within_site"]["n_cells"]==2560
    return modes,cells
def selftest():
    d=np.array([1.,-1.]);r=evaluate(d,d)
    assert r["gain"]==.5 and r["same"]==1 and r["allocation_tv"]==0
    r=evaluate(d,-d);assert r["gain"]==-.5 and r["opposite"]==1
    assert fixed_mass_baseline(d,-d,d)==1
    assert all(x is None for x in evaluate(d,np.zeros(2)).values())
    g={"s1":[(np.array([0,1]),-1.),(np.array([2,3]),1.)],"s2":[(np.array([4,5]),-1.),(np.array([6,7]),1.)]}
    W=weights(g,8,np.random.default_rng(2));assert W.shape==(512,8)
    print("One-lag sign and physical-disjoint split tests PASS")
def main():
    p=argparse.ArgumentParser()
    for k in ("root","poem","field","mancuso","fungal","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();selftest();a.root.mkdir(parents=True,exist_ok=True)
    old=harmonized(a);fm,fc=fungal(a.fungal)
    z=dict(status="HARMONIZED_ONE_LAG_RETROSPECTIVE_NOT_META",seed=SEED,n_splits=N_SPLITS,
           previous_sources=old,saine=fm)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps({**z,"previous_cells":newrows,"fungal_cells":fc},allow_nan=False)+"\n")
    print("AHIF_ONE_LAG_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
