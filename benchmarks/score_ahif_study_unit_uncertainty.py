"""Study-level effect summaries with fixed physical-unit weight perturbations.
Intervals are conditional sensitivity, never calibrated frequentist CIs or meta-analysis SEs.
"""
import argparse,csv,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from benchmarks import score_ahif_cross_source_concentration as src
from benchmarks import score_ahif_response_specificity as spec
from benchmarks import score_ahif_response_allocation as alloc
from benchmarks import score_ahif_mancuso_response as m
EPS=1e-12
N=999; SEED=2026112503; BATCH=32
METRICS=("sign_gain","same","recent_sign_gain")
SOURCES=("pond","kansas","poem","field","mancuso")
EXPECTED_CLUSTERS={"pond":20,"kansas":4,"poem":40,"field":10,"mancuso":44}
EXPECTED_PANELS={"pond":4,"kansas":1,"poem":6,"field":8,"mancuso":8}
EXPECTED_CELLS={"pond":1872,"kansas":78,"poem":1200,"field":320,"mancuso":1152}
META={
 "pond":("Narwani et al. 2019","absolute_biovolume","Narwani","Narwani publication"),
 "kansas":("Foster Kansas hay meadow, Ladouceur compilation","absolute_biomass_conditional","Foster","Original linkage unresolved"),
 "poem":("Alonso-Crespo et al. POEM","absolute_biomass","Delory_Temperton","Two cohorts, one publication"),
 "field":("van Steijn et al. 2026","unclosed_pin_hits","Related_grassland","Two sites, one publication"),
 "mancuso":("Mancuso et al. 2021 experiment 1","relative_read_composition","Mancuso","One publication; Exp2 excluded"),
}
panels=defaultdict(list); arm_arrays={}; groups={}
def keys(st,key,n):
    if st=="pond":
        b,a=key.split("-");return f"pond/{a}",f"pond/{b}",False
    if st=="kansas":return "kansas/Control","kansas/Seeds",True
    if st=="poem":
        c,b,a=key.split("_");return f"poem/{c}/{a}",f"poem/{c}/{b}",False
    if st=="field":
        s,b,a=key.split("_");return f"field/{s}/{a}",f"field/{s}/{b}",True
    if st=="mancuso":
        r,f=key.split(":");return f"mancuso/{r}:600",f"mancuso/{r}:{f}",False
    raise ValueError(st)
def register(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group=False):
    old=spec.evaluate(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group)
    ka,kb,isp=keys(st,key,a.shape[1]);assert paired==isp
    for k,x in [(ka,a),(kb,b)]:
        if k in arm_arrays: assert np.array_equal(arm_arrays[k],x),("shared arm drift",k)
        else:arm_arrays[k]=x
        # Paired block indices align across all treatments within each site.
        group=k.rsplit("/",1)[0]+"/BLOCK" if paired else k
        if group in groups:assert groups[group]==x.shape[1]
        else:groups[group]=x.shape[1]
    ga=ka.rsplit("/",1)[0]+"/BLOCK" if paired else ka
    gb=kb.rsplit("/",1)[0]+"/BLOCK" if paired else kb
    for r in old:
        r["sign_gain"]=r["signblind"]-r["signed"] if r["signblind"] is not None and r["signed"] is not None else None
        r["recent_sign_gain"]=r["signed_history"]-r["signed"] if r["signed_history"] is not None and r["signed"] is not None else None
    assert all(all(r[k] is not None for k in METRICS) for r in old)
    assert len(old)==len(list(src.splits(a.shape[1],paired)))*(len(dates)-2)
    panels[st].append(dict(key=key,a=a,b=b,ga=ga,gb=gb,paired=paired,dates=dates,old=old))
    return old
def masks(n,paired):
    out=np.zeros((4,len(list(src.splits(n,paired))),n))
    for j,indices in enumerate(src.splits(n,paired)):
        ia,ib,ha,hb=indices
        assert not set(ia)&set(ha) and not set(ib)&set(hb)
        if paired: assert ia==ib and ha==hb
        for i,v in enumerate(indices):out[i,j,list(v)]=1
    return out
def means(x,w,mask):
    z=w[:,None,:]*mask[None,:,:];s=z.sum(-1)
    assert (s>0).all()
    return np.einsum("bfn,tnk->bftk",z,x,optimize=False)/s[:,:,None,None]
def normlast(x):
    s=x.sum(-1,keepdims=True)
    return np.divide(x,s,out=np.zeros_like(x),where=s>EPS),s[...,0]>EPS
def signed(x):
    return np.concatenate((np.maximum(x,0),np.maximum(-x,0)),axis=-1)
def evaluate_panel(p,wa,wb):
    a,b=p["a"],p["b"];mask=masks(a.shape[1],p["paired"])
    train=means(b,wb,mask[1])-means(a,wa,mask[0])
    hold=means(b,wb,mask[3])-means(a,wa,mask[2])
    # exact original historical signed mean before norm (not mean of normalized profiles)
    ts=signed(train)
    t=np.arange(1,train.shape[2]-1,dtype=float)[None,None,:,None]
    base,okb=normlast(np.cumsum(ts,axis=2)[:,:,:-2,:]/t)
    lag,okl=normlast(ts[:,:,1:-1,:])
    truth,okt=normlast(signed(hold[:,:,2:,:]))
    ss=(base+lag)/2
    pos,neg=np.split(ss,2,axis=-1)
    half=(pos+neg)/2
    blind=np.concatenate((half,half),axis=-1)
    tv=lambda q:np.abs(truth-q).sum(-1)/2
    sign_gain=tv(blind)-tv(ss)
    recent_gain=tv(base)-tv(ss)
    gross=np.abs(hold[:,:,2:,:]);den=gross.sum(-1)
    agree=(np.sign(train[:,:,1:-1,:])*np.sign(hold[:,:,2:,:]))>0
    same=np.divide((gross*agree).sum(-1),den,out=np.zeros_like(den),where=den>EPS)
    ok=okb&okl&okt&(den>EPS)
    # Fold/time shape exact and fixed; any invalid cell invalidates whole draw.
    val=np.stack((sign_gain,same,recent_gain),axis=-1)
    return val.mean((1,2)),ok.all((1,2)),val
def group_layout(st):
    gs=sorted(k for k in groups if k.startswith(st+"/"))
    ids={};idx=0
    for k in gs:
        n=groups[k];ids[k]=list(range(idx,idx+n));idx+=n
    assert idx==EXPECTED_CLUSTERS[st],(st,idx,ids)
    return ids,idx
def eval_source(st,allw,ids,check=False):
    result=[];valid=np.ones(len(allw),dtype=bool)
    for p in panels[st]:
        values,ok,cells=evaluate_panel(p,allw[:,ids[p["ga"]]],allw[:,ids[p["gb"]]])
        if check:
            target=np.array([[r[k] for k in METRICS] for r in p["old"]])
            assert cells.shape[0]==1 and cells[0].reshape(-1,3).shape==target.shape
            err=np.max(np.abs(cells[0].reshape(-1,3)-target),axis=0)
            # Sign at a floating-point near-zero difference may differ; never silently waive replay.
            assert np.max(err)<1e-10,(st,p["key"],err.tolist())
        result.append(values);valid&=ok
    return np.mean(result,axis=0),valid
def chunks(st,W,ids):
    v=[];ok=[]
    for i in range(0,len(W),BATCH):
        x,q=eval_source(st,W[i:i+BATCH],ids);v.append(x);ok.append(q)
    return np.concatenate(v),np.concatenate(ok)
def summarize(st,rng,reference):
    ids,n=group_layout(st)
    point,ok=eval_source(st,np.ones((1,n)),ids,check=True)
    assert ok[0] and np.max(np.abs(point[0]-reference))<1e-10,(st,point,reference)
    W=rng.exponential(1.,(N,n))
    d,valid=chunks(st,W,ids)
    # Entire source-draw eligibility; no denominator drift.
    dv=d[valid];assert len(dv)>0
    L=np.ones((n,n));np.fill_diagonal(L,0)
    leave,lv=chunks(st,L,ids)
    lvx=leave[lv];assert len(lvx)>0
    out={}
    for j,k in enumerate(METRICS):
        a=dv[:,j];l=lvx[:,j]
        out[k]=dict(estimate=float(point[0,j]),perturbation_mean=float(a.mean()),
          perturbation_sd=float(a.std(ddof=1)),perturbation_percentiles=np.quantile(a,[.025,.975]).tolist(),
          fraction_positive=float((a>0).mean()),leave_one_min=float(l.min()),leave_one_max=float(l.max()),
          leave_one_positive=int((l>0).sum()),leave_one_valid=int(lv.sum()))
    return dict(n_clusters=n,cluster_groups={k:len(v) for k,v in ids.items()},
      n_panels=len(panels[st]),n_cells=sum(len(p["old"]) for p in panels[st]),
      n_draws=N,n_valid_draws=int(valid.sum()),n_invalid_draws=int((~valid).sum()),
      metric_order=list(METRICS),metrics=out,perturbation_covariance=np.cov(dv,rowvar=False).tolist(),
      leave_one_draws=[{"cluster":i,"valid":bool(lv[i]),"metrics":leave[i].tolist() if lv[i] else None} for i in range(n)],
      draws=[{"draw":i,"metrics":d[i].tolist() if valid[i] else None} for i in range(N)])
def selftest():
    # Fixed unit membership and shared weight: constant treatment contrast.
    x=np.ones((3,4,2));x[:,:,1]=3;y=x.copy();y[:,:,0]=3;y[:,:,1]=1
    p=dict(a=x,b=y,paired=False)
    a=np.array([[1.,2.,3.,4.],[4.,3.,2.,1.]])
    z,ok,_=evaluate_panel(p,a,a)
    assert ok.all() and np.allclose(z,[[.5,1,0],[.5,1,0]])
    assert masks(4,False).shape==(4,36,4)
    print("Physical-unit perturbation synthetic checks PASS")
def main():
    p=argparse.ArgumentParser()
    for k in ("root","poem","field","mancuso","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();selftest();a.root.mkdir(parents=True,exist_ok=True)
    src.panel=register
    rows=src.get_sources(a.root,a.poem,a.field,Path("database/ahif_v1"),Path("validation/ahif_v1"))
    g,names,named,audit=m.source(a.mancuso);mats,rosters,coverage,pairs=m.panels(g,names,named)
    assert set(pairs)==m.EXPECTED
    for r,f in pairs:
        A=mats[r,600];B=mats[r,f]
        rows+=register("mancuso",f"{r}:{f}",m.DAYS,A,B,names,np.ones((len(names),1)),["all"],False,B.mean(1)-A.mean(1))
    assert Counter(r["source"] for r in rows)==Counter(EXPECTED_CELLS)
    assert {k:len(v) for k,v in panels.items()}==EXPECTED_PANELS
    old=json.loads(Path("validation/ahif_v1/response_specificity_summary_20261125.json").read_text())
    oldm=json.loads(Path("validation/ahif_v1/mancuso2021_response_summary_20261125.json").read_text())
    rng=np.random.default_rng(SEED);results={}
    for st in SOURCES:
        ref=np.array([alloc.nested(rows,k)[st] for k in METRICS])
        archived=np.array([old["signed"]["by_source"][k][st] if st!="mancuso" else oldm["primary"]["direction"]["source_equal"][k] for k in METRICS])
        assert np.max(np.abs(ref-archived))<1e-10,(st,ref,archived)
        z=summarize(st,rng,ref);results[st]=z
        print("STUDY_UNIT_RESULT",st,json.dumps({k:v for k,v in z.items() if k not in ("draws","leave_one_draws")},allow_nan=False),flush=True)
    out={"status":"CONDITIONAL_UNIT_WEIGHT_PERTURBATION_NOT_META_SE","seed":SEED,"n_draws":N,"results":results}
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/"result.json").write_text(json.dumps(out,allow_nan=False)+"\n")
    ledger=[]
    for st in SOURCES:
        z=results[st];label,lane,lab,overlap=META[st]
        for k,v in z["metrics"].items():
            ledger.append(dict(source=st,study_label=label,measurement_lane=lane,lab_group=lab,dependence_note=overlap,
               n_clusters=z["n_clusters"],n_panels=z["n_panels"],n_dependent_cells=z["n_cells"],n_valid_draws=z["n_valid_draws"],
               endpoint=k,estimate=v["estimate"],perturbation_sd=v["perturbation_sd"],
               sensitivity_low=v["perturbation_percentiles"][0],sensitivity_high=v["perturbation_percentiles"][1],
               leave_one_min=v["leave_one_min"],leave_one_max=v["leave_one_max"],
               inverse_variance_meta_eligible="no",reason="conditional perturbation SD is not calibrated sampling SE; limited source independence and measurement comparability"))
    with (a.out/"study_effects.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(ledger[0]));w.writeheader();w.writerows(ledger)
    print("STUDY_UNIT_LEDGER",json.dumps(ledger,allow_nan=False),flush=True)
if __name__=="__main__":main()
