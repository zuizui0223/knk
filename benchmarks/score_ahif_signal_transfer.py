"""Training-only contrast reliability, source-held predictions; not causal/meta inference."""
import argparse,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from benchmarks import score_ahif_cross_source_concentration as src
from benchmarks import score_ahif_response_specificity as spec
from benchmarks import score_ahif_mancuso_response as m
from benchmarks import score_ahif_response_allocation as alloc
EPS=1e-12
ABS=("pond","kansas","poem","field")
MODELS={"intercept":(),"structure":tuple(range(5)),"strength":tuple(range(6)),"reliability":tuple(range(7)),"noise":tuple(range(6))+(7,)}
TARGETS=("sign_gain","same")
seen=[]
def features(A,B,paired):
    am=A.mean(0);bm=B.mean(0);d=bm-am;gross=np.abs(d).sum();state=(am+bm)/2;T=state.sum()
    Q=float(d@d)
    if paired:
        W=float((B-A).var(0,ddof=1).sum()/len(A))
    else:
        W=float(A.var(0,ddof=1).sum()/len(A)+B.var(0,ddof=1).sum()/len(B))
    if T<=EPS or Q+W<=EPS:return None
    f=[np.log1p(len(d)),np.log1p((state>EPS).sum()),np.log1p((np.abs(d)>EPS).sum()),
        np.max(np.abs(d))/gross if gross>EPS else 0.,abs(d.sum())/gross if gross>EPS else 0.,
        np.log1p(1000*Q/T**2),(Q-W)/(Q+W),np.log1p(1000*W/T**2)]
    assert all(np.isfinite(f)) and -1-1e-10<=f[6]<=1+1e-10
    return list(map(float,f))
def panel(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group=False):
    old=spec.evaluate(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group)
    seen.append((st,key,len(dates)))
    fs={}
    for fold,(ia,ib,ha,hb) in enumerate(src.splits(a.shape[1],paired)):
        assert not set(ia)&set(ha) and not set(ib)&set(hb)
        if paired:assert tuple(ia)==tuple(ib)
        for t in range(2,len(dates)):
            fs[(fold,str(dates[t]))]=features(a[t-1,list(ia)],b[t-1,list(ib)],paired)
    for r in old:
        r["features"]=fs[r["fold"],r["time"]]
        r["sign_gain"]=r["signblind"]-r["signed"] if r["signblind"] is not None and r["signed"] is not None else None
    return old
def weights(rows):
    count=Counter((r["source"],r["panel"],r["time"]) for r in rows)
    times=Counter((s,p) for s,p,t in count)
    panels=Counter(s for s,p in times)
    w=np.array([1/(len(panels)*panels[r["source"]]*times[r["source"],r["panel"]]*count[r["source"],r["panel"],r["time"]]) for r in rows])
    assert abs(w.sum()-1)<1e-9
    return w
def fit(train,test,cols,target):
    w=weights(train);X=np.array([r["features"] for r in train])[:,cols];Z=np.array([r["features"] for r in test])[:,cols]
    mu=w@X;sd=np.sqrt(w@((X-mu)**2));sd=np.where(sd>EPS,sd,1.)
    A=np.column_stack((np.ones(len(X)),(X-mu)/sd));B=np.column_stack((np.ones(len(Z)),(Z-mu)/sd))
    y=np.array([r[target] for r in train]);pen=.1*np.eye(A.shape[1]);pen[0,0]=0
    coef=np.linalg.solve(A.T@(w[:,None]*A)+pen,A.T@(w*y));p=B@coef
    return (np.clip(p,0,1) if target=="same" else p),coef.tolist()
def assess(train,test,target):
    w=weights(test);y=np.array([r[target] for r in test]);out={};pred={}
    for model,cols in MODELS.items():
        p,coef=fit(train,test,cols,target)
        out[model]=float(w@((p-y)**2));pred[model]=p.tolist()
    return out,pred
def aggregate(by):
    mse={k:float(np.mean([v[k] for v in by.values()])) for k in MODELS}
    gains={name:(1-mse[a]/mse[b] if mse[b]>EPS else None) for name,a,b in [
      ("strength_over_structure","strength","structure"),("reliability_over_strength","reliability","strength"),("noise_over_strength","noise","strength"),
      ("reliability_over_intercept","reliability","intercept")]}
    return dict(mse=mse,gains=gains,by_source=by)
def cv(rows,groups,target):
    by={};pred=[]
    for group in groups:
        tr=[r for r in rows if r["source"] not in group];te=[r for r in rows if r["source"] in group]
        assert tr and te
        _,p=assess(tr,te,target)
        for s in group:
            ids=[i for i,r in enumerate(te) if r["source"]==s]
            w=weights([te[i] for i in ids]);y=np.array([te[i][target] for i in ids])
            by[s]={k:float(w@((np.array(p[k])[ids]-y)**2)) for k in MODELS}
        pred.extend([dict(source=r["source"],panel=r["panel"],time=r["time"],fold=r["fold"],target=r[target],pred={k:p[k][i] for k in p}) for i,r in enumerate(te)])
    assert len(pred)==len(rows)
    return aggregate(by),pred
def corr(rows,target):
    w=weights(rows);u=np.array([r["features"][6] for r in rows]);y=np.array([r[target] for r in rows])
    x=u-w@u;z=y-w@y;den=np.sqrt((w@(x*x))*(w@(z*z)))
    return float(w@(x*z)/den) if den>EPS else None
def selftest():
    # Constant contrast vs heterogeneous training observations.
    a=np.array([[2.,1.],[2.,1.]]);b=np.array([[1.,2.],[1.,2.]])
    assert features(a,b,False)[6]==1
    a=np.array([[1.,3.],[3.,1.]]);b=a.copy()
    assert features(a,b,False)[6]==-1
    assert features(np.zeros((2,2)),np.zeros((2,2)),False) is None
    r=dict(source="a",panel="p",time="1",features=[0]*8,sign_gain=.2,same=.8)
    assert np.allclose(weights([r,r]),[.5,.5])
    p,_=fit([r,r],[r],(),"sign_gain");assert abs(p[0]-.2)<1e-12
    print("Training-only reliability synthetic checks PASS")
def main():
    p=argparse.ArgumentParser()
    for k in ("root","poem","field","mancuso","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();selftest();a.root.mkdir(parents=True,exist_ok=True)
    src.panel=panel
    rr=src.get_sources(a.root,a.poem,a.field,Path("database/ahif_v1"),Path("validation/ahif_v1"))
    assert len(seen)==19
    old=json.loads(Path("validation/ahif_v1/response_specificity_summary_20261125.json").read_text())
    for target in ("sign_gain","same"):
        by=alloc.nested(rr,target)
        ref=old["signed"]["by_source"][target]
        assert by.keys()==ref.keys() and all(abs(by[s]-ref[s])<1e-10 for s in ref),(target,by,ref)
    g,names,named,audit=m.source(a.mancuso);mats,rosters,coverage,pairs=m.panels(g,names,named)
    assert set(pairs)==m.EXPECTED
    mr=[]
    for r,f in pairs:
        A=mats[r,600];B=mats[r,f]
        mr+=panel("mancuso",f"{r}:{f}",m.DAYS,A,B,names,np.ones((len(names),1)),["all"],False,B.mean(1)-A.mean(1))
    assert len(seen)==27
    assert abs(alloc.nested(mr,"sign_gain")["mancuso"]-.1617459982758671)<1e-10
    assert abs(alloc.nested(mr,"same")["mancuso"]-.7999325258366095)<1e-10
    total=Counter(r["source"] for r in rr+mr)
    valid=[r for r in rr+mr if r["features"] is not None and all(r[k] is not None for k in TARGETS)]
    R=[r for r in valid if r["source"] in ABS];M=[r for r in valid if r["source"]=="mancuso"]
    assert set(r["source"] for r in R)==set(ABS) and M
    results={};preds={}
    for target in TARGETS:
        z={};pp={}
        for name,rows,groups in [
          ("loo",R,[(s,) for s in ABS]),
          ("grouped",R,[("pond",),("kansas",),("poem","field")]),
          ("exclude_kansas",[r for r in R if r["source"]!="kansas"],[(s,) for s in ABS if s!="kansas"])]:
            z[name],pp[name]=cv(rows,groups,target)
        for name,train,test in [
          ("microbial_transfer",R,M),
          ("microbial_depth_panel_transfer",R,[r for r in M if r["panel"] in ("0.2:4","0.4:1","0.4:4")]),
          ("microbial_transfer_excl_kansas",[r for r in R if r["source"]!="kansas"],M)]:
            assert train and test
            loss,p=assess(train,test,target);z[name]=aggregate({"mancuso":loss})
            pp[name]=[dict(source=r["source"],panel=r["panel"],time=r["time"],fold=r["fold"],target=r[target],pred={k:p[k][i] for k in p}) for i,r in enumerate(test)]
        z["descriptive_U_correlation"]={s:corr([r for r in valid if r["source"]==s],target) for s in ABS+("mancuso",)}
        results[target]=z;preds[target]=pp
    summary=dict(status="RETROSPECTIVE_CROSS_SOURCE_PREDICTION_NOT_META",panels=seen,total_cells=dict(total),eligible=dict(Counter(r["source"] for r in valid)),results=results)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps({**summary,"cells":valid,"predictions":preds},allow_nan=False)+"\n")
    print("AHIF_SIGNAL_TRANSFER_RESULT",json.dumps(summary,allow_nan=False))
if __name__=="__main__":main()
