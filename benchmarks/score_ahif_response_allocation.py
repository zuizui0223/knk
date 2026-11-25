"""Separate experimental response allocation diagnostic; source parsers pinned by parent."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
import benchmarks.score_ahif_cross_source_concentration as src
EPS=1e-12
seen=[]
def norm(x):
    x=np.asarray(x,float); s=x.sum()
    return x/s if s>EPS else None
def tv(q,p): return None if q is None or p is None else float(np.abs(q-p).sum()/2)
def overlap(a,b):
    if a.sum()<=EPS or b.sum()<=EPS:return None
    x=a==a.max(); y=b==b.max()
    return float((x&y).sum()/(x.sum()*y.sum()))
def evaluate(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group=False):
    assert a.shape==b.shape and a.ndim==3 and a.shape[0]==len(dates) and a.shape[2]==len(names)
    assert np.isfinite(a).all() and np.isfinite(b).all() and (a>=0).all() and (b>=0).all()
    assert np.all(G.sum(1)==1)
    p=b.mean(1)-a.mean(1);q=p@G if prior_is_group else p
    assert np.allclose(q,prior,rtol=1e-9,atol=1e-9),key
    seen.append((st,key,len(dates)))
    rows=[]
    for fold,(ia,ib,ha,hb) in enumerate(src.splits(a.shape[1],paired)):
        assert not(set(ia)&set(ha)) and not(set(ib)&set(hb))
        A=a[:,list(ia)];B=b[:,list(ib)]
        am=A.mean(1);bm=B.mean(1)
        hist=np.abs(bm-am)
        held=np.abs(b[:,list(hb)].mean(1)-a[:,list(ha)].mean(1))
        for t in range(3,len(dates)):
            q=norm(held[t]); lag=norm(hist[t-1])
            H=norm(hist[:t].mean(0)); M1=norm(am[:t].mean(0)); M2=norm(np.sqrt(am[:t].mean(0)))
            V=norm(np.sqrt((A[:t].var(1,ddof=1)/len(ia)+B[:t].var(1,ddof=1)/len(ib)).mean(0)))
            old=[];recent=[]
            if q is not None and lag is not None:
                for j in range(t-1):
                    oldj=norm(hist[j]); base=norm(hist[[k for k in range(t-1) if k!=j]].mean(0))
                    if oldj is None or base is None: continue
                    old.append(tv(q,(base+oldj)/2));recent.append(tv(q,(base+lag)/2))
            vals=dict(M1=tv(q,M1),M2=tv(q,M2),V=tv(q,V),H=tv(q,H),L=tv(q,lag))
            n=len(old)
            rows.append(dict(source=st,panel=key,time=str(dates[t]),fold=fold,n_candidates=n,
                old=float(np.mean(old)) if n else None,recent=float(np.mean(recent)) if n else None,
                recency=float(np.mean(np.array(old)-recent)) if n else None,
                dominant_overlap=overlap(am[t-1],held[t]),lag_overlap=overlap(hist[t-1],held[t]),
                **vals))
    return rows
def nested(rows,col):
    # Equal within panel/time, then time, panel, source.
    cells={}
    for r in rows:
        v=r[col]
        if v is not None:cells.setdefault((r["source"],r["panel"],r["time"]),[]).append(v)
    panels={}
    for (s,p,t),v in cells.items():panels.setdefault((s,p),[]).append(float(np.mean(v)))
    ss={}
    for (s,p),v in panels.items():ss.setdefault(s,[]).append(float(np.mean(v)))
    return {s:float(np.mean(v)) for s,v in ss.items()}
def selftest():
    # two disjoint halves, four states; prior j and lag swapped gives opposite recency for equal base
    q=np.array([0.,1.]);b=np.array([.5,.5]);l=q;o=np.array([1.,0.])
    assert tv(q,(b+o)/2)-tv(q,(b+l)/2)==.5
    assert norm([0,0]) is None
    r=[dict(source="s",panel="p",time="1",x=1.),dict(source="s",panel="p",time="1",x=3.),dict(source="s",panel="p",time="2",x=4.)]
    assert nested(r,"x")=={"s":3.}
    print("Synthetic checks PASS")
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--poem",type=Path,required=True);p.add_argument("--field",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    selftest();a.root.mkdir(parents=True,exist_ok=True)
    src.panel=evaluate
    rows=src.get_sources(a.root,a.poem,a.field,Path("database/ahif_v1"),Path("validation/ahif_v1"))
    assert len(seen)==19 and len(set((x,y) for x,y,z in seen))==19
    valid=[r for r in rows if r["recency"] is not None and all(r[k] is not None for k in ("H","M2","V","M1","L"))]
    for r in valid:
        r["H_vs_M2"]=r["M2"]-r["H"];r["H_vs_V"]=r["V"]-r["H"]
    by={k:nested(valid,k) for k in ("recency","old","recent","M1","M2","V","H","L","H_vs_M2","H_vs_V","dominant_overlap","lag_overlap")}
    avg={k:float(np.mean(list(v.values()))) if v else None for k,v in by.items()}
    excl={k:float(np.mean([x for s,x in v.items() if s!="kansas"])) if any(s!="kansas" for s in v) else None for k,v in by.items()}
    z=dict(status="RETROSPECTIVE_DESCRIPTIVE_NOT_META",panels=seen,total_cells=dict(Counter(r["source"] for r in rows)),eligible=dict(Counter(r["source"] for r in valid)),matched_candidates=sum(r["n_candidates"] for r in valid),by_source=by,source_equal=avg,exclude_kansas=excl)
    assert valid, "STOP no complete paired scores"
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps({**z,"cells":rows},allow_nan=False)+"\n")
    print("AHIF_RESPONSE_ALLOCATION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
