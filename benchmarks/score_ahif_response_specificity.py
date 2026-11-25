"""Fixed information-matched state/response and sign forecast, same four AHIF sources."""
import argparse,json
from collections import Counter
from pathlib import Path
import numpy as np
import benchmarks.score_ahif_cross_source_concentration as src
from benchmarks.score_ahif_response_allocation import norm,tv,nested,EPS
seen=[]
def mixed(a,b):
    return None if a is None or b is None else (a+b)/2
def signed(v):
    return np.concatenate((np.maximum(v,0),np.maximum(-v,0)),axis=-1)
def signblind(v):
    if v is None:return None
    pos,neg=np.split(v,2)
    m=(pos+neg)/2
    return np.concatenate((m,m))
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
        am=A.mean(1);bm=B.mean(1); delta=bm-am
        hold=b[:,list(hb)].mean(1)-a[:,list(ha)].mean(1)
        if paired:
            assert tuple(ia)==tuple(ib)
            noise=(B-A).var(1,ddof=1)/len(ia)
        else:
            noise=A.var(1,ddof=1)/len(ia)+B.var(1,ddof=1)/len(ib)
        for t in range(2,len(dates)):
            q=norm(np.abs(hold[t]));base=norm(np.abs(delta[:t-1]).mean(0))
            state=(am[t-1]+bm[t-1])/2
            C=dict(S=norm(state),SQ=norm(np.sqrt(state)),V=norm(np.sqrt(noise[t-1])),R=norm(np.abs(delta[t-1])))
            pred=dict(B=base,**C,**{"B"+k:mixed(base,v) for k,v in C.items()})
            scores={k:tv(q,v) for k,v in pred.items()}
            ys=norm(signed(hold[t]));bs=norm(signed(delta[:t-1]).mean(0));rs=norm(signed(delta[t-1]))
            ss=mixed(bs,rs)
            scores.update(signed=tv(ys,ss),signblind=tv(ys,signblind(ss)),signed_history=tv(ys,bs))
            gross=np.abs(hold[t]);den=gross.sum();sg=np.sign(delta[t-1])*np.sign(hold[t])
            scores.update(same=float(gross[sg>0].sum()/den) if den>EPS else None,
                         opposite=float(gross[sg<0].sum()/den) if den>EPS else None,
                         abstain=float(gross[sg==0].sum()/den) if den>EPS else None)
            rows.append(dict(source=st,panel=key,time=str(dates[t]),fold=fold,**scores))
    return rows
def summary(rows,keys,comparisons):
    valid=[r for r in rows if all(r[k] is not None for k in keys)]
    for r in valid:
        for name,k,j in comparisons:r[name]=r[k]-r[j]
    cols=keys+[n for n,_,_ in comparisons]
    by={k:nested(valid,k) for k in cols}
    avg={k:float(np.mean(list(v.values()))) if v else None for k,v in by.items()}
    excl={k:float(np.mean([x for s,x in v.items() if s!="kansas"])) if any(s!="kansas" for s in v) else None for k,v in by.items()}
    return dict(eligible=dict(Counter(r["source"] for r in valid)),by_source=by,source_equal=avg,exclude_kansas=excl)
def selftest():
    x=np.array([2.,-1.]);s=norm(signed(x));p=signblind(s)
    assert abs(s.sum()-1)<1e-12 and np.allclose(p,[1/3,1/6,1/3,1/6])
    assert abs(tv(s,p)-.5)<1e-12 and tv(s,s)==0
    assert np.allclose(mixed(np.array([0.,1.]),np.array([1.,0.])),[.5,.5])
    assert mixed(None,p) is None
    print("Specificity synthetic checks PASS")
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--poem",type=Path,required=True);p.add_argument("--field",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    selftest();a.root.mkdir(parents=True,exist_ok=True)
    src.panel=evaluate
    rows=src.get_sources(a.root,a.poem,a.field,Path("database/ahif_v1"),Path("validation/ahif_v1"))
    assert len(seen)==19 and len(set((x,y) for x,y,z in seen))==19
    keys=["B","S","SQ","V","R","BS","BSQ","BV","BR"]
    cs=[("R_over_"+k,k,"R") for k in ["S","SQ","V"]]+[("BR_over_"+k,k,"BR") for k in ["B","BS","BSQ","BV"]]
    main=summary(rows,keys,cs)
    nov=summary(rows,[k for k in keys if k not in ("V","BV")],[c for c in cs if c[1] not in ("V","BV")])
    signs=summary(rows,["signed","signblind","signed_history","same","opposite","abstain"],[("sign_gain","signblind","signed"),("recent_sign_gain","signed_history","signed")])
    z=dict(status="RETROSPECTIVE_DESCRIPTIVE_NOT_META",panels=seen,total_cells=dict(Counter(r["source"] for r in rows)),primary=main,without_variance=nov,signed=signs)
    assert len(main["eligible"])==4 and len(signs["eligible"])==4
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps({**z,"cells":rows},allow_nan=False)+"\n")
    print("AHIF_RESPONSE_SPECIFICITY_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
