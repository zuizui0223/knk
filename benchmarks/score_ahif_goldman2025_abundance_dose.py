"""Post-score diagnostic: unweighted source abundance versus nominal dose."""
import argparse,json
from pathlib import Path
import numpy as np
from sklearn.linear_model import LogisticRegression
from benchmarks import score_ahif_goldman2025_dose as base
COMBOS=base.COMBOS;PARENTS=base.PARENTS;RATIOS=base.RATIOS
def data_rows_extra(mats,pre,threshold):
    ks=sorted(set().union(*(set(v) for v in pre.values())))
    K=len(ks)
    pp={p:np.array([pre[p].get(k,0.) for k in ks]) for p in PARENTS}
    det={p:(v>0 if threshold==0 else v>=threshold) for p,v in pp.items()}
    regional=np.any(np.stack(list(det.values())),axis=0)
    result={}
    for c in COMBOS:
        x=[];y=[];w=[];fixed=[];dose_groups=[];a,b=c.split("-")
        for r in RATIOS:
            pairs=sorted((k,v) for k,v in mats.items() if k[0]==c and k[2]==r)
            assert pairs
            a0,b0=map(float,r.split(":"))
            for k,z in pairs:
                val=np.array([z.get(s,0.) for s in ks]);target=(val>0 if threshold==0 else val>=threshold)
                for focal,donor,q in ((a,b,b0/(a0+b0)),(b,a,a0/(a0+b0))):
                    risk=regional&~det[focal]
                    assert risk.any()
                    P=np.mean([det[p] for p in PARENTS if p not in (focal,donor)],axis=0)[risk]
                    D=det[donor][risk].astype(float)
                    L=np.full(len(P),np.log10(q/(1-q))/3)*D
                    A=np.log10(np.maximum(q*pp[donor][risk],1e-8))/8
                    A0=np.log10(np.maximum(pp[donor][risk],1e-8))/8
                    x.append(np.column_stack((P,D,L,A,A0)));y.append(target[risk].astype(int))
                    fixed.append(np.column_stack((P,(P+D)/2,(P+q*D)/2)))
                    w.append(np.full(len(P),1/(7*len(pairs)*2*len(P))))
                    dose_groups.append(np.full(len(P),r+"|"+focal))
        X=np.concatenate(x);Y=np.concatenate(y);W=np.concatenate(w);F=np.concatenate(fixed)
        assert abs(W.sum()-1)<1e-10
        result[c]={"X":X,"Y":Y,"W":W,"F":F}
    return result,K,int(regional.sum())

MODELS={"pool":[0],"identity":[0,1],"dose":[0,1,2],"dose_weighted_abundance":[0,1,2,3],"abundance_only":[0,1,4],"weighted_abundance_only":[0,1,3],"abundance_plus_dose":[0,1,4,2]}
REPLAY={"pool":"fit_pool","identity":"fit_identity","dose":"fit_dose","dose_weighted_abundance":"fit_dose_abundance"}
def evaluate(t,prior):
    out={}
    for c in COMBOS:
        r=prior["by_pair"][c];assert r["status"]=="SCORED"
        train=r["training_pairs"]
        assert all(not set(c.split("-"))&set(k.split("-")) for k in train)
        h=t[c];assert len(h["Y"])==r["n_risk_cells"]
        x=np.concatenate([t[k]["X"] for k in train]);y=np.concatenate([t[k]["Y"] for k in train]);w=np.concatenate([t[k]["W"]*100 for k in train])
        row={}
        for name,cols in MODELS.items():
            m=LogisticRegression(C=1,solver="lbfgs",max_iter=1000,tol=1e-9)
            m.fit(x[:,cols],y,sample_weight=w)
            assert m.n_iter_.max()<1000
            v=base.brier(h["Y"],m.predict_proba(h["X"][:,cols])[:,1],h["W"])
            if name in REPLAY:assert abs(v-r["fit"][REPLAY[name]])<1e-10,(c,name,v,r["fit"][REPLAY[name]])
            row[name]=v
        row.update(abundance_gain=row["identity"]-row["abundance_only"],nominal_dose_after_abundance=row["abundance_only"]-row["abundance_plus_dose"],weighting_at_equal_dimension=row["abundance_only"]-row["weighted_abundance_only"],weighted_vs_unweighted_with_dose=row["abundance_plus_dose"]-row["dose_weighted_abundance"])
        out[c]=row
    keys=list(next(iter(out.values())))
    return {"by_pair":out,"means":{k:float(np.mean([r[k] for r in out.values()])) for k in keys},"positive_counts":{k:sum(r[k]>0 for r in out.values()) for k in keys[-4:]}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    mats,lanes,audit=base.load()
    old=json.loads(Path("validation/ahif_v1/goldman2025_dose_summary_20261126.json").read_text())
    assert len(mats)==old["audit"]["observed_p5"]==167 and audit==old["audit"]
    result={}
    for lane,pre in lanes.items():
        assert old["lanes"][lane]["status"]=="SCORED_CONDITIONAL_SOURCE_PROXY"
        result[lane]={}
        for label,th in (("threshold_0.001",.001),("raw_detection",0)):
            table,_,_=data_rows_extra(mats,pre,th);orig,_,_=base.data_rows(mats,pre,th)
            for c in COMBOS:
                for k in ("Y","W","F"):assert np.array_equal(table[c][k],orig[c][k])
                assert np.array_equal(table[c]["X"][:,:4],orig[c]["X"])
            result[lane][label]=evaluate(table,old["lanes"][lane]["results"][label]["source_disjoint"])
    z={"status":"POST_SCORE_DIAGNOSTIC_NOT_CAUSAL","results":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    print("GOLDMAN_ABUNDANCE_DOSE_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
