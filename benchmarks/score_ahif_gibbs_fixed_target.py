"""Fixed train information, varying held replicate aggregation. No spatial grain claim."""
import argparse,itertools,json
from pathlib import Path
from collections import defaultdict
import numpy as np
from benchmarks import score_ahif_gibbs2026_spatial as b
from benchmarks import score_ahif_gibbs2026_group_mean as g
VAL=Path("validation/ahif_v1")
MOD=b.MODELS
def fold_predictions(rows,folds,adjusted):
    foc=sorted(set(r["focal"] for r in rows));pairs=["-".join(x) for x in b.PAIRS]
    groups=[(f,p) for f in foc for p in pairs];out=defaultdict(list)
    X=np.array([[1,*b.features(p)] for p in pairs]);pen=np.diag([0,1,1,1,1,1])
    for tr,te in folds:
        train=[r for r in rows if r["rep"] in tr];held=[r for r in rows if r["rep"] in te]
        T=defaultdict(list);H=defaultdict(list)
        for r in train:T[b.group(r)].append(r)
        for r in held:H[b.group(r)].append(r)
        assert all(len(T[k])>=2 and len(H[k])>=1 for k in groups)
        for j,p in enumerate(pairs):
            other=[r for r in train if r["pair"]!=p]
            sl=b.residual_slopes(other) if adjusted else np.zeros(2)
            def D(r):return r["D"]-float(sl@np.array([r["delta_total"],r["delta_composition"]]))
            def gm(rr):return float(np.mean([D(r) for r in rr]))
            for f in foc:
                prior=np.array([gm(T[f,q]) for q in pairs if q!=p])
                design=np.delete(X,j,axis=0)
                beta=np.linalg.solve(design.T@design+pen,design.T@prior)
                pred=np.array([0.,prior.mean(),gm(T[f,p]),np.mean([gm(T[k,p]) for k in foc if k!=f]),X[j]@beta])
                target=np.array([D(r) for r in H[f,p]])
                out[f,p].append((pred,target))
    assert set(out)==set(groups) and all(len(v)==len(folds) for v in out.values())
    return out
def scores(group,min_n):
    records={h:[] for h in range(1,min_n+1)};support={};maxdev=0.
    # Equal groups, then equal eligible folds in each group.
    for key,cells in group.items():
        acc={h:[] for h in records}
        kept=0
        for p,y in cells:
            k=len(y)
            if k<min_n:continue
            kept+=1;mu=y.mean();var=y.var()
            full=(p-mu)**2
            for h in records:
                sub=np.array(list(itertools.combinations(range(k),h)))
                means=y[sub].mean(1)
                err=p[None,:]-means[:,None]
                mae=np.abs(err).mean(0);mse=(err**2).mean(0)
                variance=(k-h)/(h*(k-1))*var
                check=np.max(np.abs(mse-full-variance))
                maxdev=max(maxdev,float(check));assert check<1e-10
                acc[h].append(np.stack((mae,mse,full,np.full(len(p),variance))))
        assert kept>0,(key,min_n)
        support["|".join(key)]=kept
        for h in records:records[h].append(np.mean(acc[h],axis=0))
    out={}
    for h,vals in records.items():
        ma,ms,full,var=np.mean(vals,axis=0)
        out[str(h)]={"MAE":dict(zip(MOD,ma.tolist())),"MSE":dict(zip(MOD,ms.tolist())),"MSE_full_target":dict(zip(MOD,full.tolist())),"common_target_variance":float(var[0]),"MSE_gain_vs_zero":dict(zip(MOD,(ms[0]-ms).tolist())),"MSE_gain_vs_anonymous":dict(zip(MOD,(ms[1]-ms).tolist()))}
    for m in MOD:
        v=[q["MSE_gain_vs_zero"][m] for q in out.values()]
        assert max(v)-min(v)<1e-10
    return {"h":out,"eligible_fold_groups":sum(support.values()),"possible_fold_groups":sum(len(x) for x in group.values()),"support_by_group":support,"max_squared_loss_identity_error":maxdev}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    rows,audit=b.load();old=json.loads((VAL/"gibbs2026_spatial_summary_20261127.json").read_text())
    g.check(audit,old["audit"]);g.check(b.describe(rows),old["descriptive"])
    for lane,adj in (("raw",False),("density_standardized_sensitivity",True)):
        g.check(b.summarize(b.held_predictions(rows,adj)),old["forecasts"][lane])
    converted,unit=g.units(rows,audit);rng=np.random.default_rng(g.SEED);folds=[]
    for _ in range(g.N):
        z=rng.permutation(b.REPS);folds.append((set(z[:10]),set(z[10:])))
    grouped=json.loads((VAL/"gibbs2026_group_mean_summary_20261127.json").read_text())
    result={}
    for name,data in (("author_raw",rows),("hypothetical_SA_x10",converted),("AC_UR_only",[r for r in rows if r["focal"]!="SA"])):
        result[name]={}
        for mode,adj in (("raw",False),("density_standardized",True)):
            q=fold_predictions(data,folds,adj)
            replay=np.mean([np.mean([np.abs(p-y.mean()) for p,y in v],axis=0) for v in q.values()],axis=0)
            ref=grouped["results"][name][mode]["models"]
            assert all(abs(replay[i]-ref[m]["MAE"])<1e-10 for i,m in enumerate(MOD))
            result[name][mode]={"exact_five_support":scores(q,5),"at_least_three_support":scores(q,3)}
    z={"status":"FIXED_TRAIN_TARGET_AGGREGATION_NOT_GEOGRAPHIC_SCALE_NOT_META","prior_replayed":True,"train_labels":10,"held_label_pool":5,"n_folds":g.N,"results":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("GIBBS_FIXED_TARGET_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
