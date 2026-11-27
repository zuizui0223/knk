"""Fixed grouped-target sensitivity for already-scored Gibbs experiment."""
import argparse,json,math,csv,io
from pathlib import Path
from collections import defaultdict
import numpy as np
from benchmarks import score_ahif_gibbs2026_spatial as b
SEED=2026112718;N=256
VAL=Path("validation/ahif_v1/gibbs2026_spatial_summary_20261127.json")
def check(a,z,path=""):
    if isinstance(a,dict):
        assert a.keys()==z.keys(),path
        for k in a:check(a[k],z[k],path+"/"+str(k))
    elif isinstance(a,list):
        assert len(a)==len(z),path
        for i,(x,y) in enumerate(zip(a,z)):check(x,y,path+"/"+str(i))
    elif isinstance(a,(int,float)) and not isinstance(a,bool):
        assert abs(a-z)<1e-10,(path,a,z)
    else:assert a==z,(path,a,z)
def units(rows,audit):
    raw=b.get("fecundity_data.csv")
    cal=b.get("seed_head_calibration.csv")
    cr=list(csv.reader(io.StringIO(cal.decode("utf-8-sig"))))[1:]
    cd=np.array([float(r[2]) for r in cr[1:-1]])
    by={};ds=[]
    nonbare={r["split_plot"] for r in rows}|{r["mix_plot"] for r in rows}
    for r in b.dicts(raw):
        k=int(r["Plot"])
        if k not in nonbare:continue
        h=[b.val(r[x]) for x in ("Salvia","2 Salvia","3 Salvia","4 Salvia","5 Salvia","6 Salvia")]
        h=[x for x in h if x is not None]
        by[k]=h;ds.extend(h)
    intercept=audit["calibration"]["SA_intercept"];slope=audit["calibration"]["SA_slope"]
    out=[]
    for r in rows:
        q=r.copy()
        if r["focal"]=="SA":
            for t in ("split","mix"):
                h=by[r[t+"_plot"]]
                assert h
                check(sum(intercept+slope*x*x for x in h),r[t+"_seeds"])
                q[t+"_seeds"]=sum(intercept+slope*(10*x)**2 for x in h)
            q["D"]=math.log1p(q["split_seeds"])-math.log1p(q["mix_seeds"])
        out.append(q)
    return out,{"calibration_header":"Diameter (mm)","nonbare_header_unit":"not stated",
                "calibration_diameter_quantiles_0_25_50_75_100":np.quantile(cd,[0,.25,.5,.75,1]).tolist(),
                "nonbare_diameter_quantiles_0_25_50_75_100":np.quantile(ds,[0,.25,.5,.75,1]).tolist(),
                "nonbare_heads_in_matched_plot_set":len(ds),"status":"HYPOTHETICAL_X10_NOT_VERIFIED_CORRECTION"}
def run(rows,folds,adjusted):
    foc=sorted(set(r["focal"] for r in rows));pairs=["-".join(x) for x in b.PAIRS]
    groups=[(f,p) for f in foc for p in pairs]
    models=b.MODELS;err=[];sgn=[];positive=[];nheld=[];skipped=0
    X=np.array([[1,*b.features(p)] for p in pairs]);pen=np.diag([0,1,1,1,1,1])
    for tr,te in folds:
        train=[r for r in rows if r["rep"] in tr];held=[r for r in rows if r["rep"] in te]
        T=defaultdict(list);H=defaultdict(list)
        for r in train:T[b.group(r)].append(r)
        for r in held:H[b.group(r)].append(r)
        if any(len(T[g])<2 or len(H[g])<1 for g in groups):
            skipped+=1;continue
        es=[];ss=[];ps=[];ns=[]
        for j,p in enumerate(pairs):
            other=[r for r in train if r["pair"]!=p]
            sl=b.residual_slopes(other) if adjusted else np.zeros(2)
            def D(r):return r["D"]-float(sl@np.array([r["delta_total"],r["delta_composition"]]))
            def gm(rr):return float(np.mean([D(r) for r in rr]))
            for f in foc:
                prior=np.array([gm(T[f,q]) for q in pairs if q!=p])
                design=np.delete(X,j,axis=0)
                beta=np.linalg.solve(design.T@design+pen,design.T@prior)
                target=gm(H[f,p]);known=gm(T[f,p])
                pred=np.array([0.,prior.mean(),known,np.mean([gm(T[k,p]) for k in foc if k!=f]),X[j]@beta])
                es.append(abs(target-pred))
                ss.append(np.where((target==0)|(pred==0),0,np.where(np.sign(target)==np.sign(pred),1,-1)))
                ps.append(target>0);ns.append(len(H[f,p]))
        err.append(es);sgn.append(ss);positive.extend(ps);nheld.extend(ns)
    assert err
    e=np.array(err);sg=np.array(sgn)
    by={}
    for j,m in enumerate(models):
        a=sg[:,:,j]
        by[m]={"MAE":float(e[:,:,j].mean()),"same":int((a==1).sum()),"opposite":int((a==-1).sum()),"unresolved":int((a==0).sum())}
    return {"folds":len(err),"skipped_folds":skipped,"groups":len(groups),"dependent_group_folds":len(err)*len(groups),
            "held_group_n_range":[min(nheld),max(nheld)],"positive_held_mean":int(sum(positive)),
            "models":by}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    rows,audit=b.load();prior=json.loads(VAL.read_text());check(audit,prior["audit"])
    check(b.describe(rows),prior["descriptive"])
    for lane,adj in (("raw",False),("density_standardized_sensitivity",True)):
        check(b.summarize(b.held_predictions(rows,adj)),prior["forecasts"][lane])
    converted,unit=units(rows,audit)
    rng=np.random.default_rng(SEED);folds=[]
    for _ in range(N):
        z=rng.permutation(b.REPS);folds.append((set(z[:10]),set(z[10:])))
    result={}
    for lane,v in (("author_raw",rows),("hypothetical_SA_x10",converted),("AC_UR_only",[r for r in rows if r["focal"]!="SA"])):
        result[lane]={"descriptive":b.describe(v),"raw":run(v,folds,False),"density_standardized":run(v,folds,True)}
    z={"status":"POST_SCORE_GROUP_MEAN_SENSITIVITY_NOT_META","prior_replayed":True,"seed":SEED,"unit_audit":unit,"results":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("GIBBS_GROUP_MEAN_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
