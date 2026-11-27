"""Hold two training cultures fixed and aggregate one versus two held cultures."""
import argparse,json
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_shen2023_context as s
from benchmarks import score_ahif_shen2023_baseline as bl
VAL=Path("validation/ahif_v1/shen2023_context_summary_20261127.json")
MOD=("no_change","additive","multiplicative")
def run(x,ref):
    rng=np.random.default_rng(s.SEED)
    orders=np.stack([rng.permutation(4) for _ in range(s.N*6)]).reshape(s.N,6,4)
    out={1:[],2:[]};variance=[];n=0;maxerr=0.
    for order in orders:
        for orientation in (0,1):
            tr=order[:,:2] if orientation==0 else order[:,2:]
            te=order[:,2:] if orientation==0 else order[:,:2]
            assert all(not set(a)&set(b) for a,b in zip(tr,te))
            tm=np.stack([x[:,a,ids].mean((1,2)) for a,ids in enumerate(tr)],axis=1)
            hm=np.stack([x[:,a,ids].mean(2) for a,ids in enumerate(te)],axis=1)
            if (tm[:,0]<=0).any():continue
            n+=1;v1=[];v2=[];vv=[]
            for i in range(2):
                a=tm[i,0]
                for b,d in ((1,4),(2,5)):
                    for c,background in ((b,3),(3,b)):
                        train_t=tm[i,c];h=hm[i,background];target=hm[i,d]
                        pp=np.stack([h,np.maximum(0.,h+train_t-a),h*train_t/a],axis=1)/a
                        p2=np.array([h.mean(),max(0.,h.mean()+train_t-a),h.mean()*train_t/a])/a
                        y=target/a
                        e1=pp[:,None,:]-y[None,:,None]
                        e2=p2-y.mean()
                        ma1=np.abs(e1).mean((0,1));ms1=(e1**2).mean((0,1))
                        ma2=np.abs(e2);ms2=e2**2
                        var=pp.var(0)+y.var()
                        # Linear no-change and multiplication exactly commute with averaging.
                        check=np.max(np.abs((ms1-ms2-var)[[0,2]]))
                        assert check<1e-10;maxerr=max(maxerr,float(check))
                        v1.append(np.stack([ma1,ms1]));v2.append(np.stack([ma2,ms2]));vv.append(var)
            out[1].append(v1);out[2].append(v2);variance.append(vv)
    assert n==ref["n_common_folds"] and len(out[1])==len(out[2])
    summary={}
    for h in (1,2):
        m=np.mean(out[h],axis=0)
        a=m.mean(0)
        summary[str(h)]={"MAE":dict(zip(MOD,a[0].tolist())),"MSE":dict(zip(MOD,a[1].tolist())),
          "MSE_gain_vs_zero":dict(zip(MOD,(a[1,0]-a[1]).tolist())),
          "by_panel":[{"focal":f,"competitor":c,"added_agent":j,"MAE":dict(zip(MOD,v[0].tolist())),"MSE":dict(zip(MOD,v[1].tolist()))} for (f,c,j),v in zip([(f,c,j) for f in s.FOCAL for c in ("P","S") for j in ("competitor","predator")],m)]}
    bl.eq(summary["2"]["MAE"],ref["mean_MAE"])
    for a,b in zip(summary["2"]["by_panel"],ref["by_panel"]):bl.eq(a["MAE"],b["MAE"])
    return {"folds":n,"panel_folds":n*8,"h":summary,"mean_prediction_plus_target_variance":dict(zip(MOD,np.mean(variance,axis=(0,1)).tolist())),"max_linear_squared_loss_identity_error":maxerr}
def main():
    p=argparse.ArgumentParser();p.add_argument("--csv",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    mat,common,strict,audit=s.load(a.csv);prior=json.loads(VAL.read_text());bl.eq(audit,prior["audit"])
    result={}
    for lane,ds in (("primary_common",common),("strict_complete_focals",strict),("after_second_predator_addition",[d for d in common if d>=9])):
        old=prior["results"][lane]
        bl.eq(ds,old["dates"]);bl.eq(s.forecast(mat(ds)),old["forecast"])
        result[lane]=run(mat(ds),old["forecast"])
    z={"status":"FIXED_TRAIN_HELD_CULTURE_AGGREGATION_NOT_META","prior_replayed":True,"train_per_arm":2,"results":result}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SHEN_FIXED_TARGET_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
