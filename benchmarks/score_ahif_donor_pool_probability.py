"""Fixed source-pool probability versus actual donor binary forecast."""
import argparse,json
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_villarreal2025_external as vv
from benchmarks import score_ahif_donor_mismatch as dm
SEED=2026112616;N=999
KEYS=("pool_loss","binary_loss","hybrid_loss","binary_gain","hybrid_gain","random_hard_loss")
def run(mats,prior):
    rng=np.random.default_rng(SEED);out={}
    for donor,rec in dm.GROUPS:
        src=mats[donor];dst=mats[rec]
        z=np.zeros((10,10,4,len(KEYS)))
        for i in range(10):
            for t in range(4):
                risk=dst[t,i]==0
                assert risk.any()
                y=(dst[t+1,i,risk]>0).astype(float)
                d=(src[t][:,risk]>0).astype(float)
                for j in range(10):
                    other=np.delete(d,j,axis=0)
                    p=other.mean(0);binary=d[j];h=(p+binary)/2
                    a=np.mean((p-y)**2);b=np.mean((binary-y)**2);c=np.mean((h-y)**2)
                    q=np.mean((other-y)**2)
                    z[i,j,t]=a,b,c,a-b,a-c,q
        obs=z[np.arange(10),np.arange(10)].mean(0)
        old=prior["detection"][rec]["hard_donor_brier_skill"]
        assert np.max(np.abs(obs[:,1]-old["actual"]))<1e-12
        assert np.max(np.abs(obs[:,5]-old["other_nine"]))<1e-12
        ns=np.array([z[np.arange(10),rng.permutation(10)].mean(0)[:,4] for _ in range(N)])
        post=ns[:,1:].mean(1);net=post-ns[:,0]
        ref=lambda v,x:{"observed":float(x),"reference_median":float(np.median(v)),"reference_range":np.quantile(v,[.025,.975]).tolist(),"upper_reference_tail":float((1+np.count_nonzero(v>=x-1e-12))/(N+1))}
        out[rec]={"intervals":{vv.DATES[t]+"-"+vv.DATES[t+1]:dict(zip(KEYS,obs[t].tolist())) for t in range(4)},
          "post":dict(zip(KEYS,obs[1:].mean(0).tolist())),
          "post_minus_placebo":dict(zip(KEYS,(obs[1:].mean(0)-obs[0]).tolist())),
          "hybrid_trajectory_reference":{"post":ref(post,obs[1:,4].mean()),"post_minus_placebo":ref(net,obs[1:,4].mean()-obs[0,4])},
          "pair_metrics":z[np.arange(10),np.arange(10)].tolist()}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    refs=vv.acquire(a.root);mats,audit=vv.load(a.root,refs)
    old=json.loads(Path("validation/ahif_v1/villarreal2025_donor_mismatch_summary_20261126.json").read_text())
    z={"status":"EXPLORATORY_FIXED_POOL_BRIER_NOT_CAUSAL","n_pairs":20,"seed":SEED,"n_permutations":N,"results":run(mats,old)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    brief={**z,"results":{k:{j:v for j,v in r.items() if j!="pair_metrics"} for k,r in z["results"].items()}}
    print("AHIF_DONOR_POOL_RESULT",json.dumps(brief,allow_nan=False),flush=True)
if __name__=="__main__":main()
