"""Fixed partition, conditional binary label sharp-null diagnostic; not a memory test."""
import argparse,json
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_one_lag_fungal import source_fungal,contrast,evaluate,fixed_mass_baseline
from benchmarks.score_ahif_saine2025_succession import holm
N=999; SPLIT_SEED=2026112505; PERM_SEED=2026112506
def run(root):
    prop,t,ss,sites,coverage=source_fungal(root)
    label=t!="control"; n=len(label); site_names=sorted(set(sites))
    rng=np.random.default_rng(SPLIT_SEED);part={};audit=[]
    for stratum in sorted(set(ss)):
        ids=np.flatnonzero(ss==stratum)
        halves=np.array_split(rng.permutation(ids),2)
        assert not set(halves[0])&set(halves[1])
        for h,idx in enumerate(halves):
            k=int(label[idx].sum());nc=len(idx)-k
            audit.append(dict(stratum=stratum,half=h,treated=k,control=nc))
            if min(k,nc)==0:
                return dict(status="STOP_FIXED_SPLIT_ARM_ABSENT",coverage=audit)
            part[stratum,h]=(idx,k)
    rng=np.random.default_rng(PERM_SEED)
    # Each permutation keeps all half-specific arm counts; equal stratum weights per site.
    stats=np.zeros((N+1,2));site_stats={}
    for site in site_names:
        W=np.zeros((N+1,2,n))
        for (stratum,h),(ids,k) in part.items():
            if stratum.split("|")[0]!=site:continue
            obs=ids[label[ids]];ctl=ids[~label[ids]]
            W[0,h,obs]=.5/k;W[0,h,ctl]=-.5/len(ctl)
            for j in range(1,N+1):
                p=rng.permutation(ids);W[j,h,p[:k]]=.5/k;W[j,h,p[k:]]=-.5/(len(ids)-k)
        assert np.max(np.abs(W.sum(2)))<1e-10
        assert np.max(np.abs(np.abs(W).sum(2)-2))<1e-10
        assert not np.any((W[:,0]!=0)&(W[:,1]!=0))
        active=np.flatnonzero(sites==site)
        assert not np.any(W[:,:,sites!=site])
        D=contrast(W[:,:,active].reshape(-1,len(active)),prop[active]).reshape(N+1,2,3,prop.shape[-1])
        s=np.zeros((N+1,2))
        for j in range(N+1):
            for h in (0,1):
                tr=D[j,h];te=D[j,1-h]
                gain=evaluate(tr[1],te[2])["gain"]
                pre=fixed_mass_baseline(tr[1],tr[0],te[2])
                if gain is None or pre is None:return dict(status="STOP_UNDEFINED_FIXED_DENOMINATOR",coverage=audit,site=site,permutation=j)
                s[j]+=.5*np.array([gain,pre])
        stats+=s/len(site_names)
        site_stats[site]=s[0].tolist()
    p=(1+np.count_nonzero(stats[1:]>=stats[0]-1e-12,axis=0))/(N+1)
    adj=holm(p.tolist());metrics={}
    for j,key in enumerate(("one_lag_gain","lag_vs_pre_sign_fixed_mass")):
        metrics[key]=dict(observed=float(stats[0,j]),null_median=float(np.median(stats[1:,j])),
          null_percentiles=np.quantile(stats[1:,j],[.025,.975]).tolist(),conditional_p=float(p[j]),holm_p=float(adj[j]))
    return dict(status="CONDITIONAL_SHARP_NULL_NOT_MEMORY_TEST",n_logs=n,n_sites=len(site_names),n_permutations=N,
      split_seed=SPLIT_SEED,permutation_seed=PERM_SEED,coverage=audit,by_site_observed=site_stats,metrics=metrics,null=stats[1:].tolist())
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run(a.source);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    print("SAINE_ONE_LAG_SHARP_NULL_RESULT",json.dumps({k:v for k,v in z.items() if k!="null"},allow_nan=False),flush=True)
if __name__=="__main__":main()
