"""Exploratory mismatch and paired donor detection diagnostics, source-locked."""
import argparse,json
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
from benchmarks import score_ahif_villarreal2025_external as vv
from benchmarks import score_ahif_mass_direction_transfer as md
from benchmarks import score_ahif_replication_scale as rs
SEED=2026112614
PERM_SEED=2026112615
N=256;PERMS=999
GROUPS=(("A0","Aa"),("B0","Ab"))
def base_features(d,r):
    dp=d/d.sum();rp=r/r.sum()
    return {"TV":float(np.abs(dp-rp).sum()/2),
            "Jaccard":float(1-np.count_nonzero((d>0)&(r>0))/np.count_nonzero((d>0)|(r>0))),
            "novelty_mass":float(dp[r==0].sum()),"novelty_taxa":int(((d>0)&(r==0)).sum()),
            "donor_total":int(d.sum()),"recipient_total":int(r.sum()),
            "donor_richness":int((d>0).sum()),"recipient_richness":int((r>0).sum())}
def compare(a,b,h):
    v,ok=md.factorial(a,b,h)
    return dict(zip(md.KEYS,v.mean(axis=0).tolist())) if ok.all() else None
def replay_attribution(mats,old):
    rng=np.random.default_rng(vv.SEED)
    orders={k:np.array([rng.permutation(10) for _ in range(vv.N)]) for k in ("A","B")}
    out={}
    for donor,recipient in GROUPS:
        p=("A0",recipient);g="A" if recipient=="Aa" else "B"
        d,h=vv.split(mats,orders,2,p)
        v,ok=md.factorial(d[:,:,2:3],d[:,:,3:4],h[:,:,4:5])
        d4,h4=vv.split(mats,orders,4,p)
        v4,ok4=md.factorial(d4[:,:,2:3],d4[:,:,3:4],h4[:,:,4:5])
        w,wok=rs.scores(d[:,:,3:4],d[:,:,4:5],h[:,:,4:5])
        w4,wok4=rs.scores(d4[:,:,3:4],d4[:,:,4:5],h4[:,:,4:5])
        common=ok&ok4&wok&wok4
        assert common.all() and int(common.sum())==512
        ref=old["by_panel"][recipient+"-A0"]["2"]["common"]["factorial"]
        actual=dict(zip(md.KEYS,v[common].mean(0).tolist()))
        assert all(abs(actual[k]-ref[k])<1e-10 for k in md.KEYS)
        held=np.stack((orders[g][:,-2:],orders[g][:,:2]),axis=1)
        cells=v[:,:,0,:]
        for i in range(10):
            mask=(held==i).any(-1)&common[:,:,0]
            assert mask.any()
            out[recipient,i]={"n_held_folds":int(mask.sum()),"mean":dict(zip(md.KEYS,cells[mask].mean(0).tolist()))}
    return out
def singleton(mats):
    rng=np.random.default_rng(SEED);out={}
    for donor,recipient in GROUPS:
        for i in range(10):
            r=[]
            for j in range(N):
                rem=np.delete(np.arange(10),i)
                trb=rng.choice(rem,2,replace=False)
                if recipient=="Aa":
                    ha=i;tra=trb
                else:
                    z=rng.permutation(10);ha=z[0];tra=z[1:3]
                assert i not in trb and ha not in tra
                d=mats[recipient][:,trb].mean(1)-mats["A0"][:,tra].mean(1)
                h=mats[recipient][4,i]-mats["A0"][4,ha]
                v,ok=md.factorial(d[2],d[3],h)
                if ok:r.append(v)
            assert r
            out[recipient,i]={"n_defined":len(r),"n_total":N,"mean":dict(zip(md.KEYS,np.mean(r,axis=0).tolist()))}
    return out
def corr(x,y,rank=True):
    x=np.asarray(x);y=np.asarray(y)
    if len(x)<3 or np.std(x)<1e-15 or np.std(y)<1e-15:return None
    return float(spearmanr(x,y).statistic if rank else np.corrcoef(x,y)[0,1])
def correlations(rows):
    out={}
    for lane in ("held_group","singleton"):
        out[lane]={}
        for predictor in ("TV","novelty_mass"):
            for outcome in ("mass_symmetric","mass_only_gain","sign_symmetric"):
                x=np.array([r["T2"][predictor] for r in rows])
                y=np.array([r[lane]["mean"][outcome] for r in rows])
                g=np.array([r["recipient_arm"] for r in rows])
                by={}
                xc=x.copy();yc=y.copy()
                for arm in ("Aa","Ab"):
                    k=g==arm;a=x[k];b=y[k];n=len(a);assert n==10
                    loo=[corr(np.delete(a,i),np.delete(b,i)) for i in range(n)]
                    defined=[z for z in loo if z is not None]
                    by[arm]={"spearman":corr(a,b),"delete_one_range":[min(defined),max(defined)] if defined else None,"n_pairs":n}
                    xc[k]-=a.mean();yc[k]-=b.mean()
                out[lane][predictor+"__"+outcome]={"within_group":by,"pooled_spearman":corr(x,y),"group_centered_pearson":corr(xc,yc,False)}
    return out
def mismatch(mats,old):
    att=replay_attribution(mats,old);sing=singleton(mats)
    rows=[]
    for donor,rec in GROUPS:
        for i in range(10):
            rows.append({"recipient_arm":rec,"donor_tank":vv.EXPECTED[donor][i],"recipient_tank":vv.EXPECTED[rec][i],
                         "T1":base_features(mats[donor][0,i],mats[rec][0,i]),
                         "T2":base_features(mats[donor][1,i],mats[rec][1,i]),
                         "held_group":att[rec,i],"singleton":sing[rec,i]})
    groups={}
    for rec in ("Aa","Ab"):
        rr=[r for r in rows if r["recipient_arm"]==rec]
        groups[rec]={t:{k:float(np.mean([r[t][k] for r in rr])) for k in rr[0][t]} for t in ("T1","T2")}
        groups[rec]["T2_ranges"]={k:[min(r["T2"][k] for r in rr),max(r["T2"][k] for r in rr)] for k in rr[0]["T2"]}
    return rows,groups,correlations(rows)
def detection(mats):
    # Axis: recipient pair, candidate donor, time interval.
    result={};rng=np.random.default_rng(PERM_SEED)
    for donor,rec in GROUPS:
        src=mats[donor];dst=mats[rec]
        enr=np.full((10,10,4),np.nan);loss=np.zeros((10,10,4));cover=np.full((10,10,4),np.nan)
        events=[]
        for i in range(10):
            for t in range(4):
                risk=dst[t,i]==0
                y=dst[t+1,i,risk]>0
                if len(y)==0:raise ValueError("empty risk set")
                for j in range(10):
                    d=src[t,j,risk]>0
                    loss[i,j,t]=np.mean(d!=y)
                    if y.any():cover[i,j,t]=d[y].mean()
                    if y.any() and (~y).any():enr[i,j,t]=d[y].mean()-d[~y].mean()
                events.append({"recipient_tank":vv.EXPECTED[rec][i],"interval":vv.DATES[t]+"-"+vv.DATES[t+1],
                               "n_risk":int(risk.sum()),"n_new":int(y.sum()),
                               "actual_coverage":None if not y.any() else float(cover[i,i,t]),
                               "actual_enrichment":None if not (y.any() and (~y).any()) else float(enr[i,i,t]),
                               "other_mean_enrichment":None if not (y.any() and (~y).any()) else float(np.delete(enr[i,:,t],i).mean())})
        elig=np.isfinite(enr[:,0,:])
        assert np.array_equal(np.isfinite(enr),np.broadcast_to(elig[:,None,:],enr.shape))
        def summarize(mat,perm):
            a=mat[np.arange(10),perm,:]
            means=[float(np.nanmean(a[:,t])) if np.isfinite(a[:,t]).any() else None for t in range(4)]
            if any(x is None for x in means):raise ValueError("undefined whole interval")
            return means
        def reference(mat,lower=False):
            actual=summarize(mat,np.arange(10))
            other=[]
            for t in range(4):
                v=[float(np.delete(mat[i,:,t],i).mean()) for i in range(10) if np.isfinite(mat[i,0,t])]
                assert v;other.append(float(np.mean(v)))
            adv=np.array(other)-actual if lower else np.array(actual)-other
            return {"actual":actual,"other_nine":other,"actual_advantage":adv.tolist(),"post_advantage":float(adv[1:].mean()),"post_minus_placebo_advantage":float(adv[1:].mean()-adv[0])}
        permutations=np.array([rng.permutation(10) for _ in range(PERMS)])
        observed=np.array(summarize(enr,np.arange(10)))
        null=np.array([summarize(enr,p) for p in permutations])
        post=null[:,1:].mean(1);diff=post-null[:,0]
        observed_post=observed[1:].mean();observed_diff=observed_post-observed[0]
        ref=lambda v,x:{"observed":float(x),"reference_median":float(np.median(v)),"reference_range":np.quantile(v,[.025,.975]).tolist(),
                         "upper_reference_tail":float((1+(v>=x-1e-12).sum())/(len(v)+1))}
        result[rec]={"n_pairs":10,"intervals":[vv.DATES[t]+"-"+vv.DATES[t+1] for t in range(4)],
                     "enrichment_defined_pairs":elig.sum(0).tolist(),
                     "newly_detected_counts":[sum(r["n_new"] for r in events if r["interval"]==vv.DATES[t]+"-"+vv.DATES[t+1]) for t in range(4)],
                     "enrichment":reference(enr),"hard_donor_brier_skill":reference(loss,True),"coverage":reference(cover),
                     "whole_donor_reference":{"n_permutations":PERMS,"post":ref(post,observed_post),"post_minus_placebo":ref(diff,observed_diff)},
                     "events":events}
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    refs=vv.acquire(a.root);mats,audit=vv.load(a.root,refs)
    old=json.loads(Path("validation/ahif_v1/villarreal2025_external_summary_20261126.json").read_text())
    rows,groups,corrs=mismatch(mats,old);det=detection(mats)
    z={"status":"EXPLORATORY_WITHIN_STUDY_DONOR_DIAGNOSTIC_NOT_META","source_commit_parent":"01075120396027d571517028b4cc2b25eb7c5577",
       "source_audit":audit,"group_summary":groups,"correlations":corrs,"pair_rows":rows,"detection":det,
       "limitations":["Dependent per-pair score attribution, not independent effects; singleton shares training and references",
                      "Donor identity reference is not exact randomization inference; no ancestry or establishment evidence"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    brief={k:v for k,v in z.items() if k not in ("pair_rows","source_audit","detection")}
    brief["pair_rows"]=[{k:r[k] for k in ("recipient_arm","donor_tank","recipient_tank","T2")}|{"held_group_mass":r["held_group"]["mean"]["mass_symmetric"],"singleton_mass":r["singleton"]["mean"]["mass_symmetric"]} for r in rows]
    brief["detection"]={k:{j:v for j,v in r.items() if j!="events"} for k,r in det.items()}
    print("AHIF_DONOR_MISMATCH_RESULT",json.dumps(brief,allow_nan=False),flush=True)
if __name__=="__main__":main()
