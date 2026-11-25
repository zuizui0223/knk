"""Held-well directional forecast on author 3000-read complete trajectories.
Predictive diagnostic, not founder-density causal inference."""
import argparse,csv,json
from collections import Counter
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_one_lag_fungal import evaluate
from benchmarks.score_ahif_signal_transfer import features
DAYS=(2,4,6,8); LEVELS=(1.,.1,.01,.001); SOURCES=("soil","water")
N=256; SEED=2026112507
META=("Sample_name","Day","Location","Source","Inoculum_density","Plate_position","unit")
def read(root,name):
    with (root/(name+".csv")).open() as f:
        r=csv.DictReader(f);rows=list(r);cols=r.fieldnames
    assert tuple(cols[:7])==META and len(set(cols))==len(cols)
    taxa=cols[7:];by={};units={};a={}
    for r in rows:
        s=r["Source"];d=float(r["Inoculum_density"]);day=int(r["Day"]);u=int(r["unit"])
        assert s in SOURCES and d in LEVELS and day in DAYS and u==int(r["Sample_name"].split("_")[-1])%384
        assert int(r["Plate_position"])==LEVELS.index(d)+1
        key=s,u
        label=(d,r["Location"],r["Plate_position"])
        if key in units:assert units[key]==label
        units[key]=label
        z=np.array([float(r[k]) for k in taxa])
        assert (z>=0).all() and np.isfinite(z).all() and (z==np.floor(z)).all() and z.sum()==3000
        assert (s,u,day) not in by
        by[s,u,day]=z/3000
    assert len(by)==len(rows)==len(units)*4
    ids={}
    for s in SOURCES:
        for d in LEVELS:
            us=sorted(u for (sr,u),v in units.items() if sr==s and v[0]==d)
            assert 4<=len(us)<=96
            ids[s,d]=us
            a[s,d]=np.array([[by[s,u,t] for u in us] for t in DAYS])
    return a,ids,taxa
def groups(a,ids,small):
    rng=np.random.default_rng(SEED);out=[]
    for s in SOURCES:
        for j in range(N):
            split={}
            for d in LEVELS:
                z=rng.permutation(len(ids[s,d]))
                if small:tr,te=z[:2],z[2:4]
                else:tr,te=np.array_split(z,2)
                assert len(tr)>=2 and len(te)>=2 and not set(tr)&set(te)
                split[d]=(tr,te)
            for o in (0,1):
                for d in LEVELS[1:]:
                    ia,ha=split[1.][o],split[1.][1-o];ib,hb=split[d][o],split[d][1-o]
                    A=a[s,1.];B=a[s,d]
                    td=B[:,ib].mean(1)-A[:,ia].mean(1)
                    hd=B[:,hb].mean(1)-A[:,ha].mean(1)
                    for t in range(1,4):
                        z=evaluate(td[t-1],hd[t]);c=evaluate(td[t],hd[t]);f=features(A[t-1,ia],B[t-1,ib],False)
                        q=float(td[t-1]@td[t-1])
                        w=float(A[t-1,ia].var(0,ddof=1).sum()/len(ia)+B[t-1,ib].var(0,ddof=1).sum()/len(ib))
                        out.append(dict(source=s,density=d,split=j,orientation=o,day=DAYS[t],primary=t>=2,
                            n_train_control=len(ia),n_train_treated=len(ib),n_held_control=len(ha),n_held_treated=len(hb),
                            **z,concurrent_gain=c["gain"],lag_minus_concurrent=z["gain"]-c["gain"] if z["gain"] is not None and c["gain"] is not None else None,
                            U=f[6] if f is not None else None,Q=q,W=w,train_gross=float(np.abs(td[t-1]).sum())))
    return out
KEYS=("gain","signed","blind","same","opposite","abstain","allocation_tv","concurrent_gain","lag_minus_concurrent","U","Q","W","train_gross")
def summarize(rows):
    out={}
    for name,subset in (("primary",[r for r in rows if r["primary"]]),("all_adjacent",rows)):
        b={}
        for s in SOURCES:
            for d in LEVELS[1:]:
                rr=[r for r in subset if r["source"]==s and r["density"]==d]
                zz=[r for r in rr if all(r[k] is not None for k in KEYS)]
                # Common denominator, no omitted treatment panels.
                assert rr and zz,(s,d,name)
                b[f"{s}:{d}"]={"n":len(zz),"total":len(rr),**{k:float(np.mean([r[k] for r in zz])) for k in KEYS}}
        out[name]={"by_panel":b,"equal_panel":{k:float(np.mean([z[k] for z in b.values()])) for k in KEYS},
                   "n_defined":sum(z["n"] for z in b.values()),"n_total":sum(z["total"] for z in b.values())}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    otu,ids,taxa=read(a.source,"otu99");asv,ids2,taxa2=read(a.source,"asv")
    assert ids==ids2,"ASV/OTU complete culture eligibility mismatch"
    coverage={f"{s}:{d}":len(v) for (s,d),v in ids.items()}
    r={"status":"DESCRIPTIVE_HELD_WELL_FORECAST_POSITION_ALIASED_NOT_META","source_commit":"a37af8c186a9438edccf5f729a0fdd48eed4469f","n_taxa_otu99":len(taxa),"n_taxa_asv":len(taxa2),"n_complete_wells":sum(coverage.values()),"n_complete_profiles":4*sum(coverage.values()),"coverage":coverage,"seed":SEED,"splits":N}
    cells={}
    for label,x,small in (("otu99_full_half",otu,False),("otu99_2plus2",otu,True),("asv_full_half",asv,False)):
        z=groups(x,ids,small);cells[label]=z;r[label]=summarize(z)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps({**r,"cells":cells},allow_nan=False)+"\n")
    print("HAYASHI_DIRECTION_RESULT",json.dumps(r,allow_nan=False),flush=True)
if __name__=="__main__":main()
