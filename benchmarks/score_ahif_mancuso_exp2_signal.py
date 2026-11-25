"""Separate second-experiment and signal scale diagnostic; prior score untouched."""
import argparse,hashlib,json,re
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_mancuso_response as m
from benchmarks import score_ahif_cross_source_concentration as splitsrc

RATES2=(.3,.6,.9,1.5)
FREQ2=(1,4)
EXPECT2={(.3,1),(.6,1),(.9,1),(.9,4),(1.5,1),(1.5,4)}
def source2(root):
    h1,a=m.read(root,1);h2,b=m.read(root,2)
    aa={r["SampleID"]:r for r in a};bb={r["index"]:r for r in b}
    assert len(aa)==len(a)==824 and len(bb)==len(b)==880 and aa.keys()<=bb.keys()
    assert Counter(r["Technical Issues"] for r in b)==Counter({"":824,"contamination":14,"missed dilution":42})
    assert aa.keys()=={k for k,r in bb.items() if r["Technical Issues"]==""}
    mismatch=set()
    for k,r in aa.items():
        s=bb[k]
        for col in ("Vial","Replicate","Experiment","Type","Technical Issues"):
            assert r[col]==s[col],(k,col)
        for col in ("Day","Dilution","Frequency"):
            assert float(r[col])==float(s[col]),(k,col)
        if r["SequencingType"]!=s["SequencingType"]:mismatch.add(k)
    assert len(mismatch)==42
    names=[k for k in h2 if k not in ["index"]+m.META]
    assert len(names)==103 and all(k.startswith("D_0__") or k.startswith("Unassigned") for k in names)
    named=[i for i,k in enumerate(names) if k.startswith("D_0__Bacteria;") and re.search(r";D_5__[^;]+$",k) and not re.search(r"uncultured|metagenome|unknown|unassigned|unidentified|chloroplast|mitochondri",k,re.I)]
    gr=defaultdict(dict);unit={};seq=Counter();rawsum=[];shares=[]
    for k,r in bb.items():
        if k not in aa or r["Experiment"]!="2" or r["Type"]!="timepoint":continue
        d=int(float(r["Day"]));assert d==float(r["Day"])
        if d not in m.DAYS:continue
        rate=float(r["Dilution"]);f=int(float(r["Frequency"]));rep=int(float(r["Replicate"]));v=r["Vial"]
        assert rate in RATES2 and f in FREQ2+(600,) and rep in range(1,5)
        group=(rate,f,rep)
        if v in unit:assert unit[v]==group
        unit[v]=group
        key=(rep,v)
        if key not in gr[(rate,f)]:gr[(rate,f)][key]={}
        assert d not in gr[(rate,f)][key]
        arr=np.array([float(r[n]) for n in names])
        assert np.isfinite(arr).all() and (arr>=0).all() and arr.sum()>0,(k,"invalid profile")
        rawsum.append(float(arr.sum()));arr=arr/arr.sum()
        shares.append(float(arr[named].sum()))
        gr[(rate,f)][key][d]=(arr,k in mismatch)
        seq[(str(rate),str(f),str(d),r["SequencingType"],aa[k]["SequencingType"])]+=1
    audit=dict(n_day1_6_qc_profiles=len(rawsum),row_sum_range=[min(rawsum),max(rawsum)],named_share_range=[min(shares),max(shares)],n_seq_mismatch_day1_6=sum(n for k,n in seq.items() if k[-1]!=k[-2]),sequencing_table=[[*k,v] for k,v in sorted(seq.items())])
    return gr,names,named,audit
def signal_by_panel(gr,names,named):
    mats,rosters,coverage,pairs=m.panels(gr,names,named)
    ans={}
    for r,f in pairs:
        a=mats[r,600];b=mats[r,f]
        dt=[]
        for t in range(len(m.DAYS)):
            cells=[]
            for ia,ib,ha,hb in splitsrc.splits(4,False):
                assert not set(ia)&set(ha) and not set(ib)&set(hb)
                x=b[t,list(ib)].mean(0)-a[t,list(ia)].mean(0)
                y=b[t,list(hb)].mean(0)-a[t,list(ha)].mean(0)
                cells.append([float(x@y),float(np.sum((x-y)**2)/2),float(np.sum(x*x)),float(np.sum(y*y))])
            assert len(cells)==36
            s,n,train,held=np.mean(cells,axis=0)
            full=b[t].mean(0)-a[t].mean(0)
            # Algebraic identity under complementary 2+2 splits, not independent evidence.
            assert abs((train+held)/2-s-n)<1e-10
            dt.append(dict(day=m.DAYS[t],signal=float(s),noise=float(n),train_sq=float(train),held_sq=float(held),full_tv=float(np.abs(full).sum()/2)))
        def reduce(ds):
            out={k:float(np.mean([x[k] for x in ds])) for k in ["signal","noise","train_sq","held_sq","full_tv"]}
            out["signal_noise_ratio"]=out["signal"]/out["noise"] if out["noise"]>1e-15 else None
            return out
        ans[f"{r}:{f}"]={"days":dt,"target_days3_6":reduce(dt[2:6]),"predictor_days2_5":reduce(dt[1:5])}
    return ans
def brief(z):
    return {k:v for k,v in z.items() if k not in ("cells","without_variance")}
def freq_summary(sig,direction,rates,freqs):
    out={}
    for f in freqs:
        keys=[f"{r}:{f}" for r in rates]
        assert all(k in sig and k in direction for k in keys)
        t={}
        for period in ("target_days3_6","predictor_days2_5"):
            t[period]={k:float(np.mean([sig[z][period][k] for z in keys])) for k in ("signal","noise","full_tv")}
            a=t[period];a["signal_noise_ratio"]=a["signal"]/a["noise"] if a["noise"]>1e-15 else None
        t["direction"]={k:float(np.mean([direction[z][k] for z in keys])) for k in ("same","opposite","sign_gain","recent_sign_gain")}
        out[str(f)]=t
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    m.alloc.selftest();m.spec.selftest()
    g1,n1,ng1,a1=m.source(a.root)
    s1=signal_by_panel(g1,n1,ng1)
    r1=m.run(g1,n1,ng1)
    assert r1["status"]=="SCORED_RELATIVE_COMPOSITION_ONLY"
    # Reproduce original primary exactly before loading second experiment.
    assert abs(r1["direction"]["source_equal"]["sign_gain"]-.1617459982758671)<1e-10
    assert abs(r1["recency"]["mean"]["recency"]-.08402604391726227)<1e-10
    g2,n2,ng2,a2=source2(a.root)
    assert n1==n2 and ng1==ng2
    m.RATES=RATES2;m.FREQ=FREQ2;m.EXPECTED=EXPECT2
    p2=m.run(g2,n2,ng2)
    q2=m.run(g2,n2,ng2,strict=True)
    n2r=m.run(g2,n2,ng2,only_named=True)
    s2=signal_by_panel(g2,n2,ng2)
    f1=freq_summary(s1,r1["direction_by_panel"],[.3,.4],[1,4,16])
    f2=freq_summary(s2,p2["direction_by_panel"],[.9,1.5],[1,4])
    shared={e:{"direction":r["direction_by_panel"]["0.3:1"],"absolute":r["absolute_by_panel"]["0.3:1"]} for e,r in [("experiment1",r1),("experiment2",p2)]}
    z=dict(status="SAME_PUBLICATION_SECOND_EXPERIMENT_NOT_META",source_audit_exp2=a2,primary=p2,strict_sequencing=q2,named_genera=n2r,signal_exp1=s1,signal_exp2=s2,frequency_exp1=f1,frequency_exp2=f2,shared_0_3_1=shared)
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    summary={k:v for k,v in z.items() if k not in ("primary","strict_sequencing","named_genera","signal_exp1","signal_exp2","source_audit_exp2")}
    summary["source_audit_exp2"]={k:v for k,v in a2.items() if k!="sequencing_table"}
    for label in ["primary","strict_sequencing","named_genera"]:
        r=z[label]
        if "direction" not in r:summary[label]=brief(r);continue
        summary[label]={k:r[k] for k in ["status","coverage","pairs","unique_reactors","n_total_specificity","direction_by_panel","absolute_by_panel"]}
        summary[label]["recency"]={k:v for k,v in r["recency"].items() if k!="by_panel"}
        summary[label]["direction"]=r["direction"]["source_equal"]
        summary[label]["absolute"]=r["absolute_allocation"]["source_equal"]
        summary[label]["n_signed"]=r["direction"]["eligible"]
    print("MANCUSO_EXP2_SIGNAL_RESULT",json.dumps(summary,allow_nan=False))
if __name__=="__main__":main()
