"""Fixed-depth compositional sensitivity; no redefinition of prior eligible samples."""
import argparse,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_mancuso_response as m
from benchmarks import score_ahif_mancuso_exp2_signal as ext
DEPTH=6840
SEEDS=(202611251,202611252,202611253)
def subset_depth(gr,raw,seed=None):
    rng=None if seed is None else np.random.default_rng(seed)
    out=defaultdict(dict)
    for (r,f),units in sorted(gr.items()):
        for k,days in sorted(units.items()):
            # Unit-level, complete-case, fixed depth gate.
            if set(days)!=set(m.DAYS):continue
            counts=[raw[r,f,k,d] for d in m.DAYS]
            if any(c.sum()<DEPTH for c in counts):continue
            out[r,f][k]={}
            for d,c in zip(m.DAYS,counts):
                assert np.allclose(c/c.sum(),days[d][0],rtol=0,atol=1e-12)
                x=c if rng is None else rng.multivariate_hypergeometric(c,DEPTH)
                if rng is not None:assert x.sum()==DEPTH
                out[r,f][k][d]=(x/x.sum(),days[d][1])
    return out
def summarize(res):
    if "direction" not in res:return res
    return {k:res[k] for k in ("status","coverage","pairs","unique_reactors","n_total_specificity","direction_by_panel","absolute_by_panel")}|{
       "recency":res["recency"]["mean"],"n_recency":res["recency"]["n"],
       "direction":res["direction"]["source_equal"],"absolute":res["absolute_allocation"]["source_equal"],
       "n_signed":res["direction"]["eligible"]}
def audit_exp(root,exp):
    if exp==1:
        g,names,named,audit=m.source(root)
    else:
        g,names,named,audit=ext.source2(root)
    _,b=m.read(root,2)
    raw={}
    for row in b:
        if row["Experiment"]!=str(exp) or row["Technical Issues"]!="" or row["Type"]!="timepoint":continue
        d=int(float(row["Day"]))
        if d not in m.DAYS:continue
        r=float(row["Dilution"]);f=int(float(row["Frequency"]));k=(int(float(row["Replicate"])),row["Vial"])
        assert k in g[r,f] and d in g[r,f][k]
        x=np.array([float(row[n]) for n in names])
        assert np.all(x>=0) and np.all(np.isfinite(x)) and np.all(x==np.floor(x))
        assert (r,f,k,d) not in raw
        raw[r,f,k,d]=x.astype(np.int64)
    assert len(raw)==audit["n_day1_6_qc_profiles"]
    if exp==2:
        m.RATES=ext.RATES2;m.FREQ=ext.FREQ2
    m.EXPECTED=set(m.panels(g,names,named)[3])
    original=list(sorted(m.EXPECTED))
    # Depth-based eligibility determined before calculating outcomes.
    sub=subset_depth(g,raw)
    _,_,coverage,pairs=m.panels(sub,names,named)
    m.EXPECTED=set(pairs)
    assert m.EXPECTED<=set(original)
    results={"depth_restricted_unrarefied":summarize(m.run(sub,names,named)),
       "depth_restricted_strict_seq":summarize(m.run(sub,names,named,strict=True))}
    for seed in SEEDS:
        z=subset_depth(g,raw,seed)
        assert m.panels(z,names,named)[3]==pairs
        results[str(seed)]=summarize(m.run(z,names,named))
    return {"original_pairs":original,"depth_pairs":pairs,"coverage":coverage,
      "n_profiles_above_depth":sum(x.sum()>=DEPTH for x in raw.values()),
      "n_profiles":len(raw),"results":results}
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z={"status":"DEPTH_SENSITIVITY_NOT_ORIGINAL_ASV_QC","depth":DEPTH,"seeds":list(SEEDS),
      "experiment1":audit_exp(a.root,1),"experiment2":audit_exp(a.root,2)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    for e in ["experiment1","experiment2"]:
        for r in z[e]["results"].values():
            r.pop("direction_by_panel",None);r.pop("absolute_by_panel",None)
    print("MANCUSO_DEPTH_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
