"""Mancuso independent microbial relative-composition lane; fixed PR8 diagnostics."""
import argparse,csv,hashlib,json,re
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from benchmarks import score_ahif_response_allocation as alloc
from benchmarks import score_ahif_response_specificity as spec
SHA={1:"a3f9989e3f3ef7d090b65d3e78d9e89fd208e4dcc40daaeb484d9b0ae80500e6",
     2:"6347ea67083f25634d00c08da599fe89e3c67d73f656908ebde54d2087c7198b"}
RATES=(.1,.2,.3,.4);FREQ=(1,4,16);DAYS=list(range(1,7))
EXPECTED={(.2,4),(.2,16),(.3,1),(.3,4),(.3,16),(.4,1),(.4,4),(.4,16)}
META=["Nickname","Vial","Dilution","Frequency","Replicate","Day","Generations","Biofilm","DNAconc","PrepPlate","Experiment","Type","Technical Issues","SequencingType"]
def read(root,i):
    p=root/f"data{i}.csv"; assert hashlib.sha256(p.read_bytes()).hexdigest()==SHA[i]
    with p.open(encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f);h=r.fieldnames;assert len(set(h))==len(h)
        rows=list(r)
    assert all(None not in row and all(v is not None for v in row.values()) for row in rows)
    return h,rows
def source(root):
    h1,a=read(root,1);h2,b=read(root,2)
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
    names=[k for k in h2 if k not in ["index"]+META]
    assert len(names)==103 and all(k.startswith("D_0__") or k.startswith("Unassigned") for k in names)
    named=[i for i,k in enumerate(names) if k.startswith("D_0__Bacteria;") and re.search(r";D_5__[^;]+$",k) and not re.search(r"uncultured|metagenome|unknown|unassigned|unidentified|chloroplast|mitochondri",k,re.I)]
    gr=defaultdict(dict);unit={};seq=Counter();rawsum=[];shares=[]
    for k,r in bb.items():
        if k not in aa or r["Experiment"]!="1" or r["Type"]!="timepoint":continue
        d=int(float(r["Day"]));assert d==float(r["Day"])
        if d not in DAYS:continue
        rate=float(r["Dilution"]);f=int(float(r["Frequency"]));rep=int(float(r["Replicate"]));v=r["Vial"]
        assert rate in RATES and f in FREQ+(600,) and rep in range(1,5)
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
    audit=dict(source_rows=[len(a),len(b)],n_join=len(aa),n_qc_excluded=56,n_seq_mismatch_all=42,bin_count=len(names),named_bacterial_genus_bins=len(named),n_day1_6_qc_profiles=len(rawsum),row_sum_range=[min(rawsum),max(rawsum)],named_share_range=[min(shares),max(shares)],n_seq_mismatch_day1_6=sum(n for k,n in seq.items() if k[-1]!=k[-2]),sequencing_table=[[*k,v] for k,v in sorted(seq.items())])
    return gr,names,named,audit
def panels(gr,names,named,strict=False,only_named=False):
    mats={};rosters={};coverage={}
    for rate in RATES:
        for f in FREQ+(600,):
            g=gr[(rate,f)];keys=sorted(k for k,d in g.items() if set(d)==set(DAYS) and (not strict or not any(z[1] for z in d.values())))
            coverage[f"{rate}:{f}"]=len(keys)
            if len(keys)!=4 or [k[0] for k in keys]!=[1,2,3,4]:continue
            mat=np.array([[g[k][d][0] for k in keys] for d in DAYS])
            if only_named:
                mat=mat[:,:,named];s=mat.sum(2)
                if np.any(s<=0):continue
                mat=mat/s[:,:,None]
            mats[(rate,f)]=mat;rosters[(rate,f)]=[k[1] for k in keys]
    pairs=[(r,f) for r in RATES for f in FREQ if (r,f) in mats and (r,600) in mats]
    return mats,rosters,coverage,pairs
def bypanel(rows,keys):
    out={}
    for p in sorted({r["panel"] for r in rows}):
        z=[r for r in rows if r["panel"]==p and all(r[k] is not None for k in keys)]
        out[p]={"n":len(z),**{k:float(np.mean([r[k] for r in z])) if z else None for k in keys}}
    return out
def run(gr,names,named,strict=False,only_named=False):
    mats,rosters,coverage,pairs=panels(gr,names,named,strict,only_named)
    if not strict and not only_named:assert set(pairs)==EXPECTED
    if not pairs:return {"status":"STOP_no_eligible_panels","coverage":coverage}
    taxa=[names[i] for i in named] if only_named else names
    G=np.ones((len(taxa),1));rrows=[];srows=[]
    alloc.seen.clear();spec.seen.clear()
    for r,f in pairs:
        a=mats[r,600];b=mats[r,f];prior=b.mean(1)-a.mean(1);k=f"{r}:{f}"
        rrows+=alloc.evaluate("mancuso",k,DAYS,a,b,taxa,G,["all"],False,prior)
        srows+=spec.evaluate("mancuso",k,DAYS,a,b,taxa,G,["all"],False,prior)
    assert len(rrows)==len(pairs)*36*3 and len(srows)==len(pairs)*36*4
    vr=[r for r in rrows if r["recency"] is not None and all(r[k] is not None for k in ("H","M2","V","M1","L"))]
    for r in vr:
        r["H_vs_M2"]=r["M2"]-r["H"];r["H_vs_V"]=r["V"]-r["H"]
    rkeys=["recency","old","recent","M1","M2","V","H","L","H_vs_M2","H_vs_V","dominant_overlap","lag_overlap"]
    rec={"n":len(vr),"total":len(rrows),"mean":{k:alloc.nested(vr,k).get("mancuso") for k in rkeys},"by_panel":bypanel(vr,rkeys)}
    keys=["B","S","SQ","V","R","BS","BSQ","BV","BR"]
    cs=[("R_over_"+k,k,"R") for k in ["S","SQ","V"]]+[("BR_over_"+k,k,"BR") for k in ["B","BS","BSQ","BV"]]
    absout=spec.summary(srows,keys,cs)
    nov=spec.summary(srows,[k for k in keys if k not in ("V","BV")],[c for c in cs if c[1] not in ("V","BV")])
    signkeys=["signed","signblind","signed_history","same","opposite","abstain"]
    signcs=[("sign_gain","signblind","signed"),("recent_sign_gain","signed_history","signed")]
    signs=spec.summary(srows,signkeys,signcs)
    sp=[r for r in srows if all(r[k] is not None for k in keys)]
    sg=[r for r in srows if all(r[k] is not None for k in signkeys)]
    # All comparisons use equal panel weights in spec.summary. Per-panel summaries from common denominators.
    return dict(status="SCORED_RELATIVE_COMPOSITION_ONLY",coverage=coverage,pairs=[list(k) for k in pairs],unique_reactors=len({v for k in pairs for arm in [(k[0],600),k] for v in rosters[arm]}),n_total_specificity=len(srows),
      recency=rec,absolute_allocation={k:v for k,v in absout.items() if k!="exclude_kansas"},without_variance={k:v for k,v in nov.items() if k!="exclude_kansas"},direction={k:v for k,v in signs.items() if k!="exclude_kansas"},
      absolute_by_panel=bypanel(sp,keys+[n for n,_,_ in cs]),direction_by_panel=bypanel(sg,signkeys+[n for n,_,_ in signcs]),cells={"recency":rrows,"specificity":srows})
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    alloc.selftest();spec.selftest()
    gr,names,named,audit=source(a.root)
    z={"status":"INDEPENDENT_PUBLICATION_RELATIVE_COMPOSITION_NOT_META","source_audit":audit,
       "primary":run(gr,names,named),"strict_sequencing":run(gr,names,named,strict=True),
       "named_genera":run(gr,names,named,only_named=True)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    brief={k:{j:v for j,v in q.items() if j!="cells"} if isinstance(q,dict) else q for k,q in z.items()}
    # Concise logs, full detailed aggregate and cell scores in artifact.
    brief["source_audit"]={k:v for k,v in audit.items() if k!="sequencing_table"}
    for label in ("primary","strict_sequencing","named_genera"):
        if "recency" in brief[label]:
            brief[label]["recency"].pop("by_panel",None)
    print("MANCUSO_RESPONSE_RESULT",json.dumps(brief,allow_nan=False))
if __name__=="__main__":main()
