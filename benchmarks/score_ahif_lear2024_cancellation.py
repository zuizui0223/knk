"""Resident-only aggregate response cancellation in source-locked Lear 2024 data."""
from __future__ import annotations
import argparse,csv,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from scipy.stats import t

FILES={"invader timing data.txt":"208b709aabc1b1b7d1f37cc550f5c304",
       "README.txt":"3c2e18e459a910a95b0063bd599eb474",
       "invader timing markdown.Rmd":"de0b10f9f62156c3684ee1db6789ba89"}
SPECIES=("p","o","s","a","v")
LEVELS=("A","B","C","D")
DAYS=("5","6","7","7.9","8.1","9","10","11")
PAIRS=(("7.9","8.1"),("7","9"),("6","10"),("5","11"))
SPLITS=256
SEED=20261125

def read_source(root):
    for name,md5 in FILES.items():
        if hashlib.md5((root/name).read_bytes()).hexdigest()!=md5:
            raise ValueError("source digest mismatch "+name)
    with (root/"invader timing data.txt").open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f,delimiter="\t")
        rows=list(rd)
    if set(rd.fieldnames or ())!={"day","propagule","rep",*SPECIES,"invader"}:
        raise ValueError("source schema mismatch")
    if len(rows)!=197:raise ValueError("source row count mismatch")
    g=defaultdict(list);labels=defaultdict(list)
    for x in rows:
        key=x["day"].strip(),x["propagule"].strip()
        vals=np.array([float(x[s]) for s in SPECIES])
        if not np.isfinite(vals).all() or np.any(vals<0) or np.any(vals!=np.floor(vals)) or vals.sum()<=0:
            raise ValueError("invalid resident count")
        rep=int(x["rep"])
        if rep not in range(1,7):raise ValueError("invalid replicate label")
        g[key].append(vals);labels[key].append(rep)
    expected={(d,p):5 if (d,p)==("8.1","B") else 6 for d in DAYS for p in LEVELS}
    expected["NI","NI"]=6
    if {k:len(v) for k,v in g.items()}!=expected:raise ValueError("unexpected source treatment coverage")
    bad=[]
    for k,ls in sorted(labels.items()):
        duplicates={str(a):b for a,b in Counter(ls).items() if b>1}
        if duplicates:
            bad.append({"day":k[0],"level":k[1],"duplicate_labels":duplicates,
                        "missing_labels":sorted(set(range(1,7))-set(ls))})
    if {(x["day"],x["level"]) for x in bad}!={("6","A"),("6","B")}:
        raise ValueError("unexpected label ambiguity")
    return {k:np.asarray(v) for k,v in g.items()}, {
        "source_rows":len(rows),"planned_microcosms":198,"observed_treatment_cells":33,
        "documented_missing_cell":{"day":"8.1","level":"B","missing":1},
        "ambiguous_rep_cells":bad,"rep_ids_used_for_pairing":False,
        "source_file_md5":FILES,"resident_columns":list(SPECIES)}

def metrics(d):
    d=np.asarray(d,float);gross=float(np.abs(d).sum());net=float(abs(d.sum()))
    c=gross-net
    if c<-1e-9:raise ValueError("triangle inequality")
    return {"delta":d.tolist(),"signed_total_delta":float(d.sum()),"gross":gross,"net":net,
            "cancellation":max(c,0.),"fraction":max(c,0.)/gross if gross>1e-12 else None}

def delta(a,b):
    return np.mean(b,axis=0)-np.mean(a,axis=0)

def cell_deltas(groups,before,after):
    return [(groups[before,p],groups[after,p]) for p in LEVELS]

def intervals(cells,nfamily):
    # Independent group mean contrasts, equal weighting over source dose strata.
    ds=np.mean([delta(a,b) for a,b in cells],axis=0)
    v=[];df=[]
    for a,b in cells:
        for x in (a,b):
            if len(x)<2:raise ValueError("insufficient replicate cells")
            v.append(x.var(axis=0,ddof=1)/len(x)/len(cells)**2)
            df.append(len(x)-1)
    v=np.asarray(v);df=np.asarray(df)
    vv=v.sum(0);den=np.sum(v*v/df[:,None],axis=0)
    deg=np.divide(vv*vv,den,out=np.zeros_like(vv),where=den>0)
    out=[]
    for j,s in enumerate(SPECIES):
        if vv[j]<=1e-15:
            out.append({"resident":s,"delta":float(ds[j]),"status":"STOP_ZERO_ESTIMATED_VARIANCE",
                        "simultaneous_ci95":None})
            continue
        half=float(t.ppf(1-.05/(2*nfamily),deg[j])*np.sqrt(vv[j]))
        out.append({"resident":s,"delta":float(ds[j]),"se":float(np.sqrt(vv[j])),
                    "df":float(deg[j]),"simultaneous_ci95":[float(ds[j]-half),float(ds[j]+half)],
                    "status":"APPROXIMATE_WELCH_BONFERRONI"})
    positive=[r["resident"] for r in out if r["simultaneous_ci95"] and r["simultaneous_ci95"][0]>0]
    negative=[r["resident"] for r in out if r["simultaneous_ci95"] and r["simultaneous_ci95"][1]<0]
    return {"n_family":nfamily,"rows":out,"positive":positive,"negative":negative,
            "opposing_signs_supported":bool(positive and negative)}

def validation_gap(train,test):
    return float(np.sign(train)@test-np.sign(train.sum())*test.sum())

def crossfit(cells,seed,nsplits=SPLITS):
    rng=np.random.default_rng(seed)
    pooled=[];within=[];signed_gross=[];signed_net=[]
    for _ in range(nsplits):
        train=[];held=[]
        for a,b in cells:
            pa=rng.permutation(len(a));pb=rng.permutation(len(b))
            ha=len(a)//2;hb=len(b)//2
            train.append(delta(a[pa[:ha]],b[pb[:hb]]))
            held.append(delta(a[pa[ha:]],b[pb[hb:]]))
        for tr,te in ((train,held),(held,train)):
            tr=np.array(tr);te=np.array(te)
            q=tr.mean(0);v=te.mean(0)
            pooled.append(validation_gap(q,v))
            within.append(np.mean([validation_gap(a,b) for a,b in zip(tr,te)]))
            signed_gross.append(float(np.sign(q)@v))
            signed_net.append(float(np.sign(q.sum())*v.sum()))
    return {"n_seeded_twofold_splits":nsplits,"pooled_validated_gap":float(np.mean(pooled)),
            "mean_within_level_validated_gap":float(np.mean(within)),
            "pooled_validated_signed_gross":float(np.mean(signed_gross)),
            "pooled_validated_signed_net":float(np.mean(signed_net)),
            "descriptive_only":True}

def run(root):
    g,audit=read_source(root);out=[]
    for i,(before,after) in enumerate(PAIRS):
        cells=cell_deltas(g,before,after)
        vectors=np.array([delta(a,b) for a,b in cells])
        level=[{"level":p,"n_before":len(a),"n_after":len(b),**metrics(d)}
               for p,(a,b),d in zip(LEVELS,cells,vectors)]
        q={"before":before,"after":after,"pooled":metrics(vectors.mean(0)),
           "by_level":level,"mean_within_level_cancellation":float(np.mean([x["cancellation"] for x in level])),
           "mean_within_level_gross":float(np.mean([x["gross"] for x in level])),
           "mean_within_level_net":float(np.mean([x["net"] for x in level]))}
        if before=="6":
            q.update({"status":"DESCRIPTIVE_ONLY_AMBIGUOUS_REP_LABELS","intervals":None,"crossfit":None})
        else:
            q.update({"status":"SCORED_APPROXIMATE_INDEPENDENT_GROUPS",
                     "intervals":intervals(cells,5 if i==0 else 10),
                     "crossfit":crossfit(cells,SEED+i)})
        out.append(q)
    return {"status":"SOURCE_LOCKED_RESIDENT_COUNT_CANCELLATION_NOT_INTERACTION_MEDIATION",
       "source_doi":"10.5281/zenodo.12528023","audit":audit,"contrasts":out,
       "control_resident_mean":g["NI","NI"].mean(0).tolist(),
       "limits":["Only resident count responses; no newcomer performance or growth analysis.",
          "C=G−N nonnegative by construction; simultaneous intervals and descriptive cross-fit assess sign reproducibility.",
          "Replicate labels ambiguous in two day-6 cells; no inferential scoring for affected comparison.",
          "Approximate independent-group Welch intervals; not an exact randomization test.",
          "No causal interaction coefficients, biomass inference or independent replication across contrasts."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run(a.source)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("LEAR_CANCELLATION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
