"""Source-locked pot-to-field history prediction for exact shared plants."""
from __future__ import annotations
import argparse,csv,hashlib,json
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from benchmarks.score_ahif_vansteijn2026_field import source_panels,source_csv,SITES,YEARS
from benchmarks.probe_ahif_vansteijn2026_field import MD5 as FIELD_MD5

POT_MD5={"alone.csv":"0df05c9e417947e7aaecb84913184fc5",
"simultaneous_arrival.csv":"3a6742d98e4311bf8cd64dc06ddf433b",
"second_arrival.csv":"17f876d2ba031248f6813bddde9a33be",
"traits.csv":"e7b1ffaeff71b05d79c1e9ae158aa55f",
"R_code.Rmd":"16f92a0dde59ce52e30a671dd31859ef"}
FAST={"Agrostis capillaris","Cardamine pratensis","Dianthus deltoides","Festuca rubra","Leucanthemum vulgare","Sedum acre"}
SLOW={"Antennaria dioica","Briza media","Phalaris arundinacea","Poa pratensis"}
SPECIES=sorted(FAST|SLOW)
MODELS={"trait":(0,1,2,3),"base":(0,1,2,3,4),"history":(0,1,2,3,4,5),"pre":(0,1,2,3,4,6)}
TREAT=("FtS","SFearly","SFlate","StF")
N_PERM=1999

def verify(root,manifest):
    for name,md5 in manifest.items():
        if hashlib.md5((root/name).read_bytes()).hexdigest()!=md5:raise ValueError("source hash "+name)

def read(p):
    with p.open(encoding="utf-8-sig",newline="") as f:return list(csv.DictReader(f,delimiter=";"))

def value(x):
    v=float(x)
    if not np.isfinite(v) or v<0:raise ValueError("invalid biomass")
    return v

def pot_features(root,field_taxa):
    verify(root,POT_MD5)
    taxa={r["focal"]:r["func.group"] for r in read(root/"traits.csv")}
    if Counter(taxa.values())!={"Grass":5,"Forb":5,"Legume":5}:raise ValueError("pot species roster")
    if set(taxa)&set(field_taxa)!=set(SPECIES):raise ValueError("unexpected exact species intersection")
    if {s for s in SPECIES if field_taxa[s]["Speed"]=="Fast"}!=FAST or {s for s in SPECIES if field_taxa[s]["Speed"]=="Slow"}!=SLOW:raise ValueError("field germination roster")
    alone=read(root/"alone.csv");sim=read(root/"simultaneous_arrival.csv");pre=read(root/"second_arrival.csv")
    if (len(alone),len(sim),len(pre))!=(74,206,206):raise ValueError("source pot row counts")
    controls=defaultdict(list)
    for x in alone:controls[x["focal"]].append(value(x["biomass.total.CA"]))
    if set(controls)!=set(taxa) or min(map(len,controls.values()))<4:raise ValueError("solo controls")
    def index(rows,col):
        d={}
        for x in rows:
            key=(x["focal"],x["neighbour"])
            if key in d or key[0]==key[1] or not set(key)<=set(taxa):raise ValueError("directed pair key")
            d[key]=value(x[col])
        return d
    a=index(sim,"biomass.total.CT");b=index(pre,"biomass.total.PE")
    common=set(a)&set(b)
    if len(common)!=202:raise ValueError("matched pot pairs")
    features=[];records=[]
    for s in SPECIES:
        partner=SLOW if s in FAST else FAST
        matches=sorted(i for i in partner if (s,i) in common)
        if len(matches)<max(2,len(partner)/2):raise ValueError("inadequate partner coverage "+s)
        solo=float(np.mean(controls[s]))
        if solo<=0:raise ValueError("zero solo")
        cs=np.array([(a[s,i]-solo)/(a[s,i]+solo) for i in matches])
        cp=np.array([(b[s,i]-solo)/(b[s,i]+solo) for i in matches])
        q=field_taxa[s]
        if q["Density"] not in ("High","Low") or q["Form"] not in ("Grass","Forb") or q["Form"]!=taxa[s]:raise ValueError("source trait mismatch")
        f=[float(s in FAST),float(q["Density"]=="High"),float(q["Form"]=="Grass"),float(np.log1p(solo)),float(cs.mean()),float((cp-cs).mean()),float(cp.mean())]
        features.append(f)
        records.append({"species":s,"speed":q["Speed"],"density":q["Density"],"form":q["Form"],"solo_mean":solo,
          "mean_sim_RII":f[4],"mean_history_delta_RII":f[5],"mean_pre_RII":f[6],
          "matched_opposite_partners":matches,"missing_opposite_partners":sorted(partner-set(matches))})
    return np.array(features),records

def field_arrays(panel):
    # arrays order: site, year, complete randomized block, source physical plot, focal
    x=np.empty((2,4,5,4,10));labels=np.empty((2,5,4),int)
    for i,s in enumerate(SITES):
        for j,year in enumerate(YEARS):
            p=panel[s,year]
            if not set(SPECIES)<=set(p["names"]):raise ValueError("matched source taxa absent globally")
            if j and p["plots"]!=panel[s,YEARS[0]]["plots"]:raise ValueError("plot identity drift")
            xx=p["X"][:,[p["names"].index(v) for v in SPECIES]].reshape(5,4,10)
            lab=np.array([TREAT.index(t) for t in p["history"]]).reshape(5,4)
            if any(set(row)!={0,1,2,3} for row in lab):raise ValueError("block allocation")
            if j and not np.array_equal(lab,labels[i]):raise ValueError("treatment drift")
            labels[i]=lab;x[i,j]=xx
    return x,labels

def contrasts(x,labs):
    # focal own calendar time late, compare earlier opposite group to simultaneous late.
    # lab 3=StF: fast late, slow early. lab 0=FtS: slow late, fast early. lab 2=SFlate.
    z=np.empty((2,5,10))
    for i in range(2):
        for b in range(5):
            d={k:x[i,b,np.flatnonzero(labs[i,b]==k)[0]] for k in range(4)}
            z[i,b]=np.where([s in FAST for s in SPECIES],d[3]-d[2],d[0]-d[2])
    return z.mean(1)

def design_maps(features,cols):
    # Each focal field response withheld, and entire target site withheld.
    w=np.zeros((10,10));slope=np.zeros((10,10))
    for j in range(10):
        tr=np.arange(10)!=j
        f=features[:,cols];mu=f[tr].mean(0);sd=f[tr].std(0);z=(f-mu)/np.where(sd>1e-12,sd,1.)
        a=np.column_stack((np.ones(10),z))
        p=np.array([0.]+[2.]*len(cols))
        q=np.linalg.solve(a[tr].T@a[tr]+np.diag(p),a[tr].T)
        w[j,tr]=a[j]@q
        slope[j,tr]=q[-1]
    if not np.allclose(np.diag(w),0) or not np.allclose(w.sum(1),1):raise ValueError("holdout map leakage")
    return w,slope

def predict(y,w):
    return np.stack([w@y[1],w@y[0]])

def gain(b,e):
    return (b-e)/b if b>1e-15 else 0.

def score(y,maps):
    p={k:predict(y,v[0]) for k,v in maps.items()}
    loss={k:(y-v)**2 for k,v in p.items()}
    totals={k:float(v.sum()) for k,v in loss.items()}
    b,e=loss["base"],loss["history"]
    return {"losses":totals,"history_gain":gain(totals["base"],totals["history"]),
      "sim_gain_over_trait":gain(totals["trait"],totals["base"]),
      "pre_gain_over_base":gain(totals["base"],totals["pre"]),
      "by_site":{s:gain(float(b[i].sum()),float(e[i].sum())) for i,s in enumerate(SITES)},
      "by_speed":{g:gain(float(b[:,m].sum()),float(e[:,m].sum())) for g,m in
        (("Fast",[s in FAST for s in SPECIES]),("Slow",[s in SLOW for s in SPECIES]))},
      "n_improved_cells":int(np.sum(e<b-1e-12)),
      "gain_without_focal":{s:gain(float(b[:,np.arange(10)!=j].sum()),float(e[:,np.arange(10)!=j].sum())) for j,s in enumerate(SPECIES)},
      "cells":[{"site":s,"species":sp,"observed_contrast":float(y[i,j]),
                 "predictions":{k:float(v[i,j]) for k,v in p.items()},
                 "history_training_standardized_slope":float(maps["history"][1][j]@y[1-i])}
                for i,s in enumerate(SITES) for j,sp in enumerate(SPECIES)]}

def run(pot,field,nperm=N_PERM):
    verify(field,FIELD_MD5)
    panel,audit=source_panels(field)
    taxa={x["Species"]:x for x in source_csv(field/"speciesinfo.csv")}
    f,records=pot_features(pot,taxa)
    x,lab=field_arrays(panel)
    maps={k:design_maps(f,c) for k,c in MODELS.items()}
    out={year:score(contrasts(np.log1p(x[:,j]),lab),maps) for j,year in enumerate(YEARS)}
    observed=out["2024"]["history_gain"]
    rng=np.random.default_rng(20261125);null=[]
    for _ in range(nperm):
        l=lab.copy()
        for i in range(2):
            for b in range(5):l[i,b]=rng.permutation(l[i,b])
        y=contrasts(np.log1p(x[:,-1]),l)
        b,e=(float(np.sum((y-predict(y,maps[k][0]))**2)) for k in ("base","history"))
        null.append(gain(b,e))
    null=np.array(null)
    return {"status":"SOURCE_LOCKED_POT_FIELD_TRANSFER_NOT_CAUSAL_MEDIATION","n_shared_species":10,
      "n_pot_matched_directed_pairs":202,"n_cross_guild_measured_pairs":sum(len(r["matched_opposite_partners"]) for r in records),
      "n_field_plots":40,"n_field_blocks":10,"source_species_features":records,
      "by_year_log1p_hit_contrast":out,"raw_2024_sensitivity":score(contrasts(x[:,-1],lab),maps),
      "primary_sharp_null_reference":{"n_permutations":nperm,"p_upper":float((1+np.count_nonzero(null>=observed-1e-12))/(nperm+1)),
       "null_median":float(np.median(null)),"null95":np.quantile(null,[.025,.975]).tolist()},
      "limits":["Distinct pot and field experiments but same research team and shared taxa.",
        "One pot per directed species pair; measured net biomass response not dynamic competition coefficient.",
        "Other-site and other-focal field calibration; two sites only, not independent species replication.",
        "Only 10 common sown species. Other field species remain unmeasured.",
        "Sharp null is no schedule effect, not a nonzero-effect no-moderation test. No causal mediation."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--pot",type=Path,required=True);p.add_argument("--field",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True);p.add_argument("--permutations",type=int,default=N_PERM);a=p.parse_args()
    z=run(a.pot,a.field,a.permutations)
    a.output.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("POT_FIELD_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
