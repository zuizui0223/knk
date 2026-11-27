"""Gibbs plant spatial-context transfer, source-locked, leave replicate-label/pair out."""
import argparse,csv,hashlib,io,itertools,json,math,urllib.request
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
REF="15c1840f0270792f5aad7cd1f1e22a9595a4e8d4"
BASE="https://raw.githubusercontent.com/theogibbs/SplitMixed/"+REF+"/"
FILES={"plots.csv":"57961fe2a2679a3bbe8e55c7ca1b72a249b5c94e","fecundity_data.csv":"85183b55df19da532cb80ff3e9c501a597080418","official_census.csv":"aee2e8024e166565f572da0599d7badd88f71660","seed_head_calibration.csv":"af52bf577d8feab449c832ba74a25b6250078887","FigFecundityDifferences.R":"c1911cb7d46acbbf515c71d823b1c19b0d4b5daf"}
SP=("AC","FE","PL","SA","UR");FOCAL=("AC","SA","UR")
PAIRS=list(itertools.combinations(SP,2));REPS=tuple(range(1,16))
MISSING={"","NA","NaN","nan","None"}
MODELS=("zero","other_pair_mean","known_pair_mean","other_focals_same_pair","held_pair_additive_identity")
def get(path):
    req=urllib.request.Request(BASE+path,headers={"User-Agent":"AHIF academic source audit"})
    with urllib.request.urlopen(req,timeout=60) as r:b=r.read(1_000_001)
    assert len(b)<=1_000_000 and hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==FILES[path],path
    return b
def dicts(b):return list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"))))
def val(s):
    if s is None or s.strip() in MISSING:return None
    x=float(s);assert math.isfinite(x)
    return x
def load():
    data={k:get(k) for k in FILES}
    src=data["FigFecundityDifferences.R"].decode("utf-8-sig")
    for token in ("AC_seed_head[2:33]","UR_seed_head[2:31]","SA_seed_head[-c(1, nrow(SA_seed_head)),]","Seeds = DiameterToSeeds_SA"):
        assert token in src,token
    cal=list(csv.reader(io.StringIO(data["seed_head_calibration.csv"].decode("utf-8-sig"))))[1:]
    ac=np.array([float(r[5]) for r in cal[1:33]])
    ur=np.array([float(r[0]) for r in cal[1:31]])
    sa=np.array([[float(r[2]),float(r[3])] for r in cal[1:-1]])
    assert len(ac)==32 and len(ur)==30 and np.isfinite(sa).all()
    intercept,slope=np.linalg.lstsq(np.column_stack([np.ones(len(sa)),sa[:,0]**2]),sa[:,1],rcond=None)[0]
    plots=dicts(data["plots.csv"]);assert len(plots)==315
    p={}
    for r in plots:
        k=int(r["plot"]);assert k not in p
        t=r["treatment"];rep=int(r["replicate"]);pair=tuple(sorted((r["sp1"],r["sp2"])))
        assert rep in REPS and (t=="bare" or (t in ("split","mix") and pair in PAIRS))
        p[k]=(rep,t,pair)
    assert len({(rep,t,pair) for rep,t,pair in p.values()})==315
    assert all(sum(rep==v[0] for v in p.values())==21 for rep in REPS)
    census=dicts(data["official_census.csv"]);c={}
    for r in census:
        k=int(r["plot"]);f=r["focal"];key=(k,f)
        assert key not in c and f in FOCAL
        pair=tuple(sorted((r["sp1"],r["sp2"])));assert (r["treatment"],pair)==p[k][1:]
        n1,n2=val(r["sp1_count"]),val(r["sp2_count"])
        if r["sp1"]!=pair[0]:n1,n2=n2,n1
        if n1 is not None and n2 is not None:assert min(n1,n2)>=0 and n1==int(n1) and n2==int(n2)
        c[key]=(n1,n2)
    assert len(c)==945 and set(c)=={(k,f) for k in p for f in FOCAL}
    fec=dicts(data["fecundity_data.csv"]);assert len(fec)==315
    y={};neg_heads=0;neg_totals=[];fseen=set()
    for r in fec:
        k=int(r["Plot"]);assert k not in fseen and k in p;fseen.add(k)
        if p[k][1]=="bare":continue
        for f,col,mul in (("AC","Acmispon",ac.mean()),("UR","Uropappus",ur.mean())):
            v=val(r[col]);assert v is None or v>=0
            y[k,f]=None if v is None else float(v*mul)
        h=[val(r[col]) for col in ("Salvia","2 Salvia","3 Salvia","4 Salvia","5 Salvia","6 Salvia")]
        h=[v for v in h if v is not None];assert all(v>=0 for v in h)
        vals=[float(intercept+slope*v*v) for v in h]
        neg_heads+=sum(v<0 for v in vals)
        y[k,"SA"]=None if not h else sum(vals)
        if vals and sum(vals)<0:neg_totals.append(k)
    assert fseen==set(p)
    audit={"source_commit":REF,"files":{k:{"git_blob":FILES[k],"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b)} for k,b in data.items()},
      "plots":315,"mixed_or_split_plots":300,"bare_plots":15,"replicate_labels":15,"competitor_pairs":10,"focal_species":3,
      "calibration":{"AC_n":len(ac),"AC_mean":float(ac.mean()),"UR_n":len(ur),"UR_mean":float(ur.mean()),"SA_n":len(sa),"SA_intercept":float(intercept),"SA_slope":float(slope)},
      "negative_SA_head_predictions":neg_heads,"negative_SA_focal_totals":neg_totals}
    if neg_totals:
        print("GIBBS_SOURCE_AUDIT",json.dumps(audit),flush=True);raise ValueError("STOP raw author estimate negative focal total")
    records=[];miss=Counter();potkey={v:k for k,v in p.items()}
    for f,pair,rep in itertools.product(FOCAL,PAIRS,REPS):
        split=potkey[rep,"split",pair];mix=potkey[rep,"mix",pair]
        a,b=y[split,f],y[mix,f];cs,cm=c[split,f],c[mix,f]
        for t,v,n in (("split",a,cs),("mix",b,cm)):
            if v is None:miss[f+"|"+t+"|fecundity"]+=1
            if None in n:miss[f+"|"+t+"|census"]+=1
        if a is None or b is None or None in cs or None in cm:continue
        delta=np.array(cs)-np.array(cm)
        records.append({"focal":f,"pair":"-".join(pair),"rep":rep,"split_plot":split,"mix_plot":mix,
            "split_seeds":a,"mix_seeds":b,"D":math.log1p(a)-math.log1p(b),
            "split_n1":cs[0],"split_n2":cs[1],"mix_n1":cm[0],"mix_n2":cm[1],
            "delta_total":float(delta.sum()/50),"delta_composition":float((delta[0]-delta[1])/50)})
    groups=Counter((r["focal"],r["pair"]) for r in records)
    audit.update({"missing":dict(miss),"complete_pairs":len(records),"possible_pairs":450,
       "group_counts":{"|".join(k):v for k,v in sorted(groups.items())}})
    print("GIBBS_SOURCE_AUDIT",json.dumps(audit,allow_nan=False),flush=True)
    assert len(groups)==30 and min(groups.values())>=5
    return records,audit
def features(pair):
    return [int(x in pair.split("-")) for x in SP]
def mean(x):
    assert len(x)>0
    return float(np.mean(x))
def group(r):return r["focal"],r["pair"]
def residual_slopes(rows):
    by=defaultdict(list)
    for r in rows:by[group(r)].append(r)
    xx=[];yy=[]
    for g in by.values():
        x=np.array([[r["delta_total"],r["delta_composition"]] for r in g]);y=np.array([r["D"] for r in g])
        xx.extend(x-x.mean(0));yy.extend(y-y.mean())
    x=np.asarray(xx);y=np.asarray(yy)
    return np.linalg.solve(x.T@x+np.eye(2),x.T@y)
def held_predictions(records,adjusted=False):
    result=[]
    for h in records:
        f,p,rep=h["focal"],h["pair"],h["rep"]
        train=[r for r in records if r["rep"]!=rep]
        other=[r for r in train if r["pair"]!=p]
        # Even in sensitivity, fit slope without held pair or label.
        slope=residual_slopes(other) if adjusted else np.zeros(2)
        def D(r):return float(r["D"]-slope@np.array([r["delta_total"],r["delta_composition"]]))
        def group_mean(g):return mean([D(r) for r in g])
        prior=[];X=[]
        for q in ["-".join(v) for v in PAIRS]:
            if q==p:continue
            prior.append(group_mean([r for r in other if r["focal"]==f and r["pair"]==q]))
            X.append([1,*features(q)])
        X=np.array(X);y=np.array(prior);pen=np.diag([0,1,1,1,1,1])
        beta=np.linalg.solve(X.T@X+pen,X.T@y)
        models={"zero":0.,"other_pair_mean":mean(prior),
                "known_pair_mean":group_mean([r for r in train if group(r)==(f,p)]),
                "other_focals_same_pair":mean([group_mean([r for r in train if group(r)==(i,p)]) for i in FOCAL if i!=f]),
                "held_pair_additive_identity":float(np.array([1,*features(p)])@beta)}
        target=D(h)
        result.append({"focal":f,"pair":p,"rep":rep,"target":target,**{m:abs(target-v) for m,v in models.items()}})
    return result
def summarize(rows):
    def agg(rows):
        by=defaultdict(list)
        for r in rows:by[group(r)].append(r)
        return {"n_pairs":len(rows),"n_groups":len(by),"MAE":{m:mean([mean([r[m] for r in rr]) for rr in by.values()]) for m in MODELS}}
    out=agg(rows)
    out["by_focal"]={f:agg([r for r in rows if r["focal"]==f]) for f in FOCAL}
    out["by_pair"]={p:agg([r for r in rows if r["pair"]==p]) for p in sorted(set(r["pair"] for r in rows))}
    return out
def describe(records):
    g=defaultdict(list)
    for r in records:g[group(r)].append(r)
    means=[]
    for (f,p),rr in sorted(g.items()):
        means.append({"focal":f,"pair":p,"n":len(rr),"mean_D":mean([r["D"] for r in rr]),
         "mean_split_seeds":mean([r["split_seeds"] for r in rr]),"mean_mix_seeds":mean([r["mix_seeds"] for r in rr]),
         "mean_delta_total_count":50*mean([r["delta_total"] for r in rr]),"mean_delta_composition_count":50*mean([r["delta_composition"] for r in rr])})
    signs=[]
    for p in sorted(set(r["pair"] for r in means)):
        vals=[r["mean_D"] for r in means if r["pair"]==p]
        signs.append({"pair":p,"signs":list(map(int,np.sign(vals))),"all_same_nonzero":bool(len(set(np.sign(vals)))==1 and vals[0]!=0)})
    return {"group_means":means,"positive_groups":sum(r["mean_D"]>0 for r in means),"negative_groups":sum(r["mean_D"]<0 for r in means),
        "equal_group_mean_D":mean([r["mean_D"] for r in means]),"all_focal_same_direction_pairs":sum(r["all_same_nonzero"] for r in signs),"signs_by_pair":signs}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    rows,audit=load()
    raw=held_predictions(rows);adj=held_predictions(rows,True)
    z={"status":"SPATIAL_CONTEXT_TRANSFER_NOT_SPECIES_ADDITION_NOT_META","audit":audit,"descriptive":describe(rows),
       "forecasts":{"raw":summarize(raw),"density_standardized_sensitivity":summarize(adj)}}
    a.out.mkdir(parents=True,exist_ok=True)
    for name,r in (("matched_plot_pairs",rows),("held_errors_raw",raw),("held_errors_density_sensitivity",adj)):
        with (a.out/(name+".csv")).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(r[0]));w.writeheader();w.writerows(r)
    (a.out/"summary.json").write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("GIBBS_SPATIAL_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
