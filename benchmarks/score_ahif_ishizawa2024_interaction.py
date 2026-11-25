"""Source-locked low-order assay prediction of held-out multi-strain composition and rank."""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np

FILES={"data1_rawcfu.csv":"b4e24d2348b7a90e63b0fa96a0eaa38e",
"data6_abundance_prediction_trio.csv":"c5fcddc05a46da6476cc7281230e3054",
"data7_abundance_prediction_pair.csv":"360be5f5ea1971d436ea5907a2d07ee8",
"code2_interactions.R":"047704c3b4f977b92d4ab2f6367f1768",
"code4_abundance_prediction_trio.R":"fe6986e58217d7cb110c15f66e3b1d08",
"code5_abundance_prediction_pair.R":"0ab05826f92c5c5ce6458aa333a7d473",
"functions.R":"9c2e41dbfca1d69dc625cb8e5428f684"}

def read_csv(path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f)
        return list(r),r.fieldnames

def raw(source):
    for name,h in FILES.items():
        if hashlib.md5((source/name).read_bytes()).hexdigest()!=h:
            raise ValueError("original source digest mismatch: "+name)
    rows,fields=read_csv(source/"data1_rawcfu.csv")
    print("ISHIZAWA_RAW_SCHEMA",json.dumps({"fields":fields,"n_rows":len(rows),"head":rows[:2]}))
    if not {"strain","system","nspecies","cfu"}.issubset(fields):
        raise ValueError("raw source fields differ")
    cells=defaultdict(list);members=defaultdict(set);nsize=defaultdict(set)
    for z in rows:
        s=z["system"];i=z["strain"];n=int(z["nspecies"]);v=float(z["cfu"])
        if not math.isfinite(v) or v<=0:raise ValueError("nonpositive/nonfinite original CFU")
        cells[s,i].append(math.log10(v));members[s].add(i);nsize[s].add(n)
    species=sorted({i for s,i in cells})
    if len(species)!=7 or len(members)!=127 or len(cells)!=448:
        raise ValueError("original factorial grain differs")
    systems={}
    for s,ids in members.items():
        key=frozenset(ids)
        if key in systems or nsize[s]!={len(key)}:
            raise ValueError("duplicate or inconsistent membership")
        systems[key]=s
        counts=[len(cells[s,i]) for i in ids]
        if len(set(counts))!=1 or not 5<=counts[0]<=12:
            raise ValueError("source per-strain sample count differs")
    if set(systems)!={frozenset(x) for k in range(1,8) for x in itertools.combinations(species,k)}:
        raise ValueError("not all original subsets")
    means={(S,i):float(np.mean(cells[s,i])) for S,s in systems.items() for i in S}
    audit={"fields":fields,"raw_rows":len(rows),"strain_subset_cells":len(cells),
        "source_strains":species,"subset_count_by_size":dict(Counter(map(len,systems))),
        "replicate_count_by_size":{str(k):dict(Counter(len(cells[s,next(iter(S))]) for S,s in systems.items() if len(S)==k)) for k in range(1,8)},
        "source_columns_other_than_required":sorted(set(fields)-{"strain","system","nspecies","cfu"}),
        "physical_flask_covariance_not_used":True}
    return species,systems,means,audit

def predict(means,S,k):
    ids=sorted(S);n=len(ids)
    if k not in (2,3) or n<=k:raise ValueError("training order/target invalid")
    subsets=[frozenset(x) for x in itertools.combinations(ids,k)]
    A=np.eye(n);b=np.zeros(n)
    for i,a in enumerate(ids):
        b[i]=np.mean([means[T,a] for T in subsets if a in T])
        for j,c in enumerate(ids):
            if i==j:continue
            cs=[]
            for T in subsets:
                if not {a,c}.issubset(T):continue
                den=means[T,c]
                if abs(den)<1e-12:raise ValueError("zero log denominator")
                cs.append((means[T,a]-means[T-{c},a])/den)
            A[i,j]=-(n-k)/(n-1)*np.mean(cs)
    cond=float(np.linalg.cond(A))
    if not np.isfinite(cond) or cond>1e9:raise ValueError("near-singular model")
    v=np.linalg.solve(A,b)
    if not np.isfinite(v).all():raise ValueError("nonfinite predicted log CFU")
    return v,cond

def profile(v,ids,species):
    # Stable softmax of log10 CFU, absent strains are genuinely absent by design.
    a=np.power(10.,v-np.max(v));a/=a.sum()
    out=np.zeros(len(species))
    for k,i in enumerate(ids):out[species.index(i)]=a[k]
    return out

def logtotal(v):
    m=float(np.max(v));return m+math.log10(float(np.power(10.,v-m).sum()))

def score_errors(ob,pred,observed_log,pred_log):
    named=float(np.abs(ob-pred).sum()/2)
    rank=float(np.abs(np.sort(ob)-np.sort(pred)).sum()/2)
    gap=named-rank
    if gap<-1e-12:raise ValueError("TV rearrangement inequality")
    return {"named_TV":named,"rank_TV":rank,"identity_assignment_TV":max(gap,0.),
            "simpson_abs_error":abs(float(ob@ob-pred@pred)),
            "top1_abs_error":abs(float(ob.max()-pred.max())),
            "log_sum_geometric_cfu_abs_error":abs(logtotal(observed_log)-logtotal(pred_log))}

def deposited(source,filename,systems,means):
    rows,fields=read_csv(source/filename)
    print("ISHIZAWA_AUTHOR_SCHEMA",json.dumps({"file":filename,"fields":fields,"rows":len(rows)}))
    if not {"strain","system","nspecies","obs","pred"}.issubset(fields):
        raise ValueError("published prediction columns mismatch")
    lookup={v:k for k,v in systems.items()};out={}
    for z in rows:
        S=lookup[z["system"]];i=z["strain"];key=(S,i)
        if key in out or i not in S or int(z["nspecies"])!=len(S):
            raise ValueError("duplicate or wrong published target")
        o,p=float(z["obs"]),float(z["pred"])
        if not np.isfinite([o,p]).all() or abs(o-means[key])>1e-6:
            raise ValueError("published observation mismatch")
        out[key]=p
    if set(out)!=set(means):raise ValueError("published point table coverage")
    return out

def summarize(rows):
    keys=list(rows[0]["errors"]["mono"])
    mean={m:{k:float(np.mean([r["errors"][m][k] for r in rows])) for k in keys}
          for m in ("mono","pair","trio")}
    dif={k:mean["pair"][k]-mean["trio"][k] for k in keys}
    if abs(dif["named_TV"]-dif["rank_TV"]-dif["identity_assignment_TV"])>1e-12:
        raise ValueError("improvement decomposition")
    rel={k:(dif[k]/mean["pair"][k] if mean["pair"][k]>1e-12 else None) for k in keys}
    counts={}
    for k in keys:
        v=np.asarray([r["errors"]["pair"][k]-r["errors"]["trio"][k] for r in rows])
        counts[k]={"improved":int((v>1e-12).sum()),"worsened":int((v<-1e-12).sum()),"tied":int((abs(v)<=1e-12).sum())}
    return {"n_target_subsets":len(rows),"mean_errors":mean,
            "pair_minus_trio_absolute_improvement":dif,"pair_to_trio_relative_gain":rel,
            "improvement_counts":counts}

def run(source):
    sp,systems,m,audit=raw(source)
    dp=deposited(source,"data7_abundance_prediction_pair.csv",systems,m)
    dt=deposited(source,"data6_abundance_prediction_trio.csv",systems,m)
    small={q:v for q,v in m.items() if len(q[0])<=3}
    pair={q:v for q,v in small.items() if len(q[0])<=2}
    rows=[];maxdiff={"pair":0.,"trio":0.}
    for S,s in sorted(systems.items(),key=lambda x:(len(x[0]),sorted(x[0]))):
        if len(S)<4:continue
        ids=sorted(S);o=np.array([m[S,i] for i in ids])
        mono=np.array([small[frozenset([i]),i] for i in ids])
        p,cp=predict(pair,S,2);t,ct=predict(small,S,3)
        for k,v,ref in (("pair",p,dp),("trio",t,dt)):
            maxdiff[k]=max(maxdiff[k],float(np.max(abs(v-np.array([ref[S,i] for i in ids])))))
        ob=profile(o,ids,sp)
        errors={k:score_errors(ob,profile(v,ids,sp),o,v) for k,v in (("mono",mono),("pair",p),("trio",t))}
        rows.append({"system":s,"size":len(S),"strains":ids,"errors":errors,"condition_numbers":{"pair":cp,"trio":ct}})
    if max(maxdiff.values())>1e-6:
        print("ISHIZAWA_POINT_MISMATCH",json.dumps(maxdiff))
        raise ValueError("author predictions not reproduced within log tolerance")
    if Counter(r["size"] for r in rows)!={4:35,5:21,6:7,7:1}:raise ValueError("held-target counts")
    by={str(k):summarize([r for r in rows if r["size"]==k]) for k in range(4,8)}
    eql={k:float(np.mean([b["pair_minus_trio_absolute_improvement"][k] for b in by.values()]))
         for k in rows[0]["errors"]["mono"]}
    return {"status":"SOURCE_LOCKED_LOW_ORDER_ASSAY_PREDICTION_NOT_HISTORY_MEDIATION",
            "source_version":8433832,"source_md5":FILES,"audit":audit,
            "max_abs_log_point_prediction_difference_from_authors":maxdiff,
            "training_subset_counts":{"mono":7,"pair":21,"trio":35},
            "summary":summarize(rows),"by_size":by,"equal_size_weighted_improvement":eql,
            "targets":rows,
            "limits":["All subsets share seven strains and training assays; finite combinatorial descriptive prediction only.",
                      "No main randomized arrival-order treatment; not history-specific mediation.",
                      "TV identity component is mathematical residual, not a directly observed process.",
                      "Observed profiles use geometric mean strain CFU, not average within-flask composition.",
                      "Original article already demonstrated improvement of multi-strain prediction with trio information."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run(a.source);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("ISHIZAWA_INTERACTION_RESULT",json.dumps({k:v for k,v in z.items() if k!="targets"},allow_nan=False))
    print("ISHIZAWA_TARGET_ROWS",json.dumps(z["targets"],allow_nan=False))
if __name__=="__main__":main()
