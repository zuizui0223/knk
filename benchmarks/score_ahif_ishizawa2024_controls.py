"""Ishizawa equal-information controls for low-order prediction."""
from __future__ import annotations
import argparse,itertools,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_ishizawa2024_interaction import (
    raw,deposited,predict,profile,score_errors,read_csv
)
MODELS=("mono","pair_mean","pair","trio_mean","pair_on_trio_rhs",
        "global_trio_coeff","trio","additive_log","inclusion_exclusion")
CONTRASTS=(("mono","pair"),("pair","trio"),("pair_mean","trio_mean"),
           ("trio_mean","trio"),("pair_on_trio_rhs","trio"),
           ("global_trio_coeff","trio"),("additive_log","trio"),
           ("inclusion_exclusion","trio"))

def subsets(S,k):
    return [frozenset(v) for v in itertools.combinations(sorted(S),k)]

def rhs(m,S,k):
    ss=subsets(S,k)
    return np.array([np.mean([m[T,i] for T in ss if i in T]) for i in sorted(S)])

def coefficient(m,T,i,j):
    den=m[T,j]
    if abs(den)<1e-12:raise ValueError("zero log denominator")
    return (m[T,i]-m[T-{j},i])/den

def trio_control(m,S,pool,mode):
    ids=sorted(S);n=len(S);A=np.eye(n);b=rhs(m,S,3)
    for a,i in enumerate(ids):
        for c,j in enumerate(ids):
            if a==c:continue
            if mode=="pair_on_trio_rhs":
                v=coefficient(m,frozenset([i,j]),i,j)
            elif mode=="global_trio_coeff":
                v=np.mean([coefficient(m,frozenset([i,j,k]),i,j)
                           for k in sorted(set(pool)-{i,j})])
            else:raise ValueError("unknown control")
            A[a,c]=-(n-3)/(n-1)*v
    cond=float(np.linalg.cond(A))
    if not np.isfinite(cond) or cond>1e9:raise ValueError("singular control")
    return np.linalg.solve(A,b),cond

def additive_fit(m,pool):
    fits={}
    for i in sorted(pool):
        partners=sorted(set(pool)-{i})
        cells=[(S,v) for (S,j),v in m.items() if j==i]
        if len(cells)!=22:raise ValueError("incomplete additive training set")
        X=np.array([[1.]+[float(k in S) for k in partners] for S,v in cells])
        y=np.array([v for S,v in cells])
        b,_,rank,_=np.linalg.lstsq(X,y,rcond=None)
        if rank!=7:raise ValueError("additive design rank")
        fits[i]=(partners,b)
    return fits

def additive_predict(fits,S):
    return np.array([float(np.array([1.]+[float(k in S) for k in fits[i][0]])@fits[i][1])
                     for i in sorted(S)])

def inclusion_exclusion(m,S):
    out=[]
    for i in sorted(S):
        partners=sorted(S-{i});a=m[frozenset([i]),i]
        b={j:m[frozenset([i,j]),i]-a for j in partners}
        c=[m[frozenset([i,j,k]),i]-a-b[j]-b[k]
           for j,k in itertools.combinations(partners,2)]
        out.append(a+sum(b.values())+sum(c))
    return np.array(out)

def replicate_audit(source,systems):
    rows,fields=read_csv(source/"data1_rawcfu.csv")
    lookup={v:k for k,v in systems.items()}
    groups=defaultdict(list);triples=Counter()
    for row in rows:
        s,rep,i=row["system"],row["rep"],row["strain"]
        if not rep:raise ValueError("missing replicate label")
        groups[s,rep].append(i);triples[s,rep,i]+=1
    dup=sum(n-1 for n in triples.values() if n>1)
    bad=sum(set(v)!=set(lookup[s]) or len(v)!=len(lookup[s]) for (s,rep),v in groups.items())
    return {"n_system_rep_label_groups":len(groups),"duplicate_system_rep_strain_rows":dup,
      "incomplete_or_duplicated_roster_groups":bad,"rep_label_groups_by_size":dict(Counter(len(lookup[s]) for s,rep in groups)),
      "physical_flask_pairing_independently_verified":False,"replicate_inference_used":False}

def summarize(rows):
    keys=list(rows[0]["errors"]["mono"])
    means={m:{k:float(np.mean([r["errors"][m][k] for r in rows])) for k in keys} for m in MODELS}
    comparisons=[]
    for b,e in CONTRASTS:
        dif={k:means[b][k]-means[e][k] for k in keys}
        if not np.isclose(dif["named_TV"],dif["rank_TV"]+dif["identity_assignment_TV"],atol=1e-12):
            raise ValueError("TV decomposition")
        counts={}
        for k in keys:
            d=np.array([r["errors"][b][k]-r["errors"][e][k] for r in rows])
            counts[k]={"improved":int((d>1e-12).sum()),"worsened":int((d<-1e-12).sum()),"tied":int((abs(d)<=1e-12).sum())}
        comparisons.append({"baseline":b,"enhanced":e,"absolute_improvement":dif,
           "relative_gain":{k:dif[k]/means[b][k] if means[b][k]>1e-12 else None for k in keys},
           "counts":counts})
    return {"n_targets":len(rows),"mean_errors":means,"comparisons":comparisons}

def run(source):
    sp,systems,m,audit=raw(source)
    dp=deposited(source,"data7_abundance_prediction_pair.csv",systems,m)
    dt=deposited(source,"data6_abundance_prediction_trio.csv",systems,m)
    small={q:v for q,v in m.items() if len(q[0])<=3}
    pair={q:v for q,v in small.items() if len(q[0])<=2}
    fits=additive_fit(small,sp)
    rep=replicate_audit(source,systems);rows=[];maxdiff={"pair":0.,"trio":0.}
    for S,system in sorted(systems.items(),key=lambda x:(len(x[0]),sorted(x[0]))):
        if len(S)<4:continue
        ids=sorted(S);o=np.array([m[S,i] for i in ids])
        p,cp=predict(pair,S,2);t,ct=predict(small,S,3)
        pc,cc=trio_control(small,S,sp,"pair_on_trio_rhs")
        g,cg=trio_control(small,S,sp,"global_trio_coeff")
        if len(S)==7 and not np.allclose(g,t,atol=1e-12,rtol=0):
            raise ValueError("full pool global/local coefficient mismatch")
        for k,v,ref in (("pair",p,dp),("trio",t,dt)):
            maxdiff[k]=max(maxdiff[k],float(np.max(abs(v-np.array([ref[S,i] for i in ids])))))
        preds={"mono":rhs(small,S,1),"pair_mean":rhs(small,S,2),"pair":p,
               "trio_mean":rhs(small,S,3),"pair_on_trio_rhs":pc,
               "global_trio_coeff":g,"trio":t,"additive_log":additive_predict(fits,S),
               "inclusion_exclusion":inclusion_exclusion(small,S)}
        if any(not np.isfinite(v).all() for v in preds.values()):raise ValueError("nonfinite control prediction")
        ob=profile(o,ids,sp)
        rows.append({"system":system,"size":len(S),"strains":ids,
            "errors":{k:score_errors(ob,profile(v,ids,sp),o,v) for k,v in preds.items()},
            "condition_numbers":{"pair":cp,"trio":ct,"pair_on_trio_rhs":cc,"global_trio_coeff":cg}})
    if max(maxdiff.values())>1e-6:raise ValueError("original point predictions not reproduced")
    if Counter(r["size"] for r in rows)!={4:35,5:21,6:7,7:1}:raise ValueError("heldout target counts")
    by={str(k):summarize([r for r in rows if r["size"]==k]) for k in range(4,8)}
    summary=summarize(rows)
    equal=[]
    for j,c in enumerate(summary["comparisons"]):
        equal.append({"baseline":c["baseline"],"enhanced":c["enhanced"],
            "absolute_improvement":{k:float(np.mean([z["comparisons"][j]["absolute_improvement"][k] for z in by.values()]))
                                      for k in c["absolute_improvement"]}})
    return {"status":"SOURCE_LOCKED_EQUAL_INFORMATION_DESCRIPTIVE_NOT_HISTORY_MEDIATION",
            "audit":audit,"replicate_label_audit":rep,"max_author_point_difference":maxdiff,
            "summary":summary,"by_size":by,"equal_size_weighted_comparisons":equal,
            "targets":rows,"limits":["Same original seven strains, shared assays, no independent target-level inference.",
              "Direct smaller-community means already contain context information.",
              "Controls share RHS/factor or low-order data, not identical parameterization.",
              "Geometric mean CFU profiles, not physical flask-level compositions.",
              "No randomized arrival order or causal mechanism mediation."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run(a.source);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("ISHIZAWA_CONTROLS_RESULT",json.dumps({k:v for k,v in z.items() if k!="targets"},allow_nan=False))
    print("ISHIZAWA_CONTROLS_TARGETS",json.dumps(z["targets"],allow_nan=False))
if __name__=="__main__":main()
