"""Forecast signed biomass responses across both source time and held blocks."""
from __future__ import annotations
import argparse,itertools,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from benchmarks.audit_ahif_seed_synthesis import acquire
from benchmarks.score_ahif_seed_synthesis_opposition import read

EXPERIMENT="Kansas_KUFS_LTER_Hay_Meadow_Exp_2"
MODELS=("zero","previous","cumulative","history_ridge","state_ridge")
PENALTY=10.
FIRST=5

def grid(states):
    by=defaultdict(list)
    for r in states:
        if r["mode"]!="direct":raise ValueError("not direct mode")
        by[r["age"],r["block"],r["tr"]].append(r["x"])
    ages=sorted({r["age"] for r in states});blocks=sorted({r["block"] for r in states})
    if ages!=list(range(1,16)) or len(blocks)!=4 or len(states)!=120:
        raise ValueError("unexpected direct long-term grain")
    x=np.zeros((15,4,2,100))
    for a,age in enumerate(ages):
        for b,block in enumerate(blocks):
            for k,tr in enumerate(("Control","Seeds")):
                z=by[age,block,tr]
                if len(z)!=1 or len(z[0])!=100:raise ValueError("unbalanced original block grid")
                x[a,b,k]=z[0]
    if not np.isfinite(x).all() or np.any(x<0) or np.any(x.sum(-1)<=0):
        raise ValueError("invalid absolute biomass")
    return x,blocks

def feature(delta,control,k,state):
    # Index k's response is strictly excluded from predictors.
    if k<1:raise ValueError("lag needed")
    prev=delta[k-1];avg=delta[:k].mean(0)
    cols=[prev,avg]
    if state:
        b=control[k-1];total=b.sum()
        p=b/total
        simpson=float(p@p)
        cols.extend((np.log1p(b),p,np.full(len(b),np.log1p(total)),
                     np.full(len(b),float(np.count_nonzero(b))),
                     np.full(len(b),simpson),prev*simpson))
    return np.column_stack(cols)

def ridge(x,y,new):
    mu=x.mean(0);sd=x.std(0);sd=np.where(sd>1e-12,sd,1.)
    a=np.column_stack((np.ones(len(x)),(x-mu)/sd))
    n=np.column_stack((np.ones(len(new)),(new-mu)/sd))
    pen=np.eye(a.shape[1])*PENALTY;pen[0,0]=0.
    return n@np.linalg.solve(a.T@a+pen,a.T@y)

def forecasts(control,delta,seen,k):
    result={"zero":np.zeros(delta.shape[1]),"previous":delta[k-1].copy(),
            "cumulative":delta[:k].mean(0)}
    for name,state in (("history_ridge",False),("state_ridge",True)):
        x=np.concatenate([feature(delta,control,j,state) for j in range(1,k)],axis=0)
        y=delta[1:k].reshape(-1)
        result[name]=ridge(x,y,feature(delta,control,k,state))
    for name in result:result[name]=np.where(seen,result[name],0.)
    return result

def signscore(pred,obs):
    u=float(obs[pred>0].sum());d=float(-obs[pred<0].sum())
    return {"held_increase":u,"held_decrease":d,"held_min":min(u,d),
            "both_positive":bool(u>0 and d>0)}

def summary(rows):
    if not rows:raise ValueError("empty summary")
    out={}
    for m in (*MODELS,"concurrent_sign_reference"):
        v=[r["scores"][m] for r in rows]
        U=float(np.mean([a["held_increase"] for a in v]));D=float(np.mean([a["held_decrease"] for a in v]))
        mn=float(np.mean([a["held_min"] for a in v]))
        out[m]={"mean_held_increase":U,"mean_held_decrease":D,
                "mean_held_min":mn,"fraction_both_positive":float(np.mean([a["both_positive"] for a in v])),
                "strong_descriptive_gate":bool(U>0 and D>0 and mn>0)}
        if m in MODELS:out[m]["SSE"]=float(sum(a["SSE"] for a in v))
    b=out["history_ridge"]["SSE"];s=out["state_ridge"]["SSE"]
    g=(b-s)/b if b>1e-15 else None
    return {"n_orientations":len(rows),"models":out,"state_gain_over_history":g,
            "state_gain_over_previous":1-s/out["previous"]["SSE"] if out["previous"]["SSE"]>1e-15 else None,
            "state_gain_over_zero":1-s/out["zero"]["SSE"] if out["zero"]["SSE"]>1e-15 else None}

def evaluate(x,blocks):
    n=x.shape[-1];rows=[];concurrent=[]
    for tr in itertools.combinations(range(4),2):
        tr=list(tr);te=[i for i in range(4) if i not in tr]
        a=x[:,tr].mean(1);b=x[:,te].mean(1)
        control=a[:,0];delta=a[:,1]-a[:,0];held=b[:,1]-b[:,0]
        for k in range(len(x)):
            obs=held[k]
            concurrent.append({"age":k+1,"train_blocks":[blocks[i] for i in tr],
              "test_blocks":[blocks[i] for i in te],**signscore(delta[k],obs)})
            if k<FIRST:continue
            seen=np.any(a[:k]>0,axis=(0,1))
            pred=forecasts(control,delta,seen,k)
            scores={m:{"SSE":float(np.sum((v-obs)**2)),**signscore(v,obs)} for m,v in pred.items()}
            scores["concurrent_sign_reference"]=signscore(delta[k],obs)
            p=control[k-1]/control[k-1].sum()
            rows.append({"age":k+1,"train_blocks":[blocks[i] for i in tr],
             "test_blocks":[blocks[i] for i in te],"n_seen_taxa":int(seen.sum()),
             "n_held_nonzero_unseen_taxa":int(np.sum((~seen)&(obs!=0))),
             "lag_training_control":{"total":float(control[k-1].sum()),
                 "richness":int(np.count_nonzero(p)),"simpson":float(p@p),"top1":float(p.max())},
             "scores":scores})
    by=[]
    for age in range(FIRST+1,len(x)+1):
        s=summary([r for r in rows if r["age"]==age])
        by.append({"age":age,**s})
    exact=[]
    for age in range(1,len(x)+1):
        v=[a for a in concurrent if a["age"]==age]
        U=float(np.mean([a["held_increase"] for a in v]));D=float(np.mean([a["held_decrease"] for a in v]))
        mn=float(np.mean([a["held_min"] for a in v]))
        exact.append({"age":age,"mean_held_increase":U,"mean_held_decrease":D,
                      "mean_held_min":mn,"fraction_both_positive":float(np.mean([a["both_positive"] for a in v])),
                      "strong_descriptive_gate":bool(U>0 and D>0 and mn>0)})
    return {"aggregate_forecast":summary(rows),"by_age":by,"all_forecast_cells":rows,"exact_concurrent":exact}

def run(root):
    manifest=acquire(root);states,roster,dup,drift=read(root,manifest)
    if dup[EXPERIMENT] or drift[EXPERIMENT]:raise ValueError("physical identity ambiguity")
    if len(roster[EXPERIMENT])!=100:raise ValueError("source roster")
    x,blocks=grid(states[EXPERIMENT])
    return {"status":"CURATED_DIRECT_KANSAS_TIME_BLOCK_FORECAST_DESCRIPTIVE",
      "source_manifest":manifest,"experiment":EXPERIMENT,"n_states":120,"n_taxa":100,
      "n_blocks":4,"source_age_codes":list(range(1,16)),"n_oriented_partitions":6,
      "n_forecast_times":10,"ridge_penalty":PENALTY,**evaluate(x,blocks),
      "limits":["Same selected curated experiment; no independent replication or causal mechanism.",
       "Six overlapping block partitions and ten repeated age codes not independent.",
       "No held/current data in forecasts; fixed retrospective taxon roster.",
       "Lagged control is not pretreatment state; source age codes unverified durations.",
       "No p-values, no source-specific seed identities, not a formal meta-analysis."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True);z=run(a.root)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("KANSAS_FORECAST_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
