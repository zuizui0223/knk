"""Equal-target bottom-up species versus direct group response forecasts."""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json
from pathlib import Path
import numpy as np
from benchmarks.audit_ahif_seed_synthesis import acquire
from benchmarks.score_ahif_seed_synthesis_opposition import read
from benchmarks.score_ahif_kansas_forecast import EXPERIMENT,grid,forecasts,FIRST,MODELS as SPEC_MODELS
from benchmarks.score_ahif_kansas_sown_roles import read_mapping

MODELS=("zero","previous","cumulative","bottom_history","bottom_state",
        "direct_history","direct_state","direct_history_scaled","direct_state_scaled")
MAP_BLOB="4be704e6086a18f016cf532391d91e193bc06036"
ROSTER_BLOB="d7f1e7b41035f03733e92ae65ed56f8d44d591b8"

def gitblob(path):
    b=path.read_bytes()
    return hashlib.sha1(f"blob {len(b)}\0".encode()+b).hexdigest()

def groupings(mapping,paper,names):
    if gitblob(mapping)!=MAP_BLOB or gitblob(paper)!=ROSTER_BLOB:raise ValueError("frozen membership blobs")
    s,o,u,pre,ex=read_mapping(mapping,names)
    with mapping.open(newline="",encoding="utf-8") as f:
        m={r["source_taxon"]:r for r in csv.DictReader(f)}
    with paper.open(newline="",encoding="utf-8") as f:
        pg={r["paper_species"]:r["functional_guild"] for r in csv.DictReader(f)}
    guilds=("AF","C4-PG","PLF","PNLF")
    a=np.column_stack([np.array([m[n]["role"]=="sown_list" and pg[m[n]["paper_taxon"]]==g for n in names]) for g in guilds]+[o,u])
    b=np.column_stack((s,o,u));c=np.column_stack((ex,o,u|(s&~ex)))
    out={"broad3":(["sown_list","other_named","unresolved"],b),
         "exact23":(["sown_exact23","other_named","unresolved"],c),
         "sown_guilds6":([*guilds,"other_named","unresolved"],a)}
    for n,(labels,G) in out.items():
        if not np.all(G.sum(1)==1) or np.any(G.sum(0)==0):raise ValueError("partition coverage "+n)
    if not np.array_equal(a[:,:4].sum(1),s):raise ValueError("guild-to-role alignment")
    return {n:(labels,G.astype(float)) for n,(labels,G) in out.items()}

def features(delta,control,source_control,k,state):
    prev=delta[k-1];avg=delta[:k].mean(0);cols=[prev,avg]
    if state:
        b=control[k-1];raw=source_control[k-1];tot=raw.sum()
        if tot<=0 or not np.isclose(b.sum(),tot):raise ValueError("control aggregation")
        p=raw/tot;c=float(p@p);n=len(b)
        cols.extend((np.log1p(b),b/tot,np.full(n,np.log1p(tot)),
            np.full(n,float(np.count_nonzero(raw))),np.full(n,c),prev*c))
    return np.column_stack(cols)

def fit(x,y,new,penalty):
    mu=x.mean(0);sd=x.std(0);sd=np.where(sd>1e-12,sd,1.)
    a=np.column_stack((np.ones(len(x)),(x-mu)/sd))
    b=np.column_stack((np.ones(len(new)),(new-mu)/sd))
    p=np.eye(a.shape[1])*penalty;p[0,0]=0.
    return b@np.linalg.solve(a.T@a+p,a.T@y)

def direct(a,k,G,state,penalty):
    # a time x arm x species, only past time observations may enter forecasts.
    x=a@G;control=x[:,0];delta=x[:,1]-control
    feat=np.concatenate([features(delta,control,a[:,0],j,state) for j in range(1,k)],axis=0)
    pred=fit(feat,delta[1:k].reshape(-1),features(delta,control,a[:,0],k,state),penalty)
    return np.where(np.any(x[:k]>0,axis=(0,1)),pred,0.)

def model_predictions(a,k,G,sp):
    ans={"zero":sp["zero"]@G,"previous":sp["previous"]@G,"cumulative":sp["cumulative"]@G,
         "bottom_history":sp["history_ridge"]@G,"bottom_state":sp["state_ridge"]@G}
    ag=a@G;d=ag[:,1]-ag[:,0]
    if not np.allclose(ans["previous"],d[k-1]) or not np.allclose(ans["cumulative"],d[:k].mean(0)):
        raise ValueError("simple aggregation invariance")
    for mode,state in (("history",False),("state",True)):
        ans["direct_"+mode]=direct(a,k,G,state,10.)
        ans["direct_"+mode+"_scaled"]=direct(a,k,G,state,10.*G.shape[1]/G.shape[0])
    if any(not np.isfinite(v).all() for v in ans.values()):raise ValueError("nonfinite forecast")
    return ans

def summarize(rows,labels,broad):
    obs=np.array([r["observed"] for r in rows]);z={}
    if not len(rows):raise ValueError("empty score")
    for m in MODELS:
        p=np.array([r["predictions"][m] for r in rows]);err=(p-obs)**2
        a={"sse":float(err.sum()),"sse_by_group":dict(zip(labels,err.sum(0).tolist())),
           "mae_by_group":dict(zip(labels,np.abs(p-obs).mean(0).tolist())),
           "total_sse":float(np.sum((p.sum(1)-obs.sum(1))**2))}
        if broad:
            truth=(obs[:,0]>0)&(obs[:,1]<0);pred=(p[:,0]>0)&(p[:,1]<0)
            tp=int(np.sum(truth&pred));tn=int(np.sum(~truth&~pred))
            fp=int(np.sum(~truth&pred));fn=int(np.sum(truth&~pred))
            a["opposition_forecast"]={"tp":tp,"fp":fp,"tn":tn,"fn":fn,
               "accuracy":float((tp+tn)/len(obs)),
               "balanced_accuracy":float(.5*(tp/(tp+fn)+tn/(tn+fp))) if (tp+fn)*(tn+fp)>0 else None,
               "sown_sign_accuracy":float(np.mean(np.sign(p[:,0])==np.sign(obs[:,0]))),
               "other_sign_accuracy":float(np.mean(np.sign(p[:,1])==np.sign(obs[:,1])))}
        z[m]=a
    def gain(b,e):return (z[b]["sse"]-z[e]["sse"])/z[b]["sse"] if z[b]["sse"]>1e-15 else None
    return {"n_cells":len(rows),"models":z,"gains":{
       "direct_history_over_bottom_history":gain("bottom_history","direct_history"),
       "direct_state_over_bottom_state":gain("bottom_state","direct_state"),
       "direct_history_scaled_over_bottom_history":gain("bottom_history","direct_history_scaled"),
       "direct_state_scaled_over_bottom_state":gain("bottom_state","direct_state_scaled"),
       "direct_history_over_zero":gain("zero","direct_history"),
       "direct_state_over_zero":gain("zero","direct_state"),
       "direct_state_over_direct_history":gain("direct_history","direct_state")}}

def evaluate(x,blocks,groups,prior=None):
    result={n:[] for n in (*groups,"guilds_to_broad3")}
    prioridx={}
    if prior:
        for c in prior["all_forecast_cells"]:
            key=(c["age"],tuple(tuple(b) for b in c["train_blocks"]))
            if key in prioridx:raise ValueError("duplicate prior scoring cell")
            prioridx[key]=c
    maxdiff=0.
    for tr0 in itertools.combinations(range(4),2):
        tr=list(tr0);te=[i for i in range(4) if i not in tr]
        a=x[:,tr].mean(1);b=x[:,te].mean(1)
        d=a[:,1]-a[:,0];obs=b[:,1]-b[:,0]
        for k in range(FIRST,len(x)):
            seen=np.any(a[:k]>0,axis=(0,1));sp=forecasts(a[:,0],d,seen,k)
            if prior:
                q=prioridx[k+1,tuple(tuple(blocks[i]) for i in tr)]
                if q["test_blocks"]!=[list(blocks[i]) for i in te]:raise ValueError("prior fold mismatch")
                for m in SPEC_MODELS:
                    old=q["scores"][m]["SSE"];new=float(np.sum((sp[m]-obs[k])**2))
                    maxdiff=max(maxdiff,abs(old-new))
                    if not np.isclose(old,new,rtol=1e-9,atol=1e-9):raise ValueError("prior species SSE drift")
            for name,(labels,G) in groups.items():
                pred=model_predictions(a,k,G,sp);y=obs[k]@G
                rec={"age":k+1,"train_blocks":[blocks[i] for i in tr],"test_blocks":[blocks[i] for i in te],
                     "observed":y.tolist(),"predictions":{m:v.tolist() for m,v in pred.items()}}
                result[name].append(rec)
                if name=="sown_guilds6":
                    collapse=lambda v:np.array([v[:4].sum(),v[4],v[5]])
                    result["guilds_to_broad3"].append({**rec,"observed":collapse(y).tolist(),
                                   "predictions":{m:collapse(v).tolist() for m,v in pred.items()}})
    out={}
    for name,rows in result.items():
        labels,G=groups["broad3" if name=="guilds_to_broad3" else name]
        broad=name!="sown_guilds6"
        if len(rows)!=60:raise ValueError("forecast cells")
        out[name]={"labels":labels,"group_sizes":G.sum(0).astype(int).tolist(),
          "aggregate":summarize(rows,labels,broad),
          "by_age":[{"age":k,**summarize([r for r in rows if r["age"]==k],labels,broad)} for k in range(6,16)],
          "cells":rows}
    # Comparing broad targets: bottom-up and simple predictions are identical
    # across coarse and fine mappings before direct group fitting.
    a=out["broad3"]["cells"];b=out["guilds_to_broad3"]["cells"]
    for p,q in zip(a,b):
        if p["age"]!=q["age"] or p["train_blocks"]!=q["train_blocks"] or not np.allclose(p["observed"],q["observed"]):raise ValueError("collapsed target drift")
        for m in MODELS[:5]:
            if not np.allclose(p["predictions"][m],q["predictions"][m]):raise ValueError("bottom-up collapse drift")
    return out,maxdiff

def run(root,mapping,paper,prior_path):
    manifest=acquire(root);states,roster,dup,drift=read(root,manifest)
    if dup[EXPERIMENT] or drift[EXPERIMENT]:raise ValueError("physical key ambiguity")
    x,blocks=grid(states[EXPERIMENT]);groups=groupings(mapping,paper,roster[EXPERIMENT])
    prior=json.loads(prior_path.read_text())
    if prior["source_run"]!="https://github.com/zuizui0223/knk/actions/runs/45541282653" or len(prior["all_forecast_cells"])!=60:raise ValueError("prior provenance")
    z,md=evaluate(x,blocks,groups,prior)
    return {"status":"CURATED_KANSAS_SAME_TARGET_ROLE_FORECAST_DESCRIPTIVE",
      "source_manifest":manifest,"experiment":EXPERIMENT,"n_states":120,"n_taxa":100,"n_blocks":4,
      "n_oriented_partitions":6,"n_forecast_ages":10,"prior_species_sse_max_abs_diff":md,
      "mapping_blob":MAP_BLOB,"paper_roster_blob":ROSTER_BLOB,"representations":z,
      "limits":["Probable original paper linkage, original EDI archive unavailable.",
       "Same selected curated source; no causal mechanism or independent study.",
       "Same group targets compared; no comparison of species SSE with group SSE.",
       "Other named and unresolved are membership roles, not functional guilds.",
       "Six overlapping partitions, repeated age codes; descriptive only, no p-values."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--mapping",type=Path,required=True)
    p.add_argument("--paper",type=Path,required=True);p.add_argument("--prior",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True);z=run(a.root,a.mapping,a.paper,a.prior)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("KANSAS_GROUP_FORECAST_SUMMARY",json.dumps({k:v for k,v in z.items() if k!="representations"}|{"representations":{k:{x:y for x,y in v.items() if x!="cells"} for k,v in z["representations"].items()}},allow_nan=False))
if __name__=="__main__":main()
