"""Replay persistence and label/rank representation prediction diagnostics."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_replay2025_state import (
    BLOBS,METADATA,MATRIX,gitblob,get_meta,read_matrix
)
LAMBDAS=(.1,1.,10.)
MODES={"summary":(),"summary_rank":("rank",),"summary_named":("named",),
       "summary_rank_named":("rank","named")}

def loss(y,p):
    return float(np.sum((y-p)**2))

def gain(b,e):
    return (b-e)/b if b>1e-15 else 0.

def responses(x):
    p=x/x.sum(2,keepdims=True)
    b=p[:,0]
    f=np.column_stack((np.log1p(x[:,0].sum(1)),np.log1p((x[:,0]>0).sum(1)),np.sum(b*b,1),b.max(1)))
    q=np.sort(p,axis=2)[:,:,::-1]
    k=int(np.max(np.sum(x>0,axis=2)))
    return f,{"named":(b,p[:,1:].mean(1)),
              "rank":(q[:,0,:k],q[:,1:,:k].mean(1))}

def covariates(f,tr):
    mu=f[tr].mean(0);sd=f[tr].std(0)
    z=(f-mu)/np.where(sd>1e-12,sd,1.)
    base=np.column_stack((np.ones(len(z)),z[:,0],z[:,0]**2))
    return base,np.column_stack((base,z[:,1:]))

def fitmap(a,tr):
    penalty=np.ones(a.shape[1]);penalty[0]=0.
    return np.linalg.solve(a[tr].T@a[tr]+np.diag(penalty),a[tr].T)

def persistence(a,tr,te,x,y):
    m=fitmap(a,tr)
    bx=m@x[tr];by=m@y[tr]
    rx=x[tr]-a[tr]@bx;ry=y[tr]-a[tr]@by
    den=float(np.sum(rx*rx))
    g=float(np.clip(np.sum(rx*ry)/den,0,1)) if den>1e-15 else 0.
    base=a[te]@by
    pred=base+g*(x[te]-a[te]@bx)
    return base,pred,g

def center_kernel(K,tr):
    # all rows centered relative to the training distribution ONLY
    mu=K[:,tr].mean(1);gm=float(mu[tr].mean())
    C=K-mu[:,None]-mu[None,:]+gm
    scale=float(np.diag(C)[tr].mean())
    if scale<=1e-15:return np.zeros_like(C)
    return C/scale

def build_kernels(summary,representations,tr):
    z=summary[:,1:]
    kernels={"summary":center_kernel(z@z.T,tr)}
    for k,(x,y) in representations.items():
        b=np.sqrt(x)
        kernels[k]=center_kernel(b@b.T,tr)
    return kernels

def kernel_predict(K,tr,te,y,lam):
    ym=y[tr].mean(0);z=y[tr]-ym
    # Symmetric train Gram plus fixed ridge.
    A=K[np.ix_(tr,tr)]+lam*np.eye(len(tr))
    return ym+K[np.ix_(te,tr)]@np.linalg.solve(A,z)

def evaluate(f,rep,loc):
    out={k:{"persistence_folds":[],"kernel_folds":[]} for k in rep}
    for site in np.unique(loc):
        tr=np.flatnonzero(loc!=site);te=np.flatnonzero(loc==site)
        a,b=covariates(f,tr);ks=build_kernels(b,rep,tr)
        for k,(x,y) in rep.items():
            p0,pd,gd=persistence(a,tr,te,x,y)
            p1,ps,gs=persistence(b,tr,te,x,y)
            L={"mean":loss(y[te],y[tr].mean(0)),"unchanged_initial":loss(y[te],x[te]),
               "depth":loss(y[te],p0),"summary":loss(y[te],p1),
               "persistence_depth":loss(y[te],pd),"persistence_summary":loss(y[te],ps)}
            out[k]["persistence_folds"].append({"location":str(site),"n_parents":len(te),
              "losses":L,"gamma_depth":gd,"gamma_summary":gs,
              "summary_gain_over_depth_persistence":gain(L["persistence_depth"],L["persistence_summary"])})
            for lam in LAMBDAS:
                scores={}
                for mode,extra in MODES.items():
                    K=ks["summary"].copy()
                    for label in extra:K+=ks[label]
                    scores[mode]=loss(y[te],kernel_predict(K,tr,te,y,lam))
                out[k]["kernel_folds"].append({"location":str(site),"n_parents":len(te),
                                              "lambda":lam,"losses":scores})
    for k,z in out.items():
        rows=z["persistence_folds"]
        L={m:sum(r["losses"][m] for r in rows) for m in rows[0]["losses"]}
        z["persistence_total_losses"]=L
        z["summary_gain_over_depth"]=gain(L["depth"],L["summary"])
        z["summary_gain_over_depth_persistence"]=gain(L["persistence_depth"],L["persistence_summary"])
        z["persistence_gain_over_depth"]=gain(L["depth"],L["persistence_depth"])
        z["summary_persistence_gain_over_mean"]=gain(L["mean"],L["persistence_summary"])
        z["kernel_totals"]=[]
        for lam in LAMBDAS:
            rr=[r for r in z["kernel_folds"] if r["lambda"]==lam]
            l={m:sum(r["losses"][m] for r in rr) for m in MODES}
            z["kernel_totals"].append({"lambda":lam,"losses":l,
               "rank_gain_over_summary":gain(l["summary"],l["summary_rank"]),
               "named_gain_over_summary":gain(l["summary"],l["summary_named"]),
               "named_increment_over_summary_rank":gain(l["summary_rank"],l["summary_rank_named"]),
               "rank_increment_over_summary_named":gain(l["summary_named"],l["summary_rank_named"])})
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    for name,sha in BLOBS.items():
        if gitblob((a.source/name).read_bytes())!=sha:raise ValueError("original Git blob mismatch")
    ids,loc,want=get_meta(a.source/METADATA);x,audit=read_matrix(a.source/MATRIX,want)
    f,rep=responses(x);z=evaluate(f,rep,loc)
    result={"status":"EXPLORATORY_PERSISTENCE_AND_REPRESENTATION_PREDICTION_NOT_CAUSAL",
       "n_parents":275,"n_locations":14,"n_final_replicates":1100,
       "source_matrix":audit,"source_blobs":BLOBS,"results":z,
       "limits":["Same source experiment as initial-state result; no independent new study.",
                 "All models trained on other locations; no heldout outcome tuning.",
                 "Descriptive follow-up losses, no p-values or inference of causal mechanism.",
                 "Kernel comparisons model-specific; predictions unconstrained; read profiles not biomass."]}
    a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("REPLAY_PERSISTENCE_RESULT",json.dumps(result,allow_nan=False))
if __name__=="__main__":main()
