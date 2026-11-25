"""Verify both learned directions on held blocks, not unstable total signs."""
from __future__ import annotations
import argparse,hashlib,json
from collections import defaultdict
from pathlib import Path
import numpy as np
from benchmarks.audit_ahif_seed_synthesis import acquire
from benchmarks.score_ahif_seed_synthesis_opposition import read,score,SPLITS,SEED

def signed_components(rows):
    n=len(rows)
    if n<4:raise ValueError("at least four block groups")
    x=np.array([r["delta"] for r in rows]);ids=[r["block"] for r in rows]
    if len(set(ids))!=n:raise ValueError("duplicate block group")
    rng=np.random.default_rng(SEED+int(hashlib.sha256(json.dumps(ids).encode()).hexdigest()[:8],16))
    u=[];d=[];a=[];b=[]
    for _ in range(SPLITS):
        p=rng.permutation(n)
        for tr,te in ((p[:n//2],p[n//2:]),(p[n//2:],p[:n//2])):
            train=x[tr].mean(0);test=x[te].mean(0)
            pos=train>0;neg=train<0
            u.append(float(test[pos].sum()));d.append(float(-test[neg].sum()))
            a.append(int(pos.sum()));b.append(int(neg.sum()))
    U=np.array(u);D=np.array(d);A=np.array(a);B=np.array(b)
    z={"n_blocks":n,"n_twofold_splits":SPLITS,"mean_held_increase":float(U.mean()),
       "mean_held_decrease":float(D.mean()),"mean_held_min_signed_component":float(np.minimum(U,D).mean()),
       "fraction_orientations_both_positive":float(np.mean((U>0)&(D>0))),
       "fraction_train_both_sets_nonempty":float(np.mean((A>0)&(B>0))),
       "mean_n_learned_increasing_taxa":float(A.mean()),"mean_n_learned_decreasing_taxa":float(B.mean()),
       "mean_held_signed_gross":float((U+D).mean())}
    z["strong_direction_descriptive_gate"]=bool(z["mean_held_increase"]>0 and z["mean_held_decrease"]>0 and z["mean_held_min_signed_component"]>0)
    return z

def reconstruct(ss,age,mode):
    by=defaultdict(list)
    for v in ss:
        if v["age"]==age and v["mode"]==mode:by[v["block"]].append(v)
    rs=[]
    for block,vs in sorted(by.items()):
        c=[v["x"] for v in vs if v["tr"]=="Control"];t=[v["x"] for v in vs if v["tr"]=="Seeds"]
        if c and t:rs.append({"block":block,"delta":np.mean(t,axis=0)-np.mean(c,axis=0)})
    return rs

def run(root):
    manifest=acquire(root);states,roster,dup,drift=read(root,manifest);base=score(states,roster,dup,drift)
    for e in base:
        for row in e["comparisons"]:
            if row["status"]!="SCORED_DESCRIPTIVE":continue
            cf=row["crossfit"]
            if cf["status"]!="DESCRIPTIVE_HELD_BLOCK_SIGN_SCORE":
                row["directional_validation"]={"status":cf["status"]};continue
            rs=reconstruct(states[e["experiment"]],row["age"],row["mode"])
            if len(rs)!=row["n_matched_strata"]:raise ValueError("matched blocks changed")
            z=signed_components(rs)
            if not np.isclose(z["mean_held_signed_gross"],cf["signed_gross"],rtol=1e-9,atol=1e-9):raise ValueError("split signed gross mismatch")
            row["directional_validation"]={"status":"DESCRIPTIVE_TWO_DIRECTION_HELD_BLOCKS",**z}
    return {"status":"CURATED_SYNTHESIS_TWO_DIRECTION_DIAGNOSTIC_NOT_CAUSAL","source_manifest":manifest,
      "n_experiment_labels":len(base),"n_plot_states":sum(x["n_states"] for x in base),"experiments":base,
      "limits":["Same curated source; not additional independent evidence.","Repeated split means and fractions descriptive, not p-values.",
                "Direct and allocated modes distinct; no seed identity or causal recruitment compensation."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True);z=run(a.root)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SEED_DIRECTION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
