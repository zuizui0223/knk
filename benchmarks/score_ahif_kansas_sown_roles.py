"""Conditional source-paper membership role response on curated direct biomass."""
from __future__ import annotations
import argparse,csv,itertools,json
from pathlib import Path
import numpy as np
from benchmarks.audit_ahif_seed_synthesis import acquire
from benchmarks.score_ahif_seed_synthesis_opposition import read
from benchmarks.score_ahif_kansas_forecast import EXPERIMENT,grid

def read_mapping(path,names):
    with path.open(newline="",encoding="utf-8") as f:r=list(csv.DictReader(f))
    by={x["source_taxon"]:x for x in r}
    if len(r)!=100 or len(by)!=100 or set(by)!=set(names):raise ValueError("fixed mapping source roster")
    q=[by[n] for n in names]
    s=np.array([x["role"]=="sown_list" for x in q]);o=np.array([x["role"]=="other_named" for x in q]);u=np.array([x["role"]=="unresolved" for x in q])
    pre=np.array([x["reported_present_before_sowing"]=="yes" for x in q])
    a=np.array([x["match_basis"]=="punctuation_exact" for x in q])
    if (int(s.sum()),int(o.sum()),int(u.sum()),int(pre.sum()),int(a.sum()))!=(27,62,11,3,23) or np.any(s.astype(int)+o+u!=1) or np.any(pre&~s) or np.any(a&~s):
        raise ValueError("fixed mapping counts")
    return s,o,u,pre,a

def roles(d,s,o,u,pre):
    # bounds vary unresolved taxon membership, not a CI.
    ds=float(d[s].sum());do=float(d[o].sum());du=float(d[u].sum())
    lo=float(np.minimum(d[u],0).sum());hi=float(np.maximum(d[u],0).sum())
    total=float(d.sum())
    if not np.isclose(ds+do+du,total,rtol=1e-9,atol=1e-9):raise ValueError("mass partition")
    return {"sown_list":ds,"other_named":do,"unresolved":du,"total_species_mass":total,
       "sown_preexisting_three":float(d[pre&s].sum()),"sown_other":float(d[s&~pre].sum()),
       "sown_assignment_range":[ds+lo,ds+hi],"other_assignment_range":[do+lo,do+hi],
       "opposite_named_signs":bool(ds>0 and do<0),"robust_to_unresolved":bool(ds+lo>0 and do+hi<0)}

def evaluate(d,s,o,u,pre,exact):
    if d.ndim!=2 or d.shape[0]!=4:raise ValueError("four matched blocks")
    pooled=d.mean(0);primary=roles(pooled,s,o,u,pre)
    sensitivity=roles(pooled,exact,o,u|(s&~exact),pre)
    halves=[]
    for i in itertools.combinations(range(4),2):
        r=roles(d[list(i)].mean(0),s,o,u,pre)
        halves.append({"held_indices":list(i),**r})
    up=np.array([r["sown_list"] for r in halves]);down=np.array([-r["other_named"] for r in halves])
    cv={"n_oriented_held_halves":6,"mean_held_increase":float(up.mean()),"mean_held_decrease":float(down.mean()),
        "mean_held_min":float(np.minimum(up,down).mean()),"fraction_both_positive":float(np.mean((up>0)&(down>0))),
        "fraction_robust_unresolved":float(np.mean([r["robust_to_unresolved"] for r in halves])),
        "strong_descriptive_gate":bool(up.mean()>0 and down.mean()>0 and np.minimum(up,down).mean()>0)}
    return {"pooled":primary,"exact23_sensitivity":sensitivity,"held_half_summary":cv,
        "block_role_contrasts":[roles(v,s,o,u,pre) for v in d],"held_halves":halves}

def run(root,mapping):
    m=acquire(root);states,roster,dup,drift=read(root,m)
    if dup[EXPERIMENT] or drift[EXPERIMENT]:raise ValueError("physical key ambiguity")
    x,blocks=grid(states[EXPERIMENT]);names=roster[EXPERIMENT]
    s,o,u,p,a=read_mapping(mapping,names);d=x[:,:,1]-x[:,:,0]
    rows=[{"age":k+1,**evaluate(d[k],s,o,u,p,a)} for k in range(15)]
    return {"status":"CURATED_CONDITIONAL_PUBLISHED_SOWN_ROLE_OPPOSITION_NOT_ORIGINAL_OR_CAUSAL",
       "article_doi":"10.1111/ele.13566","membership_paper_doi":"10.1890/08-0849.1","cited_original_doi":"10.6073/pasta/5703e8819c81da0beed4660f467c288e",
       "source_manifest":m,"experiment":EXPERIMENT,"n_states":120,"n_blocks":4,"n_source_ages":15,
       "mapping_counts":{"sown_list":27,"other_named":62,"unresolved":11,"preexisting_sown":3,"exact_matches":23},
       "by_age":rows,"all_age_mean":evaluate(d.mean(0),s,o,u,p,a),
       "limits":["Probable paper-to-curated lineage, original EDI archive and plot subset not independently verified.",
       "Sown-list membership not individual recruitment identity; other named provisional.",
       "Unresolved assignment bounds are not uncertainty intervals, no p-values; repeated age codes and held halves dependent.",
       "No independent causal mediation, treatment-by-fertility/hay comparison or new study."]}

def main():
    a=argparse.ArgumentParser();a.add_argument("--root",type=Path,required=True);a.add_argument("--mapping",type=Path,required=True);a.add_argument("--out",type=Path,required=True);args=a.parse_args()
    args.root.mkdir(parents=True,exist_ok=True);z=run(args.root,args.mapping);args.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("KANSAS_SOWN_ROLE_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
