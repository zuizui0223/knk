#!/usr/bin/env python3
"""Source-authenticated chemical-legacy process/output contrast (Delory et al).

Published primary source: https://doi.org/10.1111/pce.13999
Raw archive: https://github.com/BenjaminDelory/Rhizometabolome_2020/tree/v1.0.0
FROZEN before analysis:
docs/DELORY_CHEMICAL_LEGACY_EXTERNAL_PRE_SCORE_20261124.md
This is an independent source reanalysis of already published results.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

SOURCE_BLOB_SHA="d11b4b45acb13bd9c240b004332e76b09eb6c3e3"
SPECIES=("Dd","Fr")
TREATMENTS=("Soil","Forbs","Grasses")
PRIMARY=("RSD","PDW","Nuptake")
SECONDARY=("RDW","SDW","SRL","SRA","RTD","RLD","RMF","TRL")
SEED=20261124
N_BOOT=20000
N_PERM=9999


class AuditStop(ValueError):
    pass


def read_source(path:Path,exclude_sample_ids:tuple[int,...]=()):
    blob=path.read_bytes()
    blob_sha=hashlib.sha1(f"blob {len(blob)}\0".encode()+blob).hexdigest()
    if blob_sha!=SOURCE_BLOB_SHA:
        raise AuditStop(f"Source Git blob mismatch: {blob_sha}")
    text=blob.decode("utf-8-sig")
    rows=list(csv.DictReader(text.splitlines()))
    if len(rows)!=60 or len(rows[0])!=44:
        raise AuditStop(f"Unexpected sample panel {len(rows)}x{len(rows[0])}")
    needed=("SampleID","Species","Treatment","Replicate")+PRIMARY+SECONDARY
    if any(z not in rows[0] for z in needed):
        raise AuditStop("Source schema unexpectedly changed")
    keys=set()
    samples=set()
    grouped=defaultdict(list)
    for row in rows:
        sp,tr=row["Species"],row["Treatment"]
        rep=int(row["Replicate"]); sid=int(row["SampleID"])
        if sp not in SPECIES or tr not in TREATMENTS or rep not in range(1,11):
            raise AuditStop("Unexpected taxon, solution type or experimental replicate")
        key=(sp,tr,rep)
        if key in keys or sid in samples:
            raise AuditStop("Repeated biological unit")
        keys.add(key);samples.add(sid)
        if sid in exclude_sample_ids:
            continue
        vals={}
        for k in PRIMARY+SECONDARY:
            v=row[k].strip()
            if v=="":
                if k in PRIMARY:raise AuditStop(f"Missing primary response {sp} {tr} {k}")
                vals[k]=None
                continue
            number=float(v)
            if not np.isfinite(number) or number<0:
                raise AuditStop(f"Invalid response {k}")
            if k in PRIMARY and number<=0:
                raise AuditStop(f"Nonpositive primary response {k}")
            vals[k]=number
        grouped[(sp,tr)].append(vals)
    if len(grouped)!=6:
        raise AuditStop("Missing species × solution cell")
    expected=(
        {("Dd","Forbs"):9,("Dd","Grasses"):10,("Dd","Soil"):9,
         ("Fr","Forbs"):10,("Fr","Grasses"):10,("Fr","Soil"):10}
        if tuple(sorted(exclude_sample_ids))==(3,22)
        else {(sp,tr):10 for sp in SPECIES for tr in TREATMENTS}
    )
    if any(len(grouped[(sp,tr)])!=n for (sp,tr),n in expected.items()):
        raise AuditStop("Unexpected original-author sample exclusion or source panel")
    if set(samples)!=set(range(1,61)):
        raise AuditStop("Raw source SampleID panel 1..60 incomplete")
    source={
        "source_git_blob_sha":blob_sha,
        "source_sha256":hashlib.sha256(blob).hexdigest(),
        "n_plants_original":len(rows),"n_plants_analyzed":sum(len(v) for v in grouped.values()),
        "original_author_excluded_sample_ids":list(exclude_sample_ids),
        "n_columns":len(rows[0]),
        "n_species":2,"n_solution_treatments":3,
        "counts":{f"{s}|{t}":len(grouped[(s,t)]) for s in SPECIES for t in TREATMENTS},
        "primary_missing":0,"status":"SOURCE_PANEL_PASS",
        "experimental_unit":"original later-arriving plant (not metabolomic feature or plant-solution source replicate)"
    }
    return grouped,source


def outcome_bootstrap(grouped,response,n_boot=N_BOOT,seed=SEED):
    # Draw all six source groups ONCE per response; sharing the Soil controls
    # maintains the source within-species dependency across contrasts.
    rng=np.random.default_rng(seed+list(PRIMARY+SECONDARY).index(response))
    draws={}
    values={}
    for sp in SPECIES:
        for tr in TREATMENTS:
            x=np.array([r[response] for r in grouped[(sp,tr)]
                        if r[response] is not None],dtype=float)
            if len(x)<6:
                raise AuditStop("Insufficient biological plants in "+sp+" "+tr+" "+response)
            values[(sp,tr)]=x
            draws[(sp,tr)]=rng.choice(x,size=(n_boot,len(x)),replace=True).mean(1)
    out=[]
    for sp in SPECIES:
        for legacy in TREATMENTS[1:]:
            control=values[(sp,"Soil")];treat=values[(sp,legacy)]
            a=float(treat.mean()-control.mean())
            r=float(treat.mean()/control.mean())
            diff=draws[(sp,legacy)]-draws[(sp,"Soil")]
            ratio=draws[(sp,legacy)]/draws[(sp,"Soil")]
            out.append({
                "species":sp,"conditioning":legacy,"outcome":response,
                "control_mean":float(control.mean()),"treated_mean":float(treat.mean()),
                "absolute_difference":a,
                "absolute_difference_bootstrap95":np.quantile(diff,[.025,.975]).tolist(),
                "ratio_of_means":r,"log_ratio":float(np.log(r)),
                "log_ratio_bootstrap95":np.quantile(np.log(ratio),[.025,.975]).tolist(),
                "percent_change":float(100*(r-1)),
                "n_control":len(control),"n_treatment":len(treat)
            })
    did=[]
    for legacy in TREATMENTS[1:]:
        one=[r for r in out if r["conditioning"]==legacy]
        dd,ff=(r for sp in SPECIES for r in one if r["species"]==sp)
        diff=(
            np.log(draws[("Dd",legacy)]/draws[("Dd","Soil")])-
            np.log(draws[("Fr",legacy)]/draws[("Fr","Soil")])
        )
        did.append({
            "treatment":legacy,"outcome":response,
            "Dd_minus_Fr_log_ratio":float(dd["log_ratio"]-ff["log_ratio"]),
            "Dd_minus_Fr_95pct_bootstrap":np.quantile(diff,[.025,.975]).tolist(),
            "interpretation":"Difference between two species-specific log treatment/control ratios, not species interaction."
        })
    return out,did


def treatment_label_perm(values1,values0,rng,n_perm=N_PERM):
    # Plant-wise label exchangeability is a DESCRIPTIVE reference only:
    # source groups treatments by soil-solution trays and lacks tray IDs.
    obs=float(np.mean(values1)-np.mean(values0))
    both=np.r_[values1,values0]
    n1=len(values1)
    count=0
    for _ in range(n_perm):
        idx=rng.permutation(len(both))
        v=both[idx[:n1]].mean()-both[idx[n1:]].mean()
        if abs(v)>=abs(obs)-1e-12:count+=1
    return (count+1)/(n_perm+1)


def holm(raw):
    n=len(raw)
    order=np.argsort(raw)
    corr=np.zeros(n)
    last=0.
    for rank,j in enumerate(order):
        last=max(last, min(1.,raw[j]*(n-rank)))
        corr[j]=last
    return corr.tolist()


def analysis(path:Path,n_boot=N_BOOT,n_perm=N_PERM,
             exclude_sample_ids:tuple[int,...]=()):
    grouped,meta=read_source(path,exclude_sample_ids=exclude_sample_ids)
    results=[]
    interactions=[]
    rng=np.random.default_rng(SEED)
    for k in PRIMARY:
        out,did=outcome_bootstrap(grouped,k,n_boot=n_boot)
        for q in out:
            a=np.array([x[k] for x in grouped[(q["species"],q["conditioning"])]],float)
            b=np.array([x[k] for x in grouped[(q["species"],"Soil")]],float)
            q["conditional_pot_label_permutation_two_sided_p"]=treatment_label_perm(a,b,rng,n_perm)
        results.extend(out);interactions.extend(did)
    ps=[r["conditional_pot_label_permutation_two_sided_p"] for r in results]
    for r,q in zip(results,holm(ps)):
        r["holm_p_across_12_primary_exploratory_tests"]=float(q)
    secondary=[]
    for k in SECONDARY:
        o,d=outcome_bootstrap(grouped,k,n_boot=n_boot)
        secondary.extend(o)
    return {
        "status":"ORIGINAL_EXPERIMENT_REANALYSIS_CHEMICAL_LEGACY_NOT_NEW_DISCOVERY",
        "source":"Delory et al. (2021) Plant Cell and Environment DOI 10.1111/pce.13999",
        "source_data":"10.5281/zenodo.4010162 (GitHub v1.0.0 source-released version)",
        "prescore_contract":"docs/DELORY_CHEMICAL_LEGACY_EXTERNAL_PRE_SCORE_20261124.md",
        "source_audit":meta,
        "design":{"species":list(SPECIES),"solution_conditions":list(TREATMENTS),
                  "primary_outcomes":list(PRIMARY),
                  "n_per_species_condition":"9-10 after original-author exclusions" if exclude_sample_ids else 10,
                  "bootstrap_plant_unit_resamples":n_boot,
                  "bootstrap_seed":SEED,
                  "conditional_unblocked_plant_permutations":n_perm,
                  "permutation_caveat":"Not exact treatment randomization p because source used soil-solution trays and no tray or independent conditioning batch id was archived"},
        "primary_process_output":results,
        "between_species_log_response_differences":interactions,
        "secondary_ecological_states":secondary,
        "interpretive_cautions":[
            "Original publication already reported species-specific soil chemical legacy and root foraging without comparable N uptake losses.",
            "No early resident plant grows alongside later focal: illustrates environmental-conditioning channel not direct live plant competition.",
            "Soil-conditioned solution can have nutrients, metabolites and microbes; no specific chemical causal mediator is isolated by treatment label.",
            "Although individual 60 potted plants measured, all plants in a tray receive the same solution: tray/batch dependence may understate plant-level bootstrap uncertainty.",
            "Unadjusted pot permutation p and Holm p are conditional statistical diagnostics, not randomized source-level inference.",
            "Nuptake is calculated total N at harvest, not direct short-term uptake flux.",
            "Each outcome is correlated within plant and 12 endpoint contrasts are not 12 independent studies.",
            "No changes to frozen KNK v28-v38 or Ecology Letters manuscript."
        ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--csv",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--original-author-filter",action="store_true",
                    help="Exclude original authors' IDs 3 and 22 before scoring (58 plants)")
    args=ap.parse_args()
    removed=(3,22) if args.original_author_filter else ()
    data=analysis(args.csv,exclude_sample_ids=removed)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("DELORY_LEGACY_RESULT_START")
    print("SOURCE_AUDIT",json.dumps(data["source_audit"]))
    for r in data["primary_process_output"]:
        print("PROCESS_OUTPUT",json.dumps(r))
    for r in data["between_species_log_response_differences"]:
        print("DIFFERENCE_DIAN_FEST",json.dumps(r))
    print("DELORY_LEGACY_RESULT_DONE",args.output)


if __name__=="__main__":main()
