#!/usr/bin/env python3
"""Source-verified post hoc paired plant-level process vs output legacy contrast.

See docs/DELORY_PROCESS_OUTPUT_DIRECT_CONTRAST_POSTHOC_20261124.md
Source: Delory et al. 2021 Plant Cell and Environment DOI 10.1111/pce.13999.
Root foraging (RSD) and total plant mass / N uptake measured on same plants.
This is a within-source exploratory outcome comparison, not a mechanism claim.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from benchmarks.score_delory_chemical_legacy import read_source, SPECIES,TREATMENTS,SEED


OUTCOMES=("PDW","Nuptake")


def process_output_contrasts(grouped:dict,n_boot:int=20000,seed:int=SEED)->list[dict]:
    rng=np.random.default_rng(seed)
    # One bootstrap index per plant cell shared by three endpoint measurements.
    draws={}
    means={}
    n={}
    for sp in SPECIES:
        for condition in TREATMENTS:
            rows=grouped[sp,condition]
            mat=np.array([[r["RSD"],r["PDW"],r["Nuptake"]] for r in rows],float)
            if not np.isfinite(mat).all() or np.any(mat<=0):
                raise ValueError("Nonpositive or missing RSD/PDW/Nuptake")
            n[sp,condition]=len(mat)
            means[sp,condition]=mat.mean(axis=0)
            which=rng.integers(0,len(mat),size=(n_boot,len(mat)))
            draws[sp,condition]=mat[which].mean(axis=1)
    records=[]
    for sp in SPECIES:
        for legacy in TREATMENTS[1:]:
            control=means[sp,"Soil"]
            treated=means[sp,legacy]
            observed_log=np.log(treated/control)
            boot_log=np.log(draws[sp,legacy]/draws[sp,"Soil"])
            for idx,label in [(1,"PDW"),(2,"Nuptake")]:
                contrast=observed_log[0]-observed_log[idx]
                boot_delta=boot_log[:,0]-boot_log[:,idx]
                records.append({
                    "species":sp,
                    "conditioning":legacy,
                    "contrasted_output":label,
                    "n_treatment_plants":n[sp,legacy],
                    "n_control_plants":n[sp,"Soil"],
                    "RSD_log_treatment_control":float(observed_log[0]),
                    "output_log_treatment_control":float(observed_log[idx]),
                    "root_foraging_minus_output_log_response":float(contrast),
                    "conditional_plant_bootstrap95":np.quantile(boot_delta,[.025,.975]).tolist(),
                    "fraction_bootstrap_above_zero":float(np.mean(boot_delta>0)),
                    "n_bootstrap":n_boot,
                    "unit":"identical individual plant rows resampled jointly for all 3 response columns",
                    "claim_boundary":"Post hoc, tray/batch-level dependence unavailable"
                })
    return records


def audit(path:Path):
    filtered,meta58=read_source(path,exclude_sample_ids=(3,22))
    unfiltered,meta60=read_source(path,exclude_sample_ids=())
    return {
        "status":"POSTHOC_PROCESS_OUTPUT_LOG_RATIO_CONTRAST_NOT_CAUSAL",
        "source_doi":"10.5281/zenodo.4010162",
        "paper_doi":"10.1111/pce.13999",
        "source_sha256":meta58["source_sha256"],
        "analysis_protocol":"docs/DELORY_PROCESS_OUTPUT_DIRECT_CONTRAST_POSTHOC_20261124.md",
        "primary_original_author_filtered_58":{
            "source_meta":meta58,
            "eight_contrasts":process_output_contrasts(filtered)
        },
        "sensitivity_unfiltered_60":{
            "source_meta":meta60,
            "eight_contrasts":process_output_contrasts(unfiltered)
        },
        "limitations":[
            "The original author excluded SampleID 3 (measurement dosing error) and 22 (phenotypic outlier) from analysis.",
            "Ratios of endpoints are reported as log ratios within each endpoint's own measurement units, not absolute biomass-to-root ratios.",
            "Co-sampled plant endpoints share root mass and uptake-derived components; bootstrap retains this within-plant dependence.",
            "Source has no tray/batch ID, so conditioning solution treatments may have fewer independent assignment units than plant pots.",
            "Confidence intervals from 20k plant bootstrap are conditional and exploratory, not an independent experiment.",
            "Do not infer unchanged biomass/N or niche-preemption mechanism from non-significant separate endpoints.",
            "Frozen KNK v28-v38 and Ecology Letters candidate unchanged."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--csv",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=audit(args.csv)
    args.output.parent.mkdir(exist_ok=True,parents=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("DELORY_DIRECT_PROCESS_OUTPUT_START")
    for mode in ("primary_original_author_filtered_58","sensitivity_unfiltered_60"):
        for v in result[mode]["eight_contrasts"]:
            print("DIRECT_CONTRAST",mode,json.dumps(v))
    print("DELORY_DIRECT_PROCESS_OUTPUT_DONE")


if __name__=="__main__":
    main()
