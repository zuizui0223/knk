#!/usr/bin/env python3
"""Source-locked 2020/2021 POEM field identity vs rank history storage.

Run only after docs/POEM_FIELD_PRE_SCORE_CONTRACT_20261123.md is committed.
Raw ZIP is acquired ephemeral from Zenodo DOI 10.5281/zenodo.10119982;
never committed or included in the output JSON.

One experiment, two cohorts, repeated yearly measurements of persistent plots.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path
import re
import zipfile

import numpy as np

from benchmarks.audit_poem_design import (
    BIOMASS_PATTERN, SOURCE_SOWN, extract_book_rows,
)
from benchmarks.probe_poem_archive import EXPECTED_MD5

EXPECTED_SOURCE_SHA256 = "7541ce47b61244f34a6b3b7f1187a2ef16f79bc173d73e6abed206203939ef8f"
SEED = 20261123
N_PERM = 9999
FOUR = ("F", "G", "L", "S")
THREE = ("F", "G", "L")
TREATMENTS = ("F", "G", "L", "S", "B")
PLOT_PATTERN = re.compile(r"^(\d+)-([FGSLB])(\d+)$")


class SourceStop(ValueError):
    pass


def source_panels(path: Path) -> tuple[dict[tuple[str, int], dict], dict]:
    raw = path.read_bytes()
    hash_md5 = hashlib.md5(raw).hexdigest()
    hash_sha256 = hashlib.sha256(raw).hexdigest()
    if hash_md5 != EXPECTED_MD5 or hash_sha256 != EXPECTED_SOURCE_SHA256:
        raise SourceStop("Unexpected source archive digest")
    panels: dict[tuple[str,int],dict] = {}
    all_species: set[str] = set()
    by_cohort: dict[str,dict[str,str]] = defaultdict(dict)
    duplicate_cells = {}
    with zipfile.ZipFile(path) as zf:
        for name in zf.namelist():
            m = BIOMASS_PATTERN.search(name)
            if not m:
                continue
            cohort, cohort2, harvest = m.groups()
            if cohort != cohort2:
                raise SourceStop(f"Inconsistent cohort in member path: {name}")
            growing_season = int(harvest[:4]) - int(cohort) + 1
            if growing_season not in (1, 2, 3):
                raise SourceStop(f"Unexpected growing season: {name}")
            records = extract_book_rows(zf, name)
            for r in records:
                if r["year"] != cohort or r["harvest"] != harvest:
                    raise SourceStop("Harvest labels not consistent with workbook")
            plot_map: dict[str,tuple[str,str]] = {}
            unique_quad: dict[str,set[str]] = defaultdict(set)
            mass: dict[tuple[str,str,str],float] = defaultdict(float)
            key_count = Counter()
            known_species = set()
            for r in records:
                match = PLOT_PATTERN.fullmatch(r["plot"])
                if not match:
                    raise SourceStop("Unparsed experimental plot identifier " + r["plot"])
                plotid, treat, replicate = match.groups()
                if treat not in TREATMENTS or replicate not in ("1","2","3","4","5"):
                    raise SourceStop("Invalid source treatment or replicate coding")
                if r["plot"] in plot_map and plot_map[r["plot"]] != (treat,replicate):
                    raise SourceStop("Plot changes treatment or replicate within workbook")
                plot_map[r["plot"]] = (treat,replicate)
                unique_quad[r["plot"]].add(r["quadrat"])
                raw_g = r["mass"]
                if not isinstance(raw_g,(int,float)) or not np.isfinite(float(raw_g)) or float(raw_g)<0:
                    raise SourceStop("Missing/negative biomass in source")
                sp = r["taxon"]
                if not sp:
                    raise SourceStop("Blank species label")
                if sp == "Unknown":
                    continue  # original R analysis explicitly dplyr::select(-Unknown)
                known_species.add(sp)
                k=(r["plot"],r["quadrat"],sp)
                mass[k] += float(raw_g)
                key_count[k]+=1
            plot_ids=sorted(plot_map)
            if len(plot_ids)!=25 or sorted(Counter(v[0] for v in plot_map.values()).values()) != [5]*5:
                raise SourceStop(f"Unbalanced or missing plots in {name}")
            if set(plot_ids) != set(unique_quad):
                raise SourceStop("Missing quadrat for plot")
            nquad = [len(unique_quad[p]) for p in plot_ids]
            if len(set(nquad)) != 1 or nquad[0] not in (2,4):
                raise SourceStop("Quadrat count unequal within a harvest; revisit area-standardization")
            if cohort in by_cohort and by_cohort[cohort] and (
                {p:plot_map[p][0] for p in plot_ids} != by_cohort[cohort]
            ):
                raise SourceStop("Plot treatment assignment drift across harvests")
            by_cohort[cohort] = {p:plot_map[p][0] for p in plot_ids}
            plot_species_mass=defaultdict(float)
            for (p,q,sp),v in mass.items():
                plot_species_mass[(p,sp)] += v/len(unique_quad[p])
            key=(cohort,growing_season)
            if key in panels:
                raise SourceStop("Multiple workbooks for same cohort and season")
            panels[key] = {
                "plot_ids":plot_ids, "plot_map":plot_map,
                "biomass_by_plot_species":dict(plot_species_mass),
                "source_rows":len(records),
                "explicit_zeros":int(sum(r["mass"]==0 for r in records)),
                "duplicate_fragment_records":int(sum(v-1 for v in key_count.values())),
                "quadrats_per_plot":nquad[0],
                "harvest":harvest,
                "n_species_positively_recorded":len(known_species),
            }
            all_species.update(known_species)
    if len(panels)!=6 or len(all_species)!=49:
        raise SourceStop(f"Expected 6 harvests with 49 identified taxa, got {len(panels)}, {len(all_species)}")
    names=sorted(all_species)
    if not set(SOURCE_SOWN).issubset(all_species):
        raise SourceStop("Some original 14 sown species absent from all source files")
    audit={
        "status":"PASS",
        "source_zip_md5":hash_md5,"source_zip_sha256":hash_sha256,
        "n_cohorts":2,"n_unique_independent_plots_across_cohorts":50,
        "n_harvests":6,"n_species_after_excluding_unknown":len(names),
        "n_sown_taxa":len(SOURCE_SOWN),
        "total_source_rows":sum(p["source_rows"] for p in panels.values()),
        "total_duplicate_fragment_records":sum(p["duplicate_fragment_records"] for p in panels.values()),
        "total_explicit_zero_records":sum(p["explicit_zeros"] for p in panels.values()),
        "original_repeated_plot_units_per_cohort":25,
        "plot_aggregation":"sum separate entries within species×quadrat; average quadrats to plot",
        "implicit_zero_policy":"source listed species absence in otherwise observed quadrat is recorded as zero harvested biomass",
        "cohort_years":[2020,2021],
        "identified_species_names":names
    }
    for p in panels.values():
        p["identified_species_names"] = names
    return panels, audit


def matrix_for(panel:dict, treatments:tuple[str,...], restrict_sown:bool=False):
    names = (sorted(SOURCE_SOWN) if restrict_sown else panel["identified_species_names"])
    plots=sorted(p for p in panel["plot_ids"] if panel["plot_map"][p][0] in treatments)
    labels=[panel["plot_map"][p][0] for p in plots]
    counts=Counter(labels)
    if len(plots) != len(treatments)*5 or any(counts.get(k,0)!=5 for k in treatments):
        raise SourceStop(f"Unexpected treatment balance {counts}")
    x = np.array([
        [panel["biomass_by_plot_species"].get((p,sp),0.0) for sp in names]
        for p in plots
    ], dtype=float)
    return x, np.array(labels), plots, names


def bray_curtis(x:np.ndarray)->np.ndarray:
    dif=np.abs(x[:,None,:]-x[None,:,:]).sum(axis=2)
    den=(x[:,None,:]+x[None,:,:]).sum(axis=2)
    return np.divide(dif,den,out=np.zeros_like(dif),where=den>0)


def gower(d:np.ndarray)->np.ndarray:
    n=len(d)
    J=np.eye(n)-np.ones((n,n))/n
    return -.5*J@(d*d)@J


def history_r2(g:np.ndarray, y:np.ndarray)->float:
    total=float(np.trace(g))
    if total<=1e-12:
        return 0.
    groups=np.unique(y)
    if len(groups)<2:
        raise ValueError("Needs >1 treatment group")
    between=sum(float(g[np.ix_(ix,ix)].sum())/len(ix)
                for a in groups for ix in [np.flatnonzero(y==a)])
    return float(between/total)


def evaluate_one(x:np.ndarray,labels:np.ndarray,perm:np.ndarray)->dict:
    g=gower(bray_curtis(x))
    r2=history_r2(g,labels)
    null=np.array([history_r2(g,p) for p in perm],dtype=float)
    return {
        "r2":r2,
        "null_median":float(np.median(null)),
        "null_ci95":np.quantile(null,[.025,.975]).tolist(),
        "E":r2-float(np.median(null)),
        "p_upper":float((1+np.count_nonzero(null>=r2-1e-12))/(len(null)+1)),
        "n_units":len(labels),
        "n_response_columns":x.shape[1],
        "n_permutations":len(perm)
    }


def analysis(panels:dict[tuple[str,int],dict],source:dict,n_perm:int=N_PERM,
             seed:int=SEED)->dict:
    # Persistent plot IDs in each cohort are sorted identically at all ages.
    seeds={c:np.random.default_rng(seed+int(c)) for c in ("2020","2021")}
    permutation_cache={}
    primary=[]
    sensitivity=[]
    for cohort in ("2020","2021"):
        for age in (1,2,3):
            panel=panels[(cohort,age)]
            x,labels,plots,names=matrix_for(panel,FOUR)
            if ((cohort,"four") not in permutation_cache):
                permutation_cache[(cohort,"four")] = (
                    plots,labels,[seeds[cohort].permutation(labels) for _ in range(n_perm)]
                )
            original_plot_order, original_labels, perm = permutation_cache[(cohort,"four")]
            if plots!=original_plot_order or np.any(labels!=original_labels):
                raise SourceStop("Repeated-plot identity/order mismatch across harvests")
            per=np.asarray(perm)
            if np.any(x.sum(axis=1)<=0):
                raise SourceStop("Zero total plot mass; cannot form proportions")
            p=x/x.sum(axis=1,keepdims=True)
            ranked=np.sort(p,axis=1)[:,::-1]
            named=evaluate_one(p,labels,per)
            sorted_=evaluate_one(ranked,labels,per)
            total=evaluate_one(x.sum(axis=1)[:,None],labels,per)
            top1=evaluate_one(ranked[:,:1],labels,per)
            entry={
                "cohort":cohort,"age":age,"harvest":panel["harvest"],
                "n_plots":len(plots),"n_species":len(names),
                "representations":{"named_all":named,"sorted_all":sorted_},
                "delta_E_named_minus_rank":named["E"]-sorted_["E"],
                "delta_R2_named_minus_rank":named["r2"]-sorted_["r2"],
                "total_biomass":total,"largest_species_share":top1,
                "mean_total_grams_per_quadrat":float(x.sum(axis=1).mean()),
                "zeros_plot_species":int((x==0).sum()),
                "quadrat_count":panel["quadrats_per_plot"]
            }
            primary.append(entry)
            # 14 originally sown species; use same 20 full treatment assignments.
            x14,l14,plot14,name14=matrix_for(panel,FOUR,restrict_sown=True)
            assert plots==plot14 and np.array_equal(labels,l14)
            if np.any(x14.sum(axis=1)<=0):
                sown_out={"status":"STOP_ZERO_SOWN_TOTAL","n_zero_sown_plots":int(np.sum(x14.sum(axis=1)==0))}
            else:
                x14=x14/x14.sum(axis=1,keepdims=True)
                sown_out={"status":"SCORED",
                          "named":evaluate_one(x14,labels,per),
                          "rank":evaluate_one(np.sort(x14,axis=1)[:,::-1],labels,per)}
            # F/G/L-only histories (3 classes, 15 different plots).
            x3,y3,plot3,n3=matrix_for(panel,THREE)
            if (cohort,"three") not in permutation_cache:
                permutation_cache[(cohort,"three")]=(
                    plot3,y3,[seeds[cohort].permutation(y3) for _ in range(n_perm)]
                )
            p3,l3,perm3=permutation_cache[(cohort,"three")]
            if plot3!=p3 or not np.array_equal(y3,l3):
                raise SourceStop("Three-treatment plot mismatch across repeated ages")
            if np.any(x3.sum(axis=1)<=0):
                three_out={"status":"STOP_ZERO_TOTAL"}
            else:
                p3frac=x3/x3.sum(axis=1,keepdims=True)
                three_out={"status":"SCORED",
                           "named":evaluate_one(p3frac,y3,np.asarray(perm3)),
                           "rank":evaluate_one(np.sort(p3frac,axis=1)[:,::-1],y3,np.asarray(perm3))}
            sensitivity.append({"cohort":cohort,"age":age,
                                "sown14_4level":sown_out,
                                "all49_FGLonly":three_out})
    # Correct descriptive multiple comparisons for the twelve primary tests only.
    tests=[(i,r) for i,z in enumerate(primary)
           for r in ("named_all","sorted_all")]
    order=sorted(tests,key=lambda ir:primary[ir[0]]["representations"][ir[1]]["p_upper"])
    q=[min(1.,primary[i]["representations"][r]["p_upper"]*(len(order)-k))
       for k,(i,r) in enumerate(order)]
    q=list(np.maximum.accumulate(q))
    for (i,r),val in zip(order,q):
        primary[i]["representations"][r]["holm_p_12_tests"]=float(val)
    return {
        "status":"INDEPENDENT_FIELD_SOURCE_SCORED_EXPLORATORY_NOT_META_POOLED",
        "source_audit":source,
        "analysis_contract":"docs/POEM_FIELD_PRE_SCORE_CONTRACT_20261123.md",
        "design":{
            "independent_source_clusters":1,
            "cohorts":2,"persistent_plots_per_cohort":25,
            "primary_plots_per_cohort":20,
            "primary_groups":list(FOUR),"primary_n_groups":4,
            "free_succession_B_excluded":True,
            "n_permutations":n_perm,"seed":seed,
            "same_plot_permutations_across_ages":True,
            "metric":"Bray-Curtis Gower-centered one-way R2, E above permuted-null median",
            "permutation":"unblocked within cohort, whole plots only",
            "observation":"mean grams per quadrat",
            "repeated_measurement_note":"3 ages of same plot are NOT 3 independent replications"
        },
        "primary_4level_all_identified_taxa":primary,
        "predeclared_sensitivities":sensitivity,
        "claim_limits":[
            "Never pool cohort-by-age snapshots as independent studies",
            "Name vs ranked E are geometry-dependent effect sizes, not information bits",
            "Unblocked permutation assumes exchangeable plot treatments; source randomized layout should be independently checked",
            "Zero-fill applies ONLY unlisted species in known sampled quadrats",
            "No causal inference of niche pre-emption or facilitation",
            "Root image records do not identify roots by species",
            "No change to KNK frozen v28-v38",
        ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--permutations",type=int,default=N_PERM)
    args=ap.parse_args()
    panels,audit=source_panels(args.zip)
    result=analysis(panels,audit,n_perm=args.permutations)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("POEM_RESULT_START")
    print("SOURCE_AUDIT",json.dumps({k:v for k,v in audit.items()
                              if k!="identified_species_names"}))
    for v in result["primary_4level_all_identified_taxa"]:
        print("POEM_EFFECT",json.dumps({
            "cohort":v["cohort"],"age":v["age"],"harvest":v["harvest"],
            "named":v["representations"]["named_all"],
            "sorted":v["representations"]["sorted_all"],
            "delta_E":v["delta_E_named_minus_rank"],
            "total_biomass":v["total_biomass"],
            "rank1":v["largest_species_share"],
        }))
    for v in result["predeclared_sensitivities"]:
        print("POEM_SENSITIVITY",json.dumps(v))
    print("POEM_RESULT_END",args.output)

if __name__=="__main__":
    main()
