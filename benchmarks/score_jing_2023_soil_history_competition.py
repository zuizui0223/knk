#!/usr/bin/env python3
"""Zenodo 7856898 authentic aboveground biomass: soil history × drought × mixing.

Source 2023-04-23, published original DOI https://doi.org/10.5281/zenodo.7856898.
This data source might overlap Jing 2024 Dryad: do not count independent studies.
Source-locked prescore protocol:
docs/JING_ZENODO_MIXTURE_SOIL_HISTORY_PRE_SCORE_20261124.md
No source data redistributed or in-repo biological re-scoring of frozen KNK.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import re
import numpy as np

from benchmarks.probe_jing_2024_source_xlsx import xlsx_rows

SOURCE_MD5="5afb3a33ad49573aa7e255d2541ab0ed"
SOURCE_SHA256="636ad7ebbe1f8c452b40fd8fde0e03cf611d21baf1db959fcce42385bc23b899"
MONO_SPECIES=("SK","AF","AC","PT","ES")
PAIRS=("SK-AF","AF-AC","SK-AC")
ASSEMBLAGES=PAIRS+("SK-AF-AC","SK-AF-AC-PT","SK-AF-AC-ES")
SOIL=("SK","AF","AC","MIX","S")
WATER=("Wet","Dry")
B=10000
SEED=20261124

HEADER_RE=re.compile(r"^\(([A-Z-]+)\)_([A-Z]+)$")


class SourceStop(ValueError):pass


def nonnegative(value,source):
    if isinstance(value,bool):raise SourceStop("Boolean instead of biomass: "+source)
    try:v=float(value)
    except (ValueError,TypeError):
        raise SourceStop("Missing/non-numeric biomass: "+source)
    if not math.isfinite(v) or v<0:
        raise SourceStop("Negative/nonfinite biomass: "+source)
    return v


def parse_workbook(path:Path):
    contents=path.read_bytes()
    if hashlib.md5(contents).hexdigest()!=SOURCE_MD5 or \
       hashlib.sha256(contents).hexdigest()!=SOURCE_SHA256:
        raise SourceStop("Unrecognized original Zenodo source archive")
    sheets=xlsx_rows(path)
    if set(sheets)!={"monoculture_aboveground biomass","mixture_aboveground biomass"}:
        raise SourceStop("Unexpected original workbook sheets")
    mon=sheets["monoculture_aboveground biomass"]["rows"]
    mix=sheets["mixture_aboveground biomass"]["rows"]
    if len(mon)!=41 or len(mix)!=42:
        raise SourceStop("Expected 1+40 mono rows and 2+40 mix rows")
    if tuple(mon[0][:7])!=("Drought treatment","Soil treatment",*MONO_SPECIES):
        raise SourceStop("Unmatched source monoculture columns")
    if tuple(mix[1][:2])!=("Drought treatment","Soil treatment"):
        raise SourceStop("Unexpected mixture second-level headers")
    if len(mix[1])!=24:raise SourceStop("Unexpected mixture width")
    assignments={}
    for col,header in enumerate(mix[1]):
        if col<2 or header is None or str(header).strip()=="":
            continue
        matched=HEADER_RE.match(str(header).strip())
        if not matched:raise SourceStop("Unexpected source mixture header "+str(header))
        assembly,species=matched.groups()
        if assembly not in ASSEMBLAGES or species not in assembly.split("-"):
            raise SourceStop("Unexpected species/assembly pairing "+str(header))
        if (assembly,species) in assignments:
            raise SourceStop("Repeated mixture species key")
        assignments[assembly,species]=col
    expected={(m,s) for m in ASSEMBLAGES for s in m.split("-")}
    if set(assignments)!=expected:
        raise SourceStop(f"Missing mixture component columns {expected-set(assignments)}")
    # Original rows are 4 independent replicate positions per water×soil
    # within EACH planting design; a row number is not a verified block ID.
    def group_rows(raw,has_second_header=False):
        out=defaultdict(list)
        for idx,row in enumerate(raw[2:] if has_second_header else raw[1:],start=1):
            if len(row)<2:raise SourceStop("Truncated original treatment row")
            water,soil=str(row[0]).strip(),str(row[1]).strip()
            if water not in WATER or soil not in SOIL:
                raise SourceStop(f"Unrecognized source treatment {water} {soil}")
            out[water,soil].append((idx,row))
        if set(out)!={(w,s) for w in WATER for s in SOIL}:
            raise SourceStop("Not all water×soil conditions present")
        if any(len(v)!=4 for v in out.values()):
            raise SourceStop("Unequal source 4-replicate cells")
        return out
    mono_rows=group_rows(mon)
    mix_rows=group_rows(mix,has_second_header=True)
    mono={}
    mixtures={}
    for w in WATER:
        for soil in SOIL:
            mr=mono_rows[w,soil]
            xx=mix_rows[w,soil]
            for j,sp in enumerate(MONO_SPECIES,2):
                mono[w,soil,sp]=np.array(
                    [nonnegative(row[j],f"mono {w}/{soil}/{sp}/row{i}")
                     for i,row in mr],float)
            for assembly in ASSEMBLAGES:
                components=assembly.split("-")
                values=np.array([
                    [nonnegative(row[assignments[assembly,sp]],
                                 f"mixture {assembly}/{w}/{soil}/{sp}/row{i}")
                     for sp in components]
                    for i,row in xx],float)
                if values.shape!=(4,len(components)):
                    raise SourceStop("Incorrect source mixture pot×species layout")
                mixtures[w,soil,assembly]=values
    allvals=np.r_[np.concatenate(list(mono.values())),
                  np.concatenate([v.flatten() for v in mixtures.values()])]
    summary={
        "status":"SOURCE_AUTHENTICATED_AND_EXPERIMENTAL_GRAIN_VALIDATED",
        "source_doi":"10.5281/zenodo.7856898",
        "source_workbook_md5":SOURCE_MD5,"source_workbook_sha256":SOURCE_SHA256,
        "n_monoculture_rows":len(mon)-1,"n_mixture_rows":len(mix)-2,
        "n_water":len(WATER),"n_soil":len(SOIL),
        "n_source_plants":len(MONO_SPECIES),
        "n_mixture_assemblages":len(ASSEMBLAGES),
        "mixture_assemblage_names":list(ASSEMBLAGES),
        "n_monoculture_species_groupings":len(mono),
        "n_mixture_pot_groupings":len(mixtures),
        "n_mono_original_pot_units_if_column_specific":sum(len(v) for v in mono.values()),
        "n_mixture_original_pot_units_if_assembly_specific":sum(len(v) for v in mixtures.values()),
        "n_mixture_species_measurements":sum(v.size for v in mixtures.values()),
        "n_observed_zeros":int(np.count_nonzero(allvals==0)),
        "n_nonfinite":int(np.count_nonzero(~np.isfinite(allvals))),
        "minimum_biomass":float(allvals.min()),
        "max_biomass":float(allvals.max()),
        "assumed_unit_rule":"one pot for each original water×soil×assembly×replicate; replicates ordered within group but physical pot IDs not archived",
        "potential_duplicate_with_2024_dryad":"UNKNOWN, treat same source cluster pending original study audit"
    }
    return mono,mixtures,summary


def mean_log_ratio(a:np.ndarray,b:np.ndarray,about:str):
    u=float(np.mean(a));v=float(np.mean(b))
    if u<=0 or v<=0:
        raise SourceStop("Zero group mean; no pseudocount allowed "+about)
    return math.log(u/v)


def bootstrap_group_means(mono:dict,mixtures:dict,seed:int=SEED,n_boot:int=B):
    rng=np.random.default_rng(seed)
    mb={}; xb={}
    # One control draw per species×water×soil shared across every mixture
    for w in WATER:
        for s in SOIL:
            for sp in MONO_SPECIES:
                arr=mono[w,s,sp]
                which=rng.integers(0,len(arr),size=(n_boot,len(arr)))
                mb[w,s,sp]=arr[which].mean(axis=1)
            for name in ASSEMBLAGES:
                arr=mixtures[w,s,name]
                which=rng.integers(0,len(arr),size=(n_boot,len(arr)))
                # Whole mixture pot = one row; keep all species covarying
                xb[w,s,name]=arr[which].mean(axis=1)
    return mb,xb


def score(mono:dict,mixtures:dict,source:dict,n_boot:int=B,seed:int=SEED):
    mb,xb=bootstrap_group_means(mono,mixtures,seed,n_boot)
    def summarise(z):
        return np.quantile(z,[.025,.975]).tolist()
    cells=[]
    bootC={}
    for assembly in PAIRS:
        taxa=assembly.split("-")
        for w in WATER:
            for soil in SOIL:
                pop=mixtures[w,soil,assembly]
                means=pop.mean(axis=0)
                for k,sp in enumerate(taxa):
                    ref=mono[w,soil,sp]
                    observed=mean_log_ratio(pop[:,k],ref,f"{assembly},{sp},{w},{soil}")
                    bootstrap=np.log(xb[w,soil,assembly][:,k]/mb[w,soil,sp])
                    if not np.isfinite(bootstrap).all():
                        raise SourceStop("Non-finite bootstrap log ratio")
                    bootC[assembly,sp,w,soil]=bootstrap
                    cells.append({
                        "pair":assembly,"focal":sp,"other_species":taxa[1-k],
                        "water":w,"soil":soil,"n_mixture_pots":len(pop),
                        "n_monoculture_pots":len(ref),
                        "mixture_mean_g":float(means[k]),
                        "monoculture_mean_g":float(ref.mean()),
                        "solo_referenced_log_mixture_biomass":observed,
                        "bootstrap95":summarise(bootstrap),
                        "mixture_vs_mono_ratio":float(math.exp(observed)),
                        "n_resamples":n_boot
                    })
    comparisons=[]
    for assembly in PAIRS:
        for sp in assembly.split("-"):
            for w in WATER:
                for soil in SOIL:
                    if soil=="S":continue
                    current=next(r for r in cells if
                                 (r["pair"],r["focal"],r["water"],r["soil"] )==
                                 (assembly,sp,w,soil))
                    sterile=next(r for r in cells if
                                 (r["pair"],r["focal"],r["water"],r["soil"] )==
                                 (assembly,sp,w,"S"))
                    observed=current["solo_referenced_log_mixture_biomass"]-sterile["solo_referenced_log_mixture_biomass"]
                    draws=bootC[assembly,sp,w,soil]-bootC[assembly,sp,w,"S"]
                    comparisons.append({
                        "pair":assembly,"focal":sp,"water":w,"live_soil":soil,
                        "comparison":"C(live soil)-C(sterile S) within same water",
                        "difference_in_solo_referenced_mixture_log_biomass":observed,
                        "bootstrap95":summarise(draws),
                        "positive_fraction":float(np.mean(draws>0))
                    })
    stress=[]
    for assembly in PAIRS:
        for sp in assembly.split("-"):
            for soil in SOIL:
                if soil=="S":continue
                dry=next(q for q in comparisons if q["pair"]==assembly and q["focal"]==sp and q["live_soil"]==soil and q["water"]=="Dry")
                wet=next(q for q in comparisons if q["pair"]==assembly and q["focal"]==sp and q["live_soil"]==soil and q["water"]=="Wet")
                dd=(
                    bootC[assembly,sp,"Dry",soil]-bootC[assembly,sp,"Dry","S"]-
                    bootC[assembly,sp,"Wet",soil]+bootC[assembly,sp,"Wet","S"]
                )
                stress.append({
                    "pair":assembly,"focal":sp,"live_soil":soil,
                    "drought_minus_wet_soil_effect_log":dry["difference_in_solo_referenced_mixture_log_biomass"]-
                                                       wet["difference_in_solo_referenced_mixture_log_biomass"],
                    "bootstrap95":summarise(dd)
                })
    winners=[]
    ranks={}
    for pair in PAIRS:
        a,b=pair.split("-")
        signs=[]
        for w in WATER:
            for soil in SOIL:
                x=mixtures[w,soil,pair]
                which=0 if x[:,0].mean()>x[:,1].mean() else 1
                if x[:,0].mean()==x[:,1].mean():
                    raise SourceStop("Tied mean mixture biomasses; need tie handler")
                bp=xb[w,soil,pair]
                frac=float(np.mean(bp[:,0]>bp[:,1]))
                winner=a if which==0 else b
                winners.append({
                    "pair":pair,"water":w,"soil":soil,
                    "winner_by_mean_raw_mixture_biomass":winner,
                    "mean_biomass":{a:float(x[:,0].mean()),b:float(x[:,1].mean())},
                    "bootstrap_p_species_a_higher_raw_mixture_mass":frac,
                    "within_pot_species_measurements_are_correlated":True
                })
                signs.append(winner)
        ranks[pair]={
            "different_winners_observed":len(set(signs)),
            "wins_by_species":dict(Counter(signs)),
            "rank_varies_across_water_inoculum":len(set(signs))>1
        }
    mod=[]
    for pair in PAIRS:
        for sp in pair.split("-"):
            for w in WATER:
                cur=[r["solo_referenced_log_mixture_biomass"] for r in cells if
                     r["pair"]==pair and r["focal"]==sp and r["water"]==w]
                mod.append({
                    "pair":pair,"focal":sp,"water":w,
                    "range_across_five_soils_of_solo_referenced_log_mixture_biomass":max(cur)-min(cur),
                    "n_inocula":len(cur)
                })
    return {
        "status":"SOURCE_AUTHENTIC_SOLO_REFERENCED_COMPETITION_CONTEXT_SCORED_NOT_META_POOLED",
        "source_audit":source,
        "contract":"docs/JING_ZENODO_MIXTURE_SOIL_HISTORY_PRE_SCORE_20261124.md",
        "design":{"paired_source_species_taxa":list(MONO_SPECIES),
                  "primary_mixtures":list(PAIRS),"source_soils":list(SOIL),
                  "drought":list(WATER),"n_bootstrap":n_boot,"seed":seed,
                  "experimental_unit":"whole mixture pot; measured component biomass correlated within pot",
                  "no_raw_data_redistribution":True},
        "primary_solo_referenced_2_species_mixture":cells,
        "live_inoculum_vs_sterile_differences":comparisons,
        "soil_drought_difference_in_differences":stress,
        "raw_mixture_mass_winner_by_context":winners,
        "within_pair_mass_rank_stability":ranks,
        "soil_context_variation_of_log_mixture_impact":mod,
        "interpretive_boundaries":[
            "Source 2023 archive may overlap Jing 2024 Dryad 10.5061/dryad.jq2bvq8jj; do NOT count it as a new independent study cluster.",
            "Relative mixture/mono biomass is a source output contrast; seed density and mass-per-pot equivalence not fully independently audited.",
            "Potential 440 reconstructed pots but source lacks actual physical pot IDs, layout and inoculum batch IDs.",
            "Only 4 pot replicates per water×soil×assembly; all percentile intervals are exploratory conditional on observed groups.",
            "Sterilized S vs live soil also alters other chemical/nutrient conditions; causal microbial mediation cannot be identified.",
            "No direct species arrival order treatment in this dataset: soil history differs from arrival-sequence priority.",
            "Pot-level rank of raw mixture biomass is NOT a competitive fitness ranking.",
            "No formal randomization p values or false ecological interaction alpha coefficients.",
            "No changes to KNK frozen v28-v38 or main Ecology Letters manuscript."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xlsx",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=B)
    args=p.parse_args()
    solo,mix,audit=parse_workbook(args.xlsx)
    result=score(solo,mix,audit,n_boot=args.bootstrap)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("JING_SOIL_COMPETITION_RESULT_START")
    print("SOURCE_AUDIT",json.dumps(audit))
    print("RANK_STABILITY",json.dumps(result["within_pair_mass_rank_stability"]))
    for r in result["soil_context_variation_of_log_mixture_impact"]:
        print("SOIL_AMPLITUDE",json.dumps(r))
    for r in result["live_inoculum_vs_sterile_differences"]:
        print("CONTRAST",json.dumps(r))
    print("JING_SOIL_COMPETITION_RESULT_END",args.output)


if __name__=="__main__":main()
