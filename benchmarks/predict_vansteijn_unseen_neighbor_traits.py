#!/usr/bin/env python3
"""Strict leave-one-species-out transfer of history-dependent neighbor effects.

Source: van Steijn et al. (2025) Ecology e4528, doi:10.1002/ecy.4528.
Protocol committed earlier:
docs/VANSTEIJN_UNSEEN_NEIGHBOR_TRANSFER_PRE_SCORE_20261124.md

Never train with any response involving the heldout species in either role.
Predict only directed observations where the heldout species was the EARLIER
neighbor. 15 distinct heldout neighbor roles, not 202 independent taxa.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path

import numpy as np

from benchmarks.score_vansteijn_partner_effects import prep, read_tab
from benchmarks.score_vansteijn_age_matched_neighbor_effect import transform

MORPH_TRAITS=("SLA","LA","RMF","Plant height","Lateral spread","LDMC")
LOG_TRAITS={"SLA","LA","Plant height","Lateral spread","GR"}
GUILDS=("Forb","Grass","Legume")
BASE_MODELS=("intercept","guild","focal_morph","both_morph",
             "both_morph_abs_gap","guild_both_morph_abs_gap")
RIDGE=10.0
SEED=20261124
NEG_CONTROL_SEED=20261125


class AuditStop(ValueError):
    pass


def traits_for_source(source_dir:Path,with_growth_rate:bool=False):
    src=read_tab(source_dir/"traits.csv")
    fields=MORPH_TRAITS+(("GR",) if with_growth_rate else ())
    table={}
    for row in src:
        s=row["focal"]
        if s in table:raise AuditStop("Duplicated species trait row "+s)
        values=[]
        for field in fields:
            try:v=float(row[field])
            except (ValueError,TypeError) as e:
                raise AuditStop(f"Non-numeric trait {s}/{field}") from e
            if not math.isfinite(v) or v<0:
                raise AuditStop(f"Missing/negative trait {s}/{field}")
            if field in LOG_TRAITS:v=math.log1p(v)
            values.append(v)
        guild=row["func.group"]
        if guild not in GUILDS:raise AuditStop("Unrecognized guild "+guild)
        table[s]={"traits":np.asarray(values,dtype=float),"guild":guild}
    if len(table)!=15 or Counter(x["guild"] for x in table.values()) != Counter(
            {"Forb":5,"Grass":5,"Legume":5}):
        raise AuditStop("Expected 15 taxa with 5 in each source guild")
    return table, fields


def build_features(rows:list[dict],source:dict,
                   ref_species:list[str],model:str,
                   neighbor_mapping:dict[str,str]|None=None)->np.ndarray:
    if neighbor_mapping is None:
        neighbor_mapping={name:name for name in source}
    base=np.vstack([source[sp]["traits"] for sp in ref_species])
    mu=base.mean(axis=0)
    sd=base.std(axis=0)
    sd=np.where(sd<1e-8,1.0,sd)
    # Other than fixed species/guild traits no columns use heldout test response.
    def zz(species):
        return (source[species]["traits"]-mu)/sd
    guild_ref=GUILDS[1:]
    def guild_features(species):
        return [float(source[species]["guild"]==g) for g in guild_ref]
    vectors=[]
    for row in rows:
        f=row["focal"]
        n=row["neighbor"]
        zf=zz(f)
        zn=zz(neighbor_mapping[n])
        if model=="intercept":
            vec=[]
        elif model=="guild":
            vec=guild_features(f)+guild_features(n)
        elif model=="focal_morph":
            vec=zf.tolist()
        elif model=="both_morph":
            vec=np.r_[zf,zn].tolist()
        elif model=="both_morph_abs_gap":
            vec=np.r_[zf,zn,np.abs(zf-zn)].tolist()
        elif model=="guild_both_morph_abs_gap":
            vec=guild_features(f)+guild_features(n)+np.r_[zf,zn,np.abs(zf-zn)].tolist()
        else:
            raise AuditStop("Unknown fixed model "+model)
        vectors.append(vec)
    return np.asarray(vectors,dtype=float).reshape((len(rows),-1))


def ridge_prediction(train_features:np.ndarray,train_y:np.ndarray,
                     test_features:np.ndarray,penalty:float=RIDGE)->np.ndarray:
    train_x=np.column_stack([np.ones(len(train_y)),train_features])
    test_x=np.column_stack([np.ones(len(test_features)),test_features])
    diag=np.diag(np.r_[0.0,np.full(train_features.shape[1],penalty)])
    coef=np.linalg.solve(train_x.T@train_x+diag,train_x.T@train_y)
    return test_x@coef


def run_strict_loso(rows:list[dict],traits:dict,target:str,
                    models:tuple[str,...]=BASE_MODELS,
                    neighbor_mapping:dict[str,str]|None=None,
                    ridge_penalty:float=RIDGE)->dict:
    taxon_names=sorted(traits)
    if len(taxon_names)!=15:
        raise AuditStop("All 15 traits needed")
    for row in rows:
        if row["focal"] not in traits or row["neighbor"] not in traits:
            raise AuditStop("Taxon name absent from trait table")
        if not math.isfinite(float(row[target])):
            raise AuditStop(f"Non-finite scored target {target}")
    observed=np.asarray([r[target] for r in rows],dtype=float)
    preds={m:np.full(len(rows),np.nan) for m in models}
    test_count={}
    # Strict: a species is absent in BOTH roles from all training responses.
    # Evaluate only the previously unseen species acting as early NEIGHBOR.
    for withheld in taxon_names:
        train_pos=np.array([
            i for i,r in enumerate(rows)
            if r["focal"]!=withheld and r["neighbor"]!=withheld],dtype=int)
        test_pos=np.array([
            i for i,r in enumerate(rows) if r["neighbor"]==withheld],dtype=int)
        if len(test_pos)<10 or len(train_pos)<150:
            raise AuditStop(f"Unbalanced species holdout {withheld}: {len(train_pos)} {len(test_pos)}")
        train_species=[z for z in taxon_names if z!=withheld]
        test_count[withheld]=len(test_pos)
        train=[rows[i] for i in train_pos]
        test=[rows[i] for i in test_pos]
        for model in models:
            tx=build_features(train,traits,train_species,model,neighbor_mapping)
            ex=build_features(test,traits,train_species,model,neighbor_mapping)
            preds[model][test_pos]=ridge_prediction(tx,observed[train_pos],ex,ridge_penalty)
    if any(np.any(~np.isfinite(p)) for p in preds.values()):
        raise AuditStop("Unpredicted directed pair; every pair must be test exactly once")
    metrics={}
    by_species={}
    for withheld in taxon_names:
        mask=np.asarray([r["neighbor"]==withheld for r in rows])
        by_species[withheld]={"n_heldout_directed_pairs":int(mask.sum()),
            "mse":{m:float(np.mean((observed[mask]-p[mask])**2)) for m,p in preds.items()}}
    for model,p in preds.items():
        residual=observed-p
        metrics[model]={
            "overall_rmse":float(np.sqrt(np.mean(residual**2))),
            "pooled_pair_mse":float(np.mean(residual**2)),
            "species_equal_MSE":float(np.mean([by_species[s]["mse"][model] for s in taxon_names])),
            "species_equal_RMSE":float(np.sqrt(np.mean([by_species[s]["mse"][model] for s in taxon_names]))),
            "n_heldout_species":len(taxon_names),
            "n_test_pairs":len(observed),
        }
    for model in models:
        metrics[model]["overall_mse_skill_vs_intercept"]=float(
            metrics["intercept"]["pooled_pair_mse"]-metrics[model]["pooled_pair_mse"])
    return {"target":target,"models":list(models),"penalty_fixed":ridge_penalty,
            "n_directed_pairs":len(observed),
            "heldout_each_early_neighbor":test_count,
            "metrics":metrics,"by_species":by_species}


def paired_species_comparisons(report:dict,n_boot:int=20000,seed:int=SEED)->dict:
    cases=(
        ("both_morph_abs_gap","guild"),
        ("both_morph_abs_gap","focal_morph"),
        ("both_morph","focal_morph"),
        ("both_morph_abs_gap","both_morph"),
        ("guild_both_morph_abs_gap","guild")
    )
    species=sorted(report["by_species"])
    rng=np.random.default_rng(seed)
    out={}
    for new,baseline in cases:
        diff=np.asarray([
            report["by_species"][sp]["mse"][baseline]-
            report["by_species"][sp]["mse"][new] for sp in species])
        draws=rng.choice(diff,size=(n_boot,len(species)),replace=True).mean(axis=1)
        out[f"{new}_vs_{baseline}"]={
            "MSE_gain_positive_is_better":float(diff.mean()),
            "bootstrap95":np.quantile(draws,[.025,.975]).tolist(),
            "n_species_improved":int(np.sum(diff>0)),
            "n_species_worse":int(np.sum(diff<0)),
            "n_species_tied":int(np.sum(diff==0)),
            "n_species":len(diff),
            "species_equal_difference":{sp:float(d) for sp,d in zip(species,diff)}
        }
    return out


def shuffled_partner_traits(traits:dict,rng:np.random.Generator)->dict[str,str]:
    mapping={}
    for guild in GUILDS:
        within=sorted(x for x in traits if traits[x]["guild"]==guild)
        shuffled=rng.permutation(within).tolist()
        mapping.update(dict(zip(within,shuffled)))
    if set(mapping)!=set(traits):raise AuditStop("Invalid 5x3 within-guild mapping")
    return mapping


def mapping_placebo(rows:list[dict],traits:dict,target:str,
                    n_perm:int=199,seed:int=NEG_CONTROL_SEED)->dict:
    rng=np.random.default_rng(seed)
    observed=run_strict_loso(rows,traits,target,models=("intercept","both_morph_abs_gap"))
    observed_mse=observed["metrics"]["both_morph_abs_gap"]["species_equal_MSE"]
    null=[]
    for _ in range(n_perm):
        mapping=shuffled_partner_traits(traits,rng)
        # Keep all focal morphology unchanged and scramble only WHICH earlier
        # species receives each vector, respecting grass/forb/legume grouping.
        res=run_strict_loso(rows,traits,target,
                            models=("intercept","both_morph_abs_gap"),
                            neighbor_mapping=mapping)
        null.append(res["metrics"]["both_morph_abs_gap"]["species_equal_MSE"])
    a=np.array(null,dtype=float)
    return {
        "observed_species_equal_MSE":float(observed_mse),
        "mapping_null_median":float(np.median(a)),
        "mapping_null_95pct":np.quantile(a,[.025,.975]).tolist(),
        "mapping_fraction_at_least_as_predictive":float((1+np.sum(a<=observed_mse+1e-12))/(n_perm+1)),
        "n_shuffled_neighbor_trait_mappings":n_perm,
        "seed":seed,
        "within_each_5_species_guild":True,
        "interpretation":"Trait-label mapping placebo, not species-response randomized causal p value.",
    }


def analysis(source_dir:Path,shuffle_count:int=199)->dict:
    raw,source,_=prep(source_dir)
    matched=transform(raw)
    traits,fields=traits_for_source(source_dir,with_growth_rate=False)
    primary=run_strict_loso(matched,traits,"Cdelta_total")
    primary["paired_species_bootstrap"]=paired_species_comparisons(primary)
    primary["neighbor_trait_label_placebo"]=mapping_placebo(matched,traits,"Cdelta_total",n_perm=shuffle_count)

    sensitivity={}
    for name in ("RII_late_total","RII_sim_total"):
        res=run_strict_loso(matched,traits,name)
        sensitivity[name]={
            "metrics":res["metrics"],
            "comparisons":paired_species_comparisons(res,n_boot=20000,seed=SEED)
        }
    traits_growth,fields_growth=traits_for_source(source_dir,with_growth_rate=True)
    growth=run_strict_loso(matched,traits_growth,"Cdelta_total")
    sensitivity["including_growth_rate_derived_from_solo_control"]={
        "metrics":growth["metrics"],
        "comparisons":paired_species_comparisons(growth,n_boot=20000,seed=SEED),
        "flag":"POST HOC sensitivity: source GR is measured from control growth and may be mathematically coupled to Cdelta."
    }
    return {
        "status":"SOURCE_VERIFIED_SPECIES_OUT_PREDICTION_EXPLORATORY_NOT_CAUSAL",
        "source_doi":"10.5281/zenodo.10512662",
        "article_doi":"10.1002/ecy.4528",
        "source_audit":source,
        "source_method_note":"Focal growth duration 10 weeks in all modes; early neighbor had 13 weeks",
        "protocol":"docs/VANSTEIJN_UNSEEN_NEIGHBOR_TRANSFER_PRE_SCORE_20261124.md",
        "target_definition":"C_pre−C_sim, both benchmarked to same aged focal-alone",
        "main_traits":list(fields),
        "folds":"15 leave-one-species-out folds; exclude heldout species in both training roles; test its early-neighbor role only",
        "primary":primary,
        "predeclared_sensitivities":sensitivity,
        "limitations":[
            "Original source publication already found relationships between traits and priority effects.",
            "Each species appears in 14 or fewer directed-pair pots; no independent replication for each directed species combination.",
            "Species are not randomly sampled from flora and phylogenetic dependence is not modeled.",
            "Trait measurements derive from solo-grown control plants after 10 weeks, not ex-ante remote measurements.",
            "Study performs trait transfer within one 15-species panel, not an external experimental replication.",
            "Not a dynamical per-capita interaction alpha coefficient or an interspecific-network mediation test.",
            "Strict species-out folds have overlapping training sets; bootstrap is descriptive only.",
            "Current frozen KNK v28-v38 and Draft PRs #3/#4 unchanged."
        ]
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--mapping-permutations",type=int,default=199)
    args=parser.parse_args()
    out=analysis(args.source_dir,shuffle_count=args.mapping_permutations)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+"\n")
    print("UNSEEN_SPECIES_TRANSFER_START")
    print("SOURCE",json.dumps({"n_species":len(out["primary"]["heldout_each_early_neighbor"]),
                               "n_directed_pairs":out["primary"]["n_directed_pairs"],
                               "source_doi":out["source_doi"],"main_traits":out["main_traits"]}))
    print("PRIMARY_SCORES",json.dumps(out["primary"]["metrics"]))
    for key in ("both_morph_abs_gap_vs_guild","both_morph_abs_gap_vs_focal_morph",
                "both_morph_vs_focal_morph","both_morph_abs_gap_vs_both_morph"):
        print("PAIRWISE_MODEL",key,json.dumps(out["primary"]["paired_species_bootstrap"][key]))
    print("MAPPING_PLACEBO",json.dumps(out["primary"]["neighbor_trait_label_placebo"]))
    for key,val in out["predeclared_sensitivities"].items():
        print("SENSITIVITY",key,json.dumps(val["metrics"]))
    print("UNSEEN_SPECIES_TRANSFER_END",args.output)


if __name__=="__main__":main()
