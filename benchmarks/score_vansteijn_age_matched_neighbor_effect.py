#!/usr/bin/env python3
"""Age-matched net neighbor effect: simultaneous vs 3-week-preestablished neighbor.

Source method correction committed BEFORE this analysis:
docs/VANSTEIJN_FOCAL_AGE_METHOD_CORRECTION_AND_NET_NEIGHBOR_EFFECT_20261123.md

van Steijn et al. 2025 Ecology, DOI 10.1002/ecy.4528. The focal individual in
all treatments, including alone and after another plant has been established,
grew for 10 weeks. The *neighbor* in sequential pots grew for 13 weeks.

This identifies history's change in realized net neighbor effect on focal
biomass at a fixed focal age, not a dynamic per-capita alpha_ij mechanism.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from benchmarks.score_vansteijn_partner_effects import (
    prep, focal_cluster_test, neighbor_consistency,
    reciprocal_stats, guild_matrix, predictive_loo_unordered,
)
from benchmarks.audit_vansteijn_partner_specificity import (
    conditional_species_identity, guild_mapping_placebo,
    maxT_neighbor_extremes,
)


def transform(rows:list[dict]):
    out=[]
    for r in rows:
        q=dict(r)
        for organ in ("total","above","below"):
            a=r[f"biomass_alone_{organ}"]
            s=r[f"biomass_sim_{organ}"]
            p=r[f"biomass_pre_{organ}"]
            # Both compared to the same AGE-MATCHED 10-week solo control.
            if s+a<=0 or p+a<=0:
                raise ValueError("Undefined neighbor net effect")
            sim=(s-a)/(s+a)
            pre=(p-a)/(p+a)
            if not math.isclose(sim,r[f"RII_sim_{organ}"],abs_tol=1e-12):
                raise ValueError("Cannot reconcile solo matched source RII")
            q[f"C_sim_{organ}"]=sim
            q[f"C_pre_{organ}"]=pre
            q[f"Cdelta_{organ}"]=pre-sim
        q["Cdelta_above_minus_below"]=q["Cdelta_above"]-q["Cdelta_below"]
        out.append(q)
    return out


def sign_changes(rows:list[dict], organ="total"):
    x=f"C_sim_{organ}"
    y=f"C_pre_{organ}"
    count={
        "both_negative":0,"both_positive":0,
        "positive_sim_to_negative_pre":0,"negative_sim_to_positive_pre":0,
        "sim_exact_zero":0,"pre_exact_zero":0,
    }
    transitions=[]
    for r in rows:
        v=r[x];z=r[y]
        if v==0:count["sim_exact_zero"]+=1
        if z==0:count["pre_exact_zero"]+=1
        if v<0 and z<0:count["both_negative"]+=1
        if v>0 and z>0:count["both_positive"]+=1
        if v>0 and z<0:
            count["positive_sim_to_negative_pre"]+=1
            transitions.append({"focal":r["focal"],"neighbor":r["neighbor"],
                                "sim_net_RII":v,"pre_net_RII":z})
        if v<0 and z>0:
            count["negative_sim_to_positive_pre"]+=1
            transitions.append({"focal":r["focal"],"neighbor":r["neighbor"],
                                "sim_net_RII":v,"pre_net_RII":z})
    return {"total_directed_pairs":len(rows),"counts":count,
            "transitions":transitions[:30],
            "warning":"One pot per directed species pair; a sign shift is observed net biomass change, not proof of facilitation or mechanism."}


def run(source_dir:Path):
    data,audit,_=prep(source_dir)
    rows=transform(data)
    response=focal_cluster_test(rows,"Cdelta_total")
    shoot_root=focal_cluster_test(rows,"Cdelta_above_minus_below")
    return {
        "status":"SOURCE_VERIFIED_HISTORY_DEPENDENT_AGE_MATCHED_NET_NEIGHBOR_EFFECT",
        "source":"van Steijn et al. 2025, Ecology, DOI 10.1002/ecy.4528",
        "methods_correction":"docs/VANSTEIJN_FOCAL_AGE_METHOD_CORRECTION_AND_NET_NEIGHBOR_EFFECT_20261123.md",
        "focal_growing_duration_weeks":10,
        "preestablished_neighbor_growing_duration_weeks":13,
        "source_audit":audit,
        "paired_net_neighbor_effect":{
            "definition":"C_pre - C_sim, both relative to 10-week solo focal mean",
            "focal_cluster_delta":response,
            "stage_means":{
                stage:{
                    "mean_over_directed_pairs":float(np.mean([r[f"C_{stage}_total"] for r in rows])),
                    "fraction_negative":float(np.mean([r[f"C_{stage}_total"]<0 for r in rows])),
                    "fraction_positive":float(np.mean([r[f"C_{stage}_total"]>0 for r in rows])),
                } for stage in ("sim","pre")
            },
            "sign_changes":sign_changes(rows),
            "shoot_minus_root_delta":shoot_root,
        },
        "neighbor_consistency":{
            "after_focal_centering":neighbor_consistency(
                rows,key="Cdelta_total",n_perm=9999,seed=20261123),
            "after_focal_and_guild_adjustment":conditional_species_identity(
                rows,key="Cdelta_total",n_perm=9999,seed=20261124),
            "guild_specificity_random_5_5_5_partition":guild_mapping_placebo(
                rows,key="Cdelta_total",n_perm=9999,seed=20261125),
            "all_neighbor_maxT_mapping":maxT_neighbor_extremes(
                rows,key="Cdelta_total",n_perm=9999,seed=20261129),
        },
        "partner_effect_predictive":predictive_loo_unordered(rows,key="Cdelta_total",penalty=4.0),
        "reciprocal_ordered_pairs":reciprocal_stats(rows,key="Cdelta_total"),
        "guild_3x3_delta":guild_matrix(rows,key="Cdelta_total"),
        "interpretation_bounds":[
            "The focal was grown 10 weeks in solo, simultaneous and sequential treatments; earlier neighbor grew 13 weeks.",
            "C_pre−C_sim is the effect of having a pre-established neighbor vs a simultaneous neighbor on focal final biomass relative to a matched solo standard.",
            "Only one pot per directed pair; focal-species cluster uncertainty is exploratory and taxa are not a random phylogenetic sample.",
            "The total history effect combines neighbor initial size, competition, spatial preemption and conditioning; no dynamical alpha or mediation is identified.",
            "The original 2025 study already established competitive effects and traits; this is a KNK-specific decomposition.",
            "No change to the frozen KNK v28-v38 Ecology Letters letter."
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    r=run(args.source_dir)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("AGE_MATCHED_NET_NEIGHBOR_START")
    print("NET_NEIGHBOR_FOCAL_DELTA",json.dumps({
        k:v for k,v in r["paired_net_neighbor_effect"]["focal_cluster_delta"].items()
        if k!="focal_species_means"}))
    print("NET_NEIGHBOR_STAGES",json.dumps(r["paired_net_neighbor_effect"]["stage_means"]))
    print("NET_NEIGHBOR_SIGN_CHANGE",json.dumps({
        k:v for k,v in r["paired_net_neighbor_effect"]["sign_changes"].items()
        if k!="transitions"}))
    print("NET_NEIGHBOR_IDENTITY",json.dumps({
        k:({kk:vv for kk,vv in v.items() if kk not in
            ("neighbor_adjusted_means","neighbor_specific_post_selection_results")})
        for k,v in r["neighbor_consistency"].items()}))
    print("NET_NEIGHBOR_PREDICTIVE",json.dumps(r["partner_effect_predictive"]))
    print("NET_NEIGHBOR_GUILD",json.dumps(r["guild_3x3_delta"]))
    print("AGE_MATCHED_NET_NEIGHBOR_DONE",args.output)


if __name__=="__main__":main()
