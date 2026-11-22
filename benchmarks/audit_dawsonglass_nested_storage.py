#!/usr/bin/env python3
"""Exploratory hierarchical memory and 3,003 equal-size partition-placebo audit.

Source: Dawson-Glass et al. (2025), DOI 10.5061/dryad.sbcc2frh4.
This script is explicitly POST HOC and independent of frozen KNK v28–v38.
Poaceae genus membership uses external botanical taxonomy; the Dryad
species_info.csv Functional field has not been inspected here.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from benchmarks.run_dawsonglass_storage_post_submission import (
    AuditStop, between_r2, bray_curtis, centered_gower, load_panel, permutation_test,
    fitted_logistic_probability
)

SOURCE_SHA = "1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230"
GRASSES = {
    "Andropogon gerardii", "Bromus inermis", "Dactylis glomerata",
    "Sorgastrum nutans", "Tridens flavus",
}


def grass_mask(names: list[str]) -> np.ndarray:
    if len(names) != 15 or not GRASSES.issubset(set(names)):
        raise AuditStop("Source-taxonomy mismatch (includes source spelling Sorgastrum)")
    g = np.array([s in GRASSES for s in names])
    if g.sum() != 5:
        raise AuditStop("Expected exactly 5 grass and 10 forb species")
    return g


def within_group(x: np.ndarray, mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sub = x[:, mask]
    if np.any(sub.sum(axis=1) <= 0):
        raise AuditStop("A selected guild has zero measured mass in a pot")
    relative = sub / sub.sum(axis=1, keepdims=True)
    return relative, np.sort(relative, axis=1)[:, ::-1]


def observed_r2(z: np.ndarray, y: np.ndarray) -> float:
    return float(between_r2(centered_gower(bray_curtis(z)), y))


def grouping_placebo(x: np.ndarray, y: np.ndarray, names: list[str],
                     real_forb_rank_r2: float) -> dict:
    """All C(15,5) partitions, and a narrower dominant-species-conditioned audit.

    The enumerated grouping fraction is NOT a treatment-randomization p-value.
    """
    name_index = {n: i for i, n in enumerate(names)}
    b, d, r = [name_index[k] for k in
               ("Bromus inermis", "Dactylis glomerata", "Rumex crispus")]
    full_rank = []
    conditioned = []
    functional_share_r2 = []
    for chosen in itertools.combinations(range(15), 5):
        g = np.zeros(15, dtype=bool)
        g[list(chosen)] = True
        forb, sorted_forb = within_group(x, ~g)
        r2 = observed_r2(sorted_forb, y)
        full_rank.append(r2)
        if b in chosen and d in chosen and r not in chosen:
            conditioned.append(r2)
        both = np.column_stack([x[:, g].sum(axis=1), x[:, ~g].sum(axis=1)])
        both = both / both.sum(axis=1, keepdims=True)
        functional_share_r2.append(observed_r2(both, y))
    all_r2 = np.array(full_rank)
    conditional_r2 = np.array(conditioned)
    actual_g = grass_mask(names)
    actual_group = np.column_stack([
        x[:, actual_g].sum(axis=1), x[:, ~actual_g].sum(axis=1)
    ])
    actual_group /= actual_group.sum(axis=1, keepdims=True)
    actual_group_r2 = observed_r2(actual_group, y)
    return {
        "number_partitions": len(all_r2),
        "observed_real_forb_rank_R2": real_forb_rank_r2,
        "number_rank_at_least_real_forb": int(np.sum(all_r2 >= real_forb_rank_r2 - 1e-12)),
        "all_partitions_upper_fraction": float(np.mean(all_r2 >= real_forb_rank_r2 - 1e-12)),
        "conditional_rule": "Bromus and Dactylis in grass; Rumex in forb",
        "conditional_partitions": len(conditional_r2),
        "conditional_rank_at_least_real_forb": int(np.sum(conditional_r2 >= real_forb_rank_r2 - 1e-12)),
        "conditional_upper_fraction": float(np.mean(conditional_r2 >= real_forb_rank_r2 - 1e-12)),
        "real_two_functional_group_R2": actual_group_r2,
        "group_R2_all_partition_upper_fraction":
            float(np.mean(np.array(functional_share_r2) >= actual_group_r2 - 1e-12)),
        "claim_boundary": "These are exploratory mapping ranks, NOT causal or randomization p-values."
    }


def holm_adjust(p: dict[str, float]) -> dict[str, float]:
    names = sorted(p, key=p.get)
    q = np.maximum.accumulate([p[k] * (len(p) - i) for i, k in enumerate(names)])
    return {k: float(np.clip(v, 0, 1)) for k, v in zip(names, q)}


def analyze(x: np.ndarray, y: np.ndarray, names: list[str],
            n_permutations: int = 9999, seed: int = 20261122,
            mapping_placebos: bool = True) -> dict:
    if x.shape != (24, 15) or np.count_nonzero(y) != 12:
        raise AuditStop("Expected 24 whole pots by 15 species and balanced 12/12 arms")
    grass = grass_mask(names)
    forb = ~grass
    perm_rng = np.random.default_rng(seed)
    perms = np.array([perm_rng.permutation(y) for _ in range(n_permutations)])
    P = x / x.sum(axis=1, keepdims=True)
    FG = np.column_stack([P[:, grass].sum(axis=1), P[:, forb].sum(axis=1)])
    f, fr = within_group(x, forb)
    g, gr = within_group(x, grass)
    f9_mask = forb & (np.array(names) != "Rumex crispus")
    f9, f9r = within_group(x, f9_mask)
    remain = np.array(names) != "Rumex crispus"
    other = x[:, remain]
    remaining_grass = grass[remain]
    after = np.column_stack([
        other[:, remaining_grass].sum(axis=1),
        other[:, ~remaining_grass].sum(axis=1)
    ])
    after = after / after.sum(axis=1, keepdims=True)
    matrices = {
        "full15_named": P,
        "full15_sorted": np.sort(P, axis=1)[:, ::-1],
        "functional2_share": FG,
        "forb10_named_internal": f,
        "forb10_sorted_internal": fr,
        "grass5_named_internal": g,
        "grass5_sorted_internal": gr,
        "forb9_omit_Rumex_named_internal": f9,
        "forb9_omit_Rumex_sorted_internal": f9r,
        "functional2_share_omit_Rumex": after,
        "forb10_top1_within_forb": fr[:, :1],
    }
    stats = {name: permutation_test(matrix, y, perms)
             for name, matrix in matrices.items()}
    # Fixed-penalty whole-pot LOO; subgroup selection was post hoc, so this
    # is a descriptive reconstruction, NOT independent external validation.
    cv = {}
    for name in ("full15_sorted", "functional2_share",
                 "forb10_sorted_internal", "forb9_omit_Rumex_named_internal",
                 "forb9_omit_Rumex_sorted_internal", "grass5_sorted_internal"):
        matrix = matrices[name]
        pred, baseline_loss = [], []
        for i in range(len(y)):
            keep = np.arange(len(y)) != i
            pi = fitted_logistic_probability(matrix[keep], y[keep], matrix[i], penalty=4.0)
            prior = (np.count_nonzero(y[keep]) + 0.5) / (len(y[keep]) + 1)
            pred.append(pi)
            baseline_loss.append(-y[i]*np.log(prior)-(1-y[i])*np.log(1-prior))
        pred = np.asarray(pred)
        loss = -y*np.log(pred)-(1-y)*np.log(1-pred)
        cv[name] = {"n_holdouts":len(y), "accuracy":float(np.mean((pred>=0.5)==y)),
                    "correct":int(np.sum((pred>=0.5)==y)),
                    "logloss":float(np.mean(loss)),
                    "baseline_logloss":float(np.mean(baseline_loss)),
                    "logloss_gain_vs_prior":float(np.mean(baseline_loss)-np.mean(loss))}
    q = holm_adjust({name: float(v["p_upper"]) for name, v in stats.items()})
    for name, v in stats.items():
        v["posthoc_holm_adjusted_p_11_tests"] = q[name]
    a, w = np.flatnonzero(y == 0), np.flatnonzero(y == 1)
    boot_rng = np.random.default_rng(20261124)
    nb = 20000
    bdelta = (
        FG[boot_rng.choice(w, (nb, len(w)), replace=True), 1].mean(axis=1)
        - FG[boot_rng.choice(a, (nb, len(a)), replace=True), 1].mean(axis=1)
    )
    report = {
        "status": "EXPLORATORY_ONLY_NOT_FROZEN_KNK_MAIN",
        "scope": "24 same pot experimental units, no new independent replication",
        "treatment_coding": "0 Ambient Order, 1 Warm Order",
        "grass_source_species_names": sorted(GRASSES),
        "forb_source_species_names": [s for s in names if s not in GRASSES],
        "functional_mapping_caveat":
            "Poaceae genera externally verified; original Dryad species_info Functional field not audited",
        "analysis": {"n_permutations": n_permutations, "seed": seed,
                     "dissimilarity": "Bray-Curtis",
                     "grouping": "5 grasses, 10 broadleaf forbs",
                     "within_group": "within-pot renormalized to each guild's mass",
                     "multiple_testing": "Holm across 11 exploratory endpoints"},
        "results": stats,
        "leave_one_pot_out_prediction": cv,
        "means": {
            "forb_share_ambient": float(FG[a, 1].mean()),
            "forb_share_warm": float(FG[w, 1].mean()),
            "forb_share_difference": float(FG[w, 1].mean() - FG[a, 1].mean()),
            "forb_share_bootstrap_95pct": np.quantile(bdelta, [0.025, 0.975]).tolist(),
            "max_share_inside_forb_ambient": float(fr[a, 0].mean()),
            "max_share_inside_forb_warm": float(fr[w, 0].mean()),
            "remaining_forbs_share_without_Rumex_ambient": float(after[a, 1].mean()),
            "remaining_forbs_share_without_Rumex_warm": float(after[w, 1].mean()),
        },
        "limitations": [
            "11 nested endpoints tested post hoc; Holm correction does not convert exploration into preregistration",
            "R2 and E change under renormalization and are not amounts of mutual information",
            "p>0.05 does not prove an endpoint's historical signal is absent",
            "No identification of niche pre-emption, relative growth duration or functional ecosystem output",
            "Mapping placebo ranks are descriptive, not randomization p-values",
            "Keep legacy v28-v38 and existing Ecology Letters letter unchanged",
        ],
    }
    if mapping_placebos:
        report["mapping_placebo"] = grouping_placebo(
            x, y, names, stats["forb10_sorted_internal"]["R2_observed"]
        )
    return report


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assembly-csv", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    sha = hashlib.sha256(args.assembly_csv.read_bytes()).hexdigest()
    if sha != SOURCE_SHA:
        raise AuditStop("Source file bytes do not match audited source SHA256")
    x, y, audit = load_panel(args.assembly_csv)
    report = analyze(x, y, audit["species"])
    report["source_sha256"] = sha
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print(f"Post hoc nested-audit written: {args.output}")


if __name__ == "__main__":
    main()
