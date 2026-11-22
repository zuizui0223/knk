#!/usr/bin/env python3
"""Exploratory winner-identity and biomass-compensation follow-up.

Use *after* the already-frozen primary analysis. Source: Dawson-Glass et al.,
Ecology (2025), https://doi.org/10.5061/dryad.sbcc2frh4.
No raw-source redistribution. Reported p values are post hoc and do not identify
niche-preemption, direct competition, or a true no-effect for rank structure.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from math import comb
from pathlib import Path

import numpy as np

from benchmarks.run_dawsonglass_storage_post_submission import (
    DEFAULT_SEED, AuditStop, bray_curtis, between_r2, centered_gower,
    load_panel,
)

SOURCE_SHA256 = "1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230"


def g_stat(table: np.ndarray) -> float:
    v = table.astype(float)
    expected = np.outer(v.sum(axis=1), v.sum(axis=0)) / v.sum()
    valid = v > 0
    return float(2 * np.sum(v[valid] * np.log(v[valid] / expected[valid])))


def exact_winner_test(winner_indices: np.ndarray, treatment: np.ndarray) -> dict:
    """Exact conditional multivariate-hypergeometric G-test of winner identities."""
    winners, inverse = np.unique(winner_indices, return_inverse=True)
    counts = np.bincount(inverse)
    n0 = int(np.count_nonzero(treatment == 0))
    total = int(len(treatment))
    table = np.stack((
        np.bincount(inverse[treatment == 0], minlength=len(winners)),
        np.bincount(inverse[treatment == 1], minlength=len(winners)),
    ))
    observed = g_stat(table)
    if np.prod(counts + 1, dtype=np.int64) > 1000000:
        raise AuditStop("Too many conditional winner tables for exact enumeration")
    tail = 0
    denominator = comb(total, n0)
    possibilities = 0
    for alloc in itertools.product(*[range(int(c) + 1) for c in counts]):
        if sum(alloc) != n0:
            continue
        possibilities += 1
        first = np.asarray(alloc)
        candidate = np.stack((first, counts - first))
        if g_stat(candidate) >= observed - 1e-10:
            ways = 1
            for n, k in zip(counts, first):
                ways *= comb(int(n), int(k))
            tail += ways
    return {
        "winner_categories": winners.tolist(),
        "table_ambient_then_warm": table.tolist(),
        "G2": observed,
        "exact_conditional_p": float(tail / denominator),
        "conditional_tables_enumerated": possibilities,
        "sample_assignments_denominator": denominator,
    }


def bh_adjust(p: np.ndarray) -> np.ndarray:
    order = np.argsort(p)
    n = len(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.clip(ranked, 0, 1)
    return result


def bootstrap_differences(x: np.ndarray, y: np.ndarray, idx: dict) -> dict:
    """Stratified whole-mesocosm bootstrap; *descriptive*, not causal."""
    rng = np.random.default_rng(20261123)
    a, w = np.flatnonzero(y == 0), np.flatnonzero(y == 1)
    a_b = rng.choice(a, (20000, len(a)), replace=True)
    w_b = rng.choice(w, (20000, len(w)), replace=True)
    r = idx["Rumex crispus"]
    g = [idx[k] for k in ["Bromus inermis", "Dactylis glomerata"]]
    totals = x.sum(axis=1)
    ranks = (x / totals[:, None]).max(axis=1)

    def mean_diff(z: np.ndarray) -> np.ndarray:
        return z[w_b].mean(axis=1) - z[a_b].mean(axis=1)

    def summarize(z: np.ndarray) -> dict:
        return {"delta": float(z[w].mean() - z[a].mean()),
                "bootstrap95": np.quantile(mean_diff(z), [0.025, 0.975]).tolist()}

    return {
        "Rumex_grams": summarize(x[:, r]),
        "two_grasses_grams": summarize(x[:, g].sum(axis=1)),
        "total_grams": summarize(totals),
        "largest_share": summarize(ranks),
        "bootstrap_resamples": 20000,
    }


def analyze(x: np.ndarray, y: np.ndarray, species: list[str],
            seed: int = DEFAULT_SEED, n_perm: int = 9999) -> dict:
    if x.shape != (24, 15) or len(species) != 15:
        raise AuditStop("Independent validation requires exactly 24 mesocosms by 15 taxa")
    if not {"Rumex crispus", "Bromus inermis", "Dactylis glomerata"}.issubset(species):
        raise AuditStop("Expected exact original species names absent")
    if np.any(np.sum(x == x.max(axis=1, keepdims=True), axis=1) != 1):
        raise AuditStop("Non-unique largest species in one or more mesocosms")
    names = list(species)
    idx = {name: j for j, name in enumerate(names)}
    proportions = x / x.sum(axis=1, keepdims=True)
    ranks = np.sort(proportions, axis=1)[:, ::-1]
    g_l = centered_gower(bray_curtis(proportions))
    g_s = centered_gower(bray_curtis(ranks))
    r2_l, r2_s = between_r2(g_l, y), between_r2(g_s, y)
    perms = np.asarray([np.random.default_rng(seed).permutation(y)]) if n_perm == 1 else None
    rng = np.random.default_rng(seed)
    perms = np.asarray([rng.permutation(y) for _ in range(n_perm)])
    deltas = np.asarray([
        between_r2(g_l, perm) - between_r2(g_s, perm)
        for perm in perms
    ])
    observed_delta = float(r2_l - r2_s)
    winner_indices = x.argmax(axis=1)
    winner = exact_winner_test(winner_indices, y)

    d_rel = proportions[y == 1].mean(0) - proportions[y == 0].mean(0)
    d_g = x[y == 1].mean(0) - x[y == 0].mean(0)
    nul_rel = np.asarray([
        proportions[z == 1].mean(0) - proportions[z == 0].mean(0) for z in perms
    ])
    nul_g = np.asarray([x[z == 1].mean(0) - x[z == 0].mean(0) for z in perms])
    p_rel = (1 + np.count_nonzero(np.abs(nul_rel) >= np.abs(d_rel) - 1e-12, axis=0)) / (n_perm + 1)
    p_g = (1 + np.count_nonzero(np.abs(nul_g) >= np.abs(d_g) - 1e-12, axis=0)) / (n_perm + 1)
    q_rel = bh_adjust(p_rel)
    q_g = bh_adjust(p_g)
    mean_a = proportions[y == 0].mean(0)
    mean_w = proportions[y == 1].mean(0)
    mean_sort_a = ranks[y == 0].mean(0)
    mean_sort_w = ranks[y == 1].mean(0)

    def bc(a: np.ndarray, b: np.ndarray) -> float:
        return float(np.abs(a - b).sum() / (a + b).sum())

    species_effects = [
        {"species": name, "share_ambient": float(mean_a[j]),
         "share_warm": float(mean_w[j]), "share_delta": float(d_rel[j]),
         "share_permutation_p": float(p_rel[j]), "share_BH_q": float(q_rel[j]),
         "grams_ambient": float(x[y == 0, j].mean()),
         "grams_warm": float(x[y == 1, j].mean()),
         "grams_delta": float(d_g[j]),
         "grams_permutation_p": float(p_g[j]), "grams_BH_q": float(q_g[j])}
        for j, name in enumerate(names)
    ]
    return {
        "status": "EXPLORATORY_AFTER_PRIMARY_ANALYSIS",
        "design": {"unit": "whole mesocosm", "n": len(y),
                   "groups": {"ambient": int((y == 0).sum()),
                              "warm_order": int((y == 1).sum())},
                   "permutations": n_perm, "seed": seed,
                   "null": "unblocked label exchangeability",
                   "multiple_testing": "Benjamini-Hochberg within 15 species"},
        "paired_representation_contrast": {
            "labeled_R2": float(r2_l), "sorted_R2": float(r2_s),
            "delta_R2": observed_delta,
            "permutation_null_median": float(np.median(deltas)),
            "permutation_null_95pct": np.quantile(deltas, [0.025, 0.975]).tolist(),
            "permutation_p_upper_under_global_no_treatment_null":
            float((1 + np.count_nonzero(deltas >= observed_delta - 1e-12)) / (n_perm + 1)),
            "caution": "Not an equivalence test or a mechanistic identification",
        },
        "winner_identity": {
            **winner,
            "winner_species_order": [names[i] for i in winner["winner_categories"]],
        },
        "between_treatment_means": {
            "BrayCurtis_species_labeled": bc(mean_a, mean_w),
            "BrayCurtis_species_unlabeled": bc(mean_sort_a, mean_sort_w),
            "caution": "Descriptive; distances are not percentages of information",
        },
        "bootstrap_key_effects": bootstrap_differences(x, y, idx),
        "individual_species_exploratory": species_effects,
        "boundaries": [
            "The original paper already established assembly timing versus species biomass",
            "Species-level tests and winner test are post hoc",
            "Mean gram-level replacement does not imply direct competition or energetic compensation",
            "Initial height, growth duration and arrival schedule are not mechanistically separable",
            "Current frozen KNK v28-v38 manuscript is unchanged",
        ],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assembly-csv", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    raw = args.assembly_csv.read_bytes()
    sha = hashlib.sha256(raw).hexdigest()
    if sha != SOURCE_SHA256:
        raise AuditStop("Source SHA mismatch, audit new raw-data version before scoring")
    x, y, audit = load_panel(args.assembly_csv)
    result = analyze(x, y, audit["species"])
    result["source_sha256"] = sha
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
    print("Saved exploratory report to", args.output)
    print("Paired contrast:", result["paired_representation_contrast"])
    print("Winner test:", result["winner_identity"])


if __name__ == "__main__":
    main()
