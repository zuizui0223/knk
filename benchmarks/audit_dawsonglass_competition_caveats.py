#!/usr/bin/env python3
"""Exploratory source-locked audit: species selection and resource-budget caveats.

Post hoc follow-up; do NOT interpret correlations as competition. This never
modifies the frozen KNK v28-v38 results or the current submission manuscript.
Raw source: Dawson-Glass et al. (2025), 10.5061/dryad.sbcc2frh4.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np

from benchmarks.run_dawsonglass_storage_post_submission import (
    AuditStop, load_panel, permutation_test,
)

EXPECTED_SHA = "1514d4501a0bdf7a61ace6c250566c5a0b8e4a24a0f4737646d182ebad7b8230"
TOP3 = ("Rumex crispus", "Bromus inermis", "Dactylis glomerata")


def centered_by_arm(z: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Remove arm means without treating species observations as independent."""
    out = np.asarray(z, dtype=float).copy()
    for arm in (0, 1):
        idx = y == arm
        out[idx] -= out[idx].mean(axis=0)
    return out


def group_retention(x: np.ndarray, y: np.ndarray, names: list[str],
                    removed: tuple[str, ...], perms: np.ndarray) -> dict:
    keep = [j for j, n in enumerate(names) if n not in removed]
    xx = x[:, keep]
    if np.any(xx.sum(axis=1) <= 0):
        raise AuditStop("Source has zero remaining mass after selected species removal")
    pp = xx / xx.sum(axis=1, keepdims=True)
    rank = np.sort(pp, axis=1)[:, ::-1]
    return {
        "removed": list(removed),
        "remaining_taxa": len(keep),
        "mean_mass_by_arm": [float(xx[y == a].sum(axis=1).mean()) for a in (0, 1)],
        "raw_grams": permutation_test(xx, y, perms),
        "within_remaining_proportions": permutation_test(pp, y, perms),
        "sorted_rank_proportions": permutation_test(rank, y, perms),
        "caution": "Post hoc removal renormalizes remaining species; these E values are new estimands.",
    }


def audit(x: np.ndarray, y: np.ndarray, species: list[str],
          n_perm: int = 49999, seed: int = 20261124,
          retention_perm: int = 9999, retention_seed: int = 20261122) -> dict:
    if len(y) != 24 or x.shape != (24, 15) or set(np.unique(y)) != {0, 1}:
        raise AuditStop("Expect source 24-pot, 15-species, 2-arm panel")
    if not set(TOP3).issubset(species):
        raise AuditStop("Source spellings differ from expected focal species")
    index = {n: j for j, n in enumerate(species)}
    rumex = x[:, index[TOP3[0]]]
    grasses = x[:, index[TOP3[1]]] + x[:, index[TOP3[2]]]
    total = x.sum(axis=1)
    others = total - rumex - grasses

    rw = centered_by_arm(rumex, y)
    gw = centered_by_arm(grasses, y)
    tw = centered_by_arm(total, y)
    ow = centered_by_arm(others, y)
    observed_r = float(np.corrcoef(rw, gw)[0, 1])
    covariance_rg = float(np.cov(rw, gw, ddof=1)[0, 1])
    covariance_rt = float(np.cov(rw, tw, ddof=1)[0, 1])
    variance_r = float(np.var(rw, ddof=1))
    covariance_ro = float(np.cov(rw, ow, ddof=1)[0, 1])
    assert abs(covariance_rg-(covariance_rt-variance_r-covariance_ro)) < 1e-8

    other_names = [n for n in species if n != TOP3[0]]
    residuals = np.column_stack([centered_by_arm(x[:, index[n]], y)
                                 for n in other_names])
    pairs = list(itertools.combinations(range(len(other_names)), 2))
    pair_residuals = np.column_stack(
        [residuals[:, a]+residuals[:, b] for a, b in pairs])
    denom = np.linalg.norm(rw) * np.linalg.norm(pair_residuals, axis=0)
    if np.any(denom < 1e-12):
        raise AuditStop("Degenerate species pair; source-specific null needs audit")
    observed_pairs = (rw @ pair_residuals) / denom
    pair_keys = [tuple(sorted((other_names[a], other_names[b]))) for a,b in pairs]
    selected = tuple(sorted(TOP3[1:]))
    selected_index = pair_keys.index(selected)
    assert np.isclose(observed_r, observed_pairs[selected_index], atol=1e-10)

    arm0, arm1 = np.flatnonzero(y == 0), np.flatnonzero(y == 1)
    rng = np.random.default_rng(seed)
    selected_tail = 0
    pair_min_tail = 0
    # Shuffle *Rumex* whole-pot residuals within arms, keeping competitor
    # covariance and observed arm-specific distributions unchanged.
    for start in range(0, n_perm, 1000):
        n = min(1000, n_perm-start)
        rp = np.empty((n, len(y)))
        rp[:, arm0] = rw[arm0][np.argsort(
            rng.random((n, len(arm0))), axis=1)]
        rp[:, arm1] = rw[arm1][np.argsort(
            rng.random((n, len(arm1))), axis=1)]
        corr = (rp @ pair_residuals) / denom
        selected_tail += int(np.count_nonzero(
            corr[:, selected_index] <= observed_r + 1e-12))
        pair_min_tail += int(np.count_nonzero(
            corr.min(axis=1) <= observed_r + 1e-12))

    fractions = x / total[:, None]
    d = fractions[y == 1].mean(axis=0) - fractions[y == 0].mean(axis=0)
    total_L1 = float(np.abs(d).sum())
    top3_turnover = float(sum(abs(d[index[n]]) for n in TOP3) / total_L1)

    # Sensitivity-only complete and post-selected subset tests. Keep one
    # shared, source-blind treatment label stream across representations.
    r_rng = np.random.default_rng(retention_seed)
    perms = np.asarray([r_rng.permutation(y) for _ in range(retention_perm)])
    retention = {
        "full_15": group_retention(x, y, species, (), perms),
        "excluding_Rumex": group_retention(x, y, species, (TOP3[0],), perms),
        "excluding_three_largest_shifters": group_retention(x, y, species, TOP3, perms),
    }

    return {
        "status": "EXPLORATORY_AND_SELECTION_ADJUSTED_NOT_CAUSAL",
        "n_pots": len(y),
        "rumex_grass_association": {
            "within_arm_pooled_r": observed_r,
            "r_ambient": float(np.corrcoef(rumex[arm0], grasses[arm0])[0, 1]),
            "r_warm_order": float(np.corrcoef(rumex[arm1], grasses[arm1])[0, 1]),
            "r_Rumex_other12": float(np.corrcoef(rw, ow)[0, 1]),
            "one_sided_p_selected_grasses": (1+selected_tail)/(n_perm+1),
            "one_sided_maxT_familywise_p_over_91_pairs":
                (1+pair_min_tail)/(n_perm+1),
            "n_pair_candidates": len(pairs),
            "n_permutations": n_perm, "seed": seed,
            "caution": "Family correction holds Rumex fixed. Rumex was also selected post hoc.",
        },
        "biomass_budget_identity": {
            "cov_R_G": covariance_rg, "cov_R_Total": covariance_rt,
            "Var_R": variance_r, "cov_R_Other12": covariance_ro,
            "identity": "Cov(R,G)=Cov(R,Total)-Var(R)-Cov(R,Other12)",
            "note": "Exact covariance accounting is algebraic, not evidence for competition.",
        },
        "mean_composition_shift": {
            "bray_curtis_between_mean_species_proportions": total_L1/2,
            "fraction_of_absolute_species_shift_in_top3": top3_turnover,
            "source_species_mean_share_differences":
                {n: float(d[index[n]]) for n in species},
        },
        "post_hoc_retention_sensitivities": retention,
        "limitations": [
            "Two grasses and Rumex were chosen after observing outcomes.",
            "The 91-pair family correction is conditional on choosing Rumex.",
            "Total pot biomass constraints can yield negative component covariance.",
            "Source mesocosms do not isolate growth duration from biotic interactions.",
            "Initial-size measurement is post timing manipulation; do not overadjust.",
            "Removed-species subsets are not new independent experimental replications.",
        ],
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--assembly-csv", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    ap.add_argument("--permutations", type=int, default=49999)
    args = ap.parse_args()
    source_hash = hashlib.sha256(args.assembly_csv.read_bytes()).hexdigest()
    if source_hash != EXPECTED_SHA:
        raise AuditStop("Unexpected source SHA: inspect version before scoring")
    x, y, info = load_panel(args.assembly_csv)
    report = audit(x, y, info["species"], n_perm=args.permutations)
    report["source_sha256"] = source_hash
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print("Source audit complete: 24 pots, 15 species")
    print("Rumex vs grass:", report["rumex_grass_association"])
    print("Dominant share of mean turnover:", report["mean_composition_shift"]["fraction_of_absolute_species_shift_in_top3"])
    print("Saved:", args.output)


if __name__ == "__main__":
    main()
