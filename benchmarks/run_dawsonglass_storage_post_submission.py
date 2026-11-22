#!/usr/bin/env python3
"""Post-submission KNK external validation. Never alters frozen v28-v38 endpoints.

Input is the author-released Dryad CSV; raw data are not redistributed.
Provenance: https://doi.org/10.5061/dryad.sbcc2frh4
Original paper: https://doi.org/10.1002/ecy.4504

This script deliberately refuses to score an incomplete/malformed 24-pot x 15-species
panel rather than silently interpreting absent observations as zero biomass.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


EXPECTED_POTS = 24
EXPECTED_SPECIES = 15
EXPECTED_PER_ARM = 12
REQUIRED = {"Pot", "Species", "Treatment", "Biomass"}
DEFAULT_SEED = 20261122
DEFAULT_PERMUTATIONS = 9999


class AuditStop(ValueError):
    """Source panel invalid or underspecified for the frozen analysis contract."""


def load_panel(csv_path: Path) -> tuple[np.ndarray, np.ndarray, dict]:
    raw = csv_path.read_bytes()
    try:
        decoded = raw.decode("utf-8-sig")
    except UnicodeDecodeError as exc:
        raise AuditStop("CSV is not UTF-8; inspect source before recoding") from exc

    rows = list(csv.DictReader(decoded.splitlines()))
    if not rows:
        raise AuditStop("Empty input")
    if not REQUIRED.issubset(rows[0]):
        raise AuditStop(f"Missing required columns: {sorted(REQUIRED - set(rows[0]))}")
    pot_treatment: dict[str, str] = {}
    biomass: dict[tuple[str, str], float] = {}
    duplicate_keys: list[tuple[str, str]] = []
    seen_species: set[str] = set()
    for i, row in enumerate(rows, 2):
        pot, species, treatment = (row[k].strip() for k in ("Pot", "Species", "Treatment"))
        if not pot or not species or not treatment:
            raise AuditStop(f"Blank pot/species/treatment, input row {i}")
        try:
            value = float(row["Biomass"])
        except (ValueError, TypeError) as exc:
            raise AuditStop(f"Missing/non-numeric biomass, row {i}") from exc
        if not math.isfinite(value) or value < 0:
            raise AuditStop(f"Nonfinite/negative biomass, row {i}")
        if pot in pot_treatment and pot_treatment[pot] != treatment:
            raise AuditStop(f"Mixed treatments within pot {pot!r}")
        pot_treatment[pot] = treatment
        key = (pot, species)
        if key in biomass:
            duplicate_keys.append(key)
        biomass[key] = biomass.get(key, 0.0) + value
        seen_species.add(species)
    if duplicate_keys:
        raise AuditStop(
            f"Repeated pot x species records ({len(duplicate_keys)}). "
            "Source README describes individual plants; audit multiplicity before summing."
        )
    pots = sorted(pot_treatment)
    species = sorted(seen_species)
    counts = Counter(pot_treatment.values())
    if len(pots) != EXPECTED_POTS or len(species) != EXPECTED_SPECIES:
        raise AuditStop(f"Expected 24 pots x 15 species, got {len(pots)} x {len(species)}")
    if len(counts) != 2 or sorted(counts.values()) != [EXPECTED_PER_ARM] * 2:
        raise AuditStop(f"Expected 12 pots per treatment, got {dict(counts)}")
    missing = [(p, s) for p in pots for s in species if (p, s) not in biomass]
    if missing:
        raise AuditStop(
            f"Missing {len(missing)} pot x species rows (first: {missing[:3]}). "
            "Do not turn unrecorded entries into zero without a source audit."
        )
    x = np.asarray([[biomass[p, s] for s in species] for p in pots], dtype=float)
    if np.any(x.sum(axis=1) <= 0):
        raise AuditStop("A pot has zero aggregate biomass")
    arms = sorted(counts)
    y = np.asarray([int(pot_treatment[p] == arms[1]) for p in pots], dtype=np.int8)
    return x, y, {
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "source_file": csv_path.name,
        "rows": len(rows),
        "n_pots": len(pots),
        "n_species": len(species),
        "treatment_counts": dict(counts),
        "treatment_coding": {"0": arms[0], "1": arms[1]},
        "pots": pots,
        "species": species,
        "zero_biomass_cells": int(np.count_nonzero(x == 0.0)),
        "duplicate_pot_species": len(duplicate_keys),
        "missing_pot_species": len(missing),
        "status": "PASS",
    }


def bray_curtis(x: np.ndarray) -> np.ndarray:
    numerator = np.abs(x[:, None, :] - x[None, :, :]).sum(axis=2)
    denominator = (x[:, None, :] + x[None, :, :]).sum(axis=2)
    return np.divide(
        numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0
    )


def centered_gower(d: np.ndarray) -> np.ndarray:
    n = d.shape[0]
    j = np.eye(n) - np.ones((n, n)) / n
    return -0.5 * j @ (d * d) @ j


def between_r2(g: np.ndarray, y: np.ndarray) -> float:
    if not set(np.unique(y)) == {0, 1}:
        raise ValueError("Expected two nonempty groups")
    total = float(np.trace(g))
    if total <= 1e-12:
        return 0.0
    by_group = sum(
        float(g[np.ix_(idx, idx)].sum()) / len(idx)
        for level in (0, 1)
        for idx in [np.flatnonzero(y == level)]
    )
    return float(by_group / total)


def permutation_test(x: np.ndarray, y: np.ndarray,
                     permutations: np.ndarray) -> dict:
    g = centered_gower(bray_curtis(x))
    observed = between_r2(g, y)
    null = np.asarray([between_r2(g, yy) for yy in permutations], dtype=float)
    return {
        "R2_observed": observed,
        "R2_null_median": float(np.median(null)),
        "R2_null_95pct": [float(q) for q in np.quantile(null, [0.025, 0.975])],
        "E": float(observed - np.median(null)),
        "p_upper": float((1 + np.count_nonzero(null >= observed - 1e-12)) /
                         (len(null) + 1)),
        "permutations": len(null),
    }


def fitted_logistic_probability(train_x: np.ndarray, train_y: np.ndarray,
                                held_out: np.ndarray, penalty: float = 4.0) -> float:
    """Fixed-L2 logistic model, fold-local scaling, Newton steps; no model selection."""
    mu = train_x.mean(axis=0)
    sd = train_x.std(axis=0)
    sd[sd < 1e-8] = 1.0
    a = np.column_stack([np.ones(len(train_x)), (train_x - mu) / sd])
    new_a = np.r_[1.0, (held_out - mu) / sd]
    p0 = (np.count_nonzero(train_y) + 0.5) / (len(train_y) + 1)
    coef = np.zeros(a.shape[1])
    coef[0] = math.log(p0 / (1 - p0))
    ridge = np.diag([0.0] + [penalty] * train_x.shape[1])
    for _ in range(60):
        prob = 1.0 / (1.0 + np.exp(-np.clip(a @ coef, -30.0, 30.0)))
        grad = a.T @ (prob - train_y) + ridge @ coef
        weights = np.maximum(prob * (1 - prob), 1e-8)
        hessian = (a.T * weights) @ a + ridge
        step = np.linalg.solve(hessian, grad)
        coef -= step
        if np.max(np.abs(step)) < 1e-9:
            break
    return float(np.clip(1.0 / (1.0 + np.exp(-np.clip(new_a @ coef, -30, 30))),
                         1e-8, 1 - 1e-8))


def heldout_logloss(x: np.ndarray, y: np.ndarray) -> dict:
    """LOO and training-prior baseline, with the *same* folds and fixed penalty."""
    losses, baseline_losses = [], []
    for i in range(len(y)):
        keep = np.arange(len(y)) != i
        p = fitted_logistic_probability(x[keep], y[keep], x[i])
        prior = (np.count_nonzero(y[keep]) + 0.5) / (len(y[keep]) + 1)
        yi = float(y[i])
        losses.append(-yi * np.log(p) - (1 - yi) * np.log(1 - p))
        baseline_losses.append(-yi * np.log(prior) - (1 - yi) * np.log(1 - prior))
    return {
        "n_heldout": len(y),
        "logloss": float(np.mean(losses)),
        "baseline_logloss": float(np.mean(baseline_losses)),
        "logloss_gain_vs_prior": float(np.mean(baseline_losses) - np.mean(losses)),
        "note": "LOO descriptive only; overlapping folds prevent naive independent-p inference.",
    }


def run(x: np.ndarray, y: np.ndarray, permutations: int = DEFAULT_PERMUTATIONS,
        seed: int = DEFAULT_SEED) -> dict:
    if permutations <= 0:
        raise ValueError("permutations must be positive")
    rng = np.random.default_rng(seed)
    perm = np.asarray([rng.permutation(y) for _ in range(permutations)])
    totals = x.sum(axis=1)
    props = x / totals[:, None]
    rank = np.sort(props, axis=1)[:, ::-1]
    # Preplanned contrast. E values are NOT information quantities.
    results = {
        "labeled_proportions": permutation_test(props, y, perm),
        "identity_free_rank_proportions": permutation_test(rank, y, perm),
        "raw_biomass_sensitivity": permutation_test(x, y, perm),
        "total_biomass_sensitivity": permutation_test(totals[:, None], y, perm),
        "rank1_sensitivity": permutation_test(rank[:, :1], y, perm),
    }
    results["delta_E_labeled_minus_rank"] = (
        results["labeled_proportions"]["E"]
        - results["identity_free_rank_proportions"]["E"]
    )
    # Scalar evenness: entropy of proportions, not a community Shannon diversity model.
    ent = -(props * np.log(np.maximum(props, 1e-300))).sum(axis=1)
    results["shannon_sensitivity"] = permutation_test(ent[:, None], y, perm)
    results["LOO_labeled"] = heldout_logloss(props, y)
    results["LOO_rank"] = heldout_logloss(rank, y)
    results["LOO_labeled_minus_rank_logloss_gain"] = (
        results["LOO_labeled"]["logloss_gain_vs_prior"]
        - results["LOO_rank"]["logloss_gain_vs_prior"]
    )
    return {
        "design": {"permutations": permutations, "seed": seed,
                   "unit": "whole mesocosm", "history_contrast": "two order regimes",
                   "dissimilarity": "Bray-Curtis", "history_null": "unblocked labels",
                   "identity_free_map": "sort each pot's 15 proportions descending",
                   "LOO_model": "fixed L2 logistic, penalty 4, fold-local scaling"},
        "outcomes": results,
        "interpretation_boundary": [
            "E and delta_E are null-calibrated R2 diagnostics, NOT mutual information.",
            "Rank sorting cannot increase true information about history.",
            "History regimes confound arrival timing/time-for-growth with interactions.",
            "Any zeros/missingness and assignment constraints must be source-audited.",
            "This is an independent post-submission candidate, not KNK v39.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=DEFAULT_PERMUTATIONS)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    args = parser.parse_args()
    x, y, audit = load_panel(args.assembly_csv)
    result = run(x, y, permutations=args.permutations, seed=args.seed)
    out = {
        "status": "RESULT_GENERATED_UNVERIFIED",
        "data_doi": "10.5061/dryad.sbcc2frh4",
        "source_audit": audit,
        **result,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n")
    print(f"Wrote {args.output}; n={audit['n_pots']}, species={audit['n_species']}")
    print(json.dumps({k: result["outcomes"][k] for k in
                      ("labeled_proportions", "identity_free_rank_proportions",
                       "delta_E_labeled_minus_rank")}, indent=2))


if __name__ == "__main__":
    main()
