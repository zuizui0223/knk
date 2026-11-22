#!/usr/bin/env python3
"""Eligibility audit for KNK cross-study structural history synthesis.

IMPORTANT: this intentionally does NOT pool non-comparable study effect sizes.
Three scored experiments are a comparative case series, not a meta-analysis.
Sources: docs/KNK_MACRO_META_SCOPING_20261123.md and study registry.
"""
from __future__ import annotations
import argparse
import csv
import json
import math
from collections import Counter
from pathlib import Path

DEFAULT_CANDIDATES = Path("database/macro_meta/candidate_studies_v1.csv")
DEFAULT_PILOT = Path("database/macro_meta/three_experiment_pilot_v1.csv")
MIN_CLUSTERS_PER_POOL = 5  # conservative gating rule, NOT sufficiency for meta-regression
REQUIRED = (
    "experiment_id", "publication_cluster", "operation_lane", "response_granularity",
    "randomization_contract", "independent_n", "E_named", "E_sorted",
    "delta_E_descriptive", "p_named", "p_sorted", "SE_delta_E",
    "analysis_status", "source_url",
)


def load_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        r = csv.DictReader(f)
        if not r.fieldnames or any(k not in r.fieldnames for k in REQUIRED):
            raise ValueError("Pilot CSV missing required columns")
        rows = list(r)
    if len(rows) == 0:
        raise ValueError("Empty pilot")
    ids = [r["experiment_id"].strip() for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError("Repeated experiment id")
    for r in rows:
        if not r["publication_cluster"] or not r["source_url"].startswith("http"):
            raise ValueError("Missing cluster or provenance URL")
        if int(r["independent_n"]) < 2:
            raise ValueError("Invalid independent unit count")
        for k in ("E_named", "E_sorted", "delta_E_descriptive", "p_named", "p_sorted"):
            v = float(r[k])
            if not math.isfinite(v):
                raise ValueError("Nonfinite: " + k)
            if k.startswith("p_") and not 0 <= v <= 1:
                raise ValueError("Invalid p")
        if abs(float(r["E_named"]) - float(r["E_sorted"]) -
               float(r["delta_E_descriptive"])) > 1e-6:
            raise ValueError("delta_E does not reconcile to named-sorted: " + r["experiment_id"])
        if r["SE_delta_E"] and (not math.isfinite(float(r["SE_delta_E"])) or
                              float(r["SE_delta_E"]) <= 0):
            raise ValueError("Invalid paired study SE")
    return rows


def evaluate(rows: list[dict[str, str]], min_clusters: int = MIN_CLUSTERS_PER_POOL) -> dict:
    clusters = set(r["publication_cluster"] for r in rows)
    candidates = [r for r in rows if r["analysis_status"] == "frozen"]
    by_group: dict[str, list[dict[str, str]]] = {}
    for r in candidates:
        key = f'{r["operation_lane"]}|{r["response_granularity"]}'
        by_group.setdefault(key, []).append(r)
    gates = {}
    for key, data in sorted(by_group.items()):
        eligible = [r for r in data if r["SE_delta_E"] and
                    "verify" not in r["randomization_contract"]]
        n = len(set(r["publication_cluster"] for r in eligible))
        reasons = []
        if n < min_clusters:
            reasons.append(f"only {n} independent source clusters with paired SE; need >= {min_clusters}")
        codes = {r["randomization_contract"] for r in eligible}
        if len(codes) > 1:
            reasons.append("design-specific variance mapping must be harmonized")
        gates[key] = {"scored_frozen_experiments": len(data),
                      "eligible_source_clusters": n,
                      "ready_to_pool": not reasons,
                      "reasons": reasons}
    if len(gates) == 0:
        gates["none"] = {"ready_to_pool": False, "reasons": ["no frozen scored study"]}
    draft_count = sum(r["analysis_status"] != "frozen" for r in rows)
    warnings = [
        "E is null-calibrated partial R2, not mutual information or fraction information stored.",
        "Designs include genotype strata, rhizobox replicate blocks, and pot treatment: raw E differences alone cannot be naively pooled.",
        "Taxonomic species and functional-group endpoint granularity are different synthesis lanes.",
        "A review database of 268 historical articles is not 268 independent scored effect sizes.",
        "Do not treat time points, pot species rows, and field/mesocosm from same publication as independent publications.",
        "This pipeline has NO statistical pooled estimate while gate is blocked.",
    ]
    return {
        "status": "ELIGIBILITY_GATE_ONLY",
        "scored_experiment_rows": len(rows),
        "distinct_publication_clusters": len(clusters),
        "frozen_scored_rows": len(candidates),
        "unmerged_or_unfrozen_rows": draft_count,
        "all_pooled_analyses_blocked": not any(v["ready_to_pool"] for v in gates.values()),
        "gates": gates, "claim_controls": warnings,
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pilot", type=Path, default=DEFAULT_PILOT)
    p.add_argument("--output", type=Path)
    args = p.parse_args()
    report = evaluate(load_rows(args.pilot))
    payload = json.dumps(report, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    print(payload)


if __name__ == "__main__":
    main()
