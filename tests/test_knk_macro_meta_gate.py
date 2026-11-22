"""Synthetic + registry tests for conservative KNK macro/meta eligibility gate."""
import csv
from pathlib import Path

import pytest

from benchmarks.knk_macro_meta_gate import REQUIRED, evaluate, load_rows


def fixture_row(key="s1", cluster="paper1", level="taxon_species", frozen=True,
                se="0.08"):
    return {
        "experiment_id": key, "publication_cluster": cluster,
        "operation_lane": "direct_arrival", "response_granularity": level,
        "randomization_contract": "verified_blocked", "independent_n": "20",
        "E_named": "0.4", "E_sorted": "0.1", "delta_E_descriptive": "0.3",
        "p_named": "0.001", "p_sorted": "0.25", "SE_delta_E": se,
        "analysis_status": "frozen" if frozen else "draft_unmerged",
        "source_url": "https://doi.org/10.0000/test",
    }


def test_real_pilot_is_not_a_meta_analysis():
    path = Path("database/macro_meta/three_experiment_pilot_v1.csv")
    rows = load_rows(path)
    x = evaluate(rows)
    assert x["scored_experiment_rows"] == 3
    assert x["distinct_publication_clusters"] == 3
    assert x["frozen_scored_rows"] == 2
    assert x["unmerged_or_unfrozen_rows"] == 1
    assert x["all_pooled_analyses_blocked"]
    assert not any(v["ready_to_pool"] for v in x["gates"].values())


def test_five_clusters_same_lane_and_valid_se_can_pass_gate():
    rows = [fixture_row(key=str(i), cluster=f"paper{i}") for i in range(5)]
    x = evaluate(rows)
    assert not x["all_pooled_analyses_blocked"]
    assert x["gates"]["direct_arrival|taxon_species"]["eligible_source_clusters"] == 5


def test_five_rows_same_paper_do_not_pass():
    rows = [fixture_row(key=str(i), cluster="paper1") for i in range(5)]
    assert evaluate(rows)["all_pooled_analyses_blocked"]


def test_unfrozen_and_missing_se_never_pass():
    rows = [fixture_row(key=str(i), cluster=f"p{i}", frozen=i<4, se="" if i==0 else "0.1") for i in range(12)]
    assert evaluate(rows)["all_pooled_analyses_blocked"]


def test_granularity_mixing_not_allowed():
    rows = [fixture_row(key=str(i), cluster=f"p{i}", level="taxon_species" if i < 3 else "functional_group") for i in range(6)]
    assert evaluate(rows)["all_pooled_analyses_blocked"]


def test_rejects_duplicate_key_and_inconsistent_delta(tmp_path):
    p = tmp_path / "x.csv"
    with p.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=REQUIRED)
        w.writeheader()
        w.writerow(fixture_row("a", "x"))
        w.writerow(fixture_row("a", "y"))
    with pytest.raises(ValueError, match="Repeated experiment"):
        load_rows(p)
    with p.open("w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=REQUIRED)
        w.writeheader()
        r=fixture_row("a","x")
        r["delta_E_descriptive"] = "0.2"
        w.writerow(r)
    with pytest.raises(ValueError, match="delta_E"):
        load_rows(p)
