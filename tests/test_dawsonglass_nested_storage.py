"""Synthetic-only guardrails for hierarchical post-hoc group audit."""
import numpy as np
import pytest

from benchmarks.audit_dawsonglass_nested_storage import (
    GRASSES, analyze, grass_mask, grouping_placebo, holm_adjust, within_group,
)
from benchmarks.run_dawsonglass_storage_post_submission import AuditStop


def make_panel():
    names = sorted(set(GRASSES) | {
        "Rumex crispus", "Solidago canadensis", "Solidago rugose",
        "Circium arvense", "Hypericum perferatum", "Hypericum punctatum",
        "Ranunculus acris", "Ranunculus pennsylvanicus",
        "Rudbeckia hirta", "Vernonia gigantea",
    })
    assert len(names) == 15
    y = np.r_[np.zeros(12, dtype=int), np.ones(12, dtype=int)]
    x = np.full((24, 15), 0.1)
    byname = {name: i for i, name in enumerate(names)}
    x[:12, byname["Bromus inermis"]] = 12
    x[:12, byname["Dactylis glomerata"]] = 12
    x[:12, byname["Solidago canadensis"]] = 4
    x[12:, byname["Bromus inermis"]] = 6
    x[12:, byname["Dactylis glomerata"]] = 6
    x[12:, byname["Rumex crispus"]] = 19
    return x, y, names


def test_source_guild_mapping_and_missing_names():
    x, y, names = make_panel()
    g = grass_mask(names)
    assert g.sum() == 5
    assert (~g).sum() == 10
    with pytest.raises(AuditStop):
        grass_mask([s for s in names if s != "Tridens flavus"] + ["Fake grass"])


def test_within_guild_sort_is_identity_free_and_normalized():
    x, y, names = make_panel()
    rel, ranks = within_group(x, ~grass_mask(names))
    assert rel.shape == ranks.shape == (24, 10)
    assert np.allclose(rel.sum(axis=1), 1)
    assert np.allclose(ranks.sum(axis=1), 1)
    assert np.all(np.diff(ranks, axis=1) <= 0)


def test_holm_adjustment_is_monotonic_and_conservative():
    p = {"a": 0.02, "b": 0.001, "c": 0.4}
    q = holm_adjust(p)
    assert all(q[k] >= p[k] for k in p)
    assert all(q[k] <= 1 for k in p)
    assert q["b"] < q["a"] <= q["c"]


def test_nested_analysis_and_all_groupings():
    x, y, names = make_panel()
    out = analyze(x, y, names, n_permutations=99, seed=11)
    assert out["results"]["functional2_share"]["E"] > 0
    assert out["results"]["forb10_sorted_internal"]["E"] > 0
    assert out["mapping_placebo"]["number_partitions"] == 3003
    assert out["mapping_placebo"]["conditional_partitions"] == 220
    assert 0 <= out["mapping_placebo"]["conditional_upper_fraction"] <= 1
    assert out["means"]["forb_share_warm"] > out["means"]["forb_share_ambient"]


def test_zero_guild_mass_refuses():
    x, y, names = make_panel()
    x[:, ~grass_mask(names)] = 0
    with pytest.raises(AuditStop):
        analyze(x, y, names, n_permutations=9, mapping_placebos=False)
