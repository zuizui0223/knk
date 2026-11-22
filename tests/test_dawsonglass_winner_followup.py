"""Synthetic, source-free guards for the post hoc winner-identity follow-up."""
import numpy as np
import pytest

from benchmarks.run_dawsonglass_storage_post_submission import AuditStop
from benchmarks.run_dawsonglass_winner_identity_followup import (
    analyze, bh_adjust, exact_winner_test, within_arm_compensation,
)


NAMES = ["Rumex crispus", "Bromus inermis", "Dactylis glomerata"] + [
    "Filler_" + str(i) for i in range(12)
]
Y = np.r_[np.zeros(12, dtype=np.int8), np.ones(12, dtype=np.int8)]


def test_exact_winner_association_when_winner_changes():
    winners = np.r_[np.ones(12, dtype=int), np.zeros(12, dtype=int)]
    r = exact_winner_test(winners, Y)
    assert r["exact_conditional_p"] < 0.001
    assert r["sample_assignments_denominator"] == 2704156


def test_winner_swap_without_rank_shift():
    x = np.full((24, 15), 0.10 / 14)
    x[:12, 1] = 0.90
    x[12:, 0] = 0.90
    r = analyze(x, Y, NAMES, n_perm=199)
    d = r["paired_representation_contrast"]
    assert d["labeled_R2"] > 0.99
    assert abs(d["sorted_R2"]) < 1e-9
    assert d["permutation_p_upper_under_global_no_treatment_null"] <= 0.01
    assert r["winner_identity"]["exact_conditional_p"] < 0.001


def test_rank_change_without_winner_change():
    x = np.full((24, 15), 0.10 / 14)
    x[:12, 0] = 0.90
    x[12:] = 0.60 / 14
    x[12:, 0] = 0.40
    r = analyze(x, Y, NAMES, n_perm=199)
    assert r["paired_representation_contrast"]["sorted_R2"] > 0.99
    assert r["winner_identity"]["exact_conditional_p"] == 1.0


def test_tied_winners_stop():
    x = np.ones((24, 15))
    with pytest.raises(AuditStop):
        analyze(x, Y, NAMES, n_perm=20)


def test_bh_bounds_and_order():
    p = np.array([0.05, 0.001, 0.5, 0.015])
    q = bh_adjust(p)
    assert np.all(q >= p)
    assert np.all(q <= 1)
    assert np.all(np.argsort(p) == np.argsort(q))


def test_within_arm_compensation_detects_synthetic_opposite_gram_changes():
    # No observed data here. Each arm has a deliberately opposite absolute-mass trend.
    x = np.full((24, 15), 0.05)
    x[:, 0] = np.tile(np.arange(1.0, 13.0), 2)
    x[:, 1] = np.tile(np.arange(20.0, 8.0, -1.0), 2)
    x[:, 2] = 1.0
    idx = {name: j for j, name in enumerate(NAMES)}
    result = within_arm_compensation(x, Y, idx)
    assert result["Rumex_vs_two_grasses_r_within_treatment"] < -0.99
    assert result["two_sided_permutation_p"] <= 0.001
    assert result["bootstrap95"][1] < 0
