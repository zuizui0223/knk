"""No raw data needed: synthetic controls for the post hoc selection audit."""
import numpy as np
import pytest

from benchmarks.audit_dawsonglass_competition_caveats import (
    TOP3, audit, centered_by_arm,
)
from benchmarks.run_dawsonglass_storage_post_submission import AuditStop


def _fake_panel():
    rng = np.random.default_rng(55)
    names = sorted(list(TOP3) + [f"S{i:02}" for i in range(12)])
    y = np.r_[np.zeros(12, dtype=np.int8), np.ones(12, dtype=np.int8)]
    x = rng.lognormal(mean=0, sigma=0.8, size=(24, 15))
    ix = {n: i for i, n in enumerate(names)}
    x[12:, ix["Rumex crispus"]] *= 2.3
    x[12:, ix["Bromus inermis"]] *= 0.6
    x[12:, ix["Dactylis glomerata"]] *= 0.6
    return x, y, names


def test_center_by_arm_preserves_zero_means():
    x, y, names = _fake_panel()
    for vals in (x[:, 0], x[:, :3]):
        z = centered_by_arm(vals, y)
        assert np.allclose(z[y == 0].mean(axis=0), 0)
        assert np.allclose(z[y == 1].mean(axis=0), 0)


def test_pair_selection_never_improves_corrected_p():
    x, y, names = _fake_panel()
    out = audit(x, y, names, n_perm=99, retention_perm=49)
    assoc = out["rumex_grass_association"]
    assert assoc["n_pair_candidates"] == 91
    assert assoc["one_sided_maxT_familywise_p_over_91_pairs"] >= \
        assoc["one_sided_p_selected_grasses"]
    budget = out["biomass_budget_identity"]
    assert np.isclose(budget["cov_R_G"],
                      budget["cov_R_Total"] - budget["Var_R"] - budget["cov_R_Other12"])
    assert 0 <= out["mean_composition_shift"]["fraction_of_absolute_species_shift_in_top3"] <= 1
    assert out["post_hoc_retention_sensitivities"]["excluding_three_largest_shifters"]["remaining_taxa"] == 12


def test_source_species_gate():
    x, y, names = _fake_panel()
    changed = list(names)
    changed[changed.index("Rumex crispus")] = "incorrect name"
    with pytest.raises(AuditStop):
        audit(x, y, changed, n_perm=9, retention_perm=9)
