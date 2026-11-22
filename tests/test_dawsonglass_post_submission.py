"""Synthetic-only tests for the post-submission external validation.

These tests validate score behavior; they contain NO experimental results.
"""
import csv
import numpy as np
import pytest

from benchmarks.run_dawsonglass_storage_post_submission import (
    AuditStop,
    bray_curtis,
    load_panel,
    run,
)


def _write_panel(tmp_path, *, missing=False, duplicate=False, negative=False):
    fn = tmp_path / "assembly_synthetic.csv"
    species = [f"Species_{i:02d}" for i in range(15)]
    with fn.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "Pot", "Species", "Treatment", "Biomass", "init_height", "Sp_code"
        ])
        writer.writeheader()
        for pot in range(24):
            for j, sp in enumerate(species):
                if missing and pot == 0 and j == 0:
                    continue
                value = 0.9 if j == (0 if pot < 12 else 1) else 0.1 / 14
                if negative and pot == 0 and j == 0:
                    value = -1
                row = dict(Pot=str(pot), Species=sp,
                           Treatment="Ambient Order" if pot < 12 else "Warm Order",
                           Biomass=value, init_height=1, Sp_code=sp)
                writer.writerow(row)
                if duplicate and pot == 0 and j == 0:
                    writer.writerow(row)
    return fn


def test_panel_and_provenance(tmp_path):
    x, y, meta = load_panel(_write_panel(tmp_path))
    assert x.shape == (24, 15)
    assert sum(y) == 12
    assert meta["rows"] == 360
    assert meta["missing_pot_species"] == 0
    assert len(meta["source_sha256"]) == 64


@pytest.mark.parametrize("flag", ["missing", "duplicate", "negative"])
def test_stop_ambiguous_panel(tmp_path, flag):
    with pytest.raises(AuditStop):
        load_panel(_write_panel(tmp_path, **{flag: True}))


def test_identical_abundance_multiset_erases_identity_signal():
    # Exactly the same rank distribution in every pot, but a different winner.
    n, s = 24, 15
    y = np.r_[np.zeros(12, dtype=np.int8), np.ones(12, dtype=np.int8)]
    x = np.full((n, s), 0.1 / 14)
    x[:12, 0] = 0.9
    x[12:, 1] = 0.9
    result = run(x, y, permutations=199, seed=5)["outcomes"]
    assert result["labeled_proportions"]["E"] > 0.5
    assert abs(result["identity_free_rank_proportions"]["E"]) < 1e-10


def test_abundance_architecture_signal_survives_sorting():
    n, s = 24, 15
    y = np.r_[np.zeros(12, dtype=np.int8), np.ones(12, dtype=np.int8)]
    x = np.empty((n, s))
    x[:12] = [0.90] + [0.10 / 14] * 14
    x[12:] = [0.40] + [0.60 / 14] * 14
    result = run(x, y, permutations=199, seed=9)["outcomes"]
    assert result["labeled_proportions"]["p_upper"] <= 0.01
    assert result["identity_free_rank_proportions"]["p_upper"] <= 0.01


def test_bray_curtis_symmetric_and_zero():
    x = np.asarray([[1.0, 0.0], [0.0, 1.0], [1.0, 0.0]])
    d = bray_curtis(x)
    assert np.allclose(d, d.T)
    assert np.allclose(np.diag(d), 0)
    assert d[0, 1] == 1.0
    assert d[0, 2] == 0.0
