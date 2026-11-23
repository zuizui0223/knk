"""Source-free guards for partner-identity signal beyond neighbor functional guild."""
import numpy as np
from benchmarks.audit_vansteijn_partner_specificity import (
    centered_by_focal, conditional_species_identity,
    guild_mapping_placebo,
)


def fake_rows(identity_effect:bool):
    kinds=("Grass","Forb","Legume")
    names=[f"{g}_{i}" for g in kinds for i in range(5)]
    guild={n:n.split("_")[0] for n in names}
    out=[]
    rng=np.random.default_rng(22)
    for focal in names:
        for neighbor in names:
            if neighbor==focal:continue
            base=-(kinds.index(guild[focal]))*.2
            extra=(.40 if neighbor.endswith("_0") else -.08) if identity_effect else 0.
            out.append({
                "focal":focal,"neighbor":neighbor,
                "focal_guild":guild[focal],"neighbor_guild":guild[neighbor],
                "RII_late_total":base+extra+rng.normal(0,.035),
                "RII_sim_total":-.2 + rng.normal(0,.02),
            })
    return out


def test_species_id_signal_within_guild_survives():
    a=fake_rows(True)
    out=conditional_species_identity(a,"RII_late_total",n_perm=199,seed=1)
    assert out["focal_and_neighbor_guild_residual_identity_R2"]>0.7
    assert out["p_exchangeability_upper"]<.02


def test_no_species_id_signal_when_only_focal_effect():
    a=fake_rows(False)
    out=conditional_species_identity(a,"RII_late_total",n_perm=199,seed=1)
    assert out["focal_and_neighbor_guild_residual_identity_R2"]<.25


def test_balanced_random_guild_placebo():
    a=fake_rows(True)
    out=guild_mapping_placebo(a,"RII_late_total",n_perm=199)
    assert 0<=out["all_mapping_upper_fraction_MC"]<=1
    assert 0<=out["guild_R2_after_focal_centering"]<=1


def test_focal_centering_removes_focal_baselines():
    a=fake_rows(True)
    x=centered_by_focal(a,"RII_late_total")
    assert len(x)==210
    assert abs(x.mean())<1e-12
