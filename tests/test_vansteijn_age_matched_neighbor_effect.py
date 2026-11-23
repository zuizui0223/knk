"""Synthetic controls for age-matched solo-relative neighbor impact."""
import math
import numpy as np
from benchmarks.score_vansteijn_age_matched_neighbor_effect import transform, sign_changes


def toy(f="Grass_A",n="Grass_B",alone=10.,sim=12.,pre=5.):
    q={"focal":f,"neighbor":n,"focal_guild":"Grass","neighbor_guild":"Grass"}
    for organ in ("total","above","below"):
        q[f"biomass_alone_{organ}"]=alone
        q[f"biomass_sim_{organ}"]=sim
        q[f"biomass_pre_{organ}"]=pre
        q[f"RII_sim_{organ}"]=(sim-alone)/(sim+alone)
        q[f"RII_late_{organ}"]=(pre-sim)/(pre+sim)
    return q


def test_same_focal_age_net_neighbor_effect_vs_solo():
    item=transform([toy()])[0]
    assert math.isclose(item["C_sim_total"],2/22)
    assert math.isclose(item["C_pre_total"],-5/15)
    assert math.isclose(item["Cdelta_total"],-5/15-2/22)
    sign=sign_changes([item])
    assert sign["counts"]["positive_sim_to_negative_pre"]==1
    assert sign["counts"]["negative_sim_to_positive_pre"]==0


def test_no_change_under_identical_neighbor_planting():
    item=transform([toy(sim=7,pre=7)])[0]
    assert item["Cdelta_total"]==0
    assert item["C_sim_total"]<0
    assert item["C_pre_total"]<0


def test_effect_is_not_alpha_coefficient():
    # Given focal grew 10 weeks in all modes, the C contrast measures
    # net neighbor impact but does not identify resource preemption vs modification.
    # Check scale bounded and differs from raw mass ratio.
    item=transform([toy(alone=4,sim=2,pre=1)])[0]
    assert -2<item["Cdelta_total"]<0
    assert item["C_pre_total"]<item["C_sim_total"]
