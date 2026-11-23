"""Synthetic-only controls for source GR role decomposition under species-out CV."""
import numpy as np

from benchmarks.decompose_vansteijn_neighbor_growth_rate import (
    BASELINE, MAIN, MODELS, feat, score, paired_bootstrap,
    permute_GR_within_guild, gr_mapping_placebo,
)
from benchmarks.predict_vansteijn_unseen_neighbor_traits import GUILDS


def _toy():
    rng=np.random.default_rng(75)
    names=[f"{guild}_{i}" for guild in GUILDS for i in range(5)]
    traits={}
    for n in names:
        v=rng.normal(size=7)
        traits[n]={"guild":n.split("_")[0],"traits":v}
    rows=[]
    for f in names:
        for n in names:
            if f==n:continue
            fm=traits[f]["traits"]
            nm=traits[n]["traits"]
            val=-.5+.30*fm[0]-.65*nm[6]+rng.normal(0,.02)
            rows.append({"focal":f,"neighbor":n,"Cdelta_total":float(val),
                         "RII_late_total":float(val/2),
                         "RII_sim_total":float(-0.1+0.2*fm[0])})
    return rows,traits


def test_holdout_neighbor_growth_rate_predicts_synthetic():
    rows,traits=_toy()
    r=score(rows,traits)
    assert r["n_species"]==15 and r["n_directed_pairs"]==210
    assert r["model_scores"][MAIN]["species_equal_MSE"] < \
           r["model_scores"][BASELINE]["species_equal_MSE"]
    b=paired_bootstrap(r,MAIN,n_boot=99)
    assert b["n_heldout_species_improved"]>=12
    assert b["mean_MSE_gain_positive_is_better"]>0


def test_role_specific_growth_rate_features_are_separate():
    rows,traits=_toy()
    x=feat(rows[:1],traits,sorted(traits)[:-1],"focal_morph6_plus_neighbor_GR")
    y=feat(rows[:1],traits,sorted(traits)[:-1],"focal_morph6_plus_focal_GR")
    assert x.shape==y.shape==(1,7)
    assert not np.isclose(x[0,-1],y[0,-1])


def test_heldout_species_focal_response_does_not_leak_into_neighbor_test():
    rows,traits=_toy()
    a=score(rows,traits,models=(BASELINE,MAIN))
    alt=[]
    for r in rows:
        z=r.copy()
        if z["focal"]=="Grass_0":z["Cdelta_total"]+=17
        alt.append(z)
    b=score(alt,traits,models=(BASELINE,MAIN))
    for model in (BASELINE,MAIN):
        assert np.isclose(
            a["per_heldout_neighbor"]["Grass_0"]["MSE"][model],
            b["per_heldout_neighbor"]["Grass_0"]["MSE"][model])


def test_shuffle_keeps_partner_guild_without_mutating_focal_traits():
    rows,traits=_toy()
    mapping=permute_GR_within_guild(traits,np.random.default_rng(1))
    assert set(mapping)==set(traits)
    assert all(traits[a]["guild"]==traits[b]["guild"] for a,b in mapping.items())
    p=gr_mapping_placebo(rows,traits,n_perm=7)
    assert 0<=p["fraction_shuffles_as_good_or_better"]<=1


def test_source_RII_endpoint_target_plumbing():
    rows,traits=_toy()
    priority=score(rows,traits,models=(BASELINE,MAIN),target="RII_late_total")
    simultaneous=score(rows,traits,models=(BASELINE,MAIN),target="RII_sim_total")
    assert priority["target"]=="RII_late_total"
    assert simultaneous["target"]=="RII_sim_total"
    assert priority["model_scores"][MAIN]["species_equal_MSE"] < \
           priority["model_scores"][BASELINE]["species_equal_MSE"]
    audit=gr_mapping_placebo(rows,traits,n_perm=7,target="RII_late_total")
    assert audit["target"]=="RII_late_total"
