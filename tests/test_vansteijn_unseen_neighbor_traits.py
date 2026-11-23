"""Synthetic-only guards: strict leave-species-out partner morphology prediction."""
from collections import defaultdict
import numpy as np

from benchmarks.predict_vansteijn_unseen_neighbor_traits import (
    GUILDS, MORPH_TRAITS, BASE_MODELS, build_features, run_strict_loso,
    paired_species_comparisons, shuffled_partner_traits, mapping_placebo,
)


def _synthetic():
    rng=np.random.default_rng(17)
    names=[f"{guild}_{i}" for guild in GUILDS for i in range(5)]
    traits={}
    for sp in names:
        traits[sp]={"guild":sp.split("_")[0],"traits":rng.normal(0,1,6)}
    rows=[]
    for f in names:
        for n in names:
            if f==n:continue
            # True outcome depends on BOTH focal and preexisting neighbor traits.
            u=traits[f]["traits"]
            v=traits[n]["traits"]
            y=(-.5+.25*u[0]-.45*v[0]+.08*abs(u[1]-v[1])+
               rng.normal(0,.005))
            rows.append({"focal":f,"neighbor":n,"Cdelta_total":float(y),
                         "RII_late_total":float(y/2),
                         "RII_sim_total":float(-.2+y/3)})
    return rows,traits


def test_strict_species_out_all_210_directed_pairs_once():
    rows,traits=_synthetic()
    res=run_strict_loso(rows,traits,"Cdelta_total")
    assert res["n_directed_pairs"]==210
    assert res["metrics"]["intercept"]["n_heldout_species"]==15
    assert len(res["heldout_each_early_neighbor"])==15
    assert set(res["models"])==set(BASE_MODELS)
    assert all(n==14 for n in res["heldout_each_early_neighbor"].values())
    assert all(np.isfinite(v["overall_rmse"]) for v in res["metrics"].values())


def test_both_traits_predicts_new_neighbor_better_than_focal_only():
    rows,traits=_synthetic()
    res=run_strict_loso(rows,traits,"Cdelta_total")
    assert (res["metrics"]["both_morph"]["overall_rmse"] <
            res["metrics"]["focal_morph"]["overall_rmse"])
    assert (res["metrics"]["both_morph_abs_gap"]["overall_rmse"] <
            res["metrics"]["guild"]["overall_rmse"])


def test_heldout_focal_outcomes_never_train_heldout_neighbor_prediction():
    rows,traits=_synthetic()
    original=run_strict_loso(rows,traits,"Cdelta_total")
    # Modify only rows with species Grass_0 as FOCAL but not neighbor.
    # The strict Grass_0-as-neighbor fold must be unaffected.
    changed=[]
    for row in rows:
        r=dict(row)
        if r["focal"]=="Grass_0":
            r["Cdelta_total"]+=100
        changed.append(r)
    follow=run_strict_loso(changed,traits,"Cdelta_total",
                           models=("intercept","guild","focal_morph",
                                   "both_morph","both_morph_abs_gap",
                                   "guild_both_morph_abs_gap"))
    for model in BASE_MODELS:
        assert abs(original["by_species"]["Grass_0"]["mse"][model] -
                   follow["by_species"]["Grass_0"]["mse"][model]) < 1e-9


def test_neighbor_trait_randomization_keeps_guild():
    rows,traits=_synthetic()
    rng=np.random.default_rng(9)
    mapping=shuffled_partner_traits(traits,rng)
    assert set(mapping)==set(traits)
    assert all(traits[name]["guild"]==traits[altered]["guild"]
               for name,altered in mapping.items())
    res=mapping_placebo(rows,traits,"Cdelta_total",n_perm=9)
    assert 0<=res["mapping_fraction_at_least_as_predictive"]<=1


def test_paired_species_uncertainty():
    rows,traits=_synthetic()
    res=run_strict_loso(rows,traits,"Cdelta_total")
    tests=paired_species_comparisons(res,n_boot=100,seed=10)
    assert tests["both_morph_abs_gap_vs_focal_morph"]["n_species"]==15
    assert len(tests["both_morph_abs_gap_vs_guild"]["bootstrap95"])==2
