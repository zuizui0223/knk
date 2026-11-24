"""Synthetic-only whole-vial source and representation guards for Zou 2026."""
from benchmarks.score_ahif_zou2026_two_species import (
    GEN,DAYS,DENSITY,group_contrasts,coefficient_context,pearson
)
import math

def make_vials():
    out=[]
    for g in GEN:
        for day in DAYS:
            for d in DENSITY:
                for rep in range(1,4):
                    # Opposite arrival switches winner without changing
                    # maximum adult share (identity-free rank concentration).
                    t=90 if day<0 else 10 if day>0 else 50
                    f=100-t
                    out.append({"gen":g,"day":day,"startT":d,"startF":d,
                                "rep":rep,"adultT":t,"adultF":f,"vial":f"{g}|{day}|{d}|{rep}"})
    return out

def test_same_unlabeled_rank_different_named_winner():
    x=make_vials()
    means,contrasts=group_contrasts(x,boot=300,seed=20261124)
    assert len(means)==21 and len(contrasts)==9
    assert all(r["n_exact_permutations"]==400 for r in contrasts)
    assert all(math.isclose(r["delta_T_identity_share_Tfirst_minus_Ffirst"],0.8,abs_tol=1e-12)
               for r in contrasts)
    assert all(math.isclose(r["delta_unlabeled_top1_share_Tfirst_minus_Ffirst"],0,abs_tol=1e-12)
               for r in contrasts)
    assert all(r["exact_400_within_density_permutation_p_two_sided_top1"]==1
               for r in contrasts)

def test_conditional_exact_label_permutation_detects_identity_shift():
    _,cs=group_contrasts(make_vials(),boot=10)
    assert all(r["exact_400_within_density_permutation_p_two_sided_T_share"]<=.01
               for r in cs)

def test_coefficient_summary_is_not_new_independent_validation():
    x=make_vials()
    means,_=group_contrasts(x,boot=10)
    params={}
    for g in GEN:
        for day in DAYS:
            params[g,day,"alpha12"]={"median":.01}
            params[g,day,"alpha21"]={"median":.01*math.exp(-day/12)}
    result=coefficient_context(means,params)
    assert result["n_cells"]==21
    assert len(result["condition_level"])==21
    assert "SAME" in result["no_independent_prediction"]
    assert 0<=result["source_fitted_coefficient_asymmetry_vs_mean_T_share_pearson_in_sample"]<=1

def test_correlation_guard():
    assert math.isclose(pearson([1,2,3],[3,6,9]),1)
