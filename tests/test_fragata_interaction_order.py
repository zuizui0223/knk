"""Synthetic source-free guards for interaction-formation anchor (Fragata 2022)."""
import numpy as np
from benchmarks.audit_fragata_interaction_order import (
    ORDER, ORDER_LABELS, N_PERMUTATIONS_EXACT, exact_blocked_pairwise,
    bootstrap_blocked, holm, score
)


def fake_boxes():
    data=[]
    for t in ORDER:
        for block in (1,2):
            for box in range(5):
                value={"10Te10Tu":0.55,
                       "10Te10Tu48h":0.31,
                       "10Te48h10Tu":0.85}[t]
                te=int(round(value*200))
                tu=200-te
                data.append({
                    "treatment":t,"block":block,"box":box+1,
                    "Te":te,"Tu":tu,"fraction_Te":te/(te+tu),
                    "fraction_Te_leaf4":0.30+0.1*value,
                    "fraction_Te_leaf3":0.25,
                    "fraction_Te_new_leaf":0.5,
                    "fraction_Tu_new_leaf":0.5,
                    "spatial_leaf_overlap":0.6,
                })
    return data


def test_balanced_exact_10v10_is_63504_assignments():
    a=fake_boxes()
    r=exact_blocked_pairwise(a,"10Te48h10Tu","10Te10Tu48h","fraction_Te")
    assert r["n_boxes_each_group"]==10
    assert r["exact_assignments"]==252**2
    assert N_PERMUTATIONS_EXACT==252**2
    assert r["observed_mean_difference"]>0.5
    assert r["exact_two_sided_p"]<0.001


def test_holm_not_less_than_unadjusted():
    p=[0.003,0.01,0.7]
    q=holm(p)
    assert all(x>=y for x,y in zip(q,p))
    assert q[0]<=q[1]<=q[2]


def test_whole_box_scoring_and_bootstrap():
    data=fake_boxes()
    result=score(data,{"status":"SYNTHETIC_TEST_ONLY"})
    assert set(result["balanced_treatments"]) == {"Te_first","Tu_first","simultaneous"}
    assert all(v["n_independent_boxes"]==10 for v in result["balanced_treatments"].values())
    c=result["primary_10_10_arrival_contrast"]
    assert c["observed_mean_difference"]>0
    assert c["bootstrap"]["difference_percentile95"][0]>0
    assert len(result["posthoc_exact_tests"]["fraction_Te_leaf4"])==3
