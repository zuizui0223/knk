import numpy as np
import pytest
from benchmarks.score_ahif_nontarget_temporal import (
    holm_values, pair_differences, exact_joint_swap, changes, roster_data, analyse
)
from benchmarks.score_ahif_vansteijn2026_field import SourceStop, SITES, TREATMENTS

def test_holm():
    assert np.allclose(holm_values([.01,.04,.03]),[.03,.06,.06])
def test_exact_conditional_signs():
    a=np.tile([1.,0.],(5,1)); b=np.tile([0.,2.],(5,1))
    z=exact_joint_swap([a,b])
    assert z["n_exact_assignments"]==1024
    assert z["p_exact"]==4/1024  # each site's sign reversal preserves squared norm
    assert exact_joint_swap([a*0,b*0])["p_exact"]==1.
def synthetic():
    names=[f"s{i}" for i in range(76)]
    target=names[:18]+["never1","never2"]
    labels=np.tile(TREATMENTS,5)
    panel={}
    for s in SITES:
        for j,y in enumerate(["2021","2022","2023","2024"]):
            x=np.ones((20,76))
            for i in range(20):
                x[i,18:20]+=(i%4)*j
            panel[s,y]={"X":x,"names":names,"sown_target":target,
                        "plots":[(s,k) for k in range(20)],"history":labels.copy()}
    return panel
def test_design_and_signed_deltas():
    panel=synthetic()
    out=analyse(panel,n_perm=19)
    assert len(out["primary"])==2
    assert len(out["secondary"])==16
    assert len(out["order_exchange"])==3
    x,_=changes(panel,"MAL")
    assert np.allclose(x["nontarget_rank"].sum(axis=1),0)
    assert np.any(x["nontarget_rank"]<0)
    assert pair_differences(x["nontarget_rank"],panel["MAL","2021"]["history"]).shape==(5,58)
def test_stop_on_zero_roster_and_plot_mismatch():
    p=synthetic()
    p["MAL","2024"]["X"][0,18:]=0
    with pytest.raises(SourceStop): roster_data(p["MAL","2024"])
    p=synthetic(); p["MAL","2024"]["plots"].reverse()
    with pytest.raises(SourceStop): changes(p,"MAL")
