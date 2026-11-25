import numpy as np
import pytest
from benchmarks.score_ahif_routing_shape import (
    fitted_curve,shape_residual,block_contrasts,fixed_site_wild,welch_scalar,analyse
)
from benchmarks.score_ahif_vansteijn2026_field import SourceStop, SITES, TREATMENTS

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


def test_curves_preserve_richness_and_concentration():
    for p in (np.array([.6,.3,.1,0.]),np.array([.5,.5,0.,0.]),np.array([1.,0.,0.,0.])):
        for kind in ("geometric","one_dominant"):
            z=fitted_curve(p,kind)
            assert np.isclose(z.sum(),1.)
            assert np.isclose(z@z,p@p)
            assert np.count_nonzero(z)==np.count_nonzero(p)
            assert np.all(z[:-1]>=z[1:])
def test_geometric_exact_shape():
    x=np.array([4,2,1,0.])/7
    assert np.allclose(shape_residual(x[None,:],"geometric"),0)
def test_factorial_contrasts():
    y=np.tile(["FtS","StF","SFearly","SFlate"],5)
    z=np.tile([2.,0.,2.,0.],5)[:,None]
    assert np.allclose(block_contrasts(z,y),np.tile([2.,0.,0.],(5,1)))
def test_fixed_site_approximation():
    a=np.arange(5.)[:,None];b=a+10
    z=fixed_site_wild(a,b)
    assert z["difference_l2"]==10
    assert z["p_approx"]<=.05
    assert welch_scalar(a,b)["ci95_unadjusted"][1]<0
    assert fixed_site_wild(a*0,b*0)["status"]=="STOP_ZERO_VARIANCE"
def test_full_analysis_and_stop():
    z=analyse(synthetic(),n_perm=19)
    assert len(z["primary_geometric_residual"])==2
    assert len(z["secondary_one_dominant_residual"])==2
    assert len(z["approximate_fixed_site_contrasts"])==3
    assert len(z["approximate_welch_order_sensitivity"])==2
    with pytest.raises(SourceStop): fitted_curve(np.array([.2,.8]),"geometric")
