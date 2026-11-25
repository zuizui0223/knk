import numpy as np
from scipy.sparse import csr_matrix
from benchmarks.score_ahif_saine2025_transfer import cv,group_response

def test_identity_transfer_signal():
    x=np.tile(np.array([[-1.,1.],[1.,-1.]])[None,None,:,:],(5,2,1,1))
    r=cv(np.stack([x,x]))
    assert np.isclose(r["identity_gain_over_pooled"],1)
    assert all(f["alpha"]==1 for f in r["folds"])

def test_constant_species_is_pooled():
    x=np.ones((5,2,2,2))
    r=cv(np.stack([x,x]))
    assert r["identity_gain_over_pooled"]==0
    assert r["pooled_gain_over_zero"]==1
    assert all(f["alpha"]==0 for f in r["folds"])

def test_heldout_not_used_to_choose_alpha():
    rng=np.random.default_rng(2)
    a=rng.normal(size=(2,5,2,3,4))
    b=a.copy();b[:,0]*=1000
    r=cv(a)["folds"][0];s=cv(b)["folds"][0]
    assert r["alpha"]==s["alpha"]

def test_group_means_and_control_subtraction():
    w=csr_matrix(np.eye(4))
    change=csr_matrix(np.arange(8).reshape(4,2))
    controls=np.ones((2,1,2,2,2))
    # shape is (year, site, context, species, feature)
    r=group_response(w,[change,change],controls)
    assert np.allclose(r[0].ravel(),np.arange(8)-1)
