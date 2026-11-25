import numpy as np
from benchmarks.score_ahif_replay2025_state import features_and_outcomes,maps,score,restricted

def test_rank_average_is_not_rank_of_mean():
    x=np.array([[[7,3],[7,3],[3,7],[7,3],[3,7]]])
    f,y=features_and_outcomes(x)
    assert np.allclose(y["named"],[[.5,.5]])
    assert np.allclose(y["rank"],[[.7,.3]])

def test_holdout_maps():
    rng=np.random.default_rng(4);f=rng.normal(size=(30,4));g=np.repeat(np.arange(5),6)
    for h in maps(f,g):
        assert np.allclose(h[g[:,None]==g[None,:]],0)
        assert np.allclose(h.sum(1),1)

def test_predictive_gain_and_wild_reference():
    rng=np.random.default_rng(2);f=rng.normal(size=(60,4));g=np.repeat(np.arange(6),10)
    y=np.column_stack((3*f[:,1]+.1*rng.normal(size=60),f[:,2]))
    result=score(f,g,{"named":y,"rank":y},39)
    assert all(r["held_location_gain"]>.8 for r in result)
    assert all(r["p_approx"]==.025 for r in result)

def test_singleton_saturated_residual():
    rng=np.random.default_rng(3);f=rng.normal(size=(30,4));g=np.array([0]+[1]*14+[2]*15)
    y=rng.normal(size=(30,2));fit,res,sat=restricted(f,g,y)
    assert sat==1 and np.allclose(res[0],0)

def test_R_header_corner_and_strict_ids(tmp_path):
    from benchmarks.score_ahif_replay2025_state import read_matrix
    import pytest
    path=tmp_path/"matrix.tsv"
    path.write_text("ASV1\tASV2\n"+"".join(f"s{i}\t{i+1}\t2\n" for i in range(5)))
    x,audit=read_matrix(path,[f"s{i}" for i in range(5)])
    assert x.shape==(1,5,2) and audit["implicit_R_header_corner"]
    assert np.array_equal(x[0,:,0],np.arange(1,6))
    with pytest.raises(ValueError,match="orientation"):
        read_matrix(path,["missing","s1","s2","s3","s4"])
