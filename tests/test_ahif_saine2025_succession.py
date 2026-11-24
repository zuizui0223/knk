import numpy as np
from benchmarks.score_ahif_saine2025_succession import BlockPermutation,holm

def test_block_partial_r2():
    labs=np.tile([0,0,1,1],3)
    strata=np.repeat(np.arange(3),4)
    d=BlockPermutation(labs,strata,99,1)
    x=3*labs+10*strata
    r=d.score(x)
    assert np.isclose(r["R2"],1)
    assert r["p"]<=.05
    assert d.score(strata)["status"]=="CONSTANT_RESPONSE"

def test_rank_erasure():
    labs=np.tile([0,0,1,1],3)
    strata=np.repeat(np.arange(3),4)
    x=np.asarray([[.7,.3] if l==0 else [.3,.7] for l in labs])
    d=BlockPermutation(labs,strata,99,2)
    assert np.isclose(d.score(x)["R2"],1)
    assert d.score(np.sort(x,axis=1))["E"]==0

def test_holm():
    assert np.allclose(holm([.01,.04,.03]),[.03,.06,.06])
