import numpy as np
from benchmarks.score_ahif_replay2025_persistence import responses,persistence,covariates,center_kernel,kernel_predict,evaluate

def test_response_rank_before_average():
    x=np.array([[[7,3],[7,3],[3,7],[7,3],[3,7]]])
    f,r=responses(x)
    assert np.allclose(r["rank"][1],[[.7,.3]])
    assert np.allclose(r["named"][1],[[.5,.5]])

def test_perfect_initial_persistence():
    rng=np.random.default_rng(3)
    x=rng.random((20,3));x/=x.sum(1,keepdims=True)
    f=rng.normal(size=(20,4));tr=np.arange(15);te=np.arange(15,20)
    a,b=covariates(f,tr)
    p,q,g=persistence(a,tr,te,x,x)
    assert np.isclose(g,1) and np.allclose(q,x[te])

def test_kernel_training_center():
    x=np.arange(12.).reshape(6,2);tr=np.arange(4);te=np.arange(4,6)
    K=center_kernel(x@x.T,tr)
    assert np.allclose(K[:,tr].mean(1),0)
    y=np.arange(6.)[:,None]
    z=kernel_predict(K,tr,te,y,1.)
    yy=y.copy();yy[te]*=1000
    assert np.allclose(z,kernel_predict(K,tr,te,yy,1.))

def test_folds_and_penalties():
    rng=np.random.default_rng(4)
    x=rng.integers(1,100,size=(12,5,5))
    f,r=responses(x);loc=np.repeat(np.arange(3),4)
    z=evaluate(f,r,loc)
    assert all(len(v["persistence_folds"])==3 and len(v["kernel_totals"])==3 for v in z.values())
    assert all(sum(a["n_parents"] for a in v["persistence_folds"])==12 for v in z.values())
