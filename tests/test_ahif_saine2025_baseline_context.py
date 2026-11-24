import numpy as np
from benchmarks.score_ahif_saine2025_baseline_context import prepare_cache,cv,propensity,baseline_features

def synthetic():
    # Balanced exact contrasts with variable pre-treatment context.
    sites=np.repeat(np.arange(5),16)
    types=np.tile(np.repeat(["R","B"],8),5)
    f=np.tile(np.repeat([-1.,1.],4),10)
    x=np.column_stack((f,np.zeros((80,3))))
    t=np.tile([0.,1.],40);w=t-.5
    wt=np.ones(80)
    return sites,types,x,t,w,wt

def test_context_gain():
    s,t,x,lab,w,wt=synthetic()
    effect=w*x[:,0]*4
    y=np.tile(effect[:,None],(1,4))
    c=prepare_cache(x,t,s,y,wt)
    o=cv(c,w,2020)
    assert np.all(o["context_gain"]>.9)
    assert all(np.all(r["alpha"]==1) for r in o["folds"])

def test_constant_effect():
    s,t,x,lab,w,wt=synthetic()
    y=np.tile((w*4)[:,None],(1,4))
    o=cv(prepare_cache(x,t,s,y,wt),w,2020)
    assert np.all(o["pooled_gain"]>.99)
    assert np.allclose(o["context_gain"],0)

def test_holdout_no_leakage():
    s,t,x,lab,w,wt=synthetic()
    rng=np.random.default_rng(4)
    y=rng.normal(size=(80,4))
    yy=y.copy();yy[s==0]*=1000
    a=cv(prepare_cache(x,t,s,y,wt),w,2020)["folds"][0]
    b=cv(prepare_cache(x,t,s,yy,wt),w,2020)["folds"][0]
    assert np.array_equal(a["alpha"],b["alpha"])

def test_baseline_and_propensity():
    a=np.ones((4,3,5));a[:,1:]*=100
    f=baseline_features(a)
    assert f.shape==(4,4) and np.allclose(f[:,1],.2)
    assert np.allclose(propensity(np.array([0,1,0,1]),np.array([0,0,1,1])),.5)
