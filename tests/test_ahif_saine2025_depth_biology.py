import numpy as np
from benchmarks.score_ahif_saine2025_depth_biology import make_cache,cv,evaluate,LinearDesign,treatment_features
from benchmarks.score_ahif_saine2025_baseline_context import propensity

def synthetic():
    rng=np.random.default_rng(1)
    sites=np.repeat(np.arange(5),40)
    types=np.tile(np.repeat(["R","B"],20),5)
    strata=np.array([f"{s}|{t}" for s,t in zip(sites,types)])
    t=np.tile([0.,1.],100);w=t-.5
    x=rng.normal(size=(200,4));wt=np.ones(200)
    # strong, truly reproducible biological modifier
    outcome=w*(1.0+x[:,3]+3*x[:,0])+.1*rng.normal(size=len(w))
    y=np.tile(outcome[:,None],(1,4))
    return x,types,sites,strata,t,w,y,wt

def test_gain_and_no_leakage():
    f,types,sites,s,t,w,y,wt=synthetic()
    c=make_cache(f,types,sites,wt,w)
    z=cv(c,y,wt,2020)
    assert np.all(z["gain"]>.5)
    yy=y.copy();yy[sites==0]*=1000
    assert np.array_equal(z["folds"][0]["alpha"],cv(c,yy,wt,2020)["folds"][0]["alpha"])

def test_wild_joint_family():
    f,types,sites,s,t,w,y,wt=synthetic()
    z=evaluate(f,types,sites,s,t,y,wt,19)
    assert len(z["nested_prediction"])==4
    assert len(z["direct_linear_diagnostic"])==4
    assert all(r["biology_gain_over_depth"]>.5 for r in z["nested_prediction"])
    assert all(r["p_approx"]==.05 for r in z["direct_linear_diagnostic"])

def test_depth_only_is_not_biology():
    f,types,sites,s,t,w,y,wt=synthetic()
    y=np.tile((w*(1+f[:,3]))[:,None],(1,4))
    z=cv(make_cache(f,types,sites,wt,w),y,wt,2021)
    assert np.all(z["gain"]<.1)

def test_weighted_HC3_dimensions():
    rng=np.random.default_rng(2)
    x=np.column_stack((np.ones(50),rng.normal(size=(50,3))))
    d=LinearDesign(x,np.ones(50))
    y=rng.normal(size=(50,4));z,b=d.wald3(y)
    assert z.shape==(4,) and b.shape==(3,4) and np.all(z>=0)
