import numpy as np
from benchmarks.score_ahif_saine2025_buffering import tv_outcomes,weighted_contrasts,tv_score,rarefy,r2_stream,r2_arrays
from benchmarks.score_ahif_saine2025_succession import BlockPermutation

def test_gap_and_swap():
    x=np.array([[[7,3],[3,7],[7,3]],[[7,3],[9,1],[7,3]]])
    y=tv_outcomes(x)
    assert np.allclose(y[0],[.4,0,.4,0,0,0])
    assert np.allclose(y[1],[.2,.2,0,0,0,0])
    assert np.all(y[:,[2,5]]>=0)

def test_weighted_contrast():
    y=np.tile(np.arange(6),(8,1)).astype(float)
    t=np.tile([True,True,False,False],2)
    s=np.repeat(["a","b"],4)
    y[t]+=2
    assert np.allclose(weighted_contrasts(y,t,s,{"a":.6,"b":.4}),2)
    z=tv_score(y,t,s,{"a":.6,"b":.4},1,99,100)
    assert len(z)==6 and all(r["difference_treated_minus_control"]==2 for r in z)

def test_rarefaction():
    x=np.array([[[10,0,30],[20,20,0],[5,5,30]]])
    r=rarefy(x,15,np.random.default_rng(1))
    assert np.all(r.sum(2)==15)
    assert np.all(r<=x)

def test_kernel_equals_original():
    rng=np.random.default_rng(3)
    y=rng.normal(size=(12,30));s=np.repeat([0,1,2],4);t=np.tile([0,0,1,1],3)
    d=BlockPermutation(t,s,99,5)
    a=r2_stream(d,y,sparse=True);b=r2_stream(d,y,sparse=False)
    original=d.score(y)
    assert np.allclose(a,b)
    assert np.isclose(a[0],original["R2"])
    assert np.isclose(a[0]-np.median(a[1:]),original["E"])
