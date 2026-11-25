import numpy as np
from benchmarks.score_ahif_vansteijn_pot_field import SPECIES,FAST,contrasts,design_maps,predict,score,MODELS

def test_fixed_focal_date():
    x=np.zeros((2,5,4,10));x[:,:,0]=2;x[:,:,1]=99;x[:,:,2]=3;x[:,:,3]=7
    l=np.tile(np.arange(4),(2,5,1))
    y=contrasts(x,l)
    assert np.allclose(y[:,[s in FAST for s in SPECIES]],4)
    assert np.allclose(y[:,[s not in FAST for s in SPECIES]],-1)

def test_maps_no_field_focal_or_site_leakage():
    f=np.random.default_rng(1).normal(size=(10,7))
    w,_=design_maps(f,MODELS["history"])
    y=np.random.default_rng(2).normal(size=(2,10));a=predict(y,w)
    y[0]+=100;y[1,0]+=100
    b=predict(y,w)
    assert np.isclose(a[0,0],b[0,0])
    assert np.allclose(np.diag(w),0)

def test_history_signal():
    f=np.random.default_rng(3).normal(size=(10,7));f[:,5]=np.arange(10)
    y=np.tile(f[:,5],(2,1))
    maps={k:design_maps(f,c) for k,c in MODELS.items()}
    r=score(y,maps)
    assert r["history_gain"]>0
    assert len(r["cells"])==20
