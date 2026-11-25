import numpy as np
from benchmarks.score_ahif_cross_source_concentration import support,feats,splits,weights,panel,fit
def test_support_and_lag():
    y,m,e=support(np.array([3.,1.,-5.]));assert y==.75 and m==2 and np.isclose(e,1.6)
    assert support(np.array([-1.,0.]))==(None,0,None)
    a=np.ones((4,3));b=a.copy();b[0]=[4,2,1]
    x=feats(a,b,1,1,"control");a[1:]*=99;b[1:]*=55
    assert x==feats(a,b,1,1,"control")
def test_splits_and_source_weights():
    assert len(list(splits(4,True)))==6 and len(list(splits(4,False)))==36
    assert len(list(splits(5,True)))==20 and len(list(splits(5,False)))==200
    r=[dict(source=s,panel="p",time=1,role="g",direction=1) for s,n in (("a",3),("b",1)) for _ in range(n)]
    w=weights(r);assert np.isclose(w[:3].sum(),.5) and np.isclose(w[3],.5)
def test_panel_held_covariates_and_missing():
    a=np.ones((3,4,3));b=a.copy();b[:,:,0]+=2
    G=np.eye(3)
    r=panel("toy","p",[0,1,2],a,b,list("abc"),G,["g"],True,b.mean(1)-a.mean(1))
    assert len(r)==24 and sum(v["target_top1"] is None for v in r)==12
    assert all(v["target_top1"]==1 for v in r if v["direction"]==1)
    x=[v for v in r if v["target_top1"] is not None]
    p,c=fit(x,x,(),"control");assert np.allclose(p,1)
