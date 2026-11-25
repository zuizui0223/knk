import numpy as np
from benchmarks.score_ahif_role_contributors import overlap,delete,evaluate
def test_overlap_and_removal():
    assert overlap([1,0],[0,1])==0
    assert overlap([1,2],[2,4])==1
    assert overlap([0,0],[2,1]) is None
    z=delete(np.array([10.,2.,-3.]),np.array([1.,3.,0.]),["a","b","c"],1)
    assert z["selected"]==["b"] and z["residual_net"]==7
    assert np.isclose(z["capture_fraction"],1/6)
def test_turnover_and_identity_and_no_future():
    G=np.array([[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,0,1]])
    names=list("abcde");train=np.tile([10.,0.,-10.,0.,0.],(15,1));held=train.copy()
    held[5:]=[0,10,0,-10,0]
    z,r=evaluate("toy","x",list(range(15)),names,G,held,[(0,train,held)])
    s=z["all_forecast_dates"];assert s["joint"]["early"]["1"]["n_retained_among_full_opposition"]==10
    assert s["groups"][0]["deletions"]["early"]["1"]["mean_support_capture"]==0
    assert s["groups"][0]["overlap"]["lag"]==0
    z,r=evaluate("toy","x",list(range(15)),names,G,train,[(0,train,train)])
    assert z["all_forecast_dates"]["joint"]["early"]["1"]["n_retained_among_full_opposition"]==0
    assert z["all_forecast_dates"]["groups"][0]["overlap"]["lag"]==1
    changed=train.copy();changed[5:]=[0,100,0,-100,0]
    _,s=evaluate("toy","x",list(range(15)),names,G,held,[(0,changed,held)])
    for rule in ("early","previous","expanding"):
        assert r[0]["groups"][0]["deletions"][rule]["1"]["selected"]==s[0]["groups"][0]["deletions"][rule]["1"]["selected"]
