import numpy as np
from benchmarks.score_ahif_kansas_forecast import forecasts,evaluate,signscore,feature

def test_strict_time_holdout():
    rng=np.random.default_rng(7)
    x=rng.uniform(1,10,(15,4,2,8))
    z=evaluate(x,["a","b","c","d"])
    xx=x.copy();xx[5:]*=100
    # Current/future outcomes may alter loss and signs, not model fit.
    a=x[:,[0,1]].mean(1);b=xx[:,[0,1]].mean(1)
    p=forecasts(a[:,0],a[:,1]-a[:,0],np.ones(8,bool),5)
    q=forecasts(b[:,0],b[:,1]-b[:,0],np.ones(8,bool),5)
    assert all(np.allclose(p[k],q[k]) for k in p)
    assert len(z["all_forecast_cells"])==60 and len(z["exact_concurrent"])==15

def test_unseen_and_two_directions():
    c=np.ones((7,3));d=np.tile([2.,-3.,7.],(7,1))
    z=forecasts(c,d,np.array([1,1,0],bool),5)
    assert all(v[2]==0 for v in z.values())
    a=signscore(np.array([1.,-1.,0.]),np.array([4.,-2.,99.]))
    assert a["held_increase"]==4 and a["held_decrease"]==2 and a["held_min"]==2

def test_exact_partition_weights_and_constant_response():
    x=np.ones((15,4,2,3))*10
    x[:,:,1]+=np.array([2.,-3.,1.])
    z=evaluate(x,["a","b","c","d"])
    assert z["aggregate_forecast"]["models"]["previous"]["SSE"]==0
    assert z["aggregate_forecast"]["models"]["cumulative"]["SSE"]==0
    assert all(q["strong_descriptive_gate"] for q in z["exact_concurrent"])
    assert all(len(r["train_blocks"])==len(r["test_blocks"])==2 and not set(r["train_blocks"])&set(r["test_blocks"]) for r in z["all_forecast_cells"])
