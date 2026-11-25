import numpy as np
from benchmarks.score_ahif_kansas_group_forecast import model_predictions,direct,evaluate,summarize,MODELS
from benchmarks.score_ahif_kansas_forecast import forecasts

def fixture():
    rng=np.random.default_rng(22)
    x=rng.uniform(1,20,(15,4,2,6));G=np.eye(6)
    broad=np.column_stack((G[:,:4].sum(1),G[:,4],G[:,5]))
    groups={"broad3":(["s","o","u"],broad),"exact23":(["s","o","u"],broad),
            "sown_guilds6":(list("abcdef"),G)}
    return x,groups

def test_direct_no_current_future_leakage():
    x,groups=fixture();a=x[:,[0,1]].mean(1);G=groups["broad3"][1]
    z=direct(a,5,G,True,10)
    b=a.copy();b[5:]*=100
    assert np.allclose(z,direct(b,5,G,True,10))

def test_same_target_and_baselines():
    x,groups=fixture();z,_=evaluate(x,[("a","1"),("a","2"),("a","3"),("a","4")],groups)
    for name in z:
        assert z[name]["aggregate"]["n_cells"]==60
    for p,q in zip(z["broad3"]["cells"],z["guilds_to_broad3"]["cells"]):
        assert np.allclose(p["observed"],q["observed"])
        assert np.allclose(p["predictions"]["bottom_history"],q["predictions"]["bottom_history"])

def test_opposition_confusion():
    obs=[[1,-1,0],[-1,1,0],[2,-1,0],[-2,1,0]]
    pred=[[1,-1,0],[1,-1,0],[0,0,0],[0,0,0]]
    rows=[{"observed":o,"predictions":{m:p for m in MODELS}} for o,p in zip(obs,pred)]
    z=summarize(rows,["s","o","u"],True)["models"]["previous"]["opposition_forecast"]
    assert (z["tp"],z["fp"],z["tn"],z["fn"])==(1,1,1,1)
    assert z["balanced_accuracy"]==.5
