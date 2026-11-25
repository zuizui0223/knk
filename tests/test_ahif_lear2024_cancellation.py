import numpy as np
from benchmarks.score_ahif_lear2024_cancellation import metrics,validation_gap,crossfit,intervals

def test_exact_cancel_and_common_sign():
    z=metrics([10,-8,0,0,0])
    assert z["gross"]==18 and z["net"]==2 and z["cancellation"]==16
    assert metrics([1,2,3,4,5])["cancellation"]==0

def test_crossfit_reproducible_opposite_signs():
    a=np.tile([10,10,10,10,10],(6,1))
    b=np.tile([16,6,10,10,10],(6,1))
    z=crossfit([(a,b)]*4,1,16)
    assert np.isclose(z["pooled_validated_gap"],8)
    assert np.isclose(z["mean_within_level_validated_gap"],8)
    assert validation_gap(np.array([1,-1]),np.array([-1,1]))==-2

def test_simultaneous_intervals():
    rng=np.random.default_rng(2);cells=[]
    for _ in range(4):
        a=rng.normal(30,1,(6,5))
        b=a+np.array([15,-15,0,0,0])
        cells.append((a,b))
    z=intervals(cells,5)
    assert z["opposing_signs_supported"] and z["positive"]==["p"] and z["negative"]==["o"]
