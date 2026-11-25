import numpy as np
from benchmarks.score_ahif_seed_synthesis_opposition import metrics,crossfit,summarize_strata,score

def test_opposition_arithmetic_and_held_signs():
    r=[{"block":("s",str(i)),"delta":np.array([10.,-8.])} for i in range(4)]
    m,w,se=summarize_strata(r)
    assert m["cancellation"]==16 and w["cancellation"]==16
    assert crossfit(r)["gap"]==16
    assert crossfit(r[:3])["status"]=="STOP_FEWER_THAN_FOUR_BLOCK_GROUPS"

def test_distinct_measurement_and_missing_arm():
    def st(block,mode,tr,v):
        return {"age":2,"dose":5 if tr=="Seeds" else 0,"tr":tr,"block":("s",block),
                "mode":mode,"x":np.array(v),"total":sum(v),"rich":2}
    a=[st("1","direct","Control",[10,10]),st("1","direct","Seeds",[20,2]),
       st("2","direct","Seeds",[8,12]),st("2","indirect","Control",[10,10])]
    z=score({"x":a},{"x":["a","b"]},{"x":0},{"x":0})[0]["comparisons"]
    assert z[0]["n_matched_strata"]==1 and z[0]["excluded_unmatched_states"]==1
    assert z[0]["pooled_species_metrics"]["cancellation"]==16
    assert z[1]["status"]=="STOP_NO_MATCHED_ARM_STRATUM"
    assert score({"x":a},{"x":["a","b"]},{"x":1},{"x":0})[0]["status"]=="STOP_AMBIGUOUS_PHYSICAL_KEY_AGE"
