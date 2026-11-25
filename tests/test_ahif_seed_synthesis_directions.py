import numpy as np
from benchmarks.score_ahif_seed_synthesis_directions import signed_components
from benchmarks.score_ahif_seed_synthesis_opposition import crossfit

def test_constant_opposition_and_gross_identity():
    r=[{"block":("s",str(i)),"delta":np.array([10.,-8.,0.])} for i in range(4)]
    z=signed_components(r)
    assert z["mean_held_increase"]==10 and z["mean_held_decrease"]==8
    assert z["mean_held_min_signed_component"]==8 and z["strong_direction_descriptive_gate"]
    assert z["mean_held_signed_gross"]==crossfit(r)["signed_gross"]

def test_same_sign_not_opposition():
    r=[{"block":("s",str(i)),"delta":np.array([10.,8.,0.])} for i in range(4)]
    z=signed_components(r)
    assert z["mean_held_decrease"]==0 and not z["strong_direction_descriptive_gate"]

def test_reversed_signs_fail():
    r=[{"block":("s",str(i)),"delta":np.array([10.,-8.])*(1 if i<2 else -1)} for i in range(4)]
    z=signed_components(r)
    assert not z["strong_direction_descriptive_gate"] and z["mean_held_min_signed_component"]<=0
