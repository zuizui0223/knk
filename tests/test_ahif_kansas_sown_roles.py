import numpy as np
from benchmarks.score_ahif_kansas_sown_roles import roles,evaluate

def test_assignment_bounds():
    s=np.array([1,0,0,0],bool);o=np.array([0,1,0,0],bool);u=~(s|o);p=s.copy()
    z=roles(np.array([10.,-8.,3.,-2.]),s,o,u,p)
    assert z["sown_assignment_range"]==[8.,13.] and z["other_assignment_range"]==[-10.,-5.]
    assert z["robust_to_unresolved"] and z["total_species_mass"]==3.

def test_fixed_membership_not_learned_signs():
    s=np.array([1,0,0],bool);o=np.array([0,1,0],bool);u=~(s|o);p=np.zeros(3,bool)
    z=evaluate(np.tile([5.,-3.,1.],(4,1)),s,o,u,p,s)
    assert z["held_half_summary"]["n_oriented_held_halves"]==6
    assert z["held_half_summary"]["strong_descriptive_gate"]
    assert z["held_half_summary"]["fraction_robust_unresolved"]==1
