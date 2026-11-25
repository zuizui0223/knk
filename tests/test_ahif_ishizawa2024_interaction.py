import itertools
import numpy as np
from benchmarks.score_ahif_ishizawa2024_interaction import predict,profile,score_errors,summarize

def test_low_order_predict_no_target_data():
    species=["a","b","c","d"]
    base=dict(zip(species,[5.,6.,7.,8.]))
    m={(frozenset(s),i):base[i] for k in (1,2,3) for s in itertools.combinations(species,k) for i in s}
    for k in (2,3):
        v,c=predict({q:z for q,z in m.items() if len(q[0])<=k},frozenset(species),k)
        assert np.allclose(v,[5,6,7,8]) and c==1.

def test_tv_rank_and_label_decomposition():
    o=np.array([.7,.2,.1]);p=np.array([.1,.2,.7])
    e=score_errors(o,p,np.log10(o),np.log10(p))
    assert np.isclose(e["named_TV"],.6) and e["rank_TV"]==0
    assert np.isclose(e["identity_assignment_TV"],.6)

def test_profile_absent():
    p=profile(np.array([2.,3.]),["a","c"],["a","b","c","d"])
    assert np.isclose(p.sum(),1) and p[1]==p[3]==0

def test_paired_improvement_identity():
    o=np.array([.7,.2,.1]);v=[np.array([.5,.3,.2]),np.array([.2,.1,.7]),np.array([.65,.25,.1])]
    e={k:score_errors(o,p,np.log10(o),np.log10(p)) for k,p in zip(("mono","pair","trio"),v)}
    z=summarize([{"errors":e}])
    d=z["pair_minus_trio_absolute_improvement"]
    assert np.isclose(d["named_TV"],d["rank_TV"]+d["identity_assignment_TV"])
