import itertools
import numpy as np
from benchmarks.score_ahif_ishizawa2024_controls import (
    rhs,trio_control,additive_fit,additive_predict,inclusion_exclusion
)
from benchmarks.score_ahif_ishizawa2024_interaction import predict

def source(add=True):
    sp=list("abcdefg");a={i:5.+k*.1 for k,i in enumerate(sp)}
    def y(S,i):return a[i]+sum(.03*(sp.index(j)+1) for j in S-{i})*add
    m={(frozenset(S),i):y(frozenset(S),i) for k in (1,2,3) for S in itertools.combinations(sp,k) for i in S}
    return sp,m,y

def test_additive_and_inclusion_exclusion():
    sp,m,y=source();S=frozenset("abcde")
    expected=np.array([y(S,i) for i in sorted(S)])
    assert np.allclose(additive_predict(additive_fit(m,sp),S),expected)
    assert np.allclose(inclusion_exclusion(m,S),expected)

def test_no_interaction_controls():
    sp,m,y=source(False);S=frozenset("abcd")
    for mode in ("pair_on_trio_rhs","global_trio_coeff"):
        assert np.allclose(trio_control(m,S,sp,mode)[0],rhs(m,S,1))

def test_full_pool_global_equals_target():
    sp,m,y=source();S=frozenset(sp)
    assert np.allclose(trio_control(m,S,sp,"global_trio_coeff")[0],predict(m,S,3)[0])

def test_context_control_outside_target():
    sp,m,y=source();S=frozenset("abcd")
    g=trio_control(m,S,sp,"global_trio_coeff")[0]
    v=predict(m,S,3)[0]
    m[frozenset("abe"),"a"]+=1
    assert np.allclose(predict(m,S,3)[0],v)
    assert not np.allclose(trio_control(m,S,sp,"global_trio_coeff")[0],g)
