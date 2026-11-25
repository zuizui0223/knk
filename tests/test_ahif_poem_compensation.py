import numpy as np
from benchmarks.score_ahif_poem_compensation import metrics,sign_gap,role_contributions,combine,welch,crossfit,COHORTS,AGES,EARLY,GUILDS
from benchmarks.audit_poem_design import SOURCE_SOWN

def test_arithmetic():
    z=metrics([10,-8,0]);assert z["gross"]==18 and z["net"]==2 and z["cancellation"]==16
    assert sign_gap(np.array([10,-8]),np.array([8,-6]))==12

def fixture():
    names=list(SOURCE_SOWN)+["unsown"]
    x=np.full((5,15),10.)
    g={}
    for c in COHORTS:
        for a in AGES:
            g[c,a,"S"]=x.copy()
            for k in EARLY:
                q=x.copy();q[:,[n in GUILDS[k] for n in names]]+=2
                q[:,-1]-=4;g[c,a,k]=q
    return g,names

def test_shared_control_contribution():
    g,n=fixture();cs=role_contributions(g,n,2)
    assert len(cs)==8
    assert np.allclose(combine(cs),[28/3,0,-4,28/3-4])
    assert abs(sum(w for _,w in cs))<1e-12

def test_crossfit_reproduces_constant_differences():
    g,n=fixture();cf,role=crossfit(g,n,8)
    assert all(np.isclose(v,8.) for v in cf.values())
    assert all(np.isclose(v,8.) for v in role.values())

def test_welch_zero_variance():
    g,n=fixture()
    assert all(r["ci95"] is None for r in welch(role_contributions(g,n,2),4))
