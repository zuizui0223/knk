import numpy as np
from benchmarks.score_ahif_narwani_forecast import mapping,CYAN,DIAT,GREEN,OTHER,persistence,WEEKS,evaluate
def fixture():
    n=sorted(CYAN|DIAT|GREEN|OTHER)+["unknown"+str(i) for i in range(15)]
    g,h=mapping(n);return n,g,h
def test_partition():
    n,g,h=fixture();assert len(n)==66;assert tuple(g.sum(0))==(8,43,15);assert np.allclose(g.sum(1),1);assert np.allclose(h[:,1],g[:,1:].sum(1))
def test_opposition_and_unknown_gate():
    n,g,h=fixture();a=np.ones((15,4,66));b=a.copy();b[:,:,g[:,0]==1]+=1;b[:,:,g[:,1]==1]-=.5
    r=persistence(a,b,g);assert len(r)==15 and all(v["opposite"] and v["robust_unresolved"] and v["held_half_opposition_fraction"]==1 for v in r)
def test_forecast_no_future_or_held_leak():
    n,g,h=fixture();rng=np.random.default_rng(911);a=rng.uniform(.01,2,(15,4,66));b=rng.uniform(.01,2,(15,4,66))
    x={"C":a,"MD":b};y={k:v[:,:,:6].copy() for k,v in x.items()};arms={"C":["c"+str(i) for i in range(4)],"MD":["m"+str(i) for i in range(4)]}
    _,v=evaluate(x,y,arms,g,h,"C","MD")
    xx={k:v.copy() for k,v in x.items()};xx["C"][5:]*=9;xx["MD"][5:]*=4;xx["C"][:,2:]*=2;xx["MD"][:,2:]*=3
    _,w=evaluate(xx,y,arms,g,h,"C","MD")
    p=v["representations"]["conservative3"][0];q=w["representations"]["conservative3"][0]
    assert p["week"]==q["week"]==WEEKS[5] and p["predictions"]==q["predictions"] and p["observed"]!=q["observed"]
