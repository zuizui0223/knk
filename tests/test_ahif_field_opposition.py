import numpy as np
from benchmarks.score_ahif_field_opposition import rolemap,contrasts,crossfit,welch,support,metrics,SITES,YEARS

def fixture():
    f=[f"f{i}" for i in range(10)];s=[f"s{i}" for i in range(10)];n=f+s+["u"]
    groups={}
    for site in SITES:
        for y in YEARS:
            a=np.tile(np.ones(21)*10,(5,1))
            groups[site,y,"SFlate"]=a.copy()
            groups[site,y,"SFearly"]=a.copy()
            for h in ("FtS","StF"):
                b=a.copy()
                b[:,10:20] -=2 if h=="FtS" else 0
                b[:,:10] -=2 if h=="StF" else 0
                b[:,-1]+=7
                groups[site,y,h]=b
    return groups,n,f,s

def test_roles_and_shared_control():
    g,n,f,s=fixture()
    d,r=contrasts(g,n,f,s)
    assert np.allclose(r["MAL","2022","SFlate"],[0,-20,7,-13])
    assert len(d)==32 and len(r)==16

def test_whole_block_crossfit():
    g,n,f,s=fixture();d,r=contrasts(g,n,f,s);a,b=crossfit(d,r,8)
    assert all(np.isclose(x,14) for x in a.values())
    assert all(np.isclose(x,14) for x in b.values())

def test_opposing_intervals():
    rng=np.random.default_rng(5)
    x=np.tile([0,-20,7,-13],(2,5,1))+rng.normal(0,.1,(2,5,4))
    z=welch(x,4)
    assert support(z)["other_sown_negative_nonsown_positive"]
    assert not support(z)["strict_early_positive_other_or_nonsown_negative"]
    assert metrics([0,-20,7])["cancellation"]==14
