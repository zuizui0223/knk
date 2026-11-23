"""Synthetic source-free gates for Delory legacy process vs function audit."""
import numpy as np

from benchmarks.score_delory_chemical_legacy import (
    SPECIES,TREATMENTS,PRIMARY,SECONDARY,
    holm,outcome_bootstrap,treatment_label_perm
)


def fake_groups():
    rng=np.random.default_rng(17)
    d={}
    for sp in SPECIES:
        for tr in TREATMENTS:
            group=[]
            for i in range(10):
                v={"RSD":3.0,"PDW":0.8,"Nuptake":0.012,
                   **{k:1.5 for k in SECONDARY}}
                if sp=="Dd" and tr!="Soil":
                    v["RSD"]*=0.70
                if sp=="Fr":
                    v["RSD"]*=1.5
                row={k:float(val*(1+rng.normal(0,.05)))
                     for k,val in v.items()}
                group.append(row)
            d[sp,tr]=group
    return d


def test_ratio_of_means_group_and_between_species_contrast():
    data=fake_groups()
    out,did=outcome_bootstrap(data,"RSD",n_boot=199,seed=3)
    assert len(out)==4 and len(did)==2
    dd=[r for r in out if r["species"]=="Dd" and r["conditioning"]=="Forbs"][0]
    ff=[r for r in out if r["species"]=="Fr" and r["conditioning"]=="Forbs"][0]
    assert dd["ratio_of_means"]<.8
    assert .88<ff["ratio_of_means"]<1.12
    assert did[0]["Dd_minus_Fr_log_ratio"]<-.2


def test_same_response_control_shared_within_species():
    out,_=outcome_bootstrap(fake_groups(),"PDW",n_boot=100,seed=23)
    control=[r for r in out if r["species"]=="Dd"]
    assert control[0]["control_mean"]==control[1]["control_mean"]


def test_holm_is_monotone_and_never_smaller_than_raw():
    p=[.1,.01,.6,.03,.007,.9]
    q=holm(p)
    assert len(q)==6
    assert all(y>=x for x,y in zip(p,q))
    assert all(y<=1 for y in q)
    order=np.argsort(p)
    assert all(q[order[i]]<=q[order[i+1]] for i in range(len(q)-1))


def test_permutation_cannot_reject_identical_data():
    rng=np.random.default_rng(11)
    x=np.ones(10)
    y=np.ones(10)
    assert treatment_label_perm(x,y,rng,n_perm=99)==1
