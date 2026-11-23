"""Synthetic only: soil-history treatment-context reconstruction and bootstrap gates."""
import numpy as np
import pytest

from benchmarks.score_jing_2023_soil_history_competition import (
    ASSEMBLAGES,MONO_SPECIES,PAIRS,SOIL,WATER,
    SourceStop,mean_log_ratio,score
)


def artificial_panel():
    rng=np.random.default_rng(71)
    mono={}
    mix={}
    base={"SK":8.,"AF":3.,"AC":12.,"PT":4.,"ES":5.}
    for w in WATER:
        for s in SOIL:
            for name in MONO_SPECIES:
                v=base[name]*(0.65 if w=="Dry" else 1)
                mono[w,s,name]=np.maximum(0.01,rng.normal(v,v*.03,4))
            for assembly in ASSEMBLAGES:
                names=assembly.split("-")
                values=[]
                for _ in range(4):
                    v=[]
                    for name in names:
                        fac=(0.42 if name=="SK" else 0.48 if name=="AF" else 0.6)
                        if s=="SK" and name=="SK":fac*=.75
                        if s=="AF" and name=="AF":fac*=1.2
                        if s=="AC" and name=="AC":fac*=.80
                        if w=="Dry" and s=="SK" and name=="SK":fac*=.9
                        v.append(max(.01,base[name]*fac*rng.uniform(.98,1.02)))
                    values.append(v)
                mix[w,s,assembly]=np.array(values)
    return mono,mix


def test_six_species_pairs_and_48_live_soil_contrasts():
    mono,mix=artificial_panel()
    res=score(mono,mix,{"status":"synthetic"},n_boot=99)
    assert len(res["primary_solo_referenced_2_species_mixture"])==60
    assert len(res["live_inoculum_vs_sterile_differences"])==48
    assert len(res["soil_drought_difference_in_differences"])==24
    assert len(res["raw_mixture_mass_winner_by_context"])==30
    assert len(res["soil_context_variation_of_log_mixture_impact"])==12
    assert all(r["n_mixture_pots"]==4 for r in res["primary_solo_referenced_2_species_mixture"])
    assert all(len(r["bootstrap95"])==2 for r in res["live_inoculum_vs_sterile_differences"])
    assert res["source_audit"]["status"]=="synthetic"


def test_single_pot_outlier_not_pseudoreplicated_as_multiple_species():
    mono,mix=artificial_panel()
    q=mix["Wet","S","SK-AF"].copy()
    mix["Wet","S","SK-AF"]=q*2
    out=score(mono,mix,{"status":"synthetic"},n_boot=99)
    key=next(x for x in out["primary_solo_referenced_2_species_mixture"]
             if x["pair"]=="SK-AF" and x["focal"]=="SK"
             and x["water"]=="Wet" and x["soil"]=="S")
    assert np.isclose(key["mixture_mean_g"],q[:,0].mean()*2)


def test_zero_reference_denominator_causes_source_stop():
    with pytest.raises(SourceStop):
        mean_log_ratio(np.ones(4),np.zeros(4),"synthetic")
