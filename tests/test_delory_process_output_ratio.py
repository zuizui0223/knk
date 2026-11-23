"""Synthetic-only paired-plant log response contrast test, no original data."""
import numpy as np
from benchmarks.audit_delory_process_output_ratio import process_output_contrasts
from benchmarks.score_delory_chemical_legacy import SPECIES,TREATMENTS


def fake():
    rng=np.random.default_rng(17)
    out={}
    for sp in SPECIES:
        for tr in TREATMENTS:
            arr=[]
            for k in range(10 if sp=="Fr" else (9 if tr!="Grasses" else 10)):
                plant={"RSD":1.0,"PDW":3.0,"Nuptake":0.010}
                if sp=="Dd" and tr=="Forbs":
                    plant["RSD"]*=.65
                for name in plant:
                    plant[name]*=rng.uniform(.97,1.03)
                arr.append(plant)
            out[sp,tr]=arr
    return out


def test_ratio_process_vs_output_recovers_dissociation():
    a=process_output_contrasts(fake(),n_boot=2000)
    assert len(a)==8
    check=next(v for v in a if v["species"]=="Dd"
               and v["conditioning"]=="Forbs"
               and v["contrasted_output"]=="Nuptake")
    assert check["root_foraging_minus_output_log_response"]<-.4
    assert check["conditional_plant_bootstrap95"][1]<0


def test_control_plants_resampled_jointly_per_endpoint():
    a=process_output_contrasts(fake(),n_boot=100)
    assert all(len(v["conditional_plant_bootstrap95"])==2 for v in a)
    assert all(v["n_bootstrap"]==100 for v in a)
    assert all(v["n_treatment_plants"]>=9 for v in a)
