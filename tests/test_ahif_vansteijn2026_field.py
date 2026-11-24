"""Synthetic source-free tests, including blocked history and roster STOP gates."""
import numpy as np
from benchmarks.score_ahif_vansteijn2026_field import (
    TREATMENTS,SITES,YEARS,blocked_assignment_stream,
    treatment_r2,bray,gower,score_matrix,within_roster,analyse
)

def simulated(seed=19):
    names=[f"sp{i:02}" for i in range(76)]
    target=names[:20]
    fast=target[:10];slow=target[10:]
    rng=np.random.default_rng(seed)
    panels={}
    for site in SITES:
        plots=[(site,str(block),t) for block in range(1,6) for t in "abcd"]
        history=np.array(TREATMENTS*5)
        for year in YEARS:
            X=np.zeros((20,76))
            for row in range(20):
                # each treatment maps top biomass to different ORIGINAL taxon,
                # but every plot has the same strong top-rank shape.
                X[row]=(1.+rng.uniform(0,0.02,size=76))
                X[row,TREATMENTS.index(history[row])]=70
            panels[site,year]=dict(
                X=X,plots=plots,history=history,
                names=names,sown_target=target,fast=fast,slow=slow,
                site=site,year=year)
    return panels

def test_repeated_block_assignments_preserved_and_not_shuffle_across_blocks():
    labels=np.array(TREATMENTS*5)
    stream=blocked_assignment_stream(labels,n_perm=49,seed=23)
    assert stream.shape==(49,20)
    assert any(not np.array_equal(row,labels) for row in stream)
    assert all(set(row[i:i+4])==set(TREATMENTS)
               for row in stream for i in range(0,20,4))

def test_same_abundance_rank_can_hide_source_treatment_identity():
    state=simulated()["MAL","2022"]
    x=state["X"]/state["X"].sum(axis=1,keepdims=True)
    per=blocked_assignment_stream(state["history"],n_perm=199)
    named=score_matrix(x,state["history"],per)
    rank=score_matrix(np.sort(x,axis=1)[:,::-1],state["history"],per)
    assert named["p_upper"]<=.05
    assert named["E"]>rank["E"]
    assert named["R2"]>rank["R2"]

def test_full_eight_site_years_and_two_single_plot_trajectories():
    P=simulated()
    result=analyse(P,{"status":"SYNTHETIC_ONLY"},n_perm=19,n_boot=15)
    assert len(result["primary_full76_4_history_by_site_year"])==8
    assert len(result["two_site_4yr_whole_plot_trajectories"])==2
    assert all(z["n_source_plot_units"]==20 for z in result["primary_full76_4_history_by_site_year"])
    assert all(z["n_features"]==304 for z in result["two_site_4yr_whole_plot_trajectories"])
    assert all(z["representations"]["named76"]["n_permutations"]==19
               for z in result["primary_full76_4_history_by_site_year"])
    assert all(z["sown20"]["status"]=="SCORED_CONDITIONAL_ON_GUILD_MEMBERSHIP"
               for z in result["source_predeclared_secondary_roster_and_scalar_sensitivities"])

def test_zero_guild_must_STOP_without_pseudocount_or_row_removal():
    p=simulated()["UME","2021"]
    z=p.copy()
    z["X"]=p["X"].copy()
    z["X"][0,:10]=0
    per=blocked_assignment_stream(p["history"],n_perm=9)
    out=within_roster(z,p["fast"],per)
    assert out["status"]=="STOP_ZERO_SUBCOMMUNITY_TOTAL"
    assert out["n_unscorable_whole_plots"]==1

def test_equal_named_and_rank_shapes_and_nonnegative_bray():
    P=simulated()
    x=P["MAL","2021"]["X"]
    q=x/x.sum(axis=1,keepdims=True)
    assert q.shape==(20,76)
    assert np.allclose(q.sum(axis=1),1)
    D=bray(q)
    assert D.shape==(20,20)
    assert np.all(np.diag(D)==0)
    assert np.all((D>=0)&(D<=1))
    assert np.isfinite(gower(D)).all()
