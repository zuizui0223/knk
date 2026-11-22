"""Synthetic POEM scores only; no experiment data or invented biological findings."""
from collections import Counter
import numpy as np

from benchmarks.audit_poem_design import SOURCE_SOWN
from benchmarks.score_poem_field_storage import (
    FOUR, THREE, analysis, evaluate_one, gower, history_r2, matrix_for,
    bray_curtis,
)


def _panel(cohort, age):
    names = sorted(set(SOURCE_SOWN) | {f"Supplement_{i:02}" for i in range(35)})
    assert len(names) == 49
    plots = []
    plot_map = {}
    biom = {}
    rng = np.random.default_rng(int(cohort)+age)
    for t in ("F","G","L","S","B"):
        for i in range(1,6):
            p = f"{cohort}_{t}{i}"
            plots.append(p)
            plot_map[p] = (t,str(i))
            vec = rng.uniform(0.01,2.5,size=len(names))
            # Distinct selected identities but identical sampling grain.
            if t == "F":
                vec[0] += 8.0 + age
            for n,v in zip(names,vec):
                biom[(p,n)] = float(v)
    return {
        "plot_ids":sorted(plots),"plot_map":plot_map,
        "biomass_by_plot_species":biom,
        "identified_species_names":names,
        "source_rows":100,
        "explicit_zeros":0,"duplicate_fragment_records":0,
        "quadrats_per_plot":2,
        "harvest":f"{int(cohort)+age-1}-06"
    }


def test_multigroup_permutation_r2_valid():
    y=np.array(["F"]*5+["G"]*5+["L"]*5+["S"]*5)
    x=np.ones((20,3))
    x[y=="F",0]=20
    x[y=="G",1]=20
    d=bray_curtis(x)
    g=gower(d)
    obs=history_r2(g,y)
    assert 0 <= obs <= 1+1e-7
    rng=np.random.default_rng(17)
    perm=np.asarray([rng.permutation(y) for _ in range(199)])
    r=evaluate_one(x,y,perm)
    assert r["r2"]>r["null_median"]
    assert r["p_upper"]<=.05


def test_cohort_by_age_persistent_plot_scoring():
    panels={(c,age):_panel(c,age) for c in ("2020","2021") for age in (1,2,3)}
    ans=analysis(panels,{"source_md5":"synthetic-only"},n_perm=19,seed=15)
    assert ans["design"]["primary_plots_per_cohort"]==20
    assert len(ans["primary_4level_all_identified_taxa"])==6
    assert len(ans["predeclared_sensitivities"])==6
    for row in ans["primary_4level_all_identified_taxa"]:
        assert row["n_plots"]==20 and row["n_species"]==49
        assert row["representations"]["named_all"]["n_permutations"]==19
        assert row["representations"]["sorted_all"]["n_response_columns"]==49
    assert all(
        row["sown14_4level"]["status"]=="SCORED" and
        row["all49_FGLonly"]["status"]=="SCORED"
        for row in ans["predeclared_sensitivities"]
    )


def test_same_panel_dimensions_for_named_and_sorted():
    p=_panel("2020",2)
    x,y,plot,names=matrix_for(p,FOUR)
    assert x.shape==(20,49)
    x3,y3,plot3,names3=matrix_for(p,THREE)
    assert x3.shape==(15,49)
    assert all(Counter(y)[a]==5 for a in FOUR)
    assert all(Counter(y3)[a]==5 for a in THREE)
    assert np.allclose(np.sort(x/x.sum(axis=1,keepdims=True),axis=1).sum(axis=1),1)
