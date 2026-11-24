"""AHIF pilot source, numerical and causal-claim validation guards."""
from __future__ import annotations

import csv
from pathlib import Path
import sqlite3

import pytest

from benchmarks.build_ahif_db import (
    FIELDS, load, validate, build,
)

ROOT=Path("database/ahif_v1")


def fixtures():
    return {k:load(ROOT,k) for k in FIELDS}


def test_real_pilot_is_source_audited_but_not_pooled(tmp_path):
    summary=build(ROOT,tmp_path)
    assert summary["study_records"]>=20
    assert summary["experimental_panels"]>=20
    assert summary["aggregate_effect_rows"]>=40
    assert summary["source_audited_studies"]>=8
    assert summary["formal_meta_analysis_ready"] is False
    assert summary["panels_with_identified_causal_interaction_to_storage_path"]==0
    db=sqlite3.connect(tmp_path/"ahif_v1.sqlite")
    assert db.execute("PRAGMA foreign_key_check").fetchall()==[]
    assert db.execute("SELECT COUNT(*) FROM effects").fetchone()[0]==summary["aggregate_effect_rows"]
    assert db.execute("SELECT COUNT(*) FROM v_bridge_gate WHERE evidence_status='identified_joint_mechanism_to_storage'").fetchone()[0]==0
    db.close()


def test_review_268_not_automatically_curated_as_empirical_effect():
    r=fixtures()
    assert len(r["studies"])<268
    assert "Zou_2023_review_268" not in {x["study_id"] for x in r["studies"]}
    assert not any(e["contrast"].startswith("review_study_")
                   for e in r["effects"])


def test_same_original_jing_family_is_not_double_counted():
    r=fixtures()
    x={s["study_id"]:s for s in r["studies"]}
    assert x["ST_JINGZENODO2023"]["publication_cluster"]==\
           x["ST_JINGDRYAD2024"]["publication_cluster"]
    assert x["ST_JINGDRYAD2024"]["raw_source_audited"]=="no"
    assert x["ST_JINGZENODO2023"]["raw_source_audited"]=="yes"


def test_poem_cohorts_same_publication_cluster():
    r=fixtures()
    pp=[p for p in r["panels"] if p["study_id"]=="ST_POEM2025"]
    assert len(pp)==2
    assert {p["panel_id"] for p in pp}=={"P_POEM2020","P_POEM2021"}


def test_invalid_causal_claim_refused():
    r=fixtures()
    z=next(p for p in r["panels"] if p["panel_id"]=="P_VANSTEIJN")
    z["same_panel_causal_interaction_to_structure_identified"]="yes"
    with pytest.raises(ValueError,match="Unjustified"):
        validate(r)


def test_invalid_p_value_refused():
    r=fixtures()
    z=next(e for e in r["effects"] if e["effect_id"]=="E_DG_NAMED")
    z["p_value"]="1.1"
    with pytest.raises(ValueError,match="Invalid p"):
        validate(r)


def test_invalid_confidence_interval_refused():
    r=fixtures()
    z=next(e for e in r["effects"] if e["effect_id"]=="E_FRA_ORDER")
    z["ci95_low"]="0.5"  # above observed +.40021
    with pytest.raises(ValueError,match="Misaligned"):
        validate(r)


def test_uncertain_methods_not_automatically_promoted():
    r=fixtures()
    s=next(s for s in r["studies"] if s["study_id"]=="ST_JINGDRYAD2024")
    s["screen_status"]="scored_draft"
    with pytest.raises(ValueError,match="Unverified scored"):
        validate(r)


def test_null_and_non_significant_effects_retained():
    r=fixtures()
    assert any(e["p_value"] and float(e["p_value"])>.05 for e in r["effects"])
    assert any(e["estimate"] and float(e["estimate"])<0 for e in r["effects"])
    assert any(e["p_value"]=="" for e in r["effects"])


def test_orphan_effect_id_refused():
    r=fixtures()
    z=next(e for e in r["effects"] if e["effect_id"]=="E_LEO_NAMED")
    z["panel_id"]="P_NEVER_MEASURED"
    with pytest.raises(ValueError,match="Effect orphan"):
        validate(r)


def test_ten_substantive_prior_art_checks_not_independent_effect_rows(tmp_path):
    result=build(ROOT,tmp_path)
    assert result["substantively_verified_prior_art_papers_n"]==10
    assert result["source_audited_studies"]>=9
    assert result["aggregate_effect_rows"]==161
    assert result["study_records"]==33
    assert result["experimental_panels"]==35
    assert result["source_audited_studies"]==11
    assert result["panels_with_named_and_unlabeled_rank_measured"]==9
    assert result["panels_with_multispecies_or_functional_group_named_rank_n"]==8
    assert result["panels_with_two_species_top1_named_rank_n"]==1
    assert result["methodological_novelty_threats_n"]==5
    with sqlite3.connect(tmp_path/"ahif_v1.sqlite") as db:
        rows=db.execute(
            "SELECT study_id,source_raw_audited,"
            "causally_independent_interaction_to_rank_validated "
            "FROM prior_art_verified ORDER BY prior_art_id").fetchall()
        assert len(rows)==10
        assert all(r[2]=="no" for r in rows)
        row=db.execute(
            "SELECT source_raw_audited FROM prior_art_verified WHERE study_id='ST_GENERATIONS2026'"
        ).fetchone()
        assert row==("yes",)
    assert result["formal_meta_analysis_ready"] is False
    assert result["review_candidates_promoted_to_independent_AHIF_evidence"]==0


def test_methodological_precedents_are_separate_from_empirical_effects(tmp_path):
    summary=build(ROOT,tmp_path)
    assert summary["methodological_novelty_threats_n"]==5
    assert summary["substantively_verified_prior_art_papers_n"]==10
    with sqlite3.connect(tmp_path/"ahif_v1.sqlite") as db:
        r=db.execute("SELECT record_id,not_independent_effect,doi_or_unverified FROM methodological_prior_art").fetchall()
        assert len(r)==5
        assert all(z[1]=="yes" for z in r)
        assert db.execute("SELECT COUNT(*) FROM effects").fetchone()[0]==161
        assert any(z[2]=="doi_unverified" for z in r)


def test_new_2026_field_study_has_repeated_40_plot_grain_and_no_mechanism_claim():
    rows=fixtures()
    study=next(x for x in rows["studies"] if x["study_id"]=="ST_VANSTEIJNFIELD2026")
    panel=next(x for x in rows["panels"] if x["panel_id"]=="P_VANSTEIJNFIELD2026")
    assert study["raw_source_audited"]=="yes"
    assert panel["n_scored_units_or_unknown"]=="40"
    assert panel["named_structure_scored"]=="yes"
    assert panel["identity_free_rank_scored"]=="yes"
    assert panel["same_panel_causal_interaction_to_structure_identified"]=="no"
    effect=[e for e in rows["effects"] if e["panel_id"]=="P_VANSTEIJNFIELD2026"]
    assert len(effect)==75
    assert sum(e["endpoint_class"]=="structure" for e in effect)==16
    assert sum(e["endpoint_class"]=="structure_nontarget_roster" for e in effect)==16
    assert sum(e["endpoint_class"]=="structure_sown_roster" for e in effect)==16
    assert sum(e["endpoint_class"]=="structure_trajectory" for e in effect)==4
    assert sum(e["endpoint_class"]=="structure_trajectory_contrast" for e in effect)==2
    source_years={y for y in ("2021","2022","2023","2024") if
                  any(y in e["contrast"] for e in effect)}
    assert len(source_years)==4


def test_carlstrom_metadata_not_promoted_to_effect_source(tmp_path):
    build(ROOT,tmp_path)
    with sqlite3.connect(tmp_path/"ahif_v1.sqlite") as db:
        assert db.execute("SELECT raw_source_audited FROM studies WHERE study_id=?",("ST_CARLSTROM2019",)).fetchone()==("no",)
        assert db.execute("SELECT source_raw_audited,identity_free_rank_scored_in_ahif,causally_independent_interaction_to_rank_validated FROM prior_art_verified WHERE study_id=?",("ST_CARLSTROM2019",)).fetchone()==("no","no","no")
        assert db.execute("SELECT COUNT(*) FROM effects WHERE panel_id=?",("P_CARLSTROM",)).fetchone()==(0,)

def test_routing_metadata_and_descriptive_cages_separate(tmp_path):
    summary=build(ROOT,tmp_path)
    assert summary["routing_field_test_rows"]==21
    assert summary["routing_cage_descriptive_contrast_rows"]==12
    assert summary["panels_with_raw_source_audited_named_rank_n"]==8
    assert summary["panels_with_derived_matrix_descriptive_named_rank_n"]==1
    assert not summary["routing_cage_meta_eligible"]
    with sqlite3.connect(tmp_path/"ahif_v1.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM effects WHERE panel_id='P_PRASAD_DROPOUT'").fetchone()==(0,)
        assert db.execute("SELECT COUNT(*) FROM hypotheses").fetchone()==(9,)
        assert dict(db.execute("SELECT family,COUNT(*) FROM v_routing_tests GROUP BY family"))=={"primary":2,"secondary":16,"order_exchange":3}
        assert dict(db.execute("SELECT family,COUNT(*) FROM v_routing_tests WHERE p_holm<=.05 GROUP BY family"))=={"primary":1,"secondary":8,"order_exchange":2}
        assert db.execute("SELECT COUNT(*) FROM effects WHERE effect_id LIKE 'E_POEM_2020_A2_UNSOWN_%'").fetchone()==(2,)
        assert dict(db.execute("SELECT qc_lane,COUNT(*) FROM cage_descriptive_contrasts GROUP BY qc_lane"))=={"source_qc":6,"all_exported_sensitivity":6}


def test_descriptive_cage_cannot_be_promoted_to_effect():
    r=fixtures()
    r["effects"][0]["panel_id"]="P_PRASAD_DROPOUT"
    with pytest.raises(ValueError,match="Candidate presented as scored"):
        validate(r)

def test_saine_independent_experiment_preserves_nulls_and_placebo(tmp_path):
    result=build(ROOT,tmp_path)
    assert result["saine_test_rows"]==20
    assert result["saine_ecological_postintroduction_rows"]==18
    assert result["saine_baseline_placebo_rows"]==2
    assert result["saine_eligibility_tests_not_effects"]==2
    assert result["saine_complete_logs"]==446
    assert result["saine_meta_eligible"] is False
    with sqlite3.connect(tmp_path/"ahif_v1.sqlite") as db:
        assert db.execute("SELECT COUNT(*) FROM effects WHERE panel_id='P_SAINE2025'").fetchone()==(20,)
        assert db.execute("SELECT lane,year,endpoint FROM v_saine_tests WHERE p_holm<=.05").fetchall()==[("any_introduction","2021","named")]
        assert db.execute("SELECT COUNT(*) FROM v_saine_tests WHERE endpoint='rank' AND p_holm<=.05").fetchone()==(0,)
        assert db.execute("SELECT COUNT(DISTINCT publication_cluster) FROM v_saine_tests").fetchone()==(1,)
