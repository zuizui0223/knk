#!/usr/bin/env python3
"""Build AHIF-DB: source-linked experiment/effect/hypothesis SQLite and audits.

Unit of independent source evidence = publication cluster, not pots, taxa,
repeated snapshots, directed pairs or correlated endpoints.
This is a curated PILOT and does not ingest 268 independently scored effects.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
from datetime import date
import json
import math
from pathlib import Path
import shutil
import sqlite3

from benchmarks.audit_zou2023_review_source import STAGING_COLUMNS, SRC_GIT_BLOB

FILES={"studies":"studies.csv","panels":"panels.csv",
       "effects":"effects.csv","hypotheses":"hypotheses.csv"}
FIELDS={
"studies":("study_id","year","citation","domain","doi","source_data_url","screen_status","raw_source_audited","publication_cluster","lab_or_overlap_group","scope_note","evidence_anchor_url"),
"panels":("panel_id","study_id","intervention_lane","putative_mechanism_channel","n_scored_units_or_unknown","independent_unit_type","response_grain","named_structure_scored","identity_free_rank_scored","partner_or_absence_control","biological_process_measured","community_structure_scope","aggregate_output_scope","same_panel_causal_interaction_to_structure_identified","source_status","scope_caveat"),
"effects":("effect_id","panel_id","endpoint_class","estimand","contrast","estimate","ci95_low","ci95_high","p_value","unit_or_statistic","experimental_grain","analysis_status","evidence_url","limitations"),
"hypotheses":("hypothesis_id","testable_hypothesis","intervention_family","required_joint_endpoints","required_independence_grain","falsifiable_prediction","current_source_evidence","interpretation_status","relevant_panel_ids_semicolon")
}
STUDY_STATES={"scored_frozen","scored_draft","scored_draft_overlap_unresolved",
 "source_binary_blocked","candidate_unscored","paper_only","proxy_only"}
EFFECT_STATES={"frozen_knk","draft_knk","source_reanalysis",
 "posthoc_exploratory","posthoc_sensitivity","sensitivity"}
CHANNELS={"live_arrival","live_neighbor","live_competitor","conditioned_environment",
 "spatial_preemption","microbial_colonization","microbial_founder",
 "rhizobium_mutualism","resident_removal","succession"}


def load(source:Path,table:str):
    p=source/FILES[table]
    with p.open(encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f)
        if tuple(r.fieldnames or ())!=FIELDS[table]:
            raise ValueError(f"Unexpected {table} schema")
        values=list(r)
    if not values or any(None in row or any(v is None for v in row.values()) for row in values):
        raise ValueError(f"Truncated or empty {table}")
    return values


def num(z,label,required=False):
    if z.strip()=="":
        if required:raise ValueError("Missing "+label)
        return None
    try:f=float(z)
    except ValueError:raise ValueError("Not numeric "+label)
    if not math.isfinite(f):raise ValueError("Nonfinite "+label)
    return f


def validate(tables):
    studies,panels,effects,hypotheses=(tables[x] for x in FILES)
    seen=set()
    citations=set()
    for s in studies:
        k=s["study_id"]
        if k in seen or not k.startswith("ST_"):raise ValueError("Duplicate study "+k)
        seen.add(k)
        y=num(s["year"],"year",True)
        if int(y)!=y or y<1900 or y>date.today().year:raise ValueError("Invalid year "+k)
        if s["domain"] not in {"plants","microbes","fungi","arthropods"}:raise ValueError("Domain "+k)
        if s["screen_status"] not in STUDY_STATES:raise ValueError("Study status "+k)
        if s["raw_source_audited"] not in {"yes","no"}:raise ValueError("Audit status "+k)
        if s["screen_status"].startswith("scored") and s["raw_source_audited"]!="yes":
            raise ValueError("Unverified scored study "+k)
        if not s["publication_cluster"] or not s["lab_or_overlap_group"]:
            raise ValueError("Missing overlap family "+k)
        if not all(s[f].startswith("https://") for f in ("source_data_url","evidence_anchor_url")):
            raise ValueError("Missing evidence link "+k)
        if s["doi"]:
            if not s["doi"].startswith("10.") or s["doi"] in citations:
                raise ValueError("Duplicate or invalid DOI "+k)
            citations.add(s["doi"])
    ids=set()
    for p in panels:
        k=p["panel_id"]
        if k in ids or k[:2]!="P_":raise ValueError("Duplicate panel "+k)
        ids.add(k)
        if p["study_id"] not in seen:raise ValueError("Panel orphan "+k)
        if p["putative_mechanism_channel"] not in CHANNELS:raise ValueError("Channel "+k)
        if any(p[x] not in {"yes","no","unknown"} for x in ("named_structure_scored","identity_free_rank_scored")):
            raise ValueError("Evidence status "+k)
        if p["identity_free_rank_scored"]=="yes" and p["named_structure_scored"]!="yes":
            raise ValueError("Identity-free rank needs original species vector "+k)
        if p["same_panel_causal_interaction_to_structure_identified"] not in {"no","unknown"}:
            raise ValueError("Unjustified same-panel causal mediation "+k)
        if p["source_status"] in {"candidate","paper_only","binary_source_blocked","proxy_excluded"} and \
           p["same_panel_causal_interaction_to_structure_identified"]!="unknown":
            raise ValueError("Unscreened study must retain unknown causal-bridge status "+k)
        n=num(p["n_scored_units_or_unknown"],"number of units")
        if n is not None and (int(n)!=n or n<2):raise ValueError("Unit count "+k)
    eid=set()
    panelmap={p["panel_id"]:p for p in panels}
    for e in effects:
        k=e["effect_id"]
        if k in eid or not k.startswith("E_"):raise ValueError("Duplicate effect "+k)
        eid.add(k)
        if e["panel_id"] not in ids:raise ValueError("Effect orphan "+k)
        if e["analysis_status"] not in EFFECT_STATES:raise ValueError("Effect status "+k)
        v=num(e["estimate"],"estimate",True)
        lo=num(e["ci95_low"],"ci_low");hi=num(e["ci95_high"],"ci_high")
        p=num(e["p_value"],"p_value")
        if (lo is None)!=(hi is None):raise ValueError("Partial CI "+k)
        if lo is not None and (v<lo-1e-8 or v>hi+1e-8 or lo>hi):
            raise ValueError("Misaligned CI "+k)
        if p is not None and not 0<=p<=1:raise ValueError("Invalid p "+k)
        if e["estimand"]=="E_null_adj_R2" and e["unit_or_statistic"]!="E_distR2":
            raise ValueError("False equivalence of null-calibrated R2 and information "+k)
        if not e["evidence_url"].startswith("https://") or not e["limitations"]:
            raise ValueError("Effect missing source or caveat "+k)
        if panelmap[e["panel_id"]]["source_status"] in {"candidate","paper_only","binary_source_blocked"}:
            raise ValueError("Candidate presented as scored "+k)
    hids=set()
    for h in hypotheses:
        k=h["hypothesis_id"]
        if k in hids:raise ValueError("Duplicate hypothesis "+k)
        hids.add(k)
        if not h["falsifiable_prediction"] or not h["required_independence_grain"]:
            raise ValueError("Untestable hypothesis "+k)
        for ref in h["relevant_panel_ids_semicolon"].split(";"):
            if ref and ref not in ids:raise ValueError("Hypothesis panel orphan "+k)
    return {"n_studies":len(studies),"n_panels":len(panels),
            "n_effects":len(effects),"n_hypotheses":len(hypotheses)}


def schema():
    return """
PRAGMA foreign_keys=ON;
CREATE TABLE studies(
study_id TEXT PRIMARY KEY,year INTEGER NOT NULL,citation TEXT NOT NULL,
domain TEXT NOT NULL,doi TEXT,source_data_url TEXT NOT NULL,
screen_status TEXT NOT NULL,raw_source_audited TEXT NOT NULL,
publication_cluster TEXT NOT NULL,lab_or_overlap_group TEXT NOT NULL,
scope_note TEXT NOT NULL,evidence_anchor_url TEXT NOT NULL
);
CREATE TABLE panels(
panel_id TEXT PRIMARY KEY,study_id TEXT NOT NULL REFERENCES studies(study_id),
intervention_lane TEXT NOT NULL,putative_mechanism_channel TEXT NOT NULL,
n_scored_units_or_unknown INTEGER,independent_unit_type TEXT NOT NULL,
response_grain TEXT NOT NULL,named_structure_scored TEXT NOT NULL,
identity_free_rank_scored TEXT NOT NULL,partner_or_absence_control TEXT NOT NULL,
biological_process_measured TEXT NOT NULL,community_structure_scope TEXT NOT NULL,
aggregate_output_scope TEXT NOT NULL,
same_panel_causal_interaction_to_structure_identified TEXT NOT NULL,
source_status TEXT NOT NULL,scope_caveat TEXT NOT NULL
);
CREATE TABLE effects(
effect_id TEXT PRIMARY KEY,panel_id TEXT NOT NULL REFERENCES panels(panel_id),
endpoint_class TEXT NOT NULL,estimand TEXT NOT NULL,contrast TEXT NOT NULL,
estimate REAL NOT NULL,ci95_low REAL,ci95_high REAL,p_value REAL,
unit_or_statistic TEXT NOT NULL,experimental_grain TEXT NOT NULL,
analysis_status TEXT NOT NULL,evidence_url TEXT NOT NULL,
limitations TEXT NOT NULL,
CHECK(ci95_low IS NULL OR ci95_low<=estimate),
CHECK(ci95_high IS NULL OR estimate<=ci95_high),
CHECK(p_value IS NULL OR (p_value>=0 AND p_value<=1))
);
CREATE TABLE review_candidates(
review_record_id TEXT PRIMARY KEY,
source_excel_row INTEGER NOT NULL UNIQUE,
source_review_doi TEXT NOT NULL,source_git_blob TEXT NOT NULL,
authors TEXT NOT NULL,year INTEGER NOT NULL,journal TEXT,title TEXT NOT NULL,
focal_taxa TEXT,secondary_taxa TEXT,manipulated_species TEXT,
temporal_gradient TEXT,generation_time TEXT,duration_generations TEXT,
competition_model TEXT,original_review_inclusion TEXT,
ahif_independent_screen_status TEXT NOT NULL
);
CREATE INDEX idx_review_candidate_taxa ON review_candidates(focal_taxa,year);
CREATE TABLE hypotheses(
hypothesis_id TEXT PRIMARY KEY,testable_hypothesis TEXT NOT NULL,
intervention_family TEXT NOT NULL,required_joint_endpoints TEXT NOT NULL,
required_independence_grain TEXT NOT NULL,falsifiable_prediction TEXT NOT NULL,
current_source_evidence TEXT NOT NULL,interpretation_status TEXT NOT NULL,
relevant_panel_ids_semicolon TEXT NOT NULL
);
CREATE INDEX idx_study_cluster ON studies(publication_cluster);
CREATE INDEX idx_study_lab ON studies(lab_or_overlap_group);
CREATE INDEX idx_panel_treatment ON panels(intervention_lane,putative_mechanism_channel);
CREATE INDEX idx_effect_estimand ON effects(estimand,endpoint_class);
CREATE VIEW v_effect_lineage AS SELECT e.*,p.intervention_lane,p.putative_mechanism_channel,
p.independent_unit_type,s.study_id,s.publication_cluster,s.lab_or_overlap_group
FROM effects e JOIN panels p USING(panel_id) JOIN studies s USING(study_id);
CREATE VIEW v_bridge_gate AS SELECT p.panel_id,p.study_id,s.publication_cluster,
p.named_structure_scored,p.identity_free_rank_scored,
p.same_panel_causal_interaction_to_structure_identified,
CASE WHEN p.named_structure_scored='yes' AND p.identity_free_rank_scored='yes'
AND p.same_panel_causal_interaction_to_structure_identified='yes'
AND s.raw_source_audited='yes'
THEN 'identified_joint_mechanism_to_storage'
ELSE 'not_identified' END AS evidence_status
FROM panels p JOIN studies s USING(study_id);
"""


def insert_external_review_stage(connection,source_path:Path|None):
    """Citation staging ONLY: never promotes 268 review rows into scored effects."""
    if source_path is None:
        return 0
    with source_path.open(encoding="utf-8-sig",newline="") as file:
        reader=csv.DictReader(file)
        if tuple(reader.fieldnames or ())!=STAGING_COLUMNS:
            raise ValueError("External literature staging schema changed")
        records=list(reader)
    if len(records)!=268:
        raise ValueError("Expected 268 original review bibliographic rows")
    unique_ids=set()
    row_numbers=set()
    data=[]
    for row in records:
        key=row["review_record_id"]
        numrow=int(row["source_excel_row"])
        if key in unique_ids or numrow in row_numbers or numrow<2 or numrow>269:
            raise ValueError("Repeated or invalid original review row")
        unique_ids.add(key)
        row_numbers.add(numrow)
        if row["source_git_blob"]!=SRC_GIT_BLOB or \
           row["source_review_doi"]!="10.1016/j.tree.2023.08.001":
            raise ValueError("Unpinned review source")
        if row["ahif_independent_screen_status"]!="not_independently_screened":
            raise ValueError("Unverified review study promoted to curated source")
        if not row["title"] or not row["authors"]:
            raise ValueError("Empty review citation title/author")
        values=[]
        for col in STAGING_COLUMNS:
            v=row[col]
            values.append(int(v) if col in ("source_excel_row","year") else v)
        data.append(values)
    placeholders=",".join("?" for _ in STAGING_COLUMNS)
    connection.executemany(
        f"INSERT INTO review_candidates VALUES({placeholders})",data
    )
    if connection.execute("SELECT COUNT(*) FROM review_candidates").fetchone()[0]!=268:
        raise ValueError("Staging output lost review references")
    return 268


def report(tables,review_staged_n:int=0):
    s,p,e,h=(tables[k] for k in FILES)
    audited=[x for x in s if x["raw_source_audited"]=="yes"]
    architecture=[x for x in p if x["named_structure_scored"]=="yes" and x["identity_free_rank_scored"]=="yes"]
    joint=[x for x in p if x["same_panel_causal_interaction_to_structure_identified"]=="yes"]
    overlaps=defaultdict(list)
    for x in s:overlaps[x["lab_or_overlap_group"]].append(x["study_id"])
    overlaps={k:v for k,v in overlaps.items() if len(v)>1}
    sources={x["publication_cluster"] for x in s}
    pstats=Counter(x["intervention_lane"] for x in p)
    estat=Counter(x["analysis_status"] for x in e)
    ps=[num(x["p_value"],"p") for x in e]
    return {
      "result_status":"CURATED_TARGETED_PILOT_NOT_SYSTEMATIC_OR_META_POOLED",
      "db_version":"AHIF-DB v0.1",
      "generated_on":date.today().isoformat(),
      "study_records":len(s),"experimental_panels":len(p),
      "aggregate_effect_rows":len(e),"falsifiable_hypotheses":len(h),
      "source_audited_studies":len(audited),
      "candidate_publication_clusters":len(sources),
      "source_audited_publication_clusters":len({x["publication_cluster"] for x in audited}),
      "panels_with_named_and_unlabeled_rank_measured":len(architecture),
      "panels_with_identified_causal_interaction_to_storage_path":len(joint),
      "unscreened_panels_with_joint_path_unknown":sum(
          x["same_panel_causal_interaction_to_structure_identified"]=="unknown"
          for x in p),
      "audited_panels_with_joint_path_not_identified":sum(
          x["same_panel_causal_interaction_to_structure_identified"]=="no"
          for x in p),
      "possible_publication_or_lab_dependence_flags":overlaps,
      "interventions_by_panel":dict(pstats),
      "effects_by_analysis_status":dict(estat),
      "unpooled_effect_p_value_counts_not_independent_test_results":{
          "p_le_005":sum(x is not None and x<=.05 for x in ps),
          "p_gt_005":sum(x is not None and x>.05 for x in ps),
          "p_not_reported":sum(x is None for x in ps)
      },
      "review_bibliographic_candidates_staged_not_audited":review_staged_n,
      "review_candidates_promoted_to_independent_AHIF_evidence":0,
      "Zou_Rudolf_2023_prior_review_studies_not_ingested_as_effects":268,
      "formal_meta_analysis_ready":False,
      "why_not_pooled":[
         "No common study-level paired SE across >=5 independent publication clusters within one intervention/endpoint lane",
         "Some source families share labs or possible experimental overlap, including Jing Zenodo 2023 and Dryad 2024",
         "R2_E structure, solo-referenced biomass, partner-specific RII, leaf occupancy and root/soil output differ in estimand",
         "No audited source panel jointly identifies causal interaction formation and storage in named-versus-unlabeled abundance structure"
      ]
    }


def gap_report(r):
    return "\n".join([
      "# AHIF-DB v0.1 — novelty and source-coverage gate",
      "",
      "Targeted preliminary curation only: NOT a systematic survey or a random-effects meta-analysis.",
      "",
      f"- Publication records: {r['study_records']}",
      f"- Experimental panels, including cohorts and within-paper arms: {r['experimental_panels']}",
      f"- Source-audited publication records: {r['source_audited_studies']}",
      f"- Aggregate numeric evidence records: {r['aggregate_effect_rows']}",
      f"- Prespecified, falsifiable hypotheses: {r['falsifiable_hypotheses']}",
      f"- Panels with named and identity-free rank scoring: {r['panels_with_named_and_unlabeled_rank_measured']}",
      f"- Source-audited panels where the mechanism-to-structure causal path was NOT identified: {r['audited_panels_with_joint_path_not_identified']}",
      f"- Unscreened panels where the path status remains UNKNOWN: {r['unscreened_panels_with_joint_path_unknown']}",
      f"- Panels with verified causal interaction-formation to identity-vs-rank mediation: {r['panels_with_identified_causal_interaction_to_storage_path']}",
      "",
      "The count of zero applies only to THIS targeted source-audited pilot. It must never be quoted as proof that no prior paper anywhere measured both mechanisms and structure.",
      "",
      "## Proposed novelty, NOT yet validated",
      "",
      "Test, in the SAME experimental community, whether a manipulated early resident changes the later species net neighbor effect, and whether this history-specific interaction predicts which species wins versus how strongly abundances are concentrated. Add live-resident removal, soil-history replacement, and root barriers to separate mechanism channels.",
      "",
      "Prior art already establishes priority effects, competitor traits, niche pre-emption/modification and context-dependent outcomes. Our intended contribution is an auditable mechanism-to-storage crosswalk, an explicitly causal same-panel identification gate, and generalization to unseen taxa/experimental communities.",
      "",
      "## Guardrails",
      "",
      "The 268 articles reviewed in Zou & Rudolf 2023 form an external potential screening frame, not 268 re-extracted effects or 268 independent AHIF studies. Research and lab clusters, repeated field cohorts and shared source batches must be distinguished.",
      "",
      "All p-values and confidence intervals describe the original study or exploratory source reanalysis at its actual unit; the p count is NOT an independent finding tally.",
      "",
      "Meta-analysis gate: CLOSED. The 2023 Jing and 2024 Jing records may overlap. See ahif_v1_validation.json and independent source URLs for each experiment.",
      ""
    ])+"\n"


def build(source:Path,out:Path,review_screening:Path|None=None):
    tables={t:load(source,t) for t in FILES}
    counts=validate(tables)
    if counts["n_studies"]<20 or counts["n_effects"]<40:
        raise ValueError("Pilot unexpectedly lost sources or evidence")
    out.mkdir(parents=True,exist_ok=True)
    db=out/"ahif_v1.sqlite"
    if db.exists():db.unlink()
    connection=sqlite3.connect(db)
    try:
        connection.executescript(schema())
        for table in FILES:
            converted=[]
            for row in tables[table]:
                values=[]
                for col in FIELDS[table]:
                    s=row[col]
                    if s=="":
                        values.append(None)
                    elif col in ("year","n_scored_units_or_unknown"):
                        values.append(int(s))
                    elif table=="effects" and col in ("estimate","ci95_low","ci95_high","p_value"):
                        values.append(float(s))
                    else:
                        values.append(s)
                converted.append(values)
            placeholders=",".join("?" for _ in FIELDS[table])
            connection.executemany(f"INSERT INTO {table} VALUES({placeholders})",converted)
            found=connection.execute(f"SELECT COUNT(*) FROM {table}").fetchone()[0]
            if found!=len(tables[table]):raise ValueError("Source output row-count mismatch")
        n_review_candidates=insert_external_review_stage(connection,review_screening)
        if connection.execute("PRAGMA foreign_key_check").fetchall():
            raise ValueError("Foreign key audit failure")
        if connection.execute("SELECT COUNT(*) FROM v_bridge_gate WHERE evidence_status='identified_joint_mechanism_to_storage'").fetchone()[0]!=0:
            raise ValueError("Unjustified causal pathway promotion")
        connection.commit()
    finally:
        connection.close()
    result=report(tables,review_staged_n=n_review_candidates)
    (out/"ahif_v1_validation.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    (out/"ahif_v1_gap_report.md").write_text(gap_report(result))
    if review_screening is not None:
        shutil.copyfile(review_screening,out/"review_candidates_unreviewed_2023.csv")
    for t,filename in FILES.items():
        shutil.copyfile(source/filename,out/filename)
    print("AHIF_DB_BUILD_SUCCESS",json.dumps({
        "studies":result["study_records"],"panels":result["experimental_panels"],
        "effects":result["aggregate_effect_rows"],
        "joint_mechanism_storage":result["panels_with_identified_causal_interaction_to_storage_path"],
        "source_audited_studies":result["source_audited_studies"],
        "unreviewed_review_candidates":result["review_bibliographic_candidates_staged_not_audited"]
    }))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-root",type=Path,default=Path("database/ahif_v1"))
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--review-screening-csv",type=Path,
                   help="Optional unreviewed bibliography staging; no effect promotion")
    a=p.parse_args()
    build(a.source_root,a.output_dir,review_screening=a.review_screening_csv)


if __name__=="__main__":main()
