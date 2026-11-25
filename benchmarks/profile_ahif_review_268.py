#!/usr/bin/env python3
"""Reproduce original Zou & Rudolf 2023 bibliographic design coding and audit bins.

Input = hash-pinned, candidate-only 268 bibliography records staged by
benchmarks/audit_zou2023_review_source.py from SupportingData.xlsx.
No study has been independently eligible/source-scored from this file.
Ref: doi:10.1016/j.tree.2023.08.001; original LitReview.R rows 45-53.
"""
from __future__ import annotations
import argparse,csv,json,math
from collections import Counter,defaultdict
from pathlib import Path

SOURCE_BLOB="be2777d4e4bb1454db393b94a2118ada79118dfa"
SOURCE_DOI="10.1016/j.tree.2023.08.001"
EXPECTED=268
PRIOR_R_BINS=("No temporal treatment","1-2","3-10",">10")
FIXED_BINS=("No temporal treatment","1-2","3-10","11+")
GENERATIONS=("<=1","2-5",">5","Mixed","Not Applicable")

# High-relevance titles from this exact author-coded review, identified before
# the 13 source-coded model candidates were merged; review-level triage only.
AHIF_NEAREST_PRIOR_ART_A={
    "ZOU2023_R0028":"resource specialization and historical contingency",
    "ZOU2023_R0032":"assembly history links biodiversity and productivity; 18 species binned incorrectly in review R code",
    "ZOU2023_R0055":"plant-soil feedback and assembly history",
    "ZOU2023_R0075":"fungal arrival links composition and decomposition",
    "ZOU2023_R0127":"environment variation, microbial order and consumer-resource model",
    "ZOU2023_R0172":"arrival × fungivore × nitrogen factorial with fungal prevalence and wood mass loss",
    "ZOU2023_R0199":"founder effects and spatial pre-emption",
    "ZOU2023_R0212":"modern coexistence-theory interaction estimates",
    "ZOU2023_R0224":"competition model predicts grassland community composition",
    "ZOU2023_R0249":"warming × fungal arrival and decomposition",
    "ZOU2023_R0264":"priority effects on microbiome taxonomy and function",
}


def source_rows(path:Path):
    with path.open("r",encoding="utf-8-sig",newline="") as handle:
        reader=csv.DictReader(handle)
        fields=set(reader.fieldnames or [])
        needed={"review_record_id","source_excel_row","source_review_doi",
            "source_git_blob","authors","year","title","focal_taxa",
            "manipulated_species","temporal_gradient","duration_generations",
            "competition_model","ahif_independent_screen_status"}
        if not needed<=fields:raise ValueError("Review candidate columns missing")
        rows=list(reader)
    if len(rows)!=EXPECTED:raise ValueError(f"Expected 268 source review rows, got {len(rows)}")
    ids=set()
    xls_rows=set()
    for row in rows:
        key=row["review_record_id"]
        original_row=int(row["source_excel_row"])
        if key in ids or original_row in xls_rows:
            raise ValueError("Duplicate source review record/Excel row")
        ids.add(key);xls_rows.add(original_row)
        if row["source_git_blob"]!=SOURCE_BLOB or row["source_review_doi"]!=SOURCE_DOI:
            raise ValueError("Unpinned review bibliography source")
        if row["ahif_independent_screen_status"]!="not_independently_screened":
            raise ValueError("Independent AHIF review not supported")
        if row["temporal_gradient"] not in ("Yes","No"):
            raise ValueError("Unknown source temporal-gradient code")
        if row["duration_generations"] not in GENERATIONS:
            raise ValueError("Unknown source generation duration")
        if not row["authors"] or not row["title"]:
            raise ValueError("Blank source bibliographic metadata")
        if not 1980<=int(row["year"])<=2023:
            raise ValueError("Unexpected review source year")
    return rows


def authors_R_bin(code:str)->str:
    # Faithful to original R's case_when in LitReview.R lines 45-49
    if code in ("1","2"):return "1-2"
    if code in tuple(str(k) for k in range(3,11)):return "3-10"
    if code==">10":return ">10"
    return "No temporal treatment"


def numerically_corrected_bin(code:str)->str:
    if code=="NA (no temporal treatment)":return "No temporal treatment"
    if code==">10":return "11+"
    try:n=int(code)
    except ValueError:raise ValueError(f"Unresolved manipulated species coding {code!r}")
    if n<=0:raise ValueError("Nonpositive manipulated species count")
    if n<=2:return "1-2"
    if n<=10:return "3-10"
    return "11+"


def profile(rows):
    result=[]
    for r in rows:
        orig=authors_R_bin(r["manipulated_species"])
        fixed=numerically_corrected_bin(r["manipulated_species"])
        result.append({
            "review_record_id":r["review_record_id"],
            "source_excel_row":r["source_excel_row"],
            "authors":r["authors"],"year":r["year"],"title":r["title"],
            "focal_taxa_original_review":r["focal_taxa"],
            "n_manipulated_original":r["manipulated_species"],
            "original_R_bin":orig,"numeric_corrected_bin":fixed,
            "temporal_gradient":r["temporal_gradient"],
            "duration_generations":r["duration_generations"],
            "fitted_model_original_code":r["competition_model"],
            "fitted_competition_model_coded":r["competition_model"] not in ("","-"),
            "R_bin_miscoded":"yes" if orig!=fixed and not (orig==">10" and fixed=="11+") else "no",
            "ahif_source_status":"not_independently_screened"
        })
    counts=dict(Counter(r["original_R_bin"] for r in result))
    corrected=dict(Counter(r["numeric_corrected_bin"] for r in result))
    discrepancies=[r for r in result if r["R_bin_miscoded"]=="yes"]
    # Source R code classifies numeric '18' as 'No temporal treatment'
    if len(discrepancies)!=1 or discrepancies[0]["review_record_id"]!="ZOU2023_R0032" \
       or discrepancies[0]["n_manipulated_original"]!="18":
        raise ValueError("Unexpected source code coding discrepancies")
    if [counts[k] for k in PRIOR_R_BINS]!=[46,107,72,43]:
        raise ValueError("Original R reproduction does not match published design counts")
    if [corrected[k] for k in FIXED_BINS]!=[45,107,72,44]:
        raise ValueError("Corrected numeric bins do not match original source")
    modelcoded=[r for r in result if r["fitted_competition_model_coded"]]
    if len(modelcoded)!=13:
        raise ValueError("Unexpected fitted competition model code count")
    corrected_temporal=[r for r in result if r["numeric_corrected_bin"]!="No temporal treatment"]
    detailed={
        "status":"AUTHOR_REVIEW_CODE_REPRODUCED_ONE_BIN_MAPPING_ISSUE_FIXED",
        "source_review_doi":SOURCE_DOI,
        "source_workbook_git_blob":SOURCE_BLOB,
        "source_R_implementation":"LitReview.R lines 45-49",
        "publication_universe_count":len(result),
        "source_AHIF_independently_screened_count":0,
        "published_author_R_bin_counts":counts,
        "corrected_numeric_bin_counts":corrected,
        "numeric_recode_discrepancies":discrepancies,
        "review_coded_models_n":len(modelcoded),
        "review_coded_models_15_or_more_species_n":sum(
            r["numeric_corrected_bin"]=="11+" for r in modelcoded),
        "review_coded_model_studies_per_corrected_bin":{
            k:sum(r["numeric_corrected_bin"]==k for r in modelcoded)
            for k in FIXED_BINS},
        "review_coded_model_studies_per_temporal_gradient":{
            k:sum(r["temporal_gradient"]==k for r in modelcoded)
            for k in ("Yes","No")},
        "review_coded_model_studies_with_no_temporal_manipulation_n":sum(
            r["numeric_corrected_bin"]=="No temporal treatment" for r in modelcoded),
        "review_coded_model_studies_with_temporal_manipulation_n":sum(
            r["numeric_corrected_bin"]!="No temporal treatment" for r in modelcoded),
        "temporal_gradient_yes_n":sum(r["temporal_gradient"]=="Yes" for r in result),
        "duration_generation_counts":dict(Counter(
            r["duration_generations"] for r in result)),
        "top_review_taxa_counts":dict(Counter(
            r["focal_taxa_original_review"] for r in result).most_common(18)),
        "temporal_gradient_by_corrected_bin":{
            k:{flag:sum(r["numeric_corrected_bin"]==k
                       and r["temporal_gradient"]==flag for r in result)
               for flag in ("Yes","No")}
            for k in FIXED_BINS},
        "review_coded_model_candidates":[{
            "id":r["review_record_id"],"year":r["year"],"title":r["title"],
            "taxa":r["focal_taxa_original_review"],
            "model_code":r["fitted_model_original_code"],
            "corrected_bin":r["numeric_corrected_bin"]
        } for r in modelcoded],
        "caveats":[
            "Original review bins were faithfully reproduced and numeric '18' was misrouted to 'No temporal treatment' in published case_when implementation.",
            "The review's 13 fitted-model codes do not exhaust studies measuring interactions or biological mechanisms. '-' does not imply absent interactions.",
            "The 268 studies are pre-2023 reference review observations, not independently audited AHIF eligible experiments.",
            "Review-coded number of manipulated species is not the dimension of any source-verified community abundance matrix.",
            "Temporal gradient is a design covariate, not an ecological effect or meta-analytic finding.",
            "No causal mechanism-to-ranked-structure hypothesis can be established from the review metadata alone."
        ]
    }
    return detailed,result


def priority_screen_candidates(rows:list[dict])->list[dict]:
    """Prioritize published-review title or fitted-model codes, not effect data."""
    observed={r["review_record_id"]:r for r in rows}
    if not AHIF_NEAREST_PRIOR_ART_A.keys()<=observed.keys():
        raise ValueError("Original nearest-prior-art bibliography IDs missing")
    raw_ids=set(AHIF_NEAREST_PRIOR_ART_A)
    model_ids={r["review_record_id"] for r in rows if r["fitted_competition_model_coded"]}
    selected=raw_ids|model_ids
    if len(selected)!=21 or len(model_ids)!=13 or len(raw_ids)!=11:
        raise ValueError(f"Expected 21 union of original 11+13, got {len(selected)}")
    ordered=sorted(selected,key=lambda id_:(id_ not in raw_ids,int(id_.split("R")[-1])))
    priority=[]
    for id_ in ordered:
        r=observed[id_]
        priority.append({
            "review_record_id":id_,
            "priority_tier":"A_nearest_prior_art" if id_ in raw_ids else "B_review_fitted_model",
            "why_prior_art_matters":AHIF_NEAREST_PRIOR_ART_A.get(
                id_,"Source review codes a fitted competition/ecological process model"),
            "source_review_model_code":r["fitted_model_original_code"],
            "author_review_year":r["year"],
            "author_review_title":r["title"],
            "study_original_manipulated_species_code":r["n_manipulated_original"],
            "study_corrected_size_bin":r["numeric_corrected_bin"],
            "independent_AHIF_methods_screened":"no",
            "complete_community_data_source_verified":"no",
            "effect_size_extracted":"no",
            "source_review_doi":SOURCE_DOI
        })
    return priority


def write(path:Path, rows, fieldnames):
    with path.open("w",encoding="utf-8",newline="") as handle:
        wr=csv.DictWriter(handle,fieldnames=fieldnames)
        wr.writeheader();wr.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--review-screening-csv",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    info,rows=profile(source_rows(a.review_screening_csv))
    a.output_dir.mkdir(parents=True,exist_ok=True)
    (a.output_dir/"review268_design_audit.json").write_text(
        json.dumps(info,indent=2,ensure_ascii=False)+"\n")
    write(a.output_dir/"review268_design_corrected_bins.csv",
          rows,list(rows[0].keys()))
    priority=priority_screen_candidates(rows)
    write(a.output_dir/"review21_nearest_prior_art_screening.csv",
          priority,list(priority[0].keys()))
    print("AHIF_REVIEW_PROFILE",json.dumps({
        "n":info["publication_universe_count"],
        "original":info["published_author_R_bin_counts"],
        "corrected":info["corrected_numeric_bin_counts"],
        "bin_errors":len(info["numeric_recode_discrepancies"]),
        "fitted_competition_model_codes":info["review_coded_models_n"],
        "prioritized_bibliographic_candidates_n":len(priority),
        "uncertainty":"NONE of 268 independently AHIF source-screened"
    },ensure_ascii=False))


if __name__=="__main__":main()
