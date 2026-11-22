#!/usr/bin/env python3
"""POEM raw ZIP experimental-unit audit; does not score treatment effects.

Publication: Alonso-Crespo et al. (2025), Journal of Vegetation Science.
Zenodo v1: https://doi.org/10.5281/zenodo.10119982
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
from io import BytesIO
import json
from pathlib import Path
import re
import zipfile

from benchmarks.probe_poem_archive import EXPECTED_MD5

BIOMASS_PATTERN = re.compile(r"/Data/POEM(2020|2021)/Data_Biomass_POEM(2020|2021)_(\d{4}-\d{2})\.xlsx$")
CODE_PATTERN = re.compile(r"/R/Analysis_POEM_data\.Rmd$")
SOURCE_SOWN = (
    "Agrostis capillaris", "Anthoxanthum odoratum", "Bromus hordeaceus",
    "Dianthus deltoides", "Festuca ovina", "Jasione montana",
    "Lathyrus pratensis", "Lotus corniculatus", "Pilosella officinarum",
    "Pimpinella saxifraga", "Potentilla argentea", "Silene vulgaris",
    "Trifolium arvense", "Trifolium campestre"
)


def extract_book_rows(zipf, item):
    import openpyxl
    wb = openpyxl.load_workbook(BytesIO(zipf.read(item)), read_only=True, data_only=True)
    if len(wb.worksheets) != 1:
        raise ValueError(f"Unexpected number of sheets: {item}: {[s.title for s in wb]}")
    data = wb.active.iter_rows(values_only=True)
    headers = list(next(data))
    expected = ("Year", "Harvest", "Plot", "Quadrat", "Species", "SDW_g")
    if tuple(headers[:6]) != expected:
        raise ValueError(f"Unexpected columns {headers[:6]} for {item}")
    rows = []
    for ix, raw in enumerate(data, 2):
        if all(v is None for v in raw):
            continue
        if len(raw) < 6:
            raise ValueError(f"Short row {ix}: {item}")
        year, harvest, plot, quad, taxon, grams = raw[:6]
        rows.append({
            "year": str(year).strip(),
            "harvest": str(harvest).strip(),
            "plot": str(plot).strip(),
            "quadrat": str(quad).strip(),
            "taxon": str(taxon).strip(),
            "mass": grams,
            "row_number": ix
        })
    wb.close()
    return rows


def audit_zip(zippath: Path) -> dict:
    md5 = hashlib.md5(zippath.read_bytes()).hexdigest()
    if md5 != EXPECTED_MD5:
        raise ValueError("Source ZIP MD5 is not the published Zenodo v1 checksum")
    source = {"source_md5": md5, "source_sha256": hashlib.sha256(zippath.read_bytes()).hexdigest()}
    with zipfile.ZipFile(zippath) as zf:
        books = []
        by_harvest = []
        all_species = set()
        for name in zf.namelist():
            match = BIOMASS_PATTERN.search(name)
            if not match:
                continue
            c, c2, harvest_name = match.groups()
            if c != c2:
                raise ValueError(f"Cohort mismatch in filename {name}")
            recs = extract_book_rows(zf, name)
            if not recs:
                raise ValueError(f"Empty data file {name}")
            unique_keys = {(r["plot"], r["quadrat"], r["taxon"]) for r in recs}
            duplicated = len(recs) - len(unique_keys)
            cohorts = Counter(r["year"] for r in recs)
            harvest_labels = Counter(r["harvest"] for r in recs)
            vals = [r["mass"] for r in recs]
            numeric = [v for v in vals if isinstance(v, (int, float))]
            bad = sum(not isinstance(v, (int, float)) for v in vals)
            neg = sum(float(v) < 0 for v in numeric)
            plot_keys = sorted({r["plot"] for r in recs})
            quad_keys = sorted({(r["plot"], r["quadrat"]) for r in recs})
            species = sorted({r["taxon"] for r in recs})
            all_species.update(species)
            # Derive candidate treatment/replicate IDs, but NOT yet verified as randomization blocks.
            rx = re.compile(r"^(\d+)-([SFGLB])([1-9])$")
            parsed = [rx.match(p) for p in plot_keys]
            unmatched = [p for p, m in zip(plot_keys, parsed) if m is None]
            treatment_counts = Counter(m.group(2) for m in parsed if m)
            replicate_counts = Counter(m.group(3) for m in parsed if m)
            cells_per_quad = Counter((r["plot"], r["quadrat"]) for r in recs)
            replicates_per_plot = Counter(p for p, q in quad_keys)
            sample = {
                "member_name": name,
                "cohort": c, "expected_harvest": harvest_name,
                "cohort_values": dict(cohorts),
                "harvest_values": dict(harvest_labels),
                "rows": len(recs),
                "n_unique_plots": len(plot_keys),
                "n_unique_plot_quadrats": len(quad_keys),
                "n_observed_species": len(species),
                "n_sown_species_observed": len(set(species) & set(SOURCE_SOWN)),
                "n_unsown_species_observed": len(set(species) - set(SOURCE_SOWN)),
                "n_explicit_zero_rows": sum(v == 0 for v in numeric),
                "n_missing_or_non_numeric_mass": bad,
                "n_negative_mass": neg,
                "n_duplicate_plot_quadrat_species": duplicated,
                "plot_examples": plot_keys[:10],
                "unparsed_plot_examples": unmatched[:15],
                "parsed_treatment_plot_counts": dict(treatment_counts),
                "parsed_replicate_plot_counts": dict(replicate_counts),
                "quadrat_count_per_plot_distribution": dict(sorted(Counter(replicates_per_plot.values()).items())),
                "species_rows_per_plot_quadrat_distribution": dict(sorted(Counter(cells_per_quad.values()).items())),
                "all_plot_ids": plot_keys,
                "all_species": species
            }
            by_harvest.append(sample)
            books.append(name)
        if len(books) != 6:
            raise ValueError(f"Expected 6 source biomass Excel workbooks, found {len(books)}")
        source["sown_species_seen"] = sorted(set(SOURCE_SOWN) & all_species)
        source["sown_species_unseen"] = sorted(set(SOURCE_SOWN) - all_species)
        source["n_distinct_species_across_harvests"] = len(all_species)
        source["harvests"] = sorted(by_harvest, key=lambda x: (x["cohort"], x["expected_harvest"]))
        # Source paper conventions can vary; extract context without altering source.
        rmds = [n for n in zf.namelist() if CODE_PATTERN.search(n)]
        if len(rmds) == 1:
            rmd = zf.read(rmds[0]).decode("utf-8-sig", errors="replace")
            lines = rmd.splitlines()
            patterns = re.compile(
                r"POEM_shoot_biomass_data|sown|unsown|quadrat|replicate|"
                r"factor\(|contrasts|permanova|adonis|Arrival|_PFG", re.I)
            refs = [(ix, l[:240]) for ix, l in enumerate(lines, 1) if patterns.search(l)]
            source["rmd_key_lines"] = [{"line": i, "code": line} for i, line in refs[:120]]
            source["rmd_n_lines"] = len(lines)
        readmes = [n for n in zf.namelist() if n.endswith("Data/README.txt")]
        if len(readmes) == 1:
            readme = zf.read(readmes[0]).decode("utf-8-sig", errors="replace")
            source["data_readme_text"] = readme
    source["status"] = "RAW_SOURCE_DESIGN_METADATA_ONLY_NOT_EFFECT_SCORING"
    source["claim_boundaries"] = [
        "Missing species entries must not be zero-filled without source convention audit",
        "Plot replicate coding does not necessarily imply restricted treatment randomization",
        "Cohort/year repeated observations do not increase number of independent plots",
        "This pilot remains separate from frozen KNK v28-v38 and Draft PR #3",
    ]
    return source


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--zip", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    report = audit_zip(args.zip)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n")
    print("POEM_DESIGN_AUDIT_START")
    print("ARCHIVE", report["source_sha256"], "species_total", report["n_distinct_species_across_harvests"])
    for r in report["harvests"]:
        # Compact audit for tool-visible GitHub job logs.
        print("HARVEST", json.dumps({k:v for k,v in r.items()
                                   if k not in ("all_species","all_plot_ids","member_name")},ensure_ascii=False)[:4000])
    print("README_EXTRACT", report.get("data_readme_text","")[:8000])
    print("RMD_REFERENCES", json.dumps(report.get("rmd_key_lines",[])[:55],ensure_ascii=False)[:12000])
    print("POEM_DESIGN_AUDIT_DONE",args.output)


if __name__ == "__main__":
    main()
