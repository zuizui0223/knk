#!/usr/bin/env python3
"""Audit authentic van Steijn et al 2026 field data JOIN and plot-year grain.

Uses original code's filter Site != GAL AND Type == Transect.
Only data quality/provenance, no treatment-effect analysis.
Source files are hash-checked by probe_ahif_vansteijn2026_field.py.
"""
from __future__ import annotations
from collections import Counter,defaultdict
import csv,json
from pathlib import Path

def read(path:Path):
    with path.open(encoding="utf-8-sig",newline="") as f:
        lines=csv.DictReader(f,delimiter=";")
        if not lines.fieldnames:raise ValueError("Missing CSV header")
        return list(lines)

def audit(directory:Path):
    plot=read(directory/"plotdata.csv")
    records=read(directory/"Pinpoint_all_years.csv")
    sp=read(directory/"speciesinfo.csv")
    if len(plot)!=40 or len(records)!=7070 or len(sp)!=20:
        raise ValueError("Original 40/7070/20 source counts unexpectedly changed")
    if tuple(plot[0])!=("Trmt","Site","Block","Plot"):
        raise ValueError("Original plot design fields changed")
    need=("Year","Site","Block","Plot","OTC","Type","Trmt","Species","Abundance")
    if tuple(records[0])!=need:
        raise ValueError("Original pinpoint columns differ")
    design={}
    for p in plot:
        key=(p["Site"],p["Block"],p["Plot"])
        if key in design:raise ValueError("Nonunique 40 source plot identities")
        design[key]=p["Trmt"]
    targets={r["Species"] for r in sp}
    if len(targets)!=20:raise ValueError("Missing target source species")
    selected=[r for r in records if r["Site"]!="GAL" and r["Type"]=="Transect"]
    tab=Counter((r["Site"],r["Year"]) for r in selected)
    groups=defaultdict(list)
    species_per_plot=Counter()
    for r in selected:
        name=(r["Site"],r["Block"],r["Plot"],r["Year"])
        groups[name].append(r)
        species_per_plot[name,r["Species"]]+=1
        if r["Plot"] in ("","NA"):raise ValueError("No independent plot identifier in sampled transect")
        try: value=float(r["Abundance"])
        except ValueError as e:raise ValueError("Non-numeric source pinpoint abundance") from e
        if value<0:raise ValueError("Negative pinpoint count")
    out={
        "status":"ORIGINAL_TRANSECT_PLOT_JOIN_AUDIT_NOT_AN_EFFECT_SCORE",
        "source_filter":"Site != GAL & Type == Transect (original Abundance_code.Rmd line 37)",
        "source_original_rows":len(records),"filter_rows":len(selected),
        "plot_assignment_rows":len(plot),
        "original_species_info_n":len(targets),
        "site_year_rows":{f"{site}|{year}":n for (site,year),n in sorted(tab.items())},
        "site_year_independent_plot_count":{
            f"{site}|{year}":len({(block,plot) for s,block,plot,y in groups if s==site and y==year})
            for site,year in sorted(tab)},
        "site_year_block_distribution":{
            f"{site}|{year}":dict(Counter(block for s,block,plot,y in groups if s==site and y==year))
            for site,year in sorted(tab)},
        "site_year_treatment_distribution":{
            f"{site}|{year}":dict(Counter(
                groups[k][0]["Trmt"] for k in groups if k[0]==site and k[3]==year))
            for site,year in sorted(tab)},
        "unique_plot_years":len(groups),
        "site_block_plot_keys_in_survey_not_in_assignment":[
            list(k[:3]) for k in sorted(groups) if k[:3] not in design][:70],
        "assigned_plot_keys_not_in_any_survey":[list(k) for k in design if
            not any(g[:3]==k for g in groups)],
        "source_trmt_disagree_with_assignment":[
            {"plot_year":list(k),"assigned":design[k[:3]],"observed":list({r["Trmt"] for r in v})}
            for k,v in groups.items() if k[:3] in design and
            ({r["Trmt"] for r in v}!={design[k[:3]]})][:70],
        "n_duplicate_same_species_per_plot_year":sum(n-1 for n in species_per_plot.values() if n>1),
        "n_target_taxon_rows":sum(r["Species"] in targets for r in selected),
        "n_non_target_taxon_rows":sum(r["Species"] not in targets for r in selected),
        "n_selected_species":len({r["Species"] for r in selected}),
        "n_zeros_explicit":sum(float(r["Abundance"])==0 for r in selected),
        "source_abundance_min":min(float(r["Abundance"]) for r in selected),
        "source_abundance_max":max(float(r["Abundance"]) for r in selected),
        "note":"Counts are recorded plant pinpoint hits; point/hit counts, not biomass or physical individuals. Plant individuals/annual repeated rows NOT independent."
    }
    return out

def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    r=audit(a.source_dir)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
    print("FIELD_PLOT_GRAIN_AUDIT",json.dumps(r,ensure_ascii=False))
if __name__=="__main__":main()
