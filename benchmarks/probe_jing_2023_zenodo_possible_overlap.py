#!/usr/bin/env python3
"""Public Zenodo original workbook inventory, possible duplicate of Jing 2024.

DO NOT assume that Zenodo 7856898 and Dryad jq2bvq8jj are independent
experiments.  This step only exposes original row structure and source headers.
"""
from __future__ import annotations
import argparse
import hashlib,json
from pathlib import Path

from benchmarks.probe_jing_2024_source_xlsx import xlsx_rows

MD5_WORKBOOK="5afb3a33ad49573aa7e255d2541ab0ed"
MD5_README="ae5c5b22981990bdf3a701cd5c7d7393"


def audit(original_workbook:Path,readme:Path):
    raw=original_workbook.read_bytes()
    doc=readme.read_bytes()
    if hashlib.md5(raw).hexdigest()!=MD5_WORKBOOK:
        raise ValueError("Zenodo 7856898 original XLSX MD5 mismatch")
    if hashlib.md5(doc).hexdigest()!=MD5_README:
        raise ValueError("Zenodo 7856898 README MD5 mismatch")
    readme_content=doc.decode("utf-8-sig")
    tabs=xlsx_rows(original_workbook)
    sheets=[]
    for name,v in tabs.items():
        rows=v["rows"]
        if not rows:continue
        headers=[str(x) if x is not None else "" for x in rows[0]]
        samples=[r[:20] for r in rows[1:4]]
        sheets.append({
            "sheet":name,"data_rows":len(rows)-1,
            "n_cols":max(map(len,rows)),
            "headers":headers[:60],
            "sample_rows_first20cols":samples,
            "last_6_rows_first20cols":[row[:20] for row in rows[-6:]],
            "water_soil_row_counts":{
                key:sum(str(row[0]).strip()+"|"+str(row[1]).strip()==key for row in rows)
                for key in sorted(set(str(row[0]).strip()+"|"+str(row[1]).strip()
                                      for row in rows[1:] if len(row)>=2))
            },
            "n_formula_cells":v["formula_cells"],
            "n_complete_nonblank_rows":sum(all(y is not None and str(y).strip()!="" for y in row)
                                         for row in rows[1:])
        })
    return {
        "status":"SCHEMA_ONLY_POT_UNIT_NOT_YET_AUDITED",
        "source_doi":"10.5281/zenodo.7856898",
        "source_workbook_md5":MD5_WORKBOOK,
        "source_workbook_sha256":hashlib.sha256(raw).hexdigest(),
        "source_readme_md5":MD5_README,
        "source_readme_verbatim":readme_content,
        "n_sheets":len(sheets),
        "sheets":sheets,
        "important_caveats":[
            "May overlap Jing et al 2024 Dryad source from the same research group; do not count as independent.",
            "This is raw archive schema, NOT an effect size or source-replicated biological result.",
            "Column structure needs actual pot/repeated-measures audit before a competition-effect estimator",
            "No redistribution of source XLSX to KNK"
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--workbook",type=Path,required=True)
    p.add_argument("--readme",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    j=audit(args.workbook,args.readme)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(j,ensure_ascii=False,indent=2)+"\n")
    print("JING_ZENODO_2023_AUDIT_START")
    print("SOURCE_META",json.dumps({k:v for k,v in j.items() if k not in ("sheets","source_readme_verbatim")}))
    print("SOURCE_README",j["source_readme_verbatim"][:10000])
    for item in j["sheets"]:
        print("SOURCE_SHEET",json.dumps(item,ensure_ascii=False)[:10000])
    print("JING_ZENODO_2023_AUDIT_DONE")


if __name__=="__main__":main()
