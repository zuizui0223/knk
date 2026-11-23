#!/usr/bin/env python3
"""Jing et al. 2024 Dryad source-matrix inventory, standard-library XLSX reader.

No outcome scoring or cell-level source redistribution. Outputs an inspectable
schema and sample of the ORIGINAL workbook as proof of empirical grain.
Original DOI https://doi.org/10.5061/dryad.jq2bvq8jj
The 11 tabs mix raw and author-aggregated derivatives: do NOT pseudo-replicate.
"""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import posixpath
import re
from xml.etree import ElementTree as ET
import zipfile

MAIN="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
OFF="http://schemas.openxmlformats.org/officeDocument/2006/relationships"
PACK="http://schemas.openxmlformats.org/package/2006/relationships"
COL_RE=re.compile(r"^([A-Z]+)(\d+)$")


def excel_col_to_index(cell_address:str)->int:
    m=COL_RE.match(cell_address)
    if not m:raise ValueError("Invalid Excel cell address: "+cell_address)
    n=0
    for letter in m.group(1):n=n*26+ord(letter)-ord('A')+1
    return n-1


def as_number(s:str):
    try:
        x=float(s)
        return int(x) if x.is_integer() else x
    except (TypeError,ValueError,OverflowError):
        return s


def xlsx_rows(path:Path):
    with zipfile.ZipFile(path) as zin:
        if len(zin.namelist())>2000:raise ValueError("Oversize unexpected XLSX archive")
        wb=ET.fromstring(zin.read("xl/workbook.xml"))
        relationships=ET.fromstring(zin.read("xl/_rels/workbook.xml.rels"))
        relmap={x.attrib["Id"]:x.attrib["Target"] for x in relationships}
        strings=[]
        if "xl/sharedStrings.xml" in zin.namelist():
            sr=ET.fromstring(zin.read("xl/sharedStrings.xml"))
            for node in sr.findall(f"{{{MAIN}}}si"):
                strings.append("".join((e.text or "") for e in node.iter(f"{{{MAIN}}}t")))
        results={}
        for sh in wb.findall(f"{{{MAIN}}}sheets/{{{MAIN}}}sheet"):
            rid=sh.attrib.get(f"{{{OFF}}}id")
            target=relmap[rid]
            if target.startswith("/"):
                member=target.lstrip("/")
            else:
                member=posixpath.normpath(posixpath.join("xl",target))
            if member not in zin.namelist():
                raise ValueError(f"Missing worksheet entry {member}")
            root=ET.fromstring(zin.read(member))
            rows=[]
            formula_cells=0
            for row in root.findall(f".//{{{MAIN}}}sheetData/{{{MAIN}}}row"):
                values={}
                for c in row.findall(f"{{{MAIN}}}c"):
                    a=c.attrib.get("r","")
                    j=excel_col_to_index(a)
                    if c.find(f"{{{MAIN}}}f") is not None:
                        formula_cells+=1
                    celltype=c.attrib.get("t")
                    if celltype=="inlineStr":
                        content="".join(x.text or "" for x in c.iter(f"{{{MAIN}}}t"))
                    else:
                        el=c.find(f"{{{MAIN}}}v")
                        val=el.text if el is not None else ""
                        if celltype=="s":
                            content=strings[int(val)] if val else ""
                        elif celltype=="b":
                            content=(val=="1")
                        elif celltype in ("str","e"):
                            content=val
                        else:
                            content=as_number(val) if val!="" else None
                    values[j]=content
                if values:
                    size=max(values)+1
                    rows.append([values.get(j) for j in range(size)])
            results[sh.attrib["name"]]={"rows":rows,"formula_cells":formula_cells}
        return results


def audit(path:Path, expected_sha256:str|None=None):
    raw=path.read_bytes()
    sha=hashlib.sha256(raw).hexdigest()
    if expected_sha256 and sha!=expected_sha256:
        raise ValueError("Source checksum mismatch")
    if len(raw)>2_000_000:
        raise ValueError("Original source unexpectedly large")
    tabs=xlsx_rows(path)
    sheets=[]
    for name,info in tabs.items():
        data=info["rows"]
        if not data:continue
        headers=[str(x).strip() if x is not None else "" for x in data[0]]
        ncols=max(map(len,data))
        preview=[r[:min(ncols,16)] for r in data[1:4]]
        candidate_cols={}
        for ci,title in enumerate(headers):
            if re.search(r"species|arrange|inocul|drought|replic|pot|soil|mix|RY|biomass",title,re.I):
                values=[str(r[ci]) if ci<len(r) and r[ci] is not None else "(missing)"
                        for r in data[1:]]
                levels=Counter(values)
                candidate_cols[title[:80]]={
                    "distinct_values":len(levels),
                    "top_10":levels.most_common(10),
                    "n_missing":levels.get("(missing)",0)
                }
        sheets.append({
            "sheet_name":name,"rows_excluding_header":len(data)-1,
            "total_columns_max":ncols,"first_row_headers":headers[:60],
            "first_three_data_rows_truncated":preview,
            "category_cardinalities":candidate_cols,
            "cached_formula_cell_count":info["formula_cells"]
        })
    return {
        "status":"SOURCE_SCHEMA_AND_GRAIN_ONLY_NO_BIOLOGICAL_SCORING",
        "doi":"10.5061/dryad.jq2bvq8jj",
        "original_filename":path.name,
        "source_sha256":sha,"original_source_bytes":len(raw),
        "n_worksheets":len(sheets),
        "worksheets":sheets,
        "notes":[
            "Author workbook includes per-replicate and derived summary tabs, not all rows are independent pots",
            "Sheet-specific plot/replicate identifiers and raw biomass must be inspected before calculating an effect",
            "Exactly one published original experiment regardless of worksheet count",
            "Original source workbook never committed to KNK",
            "No edits to frozen KNK v28-v38 or manuscript"
        ]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--xlsx",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--expected-sha256")
    args=p.parse_args()
    report=audit(args.xlsx,args.expected_sha256)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("JING_SOURCE_SCHEMA_START")
    print("SOURCE",json.dumps({k:v for k,v in report.items() if k!="worksheets"}))
    for sh in report["worksheets"]:
        print("JING_SHEET",json.dumps(sh,ensure_ascii=False)[:8500])
    print("JING_SOURCE_SCHEMA_END",args.output)


if __name__=="__main__":main()
