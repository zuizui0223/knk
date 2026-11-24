#!/usr/bin/env python3
"""Public Zou & Rudolf 2023 review workbook source schema audit, no raw republication.

Source https://github.com/hengxingzou/Zou2023TrendsEcolEvol/
Paper DOI 10.1016/j.tree.2023.08.001. The reviewed studies are NOT
automatically AHIF-DB source-audited records/effects.
"""
from __future__ import annotations
import argparse
import csv
from collections import Counter
from pathlib import Path
import hashlib
import json
import re
import zipfile
import xml.etree.ElementTree as ET

SRC_GIT_BLOB="be2777d4e4bb1454db393b94a2118ada79118dfa"
S="http://schemas.openxmlformats.org/spreadsheetml/2006/main"
R="http://schemas.openxmlformats.org/officeDocument/2006/relationships"


def blobhash(b:bytes):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()


def sheet_rows(zf, member, strings):
    root=ET.fromstring(zf.read(member))
    output=[]
    for row in root.findall(f".//{{{S}}}sheetData/{{{S}}}row"):
        cells={}
        for cell in row.findall(f"{{{S}}}c"):
            ref=cell.attrib["r"]
            match=re.match(r"^([A-Z]+)\d+$",ref)
            if not match:raise ValueError("Unknown cell ref "+ref)
            ix=0
            for ch in match.group(1):ix=26*ix+ord(ch)-64
            ix-=1
            ty=cell.attrib.get("t")
            value=cell.find(f"{{{S}}}v")
            if ty=="inlineStr":
                v="".join((x.text or "") for x in cell.iter(f"{{{S}}}t"))
            elif value is None:v=None
            elif ty=="s":
                v=strings[int(value.text)]
            else:v=value.text
            cells[ix]=v
        if cells:
            output.append([cells.get(i) for i in range(max(cells)+1)])
    return output


def audit(path:Path):
    raw=path.read_bytes()
    sha=blobhash(raw)
    if sha!=SRC_GIT_BLOB:raise ValueError("External systematic review source blob changed")
    with zipfile.ZipFile(path) as z:
        strings=[]
        if "xl/sharedStrings.xml" in z.namelist():
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            for el in root.findall(f"{{{S}}}si"):
                strings.append("".join((x.text or "") for x in el.iter(f"{{{S}}}t")))
        book=ET.fromstring(z.read("xl/workbook.xml"))
        rels=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap={r.attrib["Id"]:r.attrib["Target"] for r in rels}
        results=[]
        for sheet in book.findall(f".//{{{S}}}sheets/{{{S}}}sheet"):
            name=sheet.attrib["name"]
            rel=sheet.attrib[f"{{{R}}}id"]
            t=relmap[rel]
            member=t.lstrip("/") if t.startswith("/") else "xl/"+t.lstrip("/")
            records=sheet_rows(z,member,strings)
            headers=[str(x) if x is not None else "" for x in records[0]]
            results.append({
                "sheet":name,"number_nonblank_rows":len(records),
                "n_data_rows_naive":len(records)-1,
                "header_fields":headers,
                "first_two_nonheader_rows":[x[:25] for x in records[1:3]],
                "n_data_cols":max(map(len,records)),
                "source_data_year_not_updated_after_2023":True
            })
    return {"status":"PUB_REVIEW_SCHEMA_ONLY_NOT_AHIF_CURATION",
            "review_doi":"10.1016/j.tree.2023.08.001",
            "source_workbook_git_blob":sha,
            "source_sha256":hashlib.sha256(raw).hexdigest(),
            "source_bytes":len(raw),"sheets":results,
            "claim_boundary":"External 2023 review studies not independently AHIF-screened or effects scored."}


STAGING_COLUMNS=(
 "review_record_id","source_excel_row","source_review_doi","source_git_blob",
 "authors","year","journal","title","focal_taxa","secondary_taxa",
 "manipulated_species","temporal_gradient","generation_time",
 "duration_generations","competition_model","original_review_inclusion",
 "ahif_independent_screen_status"
)


def stage_review_candidates(original:Path, output:Path):
    """Extract 268 bibliographic screening CANDIDATES, never empirical effect sizes.

    Original review experimental data classifications are author-owned and must
    be attributed. AHIF status for all 268 begins as NOT independently screened.
    No unlicensed full text, article bodies or source workbook are reproduced.
    """
    raw=original.read_bytes()
    if blobhash(raw)!=SRC_GIT_BLOB:
        raise ValueError("Review workbook source changed")
    with zipfile.ZipFile(original) as zf:
        strings=[]
        if "xl/sharedStrings.xml" in zf.namelist():
            root=ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for element in root.findall(f"{{{S}}}si"):
                strings.append("".join((x.text or "") for x in element.iter(f"{{{S}}}t")))
        book=ET.fromstring(zf.read("xl/workbook.xml"))
        rels=ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        mapping={node.attrib["Id"]:node.attrib["Target"] for node in rels}
        sheets=book.findall(f".//{{{S}}}sheets/{{{S}}}sheet")
        found=[q for q in sheets if q.attrib.get("name")=="Experiments for Analysis"]
        if len(found)!=1:
            raise ValueError("Published study screen sheet missing")
        rel=found[0].attrib[f"{{{R}}}id"]
        target=mapping[rel]
        member=target.lstrip("/") if target.startswith("/") else "xl/"+target.lstrip("/")
        allrows=sheet_rows(zf,member,strings)
    if len(allrows)!=269:
        raise ValueError(f"Source review expected 268 study candidate rows, got {len(allrows)-1}")
    headers=allrows[0]
    expected=["Authors","Year","Journal","Focal Taxa","Secondary Taxa",
              "Generation Time","# Manipulated Species","Temporal Gradient?",
              "Max duration","Duration (Generations)","Fitted Competition Model",
              "Include in analysis?","Notes","Title"]
    if headers!=expected:
        raise ValueError("Unexpected source review sheet header")
    mapped=[]
    for ix,vals in enumerate(allrows[1:],2):
        row=vals+[None]*(len(headers)-len(vals))
        fields=dict(zip(headers,row))
        citation=fields["Title"]
        authors=fields["Authors"]
        if not citation or not authors:
            raise ValueError(f"Blank review bibliographic row {ix}")
        year=fields["Year"]
        if not str(year or "").isdigit():
            raise ValueError(f"Bad year at source row {ix}")
        out={
            "review_record_id":f"ZOU2023_R{ix-1:04d}",
            "source_excel_row":str(ix),
            "source_review_doi":"10.1016/j.tree.2023.08.001",
            "source_git_blob":SRC_GIT_BLOB,
            "authors":authors or "",
            "year":str(year),
            "journal":fields["Journal"] or "",
            "title":citation,
            "focal_taxa":fields["Focal Taxa"] or "",
            "secondary_taxa":fields["Secondary Taxa"] or "",
            "manipulated_species":fields["# Manipulated Species"] or "",
            "temporal_gradient":fields["Temporal Gradient?"] or "",
            "generation_time":fields["Generation Time"] or "",
            "duration_generations":fields["Duration (Generations)"] or "",
            "competition_model":fields["Fitted Competition Model"] or "",
            "original_review_inclusion":fields["Include in analysis?"] or "",
            "ahif_independent_screen_status":"not_independently_screened"
        }
        mapped.append(out)
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open("w",newline="",encoding="utf-8") as file:
        writer=csv.DictWriter(file,fieldnames=STAGING_COLUMNS)
        writer.writeheader()
        writer.writerows(mapped)
    return {
        "source_review_rows":len(mapped),
        "candidate_status":"NOT independently screened, no derived effect sizes",
        "independent_AHIF_source_audits_from_review":0,
        "source_blob":SRC_GIT_BLOB,
        "staging_path":str(output)
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--screening-csv",type=Path,
                    help="Write candidate-only bibliographic metadata (268 unreviewed studies)")
    a=ap.parse_args()
    data=audit(a.xlsx)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("ZOU_SOURCE_MANIFEST",json.dumps(data,ensure_ascii=False))
    if a.screening_csv:
        c=stage_review_candidates(a.xlsx,a.screening_csv)
        print("ZOU_CANDIDATE_STAGING",json.dumps(c,ensure_ascii=False))


if __name__=="__main__":main()
