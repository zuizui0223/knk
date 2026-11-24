#!/usr/bin/env python3
"""Public Zou & Rudolf 2023 review workbook source schema audit, no raw republication.

Source https://github.com/hengxingzou/Zou2023TrendsEcolEvol/
Paper DOI 10.1016/j.tree.2023.08.001. The reviewed studies are NOT
automatically AHIF-DB source-audited records/effects.
"""
from __future__ import annotations
import argparse
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


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--xlsx",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    data=audit(a.xlsx)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(data,indent=2,ensure_ascii=False)+"\n")
    print("ZOU_SOURCE_MANIFEST",json.dumps(data,ensure_ascii=False))


if __name__=="__main__":main()
