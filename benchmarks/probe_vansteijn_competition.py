#!/usr/bin/env python3
"""Source-locked schema and methods inspection only, no ecological scoring.

van Steijn et al. (2025), Ecology doi:10.1002/ecy.4528
Zenodo v2 doi:10.5281/zenodo.10512662
"""
from __future__ import annotations
import argparse, csv, hashlib, json, re
from collections import Counter
from pathlib import Path

EXPECTED_MD5={
    "alone.csv":"0df05c9e417947e7aaecb84913184fc5",
    "simultaneous_arrival.csv":"3a6742d98e4311bf8cd64dc06ddf433b",
    "second_arrival.csv":"17f876d2ba031248f6813bddde9a33be",
    "traits.csv":"e7b1ffaeff71b05d79c1e9ae158aa55f",
    "R_code.Rmd":"16f92a0dde59ce52e30a671dd31859ef",
}
TEXT_RE=re.compile(r"RII|competitive|second_arrival|simultaneous_arrival|alone.csv|trait|pair|late|species|plant", re.I)

def audit(directory:Path)->dict:
    result={"status":"SOURCE_SCHEMA_ONLY_NOT_EFFECT_SCORING",
            "dataset_doi":"10.5281/zenodo.10512662",
            "source_article_doi":"10.1002/ecy.4528",
            "files":{}}
    for filename,md5 in EXPECTED_MD5.items():
        path=directory/filename
        b=path.read_bytes()
        got=hashlib.md5(b).hexdigest()
        if got!=md5:
            raise ValueError(f"Source MD5 mismatch {filename}: {got} != {md5}")
        info={"filename":filename,"bytes":len(b),"md5":got,
              "sha256":hashlib.sha256(b).hexdigest()}
        if filename.endswith(".csv"):
            src=b.decode("utf-8-sig",errors="strict")
            sniff=csv.Sniffer().sniff(src[:5000],delimiters=",;\t")
            reader=csv.DictReader(src.splitlines(),delimiter=sniff.delimiter)
            rows=list(reader)
            hdr=reader.fieldnames or []
            info.update({
                "rows":len(rows),"column_names":hdr,"delimiter":sniff.delimiter,
                "preview_first_3_rows":rows[:3],
                "field_cardinality_and_values":{
                    k:{"n_levels":len(Counter(row.get(k,"") for row in rows)),
                       "top_levels":Counter(row.get(k,"") for row in rows).most_common(10)}
                    for k in hdr if re.search(r"species|neighbor|treat|group|pair|arriv|replicate|name|ID|sp[12]",k,re.I)
                },
                "all_na_counts":{k:sum(str(row.get(k,"")).strip().lower() in ("","na","nan") for row in rows) for k in hdr},
            })
        else:
            text=b.decode("utf-8-sig",errors="replace")
            matches=[{"line":i,"text":line[:260]}
                     for i,line in enumerate(text.splitlines(),1)
                     if TEXT_RE.search(line)]
            info.update({
                "n_lines":len(text.splitlines()),
                "relevant_r_code_lines":matches[:135]
            })
        result["files"][filename]=info
    return result

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument("--source-dir",type=Path,required=True)
    a.add_argument("--out",type=Path,required=True)
    args=a.parse_args()
    result=audit(args.source_dir)
    args.out.parent.mkdir(exist_ok=True,parents=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("VAN_STEIJN_SCHEMA_START")
    for k,v in result["files"].items():
        print("FILE_SCHEMA",k,json.dumps({p:q for p,q in v.items() if p!="relevant_r_code_lines"},ensure_ascii=False)[:8000])
    print("SOURCE_R_LINES",json.dumps(result["files"]["R_code.Rmd"]["relevant_r_code_lines"],ensure_ascii=False)[:15500])
    print("VAN_STEIJN_SCHEMA_DONE")

if __name__=="__main__":
    main()
