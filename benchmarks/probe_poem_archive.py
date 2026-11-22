#!/usr/bin/env python3
"""Inspect the published POEM Zenodo ZIP without redistributing raw data.

Source DOI: https://doi.org/10.5281/zenodo.10119982
Expected release (v1.0) ZIP MD5: 3faa5f7d1619a3beae5941ba763b7a26
Print file names, shapes, field names and small snippets of analysis code.
No ecological results are scored by this probe.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import re
import zipfile

EXPECTED_MD5 = "3faa5f7d1619a3beae5941ba763b7a26"
DATA_SUFFIXES = (".csv", ".tsv", ".txt", ".xlsx", ".xls", ".rds", ".rdata", ".rda")
CODE_SUFFIXES = (".r", ".rmd")
INTERESTING_PAT = re.compile(
    r"read[_ .]*(csv|excel|table|xlsx|rds)|load\(|shoot|biomass|richness|plot_id|"
    r"minirhizotron|species|root|treatment|poem202[01]", re.IGNORECASE
)


def preview_text(data: bytes, file_name: str) -> dict:
    try:
        text = data.decode("utf-8-sig")
        encoding = "utf-8-sig"
    except UnicodeDecodeError:
        text = data.decode("latin-1")
        encoding = "latin-1"
    sample = text[:15000]
    if not sample.strip():
        return {"encoding": encoding, "status": "empty"}
    try:
        dialect = csv.Sniffer().sniff(sample[:4000], delimiters=",;\t|")
        delim = dialect.delimiter
    except csv.Error:
        delim = "\t" if file_name.lower().endswith((".tsv", ".txt")) else ","
    reader = csv.reader(io.StringIO(sample), delimiter=delim)
    first = []
    for i, row in enumerate(reader):
        first.append(row[:35])
        if i >= 2:
            break
    return {
        "encoding": encoding, "delimiter": delim,
        "header_columns_sample": first[0] if first else [],
        "first_2_data_rows_sample": first[1:],
        "preview_only": True
    }


def preview_xlsx(data: bytes) -> dict:
    try:
        import openpyxl
    except ImportError:
        return {"status": "openpyxl_unavailable"}
    wb = openpyxl.load_workbook(io.BytesIO(data), read_only=True, data_only=True)
    tabs = []
    for sheet in wb.worksheets:
        r = sheet.iter_rows(min_row=1, max_row=3, values_only=True)
        first = [[str(v)[:130] if v is not None else None for v in row[:35]]
                 for row in r]
        tabs.append({"sheet": sheet.title, "max_row": sheet.max_row,
                     "max_column": sheet.max_column,
                     "first_3_rows_sample": first})
    wb.close()
    return {"sheets": tabs}


def inspect_archive(path: Path, strict_md5: bool = True) -> dict:
    md5 = hashlib.md5(path.read_bytes()).hexdigest()
    if strict_md5 and md5 != EXPECTED_MD5:
        raise ValueError(
            f"POEM source archive hash mismatch: {md5} != {EXPECTED_MD5}"
        )
    records = []
    with zipfile.ZipFile(path) as zipf:
        infos = zipf.infolist()
        if len(infos) > 2000:
            raise ValueError(f"Unexpected archive member count: {len(infos)}")
        for entry in infos:
            name = entry.filename
            posix = PurePosixPath(name)
            if posix.is_absolute() or ".." in posix.parts:
                raise ValueError(f"Unsafe ZIP entry path {name!r}")
            if entry.is_dir():
                continue
            suffix = posix.suffix.lower()
            rec = {"path": name, "uncompressed_bytes": entry.file_size,
                   "suffix": suffix}
            if suffix in DATA_SUFFIXES or suffix in CODE_SUFFIXES:
                raw = zipf.read(entry)
                rec["sha256"] = hashlib.sha256(raw).hexdigest()
                if suffix in (".csv", ".tsv", ".txt") and len(raw) < 8_000_000:
                    rec["peek"] = preview_text(raw, name)
                elif suffix == ".xlsx" and len(raw) < 8_000_000:
                    rec["peek"] = preview_xlsx(raw)
                elif suffix in CODE_SUFFIXES:
                    code = raw.decode("utf-8-sig", errors="replace")
                    matches = [(i, line[:180]) for i, line
                               in enumerate(code.splitlines(), 1)
                               if INTERESTING_PAT.search(line)]
                    rec["analysis_code_matches"] = [
                        {"line": i, "text": line} for i, line in matches[:35]
                    ]
            records.append(rec)
    summary = {
        "status": "SOURCE_MANIFEST_ONLY_NO_BIOLOGICAL_SCORING",
        "source_doi": "10.5281/zenodo.10119982",
        "archive_md5": md5,
        "archive_bytes": path.stat().st_size,
        "members": len(records),
        "data_candidate_files": sum(x["suffix"] in DATA_SUFFIXES for x in records),
        "paths": records,
        "limitations": [
            "No raw zip members are committed to KNK",
            "Row-level experimental design must be audited before scoring",
            "Archive publication date 2023; final paper 2025",
            "Both 2020/2021 cohorts are one publication cluster, not two meta studies"
        ],
    }
    return summary


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--zip", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    result = inspect_archive(args.zip)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    print("POEM_SOURCE_AUDIT_START")
    print("MD5", result["archive_md5"], "members", result["members"],
          "potential_data_files", result["data_candidate_files"])
    for r in result["paths"]:
        print("MEMBER", r["path"], "BYTES", r["uncompressed_bytes"])
        if "peek" in r:
            print("SCHEMA", json.dumps(r["peek"], ensure_ascii=False)[:1400])
        if "analysis_code_matches" in r:
            print("CODE_IMPORTS", json.dumps(r["analysis_code_matches"][:8],
                                            ensure_ascii=False)[:1800])
    print("POEM_SOURCE_AUDIT_DONE", args.output)


if __name__ == "__main__":
    main()
