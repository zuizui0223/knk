#!/usr/bin/env python3
"""Original 2026 Sweden 20-species field community raw-data schema audit.

Source: van Steijn et al. (2026), Ecology DOI 10.1002/ecy.70291
Zenodo original DOI 10.5281/zenodo.15363824 (2025-05-08).
All SEVEN original files are hash-verified against public Zenodo metadata,
then read ephemerally. This script does NOT compute ecological treatment scores.
No original author dataset or code is redistributed in DB artifacts.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,urllib.request
from collections import Counter
from pathlib import Path
from urllib.parse import quote

SOURCE_DOI="10.5281/zenodo.15363824"
RECORD=15363824
MD5={
    "plotdata.csv":"82905dc2f9db8740d31f27f87c507ae5",
    "Germination_indoors.csv":"cac851c56c9c774faf0a853f0db6fc17",
    "speciesinfo.csv":"bbe25d9390d87ac297691dae120b0891",
    "Pinpoint_all_years.csv":"709a1d1e453fbe5c4319294e94e2b21c",
    "Germination_outdoors.csv":"998401b75fd2961dd78b00edd8fd2c1d",
    "Germination_code.Rmd":"81d07798f7f9a08e9de52ec80e53cc07",
    "Abundance_code.Rmd":"0589cbc728d1f77053f299fc0f698313",
}
CATEGORICAL=("plot","year","site","treatment","group","species",
             "density","sow","sowing","order","replicate","sample","pin","block","time")
R_TERMS=("Pinpoint_all_years","plotdata","speciesinfo",
         "left_join","group_by","pivot_longer","pivot_wider","NMDS",
         "random","lmer","pinpoint","Year","Sowing","Plot","Site")


def fetch_original(output:Path):
    output.mkdir(parents=True,exist_ok=True)
    req=urllib.request.Request(
        f"https://zenodo.org/api/records/{RECORD}",
        headers={"Accept":"application/json","User-Agent":"AHIF-DB/0.2 academic source reproducibility"})
    with urllib.request.urlopen(req,timeout=90) as r:
        manifest=json.load(r)
    if manifest.get("id")!=RECORD or manifest.get("doi")!=SOURCE_DOI:
        raise ValueError("Zenodo record not original pinned paper")
    files={f["key"]:f for f in manifest.get("files",[])}
    if set(files)!=set(MD5):raise ValueError("Original release file list changed")
    records={}
    for name,md5 in MD5.items():
        f=files[name]
        public_md5=f.get("checksum","")
        if public_md5!=f"md5:{md5}":
            raise ValueError(f"Zenodo checksum mismatch in metadata for {name}")
        url=f.get("links",{}).get("self")
        prefix=f"https://zenodo.org/api/records/{RECORD}/files/"
        if not isinstance(url,str) or not url.startswith(prefix) or not url.endswith("/content"):
            raise ValueError("Unexpected remote URL, allowlist source files only")
        # Zenodo's public API content endpoint can return 406 for Accept
        # application/octet-stream. Prefer its ordinary documented public file
        # link; never bypass access controls or use untrusted third-party mirrors.
        download_url=f"https://zenodo.org/records/{RECORD}/files/{quote(name)}?download=1"
        req=urllib.request.Request(
            download_url,headers={"Accept":"*/*",
                                  "User-Agent":"Mozilla/5.0 (AHIF academic reproducibility)"})
        try:
            with urllib.request.urlopen(req,timeout=120) as stream:
                body=stream.read(1_000_001)
        except Exception as exc:
            raise ValueError(f"Official Zenodo public file {name} unavailable: {exc}") from exc
        if len(body)>1_000_000 or hashlib.md5(body).hexdigest()!=md5:
            raise ValueError("Original file bytes failed source check for "+name)
        (output/name).write_bytes(body)
        records[name]={"size_bytes":len(body),"md5":md5,
                       "sha256":hashlib.sha256(body).hexdigest()}
    return records

def tabular(path:Path):
    original=path.read_text(encoding="utf-8-sig")
    sn=csv.Sniffer()
    try:delim=sn.sniff(original[:8000],delimiters=",;\t").delimiter
    except csv.Error:delim=","
    # This is only a schema inventory. Some original germination source
    # worksheets have duplicate CSV header labels. Preserve and report them,
    # assign DISAMBIGUATED audit-only keys without changing source content.
    matrix=csv.reader(io.StringIO(original),delimiter=delim)
    raw_headers=next(matrix,None) or []
    duplicate_headers={k:n for k,n in Counter(raw_headers).items() if n>1}
    seen=Counter()
    headers=[]
    for i,k in enumerate(raw_headers):
        seen[k]+=1
        headers.append(k if seen[k]==1 and k!="" else
                       f"{k or 'blank_col'}__source_occurrence_{seen[k]}_{i}")
    if len(set(headers))!=len(headers):
        raise ValueError("Audit disambiguation produced repeated field keys")
    rows=[]
    for i,values in enumerate(matrix,2):
        if not any(str(x).strip() for x in values):
            continue
        if len(values)>len(headers):
            raise ValueError("Source has more cells than headers at row "+str(i))
        rows.append(dict(zip(headers,values+[""]*(len(headers)-len(values)))))
    if not rows:raise ValueError("Empty source "+path.name)
    cardinality={}
    missing={}
    for col in headers:
        values=[r.get(col,"") for r in rows]
        unique=Counter(values)
        missing[col]=sum(v is None or str(v).strip() in ("","NA","NaN") for v in values)
        if any(key in col.lower() for key in CATEGORICAL) or len(unique)<=20:
            cardinality[col]={
                "n_distinct":len(unique),
                "top_25":unique.most_common(25),
                "n_missing":missing[col],
            }
    return {
        "n_source_rows":len(rows),
        "n_columns":len(headers),
        "delim":delim,
        "headers":headers,
        "original_source_headers":raw_headers,
        "duplicate_original_headers":duplicate_headers,
        "first_5_source_rows":rows[:5],
        "last_2_source_rows":rows[-2:],
        "categorical_cardinality":cardinality,
        "missing_by_field":missing,
    }


def script_audit(path:Path):
    content=path.read_text(encoding="utf-8-sig").splitlines()
    hits=[{"line":i+1,"source":line[:280]} for i,line in enumerate(content)
          if any(term in line for term in R_TERMS)]
    return {"source_line_count":len(content),
            "selected_source_R_lines":hits[:160]}


def audit(folder:Path)->dict:
    files=fetch_original(folder)
    report={
        "status":"ORIGINAL_DATA_HASH_AND_SCHEMA_ONLY_NO_TREATMENT_EFFECTS",
        "article_doi":"10.1002/ecy.70291",
        "record_doi":SOURCE_DOI,
        "files":files,
        "tables":{},
        "original_analysis":{},
        "claim_limit":"Do not score without exact plot, site, year, repeated plot and species vs nonsown assay audit.",
    }
    for name in MD5:
        if name.endswith(".csv"):report["tables"][name]=tabular(folder/name)
        elif name.endswith(".Rmd"):report["original_analysis"][name]=script_audit(folder/name)
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir",type=Path,required=True)
    args=p.parse_args()
    report=audit(args.outdir)
    (args.outdir/"original_field_schema.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print("FIELD_SOURCE_AUDIT_START")
    for k,v in report["tables"].items():
        print("FIELD_TABLE",k,json.dumps(v,ensure_ascii=False)[:23000])
    for k,v in report["original_analysis"].items():
        print("FIELD_ORIGINAL_R",k,json.dumps(v,ensure_ascii=False)[:22000])
    print("FIELD_SOURCE_AUDIT_DONE")
if __name__=="__main__":main()
