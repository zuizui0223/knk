#!/usr/bin/env python3
"""Retrieve one public, pinned Dryad dataset file with source checksum audit.

No authentication, no scraping private datasets; public API metadata and
public download endpoints only. Raw XLSX kept in ephemeral GitHub Actions job.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from urllib.request import Request,urlopen
from urllib.parse import quote
import re

DOI="doi:10.5061/dryad.jq2bvq8jj"
API="https://datadryad.org/api/v2"


def get_json(url:str)->dict:
    req=Request(url,headers={"User-Agent":"KNK-open-data-source-audit/1.0 (academic source verification)",
                             "Accept":"application/json"})
    with urlopen(req,timeout=75) as r:
        data=r.read(4_000_000)
    parsed=json.loads(data)
    if not isinstance(parsed,dict):
        raise ValueError("Unexpected Dryad JSON response")
    return parsed


def pick_embedded(obj:dict,keyword:str)->list:
    if "_embedded" in obj:
        for k,v in obj["_embedded"].items():
            if keyword in k.lower() and isinstance(v,list):return v
    for k,v in obj.items():
        if keyword in k.lower() and isinstance(v,list):return v
    return []


def download(url:str,outpath:Path,limit:int=2_000_000):
    req=Request(url,headers={"User-Agent":"Mozilla/5.0 (academic reproducibility; original published public dataset)",
                            "Accept":"application/octet-stream"})
    with urlopen(req,timeout=120) as r:
        kind=r.headers.get("Content-Type","")
        block=r.read(limit+1)
    if len(block)>limit:
        raise ValueError("Dryad public download exceeds expected file size")
    if not block.startswith(b"PK\x03\x04"):
        raise ValueError("Dryad download did not return an XLSX workbook; MIME="+kind+" header="+repr(block[:100]))
    outpath.write_bytes(block)
    return kind


def fetch_original(output_dir:Path)->dict:
    output_dir.mkdir(parents=True,exist_ok=True)
    root=f"{API}/datasets/{quote(DOI,safe='')}"
    versions=get_json(root+"/versions")
    ver=pick_embedded(versions,"versions")
    if not ver and isinstance(versions.get("versions"),list):ver=versions["versions"]
    if not ver:raise ValueError("Dryad dataset versions missing; cannot ground source version")
    target=max(ver,key=lambda x:x.get("versionNumber",x.get("id",0)))
    version_id=target.get("id")
    if not isinstance(version_id,int):
        # Public Dryad versions often expose their numeric ID only inside a
        # version self-link, not as a top-level ID field.
        link=(target.get("_links") or {}).get("self") or {}
        href=link.get("href","") if isinstance(link,dict) else str(link)
        match_id=(href.rsplit("/",1)[-1] if href.startswith("/api/v2/versions/") else "")
        if not match_id.isdigit():
            raise ValueError("Invalid Dryad version self href: "+repr(href))
        version_id=int(match_id)
    files_json=get_json(f"{API}/versions/{version_id}/files?per_page=100")
    files=pick_embedded(files_json,"files")
    if not files and isinstance(files_json.get("files"),list):files=files_json["files"]
    file_targets=[x for x in files if x.get("path","").split("/")[-1].lower()=="fulldata.xlsx"]
    if len(file_targets)!=1:
        raise ValueError(f"Expected one Fulldata.xlsx, found {len(file_targets)}; files={[x.get('path') for x in files]}")
    file=file_targets[0]
    file_id=file.get("id")
    if not isinstance(file_id,int):raise ValueError("Expected public file ID")
    links=file.get("_links") or {}
    link=links.get("stash:download",{})
    if isinstance(link,str):download_url=link
    else:download_url=link.get("href") or ""
    if not download_url:
        download_url=f"https://datadryad.org/downloads/file_stream/{file_id}"
    if download_url.startswith("/"):download_url="https://datadryad.org"+download_url
    if not re.match(r"^https://datadryad\.org/(downloads/file_stream|api/v2/files)/",download_url):
        raise ValueError("Unexpected outbound file download host or path")
    dest=output_dir/"Fulldata.xlsx"
    try:
        mime=download(download_url,dest)
    except Exception as e:
        fallback=f"https://datadryad.org/downloads/file_stream/{file_id}"
        if fallback==download_url:raise
        mime=download(fallback,dest)
        download_url=fallback
    sha=hashlib.sha256(dest.read_bytes()).hexdigest()
    declared=file.get("digest") or file.get("checksum") or ""
    digest_type=(file.get("digestType") or file.get("checksumType") or "").lower()
    if declared and digest_type in ("sha-256","sha256") and sha.lower()!=str(declared).lower():
        raise ValueError("Dryad published SHA256 mismatch")
    if declared and digest_type=="md5" and hashlib.md5(dest.read_bytes()).hexdigest()!=str(declared).lower():
        raise ValueError("Dryad published MD5 mismatch")
    return {
        "status":"PUBLIC_API_DRYAD_FILE_SOURCE_HASH_VERIFIED",
        "doi":DOI,
        "version_id":version_id,"version_number":target.get("versionNumber"),
        "file_id":file_id,"path":file["path"],
        "declared_digest":declared,
        "declared_digest_type":digest_type,
        "download_url":download_url,
        "file_sha256":sha,
        "bytes":dest.stat().st_size,
        "http_mime":mime,
        "other_source_file_names":[x.get("path") for x in files]
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    args=p.parse_args()
    result=fetch_original(args.output_dir)
    (args.output_dir/"api_source_manifest.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("JING_DRYAD_SOURCE_MANIFEST",json.dumps(result,ensure_ascii=False))


if __name__=="__main__":main()
