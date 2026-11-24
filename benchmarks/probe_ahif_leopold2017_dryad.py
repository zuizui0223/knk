#!/usr/bin/env python3
"""Public Dryad DOI 10.5061/dryad.7p2cv raw data ACCESS/SCHEMA probe only.

Leopold et al. 2017 Ecology Letters doi:10.1111/ele.12803.
Files listed by Dryad: sample.data.csv, species.prevalence.csv, collembola.csv.

The paper already demonstrated fungal arrival × fungivory × nitrogen effects on
species composition and decomposition. A matched source-pot survey is required
BEFORE any new identity-versus-sorted-rank ecological scoring.

If official public Dryad source access is blocked, report 'access_blocked'.
Do NOT guess bytes, use mirrors of unknown provenance or infer outcomes.
"""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re
from pathlib import Path
from urllib.parse import quote,urljoin
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

DOI="doi:10.5061/dryad.7p2cv"
API="https://datadryad.org/api/v2"
TARGET={"sample.data.csv","species.prevalence.csv","collembola.csv"}
MAX_FILE=500_000


def get(url:str,accept:str="application/json",limit:int=1_000_000):
    req=Request(url,headers={"Accept":accept,
         "User-Agent":"AHIF-DB/0.1 publicly licensed scholarly data reproducibility"})
    with urlopen(req,timeout=60) as r:
        payload=r.read(limit+1)
        if len(payload)>limit:raise ValueError("Source data exceeds audited source limit")
        return payload,r.headers.get("Content-Type",""),r.geturl()


def jget(url:str):
    data,_,_=get(url)
    x=json.loads(data)
    if not isinstance(x,dict):raise ValueError("Expected Dryad API JSON object")
    return x


def items(obj,word:str):
    embedded=obj.get("_embedded") or {}
    for k,v in embedded.items():
        if word in k.lower() and isinstance(v,list):return v
    for k,v in obj.items():
        if word in k.lower() and isinstance(v,list):return v
    return []


def api_id(meta:dict,label:str)->int:
    raw=meta.get("id")
    if isinstance(raw,int):return raw
    link=(meta.get("_links") or {}).get("self") or {}
    href=link.get("href","") if isinstance(link,dict) else str(link)
    if not href:
        raise ValueError("Missing "+label+" metadata source self link")
    source_id=href.rsplit("/",1)[-1]
    if not source_id.isdigit():raise ValueError("Invalid "+label+" id from "+href)
    return int(source_id)


def data_links(file:dict,fileid:int):
    links=file.get("_links") or {}
    primary=links.get("stash:download") or {}
    if isinstance(primary,str):href=primary
    else:href=primary.get("href") or ""
    url=[]
    if href:
        u=urljoin("https://datadryad.org",href)
        if u.startswith(("https://datadryad.org/api/v2/","https://datadryad.org/downloads/")):
            url.append(u)
    for path in (f"https://datadryad.org/api/v2/files/{fileid}/download",
                 f"https://datadryad.org/downloads/file_stream/{fileid}"):
        if path not in url:url.append(path)
    return url


def probe(out:Path):
    out.mkdir(parents=True,exist_ok=True)
    root=API+"/datasets/"+quote(DOI,safe="")
    versions=jget(root+"/versions")
    v=items(versions,"versions")
    if not v:raise ValueError("Source Dryad DOI has no public version list")
    vID=api_id(max(v,key=lambda x:x.get("versionNumber",0)),"version")
    payload=jget(f"{API}/versions/{vID}/files?per_page=100")
    available=items(payload,"files")
    if not available:raise ValueError("No public original files visible")
    names={x.get("path","").split("/")[-1]:x for x in available}
    if not TARGET<=set(names):
        raise ValueError("Missing expected original files "+str(TARGET-set(names)))
    report={"doi":DOI,"source_paper_doi":"10.1111/ele.12803",
            "version_id":vID,
            "original_files":[],
            "status":"source_metadata_verified_raw_bytes_pending",
            "claim_bounds":["These files have NOT been source scored into identity-vs-rank history effects.",
               "Dryad access errors must not be bypassed or replaced by invented outcomes.",
               "Original authors already analyzed community composition and wood decomposition.",
               "Group-level fungal composition includes 9 wood-subsample prevalence replicates per wood disk, not independent jars."]}
    success=0
    for basename in sorted(TARGET):
        meta=names[basename];id_=api_id(meta,"file")
        results={"name":basename,"file_id":id_,"declared_size":meta.get("size"),
                 "declared_digest":meta.get("digest") or meta.get("checksum"),
                 "attempts":[],"downloaded":False}
        for url in data_links(meta,id_):
            try:
                body,mime,final=get(url,accept="text/csv,*/*",limit=MAX_FILE)
                if b"<html" in body[:200].lower() or b"<!doctype html" in body[:200].lower():
                    results["attempts"].append({"url":url,"status":"html_validation_not_csv"})
                    continue
                parsed=list(csv.reader(io.StringIO(body.decode("utf-8-sig"))))
                if len(parsed)<2 or len(parsed[0])<2:
                    results["attempts"].append({"url":url,"status":"invalid_CSV"})
                    continue
                original=out/basename
                original.write_bytes(body)
                results.update({"downloaded":True,"sha256":hashlib.sha256(body).hexdigest(),
                     "byte_length":len(body),"rows_excluding_header":len(parsed)-1,
                     "column_header":parsed[0],
                     "source_url":url,"mime_type":mime})
                success+=1;break
            except (HTTPError,URLError,TimeoutError,ValueError,UnicodeDecodeError) as exc:
                results["attempts"].append({"url":url,"status":str(exc)[:250]})
        report["original_files"].append(results)
    report["status"]="source_files_retrieved" if success==len(TARGET) else (
        "source_partial_or_access_blocked")
    report["n_source_files_retrieved"]=success
    (out/"leopold2017_source_manifest.json").write_text(
        json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("LEOPOLD2017_SOURCE_MANIFEST",json.dumps({
       "status":report["status"],"version_id":vID,"retrieved":success,
       "files":[{k:x.get(k) for k in ("name","file_id","downloaded","rows_excluding_header","column_header","attempts")}
                for x in report["original_files"]]}))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    probe(a.output_dir)


if __name__=="__main__":main()
