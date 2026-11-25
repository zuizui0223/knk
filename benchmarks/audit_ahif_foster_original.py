"""Official EDI version-locked original-source access and grain inspection."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,urllib.request,urllib.error,xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from urllib.parse import urlparse

BASE="https://pasta.edirepository.org"
DOI="10.6073/pasta/5703e8819c81da0beed4660f467c288e"
ALLOWED={"pasta.edirepository.org","pasta.lternet.edu"}
CAP=25_000_000
TOTAL=70_000_000

def request(url,limit):
    u=urlparse(url)
    if u.scheme!="https" or u.hostname not in ALLOWED:
        raise ValueError("unexpected official host")
    req=urllib.request.Request(url,headers={"User-Agent":"AHIF source reproducibility","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=35) as response:
        v=urlparse(response.geturl())
        if v.scheme!="https" or v.hostname not in ALLOWED:raise ValueError("unexpected final host")
        b=response.read(limit+1)
    if len(b)>limit:raise ValueError("source byte limit")
    return b

def resolve(text):
    meta=[];data=[]
    for line in text.splitlines():
        if not line.strip():continue
        u=urlparse(line.strip())
        if u.scheme!="https" or u.hostname not in ALLOWED:raise ValueError("unexpected resource-map URL")
        q=u.path.strip("/").split("/")
        if len(q)==6 and q[:3]==["package","metadata","eml"]:meta.append(tuple(q[3:]))
        if len(q)==7 and q[:3]==["package","data","eml"]:data.append((tuple(q[3:6]),q[6]))
    if len(meta)!=1 or not data or any(k!=meta[0] for k,id in data) or len(set(id for k,id in data))!=len(data):
        raise ValueError("resource-map version/coverage mismatch")
    if not re.fullmatch(r"[A-Za-z0-9_.-]+",meta[0][0]) or not all(x.isdigit() for x in meta[0][1:]):
        raise ValueError("unexpected EDI identifier")
    if any(not re.fullmatch(r"[0-9a-fA-F]{32}",i) for k,i in data):
        raise ValueError("invalid entity identifier")
    return "/".join(meta[0]),[i for k,i in data]

def textnode(e):
    return " ".join(" ".join(e.itertext()).split()) if e is not None else ""

def parse_eml(raw,package):
    if b"<!DOCTYPE" in raw.upper():raise ValueError("EML DOCTYPE not accepted")
    root=ET.fromstring(raw)
    for n in root.iter():n.tag=n.tag.split("}")[-1]
    if root.attrib.get("packageId")!=package.replace("/","."):
        raise ValueError("EML packageId differs")
    ds=root.find("dataset")
    if ds is None:raise ValueError("EML dataset absent")
    title=textnode(ds.find("title"))
    if "hay meadow" not in title.lower() and "hayfield" not in title.lower():
        raise ValueError("EML title does not identify hay meadow")
    tables=[]
    for tab in ds:
        if tab.tag not in ("dataTable","otherEntity","spatialRaster","spatialVector"):continue
        attrs=[]
        for a in tab.findall("attributeList/attribute"):
            attrs.append({"name":textnode(a.find("attributeName")),"label":textnode(a.find("attributeLabel")),
               "definition":textnode(a.find("attributeDefinition"))[:500],
               "scale":textnode(a.find("measurementScale"))[:1200],
               "missing_codes":[textnode(v.find("code")) for v in a.findall("missingValueCode")]})
        tables.append({"kind":tab.tag,"entity_name":textnode(tab.find("entityName")),
           "description":textnode(tab.find("entityDescription"))[:3000],
           "objects":[textnode(x) for x in tab.findall("physical/objectName")],
           "field_delimiters":[textnode(x) for x in tab.findall(".//fieldDelimiter")],
           "record_count":textnode(tab.find("numberOfRecords")),
           "attributes":attrs})
    return {"title":title,"package_id":root.attrib["packageId"],"abstract":textnode(ds.find("abstract"))[:9000],
       "methods":textnode(ds.find("methods"))[:18000],"temporal_coverage":textnode(ds.find("coverage/temporalCoverage"))[:1000],
       "intellectual_rights":textnode(ds.find("intellectualRights"))[:2500],
       "entities":tables}

def decode(b):
    for enc in ("utf-8-sig","cp1252","latin-1"):
        try:return b.decode(enc),enc
        except UnicodeError:pass
    raise ValueError("unable to decode")

def schema_csv(b):
    s,enc=decode(b)
    try:sep=csv.Sniffer().sniff(s[:10000],delimiters=",;\t|").delimiter
    except csv.Error:return {"status":"TEXT_NOT_PARSED_AS_TABLE","encoding":enc,
                           "n_lines":len(s.splitlines())}
    rd=csv.reader(io.StringIO(s,newline=""),delimiter=sep)
    header=next(rd,[])
    if len(header)<2 or len(set(header))!=len(header):return {"status":"AMBIGUOUS_HEADER","header":header[:150]}
    cnt=[Counter() for _ in header];n=bad=0
    for row in rd:
        if not row or not any(x.strip() for x in row):continue
        n+=1
        if len(row)!=len(header):bad+=1;continue
        for c,x in zip(cnt,row):c[x.strip()]+=1
    fields=[]
    for h,c in zip(header,cnt):
        f={"name":h,"distinct":len(c),"missing":sum(v for k,v in c.items() if k in ("","NA","NaN","-999","NULL"))}
        if len(c)<=40 and re.search(r"seed|treat|block|year|plot|site|disturb|fert|hay|sow|group|rep|species|taxon",h,re.I):
            f["levels"]=dict(sorted(c.items()))
        fields.append(f)
    return {"status":"TABLE_SCHEMA_ONLY","encoding":enc,"separator":repr(sep),"n_rows":n,"nonrectangular_rows":bad,"fields":fields}

def inventory(root):
    root.mkdir(parents=True,exist_ok=True)
    z={"status":"STOP_PENDING","source_doi":DOI,"source_api":BASE}
    try:
        u=BASE+"/package/doi/doi:"+DOI
        mapping=request(u,2_000_000).decode("utf-8")
        pkg,ids=resolve(mapping)
        z["package_path"]=pkg;z["entity_ids"]=ids
        declared=request(BASE+"/package/doi/eml/"+pkg,2000).decode().strip().removeprefix("doi:")
        if declared.lower()!=DOI:raise ValueError("resolved DOI differs")
        eml=request(BASE+"/package/metadata/eml/"+pkg,4_000_000)
        z["eml_sha256"]=hashlib.sha256(eml).hexdigest()
        z["metadata"]=parse_eml(eml,pkg);z["entities"]=[]
        used=0
        for id in ids:
            q=pkg+"/"+id
            e={"entity_id":id,"status":"STOP_PENDING"}
            z["entities"].append(e)
            try:
                name=request(BASE+"/package/name/eml/"+q,4000).decode().strip()
                size=int(request(BASE+"/package/data/size/eml/"+q,1000).decode().strip())
                sha=request(BASE+"/package/data/checksum/eml/"+q,2000).decode().strip().lower()
                e.update(name=name,declared_bytes=size,declared_sha1=sha)
                if not re.fullmatch("[0-9a-f]{40}",sha) or size<0:raise ValueError("invalid declared checksum/size")
                if size>CAP or used+size>TOTAL:
                    e["status"]="STOP_SIZE_CAP";continue
                b=request(BASE+"/package/data/eml/"+q,CAP)
                if len(b)!=size or hashlib.sha1(b).hexdigest()!=sha:raise ValueError("original data byte mismatch")
                used+=len(b);e["sha256"]=hashlib.sha256(b).hexdigest();e["status"]="SOURCE_BYTES_VERIFIED"
                (root/id).write_bytes(b)
                # EML names may not have extensions. Parse text only when
                # it decodes and has no binary NUL; do not execute content.
                if b"\0" not in b[:1000] and not b.startswith((b"PK",b"\xd0\xcf\x11\xe0",b"%PDF")):
                    e["schema"]=schema_csv(b)
                else:e["schema"]={"status":"BINARY_OR_DOCUMENT_NOT_PARSED"}
            except urllib.error.HTTPError as ex:
                e.update(status="STOP_ENTITY_HTTP",http_status=ex.code)
            except (ValueError,OSError,UnicodeError) as ex:
                e.update(status="STOP_ENTITY_GATE",reason=str(ex)[:500])
        z["status"]="SOURCE_METADATA_VERIFIED_ENTITY_AUDIT_NOT_SCORED"
        z["n_verified_entities"]=sum(e["status"]=="SOURCE_BYTES_VERIFIED" for e in z["entities"])
        z["verified_bytes"]=used
    except urllib.error.HTTPError as ex:
        z.update(status="STOP_PUBLIC_SOURCE_HTTP",http_status=ex.code,source=ex.url)
    except (ValueError,OSError,ET.ParseError) as ex:
        z.update(status="STOP_SOURCE_GATE",reason=str(ex)[:500])
    return z

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=inventory(a.root);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("FOSTER_ORIGINAL_AUDIT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
