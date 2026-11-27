"""Bounded official source recruitment; no scoring of outcome values."""
import argparse,csv,hashlib,io,json,re,urllib.request,urllib.error,zipfile
from pathlib import Path,PurePosixPath
from urllib.parse import quote,urlparse
MAX=30_000_000
def get(url,limit):
    req=urllib.request.Request(url,headers={"User-Agent":"AHIF academic reproducibility source audit","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=45) as r:b=r.read(limit+1);final=r.url
    assert len(b)<=limit,"source size bound"
    return b,final
def manifest(url):
    b,final=get(url,4_000_000)
    return json.loads(b),final
def inventory(b,name):
    out=[]
    if name.lower().endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            assert len(z.infolist())<=5000 and sum(i.file_size for i in z.infolist())<=300_000_000
            for f in z.infolist():
                p=PurePosixPath(f.filename);assert not p.is_absolute() and ".." not in p.parts
                if f.is_dir():continue
                v={"name":f.filename,"bytes":f.file_size}
                if f.file_size<=10_000_000 and p.suffix.lower() in (".csv",".tsv",".txt",".r",".rmd",".md"):
                    v.update(describe(z.read(f),f.filename))
                out.append(v)
    else:out=[{"name":name,"bytes":len(b),**describe(b,name)}]
    return out
def describe(b,name):
    if name.lower().endswith(".xlsx"):
        import openpyxl
        w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)
        return {"sheets":[{"name":s.title,"rows":s.max_row,"cols":s.max_column,"first_row":[str(x)[:100] if x is not None else None for x in next(s.values,())][:60]} for s in w]}
    s=b.decode("utf-8-sig",errors="replace");ext=Path(name).suffix.lower()
    if ext in (".csv",".tsv",".txt"):
        try:
            d=csv.Sniffer().sniff(s[:12000],delimiters=",;\t")
            rows=csv.reader(io.StringIO(s),d);h=next(rows,[])
            return {"header":h[:70],"n_columns":len(h),"n_rows":sum(1 for _ in rows),"delimiter":d.delimiter}
        except csv.Error:pass
    if ext in (".r",".rmd",".md",".txt"):
        return {"mapping_lines":[[i,l[:240]] for i,l in enumerate(s.splitlines(),1) if re.search(r"read[._]|load\(|readRDS|monocult|treat|species|replic|count|densit|biomass|time",l,re.I)][:40]}
    return {}
def shen():
    z={"source":"shen2023","data_doi":"10.5281/zenodo.7896204","status":"SOURCE_ONLY"}
    m,_=manifest("https://zenodo.org/api/records/7896204")
    title=m["metadata"]["title"];assert "interaction" in title.lower()
    z.update(title=title,record=m["id"],doi=m.get("doi"),related_identifiers=m["metadata"].get("related_identifiers"),files=[])
    for f in m["files"]:
        k=f["key"];size=f["size"];checksum=f["checksum"]
        a={"name":k,"size":size,"checksum":checksum,"status":"MANIFEST_ONLY"}
        if size<=MAX and Path(k).suffix.lower() in (".zip",".csv",".tsv",".xlsx",".r",".rmd",".txt"):
            try:
                url=f"https://zenodo.org/records/{m['id']}/files/"+quote(k,safe="")+"?download=1"
                b,_=get(url,size+100);assert len(b)==size
                alg,digest=checksum.split(":",1);assert hashlib.new(alg,b).hexdigest()==digest
                a.update(status="VERIFIED",sha256=hashlib.sha256(b).hexdigest(),inventory=inventory(b,k))
            except Exception as e:a.update(status="STOP_FILE",error=type(e).__name__+": "+str(e))
        z["files"].append(a)
    return z
def mahaut():
    z={"source":"mahaut2020","data_doi":"10.5061/dryad.1jwstqjqh","status":"SOURCE_ONLY"}
    m,_=manifest("https://datadryad.org/api/v2/datasets/doi%3A10.5061%2Fdryad.1jwstqjqh")
    assert "diversity effects" in m["title"].lower()
    z["title"]=m["title"];z["license"]=m.get("license");links=m["_links"]
    def uri(s):return "https://datadryad.org"+s if s.startswith("/") else s
    v,_=manifest(uri(links["stash:version"]["href"]))
    z["version"]=v.get("versionNumber");z["version_links"]=v["_links"]
    fs,_=manifest(uri(v["_links"]["stash:files"]["href"]))
    z["files"]=[]
    for f in fs["_embedded"]["stash:files"]:
        a={k:f.get(k) for k in ("path","size","digest","digestType")};a["status"]="MANIFEST_ONLY"
        if f["path"].lower().endswith(".xlsx") and f["size"]<=MAX:
            try:
                url=uri(f["_links"]["stash:download"]["href"])
                b,_=get(url,f["size"]+100);assert len(b)==f["size"]
                assert hashlib.new(f["digestType"].lower().replace("-",""),b).hexdigest()==f["digest"]
                a.update(status="VERIFIED",sha256=hashlib.sha256(b).hexdigest(),inventory=inventory(b,f["path"]))
            except Exception as e:a.update(status="STOP_FILE",error=type(e).__name__+": "+str(e))
        z["files"].append(a)
    return z
def mickalide():
    # Official repository landing page only, no guessed direct download URLs.
    b,u=get("https://doi.org/10.13012/B2IDB-0946028_V1",3_000_000)
    h=b.decode("utf-8",errors="replace")
    assert "Mickalide" in h and "Kuehn" in h
    links=[x for x in re.findall(r'href=["\']([^"\']+)["\']',h) if re.search(r"download|\.zip|\.csv|\.mat|\.xlsx|IDB-0946028",x,re.I)]
    return {"source":"mickalide2019","status":"LANDING_PAGE_ONLY","data_doi":"10.13012/B2IDB-0946028_V1","resolved_url":u,"links":links[:80]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);p.add_argument("--source",choices=("shen","mahaut","mickalide"),required=True);a=p.parse_args()
    try:z={"shen":shen,"mahaut":mahaut,"mickalide":mickalide}[a.source]()
    except Exception as e:z={"source":a.source,"status":"STOP_SOURCE_ACCESS_OR_SCHEMA","error":type(e).__name__+": "+str(e)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,ensure_ascii=False)+"\n")
    # Full metadata in derived audit, concise logs.
    print("CROSS_TAXA_SOURCE_RESULT",json.dumps(z,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
