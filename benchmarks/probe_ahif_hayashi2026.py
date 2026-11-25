"""Source-only metadata audit, no biological response scoring."""
import argparse,csv,hashlib,io,json,re,tarfile,urllib.request,zipfile
from pathlib import Path,PurePosixPath
from urllib.parse import quote
DOI="10.5281/zenodo.19396153";REC=19396153
MAX_FILE=45_000_000;MAX_TOTAL=90_000_000
META=re.compile(r"sample|metadata|plate|replicate|timepoint|day|dilution|culture|otu|asv|read|taxon|community|filter|qc",re.I)
def get(url,limit,accept="*/*"):
    req=urllib.request.Request(url,headers={"User-Agent":"AHIF public source audit","Accept":accept})
    with urllib.request.urlopen(req,timeout=90) as r:b=r.read(limit+1)
    if len(b)>limit:raise ValueError("bounded download exceeded")
    return b
def schema(name,b):
    x={"name":name,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()}
    ext=Path(name).suffix.lower()
    if ext in (".r",".rmd",".md",".txt",".py",".qmd"):
        s=b.decode("utf-8-sig",errors="replace")
        # only lines describing data, mapping and QC; no biological scores
        x["n_lines"]=len(s.splitlines())
        x["mapping_lines"]=[{"line":i,"text":l[:220]} for i,l in enumerate(s.splitlines(),1) if META.search(l)][:50]
    elif ext in (".csv",".tsv"):
        s=b.decode("utf-8-sig",errors="replace")
        try:delim=csv.Sniffer().sniff(s[:12000],delimiters=",;\t").delimiter
        except csv.Error:delim="\t" if ext==".tsv" else ","
        r=csv.reader(io.StringIO(s),delimiter=delim);h=next(r,[])
        n=0; widths={}
        for line in r:
            n+=1;widths[len(line)]=widths.get(len(line),0)+1
        x.update(rows=n,columns=len(h),first_headers=h[:24],last_headers=h[-5:],row_widths=widths)
    return x
def inventory(name,b):
    if name.lower().endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            inf=z.infolist()
            if len(inf)>3000 or sum(x.file_size for x in inf)>180_000_000:raise ValueError("archive expansion limit")
            ans=[]
            for i in inf:
                if i.is_dir():continue
                p=PurePosixPath(i.filename)
                assert not p.is_absolute() and ".." not in p.parts
                ans.append(schema(i.filename,z.read(i)) if i.file_size<=MAX_FILE else {"name":i.filename,"bytes":i.file_size,"status":"size_gate"})
            return ans
    return [schema(name,b)]
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True);out={"status":"SOURCE_METADATA_ONLY","cited_doi":DOI,"files":[]}
    try:
        r=json.loads(get(f"https://zenodo.org/api/records/{REC}",2_000_000,"application/json"))
        assert DOI in (r.get("doi"),r.get("conceptdoi"))
        vid=int(r["id"]);assert r["doi"]==f"10.5281/zenodo.{vid}"
        out["record"]={k:r.get(k) for k in ("id","doi","conceptdoi","conceptrecid")}
        out["title"]=r["metadata"]["title"];out["license"]=r["metadata"].get("license")
        used=0
        for f in sorted(r.get("files",[]),key=lambda x:x["size"]):
            name=f["key"];size=int(f["size"]);chk=f["checksum"]
            assert PurePosixPath(name).name==name
            d={"name":name,"size":size,"checksum":chk,"status":"size_gate"}
            if size<=MAX_FILE and used+size<=MAX_TOTAL:
                try:
                    b=get(f"https://zenodo.org/records/{vid}/files/{quote(name)}?download=1",size+100)
                    assert len(b)==size and chk=="md5:"+hashlib.md5(b).hexdigest()
                    used+=size
                    d.update(status="checksum_verified",sha256=hashlib.sha256(b).hexdigest(),inventory=inventory(name,b))
                except Exception as e:d.update(status="file_access_or_schema_stop",error=str(e))
            out["files"].append(d)
        out["status"]="SOURCE_MANIFEST_VERIFIED_NO_SCORE" if out["files"] else "STOP_NO_FILE_MANIFEST"
    except Exception as e:out.update(status="STOP_SOURCE_ACCESS_OR_IDENTITY",error=str(e))
    (a.out/"audit.json").write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
    # Keep metadata log compact; full inventory in audit artifact.
    brief={k:v for k,v in out.items() if k!="files"}
    brief["files"]=[]
    for f in out["files"]:
        d={k:v for k,v in f.items() if k!="inventory"}
        d["inventory"]=[{k:v for k,v in z.items() if k!="mapping_lines"} for z in f.get("inventory",[])]
        d["mapping_excerpt"]=[{"name":z["name"],"lines":z.get("mapping_lines",[])[:12]} for z in f.get("inventory",[]) if z.get("mapping_lines")]
        brief["files"].append(d)
    print("HAYASHI_SOURCE_RESULT",json.dumps(brief,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
