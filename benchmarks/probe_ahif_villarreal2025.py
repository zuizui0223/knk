"""Acquire public Figshare v1 metadata and bounded source files, schema only."""
import argparse,csv,hashlib,io,json,re,urllib.request,urllib.error,zipfile
from pathlib import Path,PurePosixPath
from urllib.parse import urlparse
ARTICLE=28651985
DOI="10.6084/m9.figshare.28651985.v1"
MAX=20_000_000
TOTAL=50_000_000
TERMS=re.compile(r"read[._]|load\(|readRDS|tank|pair|replic|treat|species|sample|abund|day|time|A0|B0|Aa|Ab",re.I)
def get(url,limit):
    req=urllib.request.Request(url,headers={"User-Agent":"AHIF academic reproducibility","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=75) as r:b=r.read(limit+1)
    if len(b)>limit:raise ValueError("bounded source exceeded")
    return b
def schema(name,b):
    ext=Path(name).suffix.lower()
    z={"name":name,"size":len(b),"sha256":hashlib.sha256(b).hexdigest()}
    if ext in (".csv",".tsv",".txt",".dat"):
        s=b.decode("utf-8-sig",errors="replace")
        try:d=csv.Sniffer().sniff(s[:16000],delimiters=",;\t")
        except csv.Error:d=None
        if d:
            rr=csv.reader(io.StringIO(s),dialect=d);h=next(rr,[])
            widths={};n=0
            for row in rr:n+=1;widths[len(row)]=widths.get(len(row),0)+1
            z.update(header_n=len(h),header_first=h[:25],header_last=h[-6:],rows=n,row_widths=widths,delimiter=d.delimiter)
        else:z["n_lines"]=len(s.splitlines())
    if ext in (".r",".rmd",".qmd",".md",".txt"):
        s=b.decode("utf-8-sig",errors="replace")
        z["mapping_lines"]=[{"line":i,"text":l[:240]} for i,l in enumerate(s.splitlines(),1) if TERMS.search(l)][:70]
    return z
def inventory(name,b):
    if name.lower().endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(b)) as z:
            files=z.infolist();assert len(files)<=500 and sum(f.file_size for f in files)<=100_000_000
            out=[]
            for f in files:
                p=PurePosixPath(f.filename);assert not p.is_absolute() and ".." not in p.parts
                if f.is_dir():continue
                if f.file_size<=MAX:out.append(schema(f.filename,z.read(f)))
                else:out.append({"name":f.filename,"size":f.file_size,"status":"size_gate"})
            return out
    return [schema(name,b)]
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True);z={"status":"SOURCE_GATE_ONLY","article":ARTICLE,"doi":DOI,"files":[]}
    try:
        m=json.loads(get(f"https://api.figshare.com/v2/articles/{ARTICLE}/versions/1",2_000_000))
        z["observed_identity"]={k:m.get(k) for k in ("id","version","doi","title","url_public_html","published_date")}
        z["observed_file_names"]=[f.get("name") for f in m.get("files",[])]
        assert int(m["id"])==ARTICLE, "article ID mismatch"
        assert int(m["version"])==1, "version mismatch"
        assert m["doi"] in (DOI,DOI[:-3]), "DOI mismatch"
        # Title is descriptive and may not match the article exactly; require the published DOI above.
        # Source content mapping remains a separate gate before scoring.
        z.update(title=m["title"],license=m.get("license"),description=m.get("description","")[:8000])
        used=0
        for f in m["files"]:
            name=f["name"];size=int(f["size"]);md5=f.get("computed_md5") or f.get("supplied_md5")
            d={"name":name,"file_id":f["id"],"size":size,"md5":md5,"status":"size_gate"}
            assert PurePosixPath(name).name==name
            if size<=MAX and used+size<=TOTAL:
                try:
                    url=f.get("downloadUrl") or f.get("download_url")
                    assert isinstance(url,str),"published download URL missing: "+str(sorted(f))
                    assert urlparse(url).scheme=="https" and urlparse(url).hostname in ("figshare.com","api.figshare.com","ndownloader.figshare.com","s3-eu-west-1.amazonaws.com","s3.amazonaws.com")
                    assert not f.get("is_link_only") and re.fullmatch("[0-9a-f]{32}",md5)
                    b=get(url,size+100)
                    assert len(b)==size and hashlib.md5(b).hexdigest()==md5
                    used+=size;(a.out/name).write_bytes(b)
                    d.update(status="checksum_verified",sha256=hashlib.sha256(b).hexdigest(),inventory=inventory(name,b))
                except Exception as e:d.update(status="STOP_FILE",error=type(e).__name__+": "+str(e))
            z["files"].append(d)
        z["n_checksum_verified"]=sum(f["status"]=="checksum_verified" for f in z["files"])
        z["status"]="MANIFEST_AND_SOURCE_AUDITED_NO_ECOLOGICAL_SCORES" if z["files"] and z["n_checksum_verified"]==len(z["files"]) else "STOP_SOURCE_FILES_NOT_FULLY_VERIFIED"
    except Exception as e:z.update(status="STOP_SOURCE_IDENTITY_OR_ACCESS",error=type(e).__name__+": "+str(e))
    (a.out/"audit.json").write_text(json.dumps(z,ensure_ascii=False,indent=2)+"\n")
    print("VILLARREAL_SOURCE_RESULT",json.dumps(z,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
