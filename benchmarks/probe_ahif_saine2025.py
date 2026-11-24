"""Public source manifest and safe ephemeral schema audit; no outcomes scored."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,os,re,tarfile,urllib.request,zipfile
from pathlib import Path,PurePosixPath
from urllib.parse import quote
RECORD=10959957
DOI="10.5281/zenodo.10959957"
MAX=60_000_000
TOTAL=120_000_000
TERMS=re.compile(r"read[._]*(csv|table|delim|RDS|xlsx)|load\(|save\(|readCount|read_count|OTU|inocul|sample|treatment|site|log.id|logID|raref|10000|10[ ,]000|4646",re.I)
def get(url,limit,accept="*/*"):
    req=urllib.request.Request(url,headers={"Accept":accept,"User-Agent":"AHIF-DB source audit"})
    with urllib.request.urlopen(req,timeout=90) as r:
        b=r.read(limit+1)
    if len(b)>limit: raise ValueError("file exceeds safe acquisition size")
    return b
def text_schema(name,b):
    rec={"name":name,"size":len(b),"sha256":hashlib.sha256(b).hexdigest()}
    suffix=Path(name).suffix.lower()
    if suffix in (".r",".rmd",".qmd",".txt",".md",".py"):
        s=b.decode("utf-8-sig",errors="replace")
        rec["source_lines_n"]=len(s.splitlines())
        rec["mapping_lines"]=[{"line":i,"text":l[:240]} for i,l in enumerate(s.splitlines(),1) if TERMS.search(l)][:80]
    elif suffix in (".csv",".tsv"):
        s=b.decode("utf-8-sig",errors="replace")
        try: delim=csv.Sniffer().sniff(s[:16000],delimiters=",;\t").delimiter
        except csv.Error: delim="\t" if suffix==".tsv" else ","
        rows=csv.reader(io.StringIO(s),delimiter=delim)
        h=next(rows,[])
        n=0; widths={}
        for row in rows:
            n+=1; widths[len(row)]=widths.get(len(row),0)+1
        rec.update({"columns":len(h),"first_headers":h[:45],"last_headers":h[-5:],"rows":n,"row_width_counts":widths,"delimiter":delim})
    return rec
def audit_archive(name,data):
    files=[]
    if name.lower().endswith(".zip"):
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            members=z.infolist()
            if len(members)>2000 or sum(i.file_size for i in members)>TOTAL:
                raise ValueError("oversized archive")
            for i in members:
                p=PurePosixPath(i.filename)
                if p.is_absolute() or ".." in p.parts:raise ValueError("unsafe archive path")
                if i.is_dir():continue
                if i.file_size>MAX:
                    files.append({"name":i.filename,"size":i.file_size,"status":"metadata_only_size"})
                else:
                    files.append(text_schema(i.filename,z.read(i)))
    else: files.append(text_schema(name,data))
    return files
def audit(out):
    out.mkdir(parents=True,exist_ok=True)
    manifest=json.loads(get(f"https://zenodo.org/api/records/{RECORD}",2_000_000,"application/json"))
    if manifest.get("id")!=RECORD or manifest.get("doi")!=DOI: raise ValueError("record DOI mismatch")
    files=manifest.get("files",[])
    if not files: raise ValueError("no source files")
    result={"record":RECORD,"doi":DOI,"title":manifest["metadata"]["title"],"license":manifest["metadata"].get("license"),"status":"SOURCE_SCHEMA_ONLY_NO_ECOLOGICAL_SCORES","files":[]}
    used=0
    for f in files:
        name=f["key"]; size=int(f["size"]); chk=f["checksum"]
        if PurePosixPath(name).name!=name:raise ValueError("unexpected name")
        rec={"name":name,"size":size,"public_checksum":chk,"status":"not_acquired_size_gate"}
        if size<=MAX and used+size<=TOTAL:
            try:
                body=get(f"https://zenodo.org/records/{RECORD}/files/{quote(name)}?download=1",min(MAX,size+1000))
                if len(body)!=size or chk!="md5:"+hashlib.md5(body).hexdigest():raise ValueError("source byte digest/size mismatch")
                used+=size
                (out/name).write_bytes(body)
                rec.update({"status":"acquired_checksum_verified","sha256":hashlib.sha256(body).hexdigest(),"inventory":audit_archive(name,body)})
            except Exception as e:rec.update({"status":"acquisition_or_schema_failed","error":str(e)})
        result["files"].append(rec)
        print("SAINE_SOURCE_FILE",json.dumps(rec,ensure_ascii=False))
    return result
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    z=audit(a.out/"ephemeral_original")
    (a.out/"saine_manifest.json").write_text(json.dumps(z,ensure_ascii=False,indent=2)+"\n")
    print("SAINE_MANIFEST_RESULT",json.dumps(z,ensure_ascii=False))
if __name__=="__main__":main()
