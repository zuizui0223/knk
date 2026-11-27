"""Plant context recruitment: official public manifests, checksum-verified table schemas."""
import argparse,hashlib,json,io
from pathlib import Path
from urllib.parse import quote,urlparse
from benchmarks.probe_ahif_context_cross_taxa import get,manifest,inventory
LIMIT=30_000_000
def reynaert():
    m,_=manifest("https://zenodo.org/api/records/10167759")
    assert int(m["id"])==10167759
    title=m["metadata"]["title"]
    assert "weather" in title.lower() and "interaction" in title.lower()
    out={"source":"reynaert2023","record":m["id"],"doi":m.get("doi"),"title":title,"files":[]}
    for f in m["files"]:
        name=f["key"];size=f["size"];h=f["checksum"]
        z={"name":name,"size":size,"checksum":h,"status":"MANIFEST_ONLY"}
        if size<=LIMIT and Path(name).suffix.lower() in (".xlsx",".csv",".zip",".txt",".r",".rmd"):
            try:
                b,_=get("https://zenodo.org/records/10167759/files/"+quote(name,safe="")+"?download=1",size+100)
                alg,v=h.split(":");assert len(b)==size and hashlib.new(alg,b).hexdigest()==v
                z.update(status="VERIFIED",sha256=hashlib.sha256(b).hexdigest(),inventory=inventory(b,name))
            except Exception as e:z.update(status="STOP_FILE",error=type(e).__name__+": "+str(e))
        out["files"].append(z)
    return out
def yang():
    m,_=manifest("https://api.figshare.com/v2/articles/24174558/versions/5")
    assert int(m["id"])==24174558 and int(m["version"])==5 and m["doi"]=="10.6084/m9.figshare.24174558.v5"
    assert "competitiveness" in m["title"].lower()
    out={"source":"yang2024","doi":m["doi"],"title":m["title"],"description":m.get("description","")[:5000],"files":[]}
    for f in m["files"]:
        name=f["name"];size=f["size"];h=f.get("computed_md5") or f.get("supplied_md5")
        z={"name":name,"id":f["id"],"size":size,"md5":h,"status":"MANIFEST_ONLY"}
        if size<=LIMIT and Path(name).suffix.lower() in (".xlsx",".csv",".zip",".txt",".r",".rmd"):
            try:
                u=f.get("download_url") or f.get("downloadUrl")
                assert urlparse(u).scheme=="https"
                b,_=get(u,size+100)
                assert len(b)==size and hashlib.md5(b).hexdigest()==h
                z.update(status="VERIFIED",sha256=hashlib.sha256(b).hexdigest(),inventory=inventory(b,name))
            except Exception as e:z.update(status="STOP_FILE",error=type(e).__name__+": "+str(e))
        out["files"].append(z)
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--source",choices=("reynaert","yang"),required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    try:z={"status":"SOURCE_ONLY_NO_OUTCOME",**{"reynaert":reynaert,"yang":yang}[a.source]()}
    except Exception as e:z={"source":a.source,"status":"STOP_SOURCE_ACCESS_OR_SCHEMA","error":type(e).__name__+": "+str(e)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,ensure_ascii=False)+"\n")
    print("PLANT_CONTEXT_SOURCE_RESULT",json.dumps(z,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
