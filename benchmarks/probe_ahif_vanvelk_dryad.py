"""Audit official public Dryad file endpoint; stop on access or identity errors."""
import argparse,csv,hashlib,io,json,re,urllib.request
from pathlib import Path
from urllib.parse import quote,urljoin,urlparse
BASE="https://datadryad.org"
DOI="10.5061/dryad.x95x69pss"
def get(u,limit,accept="application/json"):
    req=urllib.request.Request(u,headers={"Accept":accept,"User-Agent":"AHIF public metadata/source audit"})
    with urllib.request.urlopen(req,timeout=60) as r:
        b=r.read(limit+1)
    if len(b)>limit:raise ValueError("size gate")
    return b
def link(obj,key):
    v=obj.get("_links",{}).get(key)
    if not isinstance(v,dict) or "href" not in v:raise ValueError("missing public API link "+key)
    u=urljoin(BASE,v["href"])
    if urlparse(u).netloc!="datadryad.org":raise ValueError("unrecognized API host")
    return u
def read(u):return json.loads(get(u,2_000_000))
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    z=dict(status="SOURCE_GATE_ONLY",doi=DOI,files=[])
    try:
        d=read(BASE+"/api/v2/datasets/"+quote("doi:"+DOI,safe=""))
        z["dataset"]={k:d.get(k) for k in ("identifier","title","versionNumber","id","_links")}
        assert d.get("identifier") in ("doi:"+DOI,DOI),"dataset DOI mismatch"
        assert "eco-evolutionary interactions" in d["title"].lower(),"title mismatch"
        v=read(link(d,"stash:version"))
        z["version"]={k:v.get(k) for k in ("versionNumber","id","_links")}
        f=read(link(v,"stash:files"))
        fs=f.get("_embedded",{}).get("stash:files",[])
        if not fs:raise ValueError("no public file records")
        for r in fs:
            name=r.get("path") or r.get("name");size=r.get("size")
            q=dict(name=name,size=size,checksum=r.get("digest"),digest_type=r.get("digestType"),status="metadata_only",links=r.get("_links"))
            z["files"].append(q)
            if name and name.lower() in ("cc-b.csv","allbiomass.csv","readme.md"):
                assert isinstance(size,int) and 0<size<2_000_000
                # Only the official API's advertised public download link.
                b=get(link(r,"stash:download"),size+100,"application/octet-stream")
                assert len(b)==size,"size mismatch"
                digest=r.get("digest");typ=str(r.get("digestType","")).replace("-","").lower()
                if not digest or typ not in ("md5","sha1","sha256","sha512"):raise ValueError("missing published digest")
                assert hashlib.new(typ,b).hexdigest().lower()==digest.lower(),"published digest mismatch"
                q.update(status="checksum_verified",sha256=hashlib.sha256(b).hexdigest())
                (a.out/name).write_bytes(b)
                if name.lower().endswith(".csv"):
                    txt=b.decode("utf-8-sig");dial=csv.Sniffer().sniff(txt[:10000],delimiters=",;\t")
                    rows=list(csv.DictReader(io.StringIO(txt),dialect=dial));assert rows
                    q.update(header=list(rows[0]),rows=len(rows),delimiter=dial.delimiter)
        got={x["name"].lower() for x in z["files"] if x["status"]=="checksum_verified"}
        z["status"]="VERIFIED_SOURCE_NO_SCORES" if {"cc-b.csv","allbiomass.csv"}<=got else "STOP_REQUIRED_FILES_ABSENT"
    except Exception as e:z.update(status="STOP_PUBLIC_SOURCE_ACCESS_OR_SCHEMA",error=type(e).__name__+": "+str(e))
    (a.out/"audit.json").write_text(json.dumps(z,indent=2)+"\n")
    print("VANVELK_DRYAD_RESULT",json.dumps(z),flush=True)
if __name__=="__main__":main()
