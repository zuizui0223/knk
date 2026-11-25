"""Public, cross-linked author archive: metadata and checksum gate only."""
import argparse,csv,hashlib,io,json,re,urllib.request
from pathlib import Path,PurePosixPath
from collections import Counter
from urllib.parse import quote
DOI="10.5281/zenodo.10964334";REC=10964334
def get(url,max_bytes,accept="*/*"):
    req=urllib.request.Request(url,headers={"Accept":accept,"User-Agent":"AHIF public source audit"})
    with urllib.request.urlopen(req,timeout=90) as r:b=r.read(max_bytes+1)
    if len(b)>max_bytes:raise ValueError("size gate")
    return b
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    z=dict(status="SOURCE_GATE_ONLY",doi=DOI,files=[])
    try:
        r=json.loads(get(f"https://zenodo.org/api/records/{REC}",2_000_000,"application/json"))
        assert r.get("doi")==DOI and r.get("id")==REC,"record identity mismatch"
        z["record"]={k:r.get(k) for k in ("id","doi","conceptdoi")}
        z["title"]=r["metadata"]["title"];z["license"]=r["metadata"].get("license")
        assert "eco-evolutionary interactions" in z["title"].lower(),"title mismatch"
        for f in r.get("files",[]):
            name=f["key"];size=int(f["size"]);md5=f["checksum"]
            d=dict(name=name,size=size,checksum=md5,status="metadata_only")
            if name.lower() in ("cc-b.csv","allbiomass.csv","readme.md"):
                assert PurePosixPath(name).name==name and size<2_000_000
                b=get(f"https://zenodo.org/records/{REC}/files/{quote(name)}?download=1",size+100)
                assert len(b)==size and md5=="md5:"+hashlib.md5(b).hexdigest(),"source byte mismatch"
                (a.out/name).write_bytes(b)
                d.update(status="checksum_verified",sha256=hashlib.sha256(b).hexdigest())
                if name.lower().endswith(".csv"):
                    text=b.decode("utf-8-sig");dialect=csv.Sniffer().sniff(text[:10000],delimiters=",;\t")
                    rows=list(csv.DictReader(io.StringIO(text),dialect=dialect))
                    assert rows and all(None not in x for x in rows)
                    d.update(header=list(rows[0]),rows=len(rows),delimiter=dialect.delimiter)
                    for c in ("Time","Treatment","Replicate","Replica","Name"):
                        if c in rows[0]:
                            v=Counter(x[c] for x in rows);d[c]={"n_distinct":len(v),"counts":dict(v) if len(v)<=40 else None}
                    d["missing_cells_by_column"]={k:sum(str(x[k]).strip().lower() in ("","na","nan","null","none") for x in rows) for k in rows[0] if k not in ("Treatment","Name","Replica")}
            z["files"].append(d)
        got={f["name"].lower() for f in z["files"] if f["status"]=="checksum_verified"}
        z["status"]="VERIFIED_SOURCE_NO_SCORES" if {"cc-b.csv","allbiomass.csv"}<=got else "STOP_REQUIRED_FILES_ABSENT"
    except Exception as e:z.update(status="STOP_SOURCE_ACCESS_OR_IDENTITY",error=f"{type(e).__name__}: {e}")
    (a.out/"audit.json").write_text(json.dumps(z,indent=2)+"\n")
    print("VANVELK_SOURCE_RESULT",json.dumps(z),flush=True)
if __name__=="__main__":main()
