"""Pinned author dataset schema audit only; no biological outcome scoring."""
import argparse,csv,hashlib,io,json,urllib.request,zipfile
from pathlib import Path,PurePosixPath
from collections import Counter,defaultdict
REF="b07cf57f837b8b6041031f76625792dce4e0ed64"
BASE="https://raw.githubusercontent.com/DoranG1/dose-dependent-colonization/"+REF+"/"
FILES={
"sequencing/communityCoalescence/data/intermediates/mixtureDataframe.zip":"eabfd7ff302e19607237649c73eccb7f6a698a25",
"sequencing/communityCoalescence/data/210819-e0015-e0020-16S-samplesheet.txt":"bc61d1bcc0ea65d5e68ccbec909744a210ad20a5",
"sequencing/parentCommunities/data/intermediates/e0012-mixtureDataframe.txt":"d57e4efb9c49b5f2127f14c61011394e005f60a2",
"sequencing/parentCommunities/data/20210201-e0012-samplesheet.txt":"d79a9c29dd65adaed52d2a304d0749dd82fa8369"
}
def get(path):
    with urllib.request.urlopen(urllib.request.Request(BASE+path,headers={"User-Agent":"AHIF source audit"}),timeout=90) as r:b=r.read(20_000_001)
    assert len(b)<=20_000_000 and hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==FILES[path],path
    return b
def schema(b):
    s=b.decode("utf-8-sig",errors="strict")
    reader=csv.reader(io.StringIO(s),delimiter="\t")
    h=next(reader);n=0; widths=Counter(); profiles=defaultdict(set); count=Counter()
    # Source-only metadata fields. No abundance values emitted.
    selected=[k for k in ("mixtureType","passage","ratio","combo","line","group","biosample","biosample1","biosample2","media","inoculationReplicate","inoculationreplicate","extractionReplicate","PCRreplicate","sample","filename","label") if k in h]
    combos=set();taxa=set();taxcol=next((k for k in ("OTU","ASV") if k in h),None)
    for row in reader:
        n+=1;widths[len(row)]+=1
        if len(row)!=len(h):continue
        for k in selected:profiles[k].add(row[h.index(k)])
        if taxcol:taxa.add(row[h.index(taxcol)])
        k=tuple(row[h.index(k)] for k in selected)
        combos.add(k)
        if "mixtureType" in h and "passage" in h:count[(row[h.index("mixtureType")],row[h.index("passage")])]+=1
    return {"header":h,"rows":n,"widths":dict(widths),"metadata_values":{k:sorted(v) if len(v)<=60 else {"distinct":len(v)} for k,v in profiles.items()},"metadata_unique_combinations":len(combos),"taxon_id_count":len(taxa),"taxon_id_character_lengths":dict(Counter(map(len,taxa))),"rows_by_type_passage":{"|".join(k):v for k,v in count.items()}}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    z={"status":"GOLDMAN_SOURCE_SCHEMA_ONLY","commit":REF,"publication":"10.1073/pnas.2322440122","files":[]}
    for name,sha in FILES.items():
        f={"path":name,"git_blob":sha}
        try:
            b=get(name);f.update(bytes=len(b),sha256=hashlib.sha256(b).hexdigest())
            if name.endswith(".zip"):
                with zipfile.ZipFile(io.BytesIO(b)) as q:
                    assert sum(i.file_size for i in q.infolist())<=150_000_000 and len(q.infolist())<=20
                    f["members"]=[]
                    for i in q.infolist():
                        p=PurePosixPath(i.filename);assert not p.is_absolute() and ".." not in p.parts
                        if not i.is_dir():f["members"].append({"name":i.filename,"schema":schema(q.read(i))})
            else:f["schema"]=schema(b)
            f["status"]="VERIFIED"
        except Exception as e:f.update(status="STOP",error=type(e).__name__+": "+str(e))
        z["files"].append(f)
    z["n_verified"]=sum(f["status"]=="VERIFIED" for f in z["files"])
    (a.out/"audit.json").write_text(json.dumps(z,indent=2)+"\n")
    print("GOLDMAN_SOURCE_RESULT",json.dumps(z),flush=True)
if __name__=="__main__":main()
