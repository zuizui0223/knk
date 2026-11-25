"""Narwani original source metadata and completeness; no treatment scores."""
from __future__ import annotations
import argparse,hashlib,io,json,math
from collections import Counter,defaultdict
from pathlib import Path
import openpyxl
from benchmarks.audit_ahif_independent_role_sources import fetch
URL="https://zenodo.org/records/4941856/files/Narwani_et_al_2019_Proceedings_B_data.xlsx?download=1"
MD5="c5f52c568a1119429093d3db1603e59c"
SHA="6242f9c3a92cde7f451066ceebf6ae6c9cc6f924d5ed9bc0a7f11de7b9beb3c3"
def load_source():
    b=fetch(URL,2_000_000)
    if hashlib.md5(b).hexdigest()!=MD5 or hashlib.sha256(b).hexdigest()!=SHA:raise ValueError("source hash")
    w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    raw=list(w.worksheets[0].values);m=list(w["Metadata"].values);w.close()
    names=list(raw[0]);r=[dict(zip(names,x)) for x in raw[1:]]
    meta=[dict(zip(m[0],x)) for x in m[1:]]
    if len(r)!=660 or len(names)!=170 or len(set(names))!=170 or len(meta)!=170:raise ValueError("source grain")
    return names,r,meta
def isnum(v):return isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v)
def run():
    names,r,meta=load_source()
    # schema only; distinguish no sampling from true zeros.
    cols=[]
    for name in names:
        vals=[x[name] for x in r]
        m=next((v for v in meta if v.get("Column")==name),None)
        cols.append({"name":name,"metadata":m,"n_numeric":sum(isnum(x) for x in vals),
          "n_blank":sum(x is None or x=="" for x in vals),"n_non_numeric_nonblank":sum(not isnum(x) and x is not None and x!="" for x in vals)})
    designs={}
    for x in r:
        k=x["Pond"];v=tuple(x[t] for t in ("Block","Rep","Dreisse","Macrophytes","Nutrients","TreatNR","Treatment"))
        if k in designs and designs[k]!=v:raise ValueError("allocation drift")
        designs[k]=v
    if len(designs)!=20 or len({(x["Pond"],x["Week"]) for x in r})!=660:raise ValueError("source design")
    pref=defaultdict(list)
    for name in names:
        if "_" in name:pref[name.split("_")[0]].append(name)
    cover={}
    for wk in sorted({x["Week"] for x in r}):
        sub=[x for x in r if x["Week"]==wk]
        cover[str(wk)]={p:{"n_cols":len(c),"complete_numeric_rows":sum(all(isnum(x[v]) for v in c) for x in sub),"any_numeric_rows":sum(any(isnum(x[v]) for v in c) for x in sub)} for p,c in pref.items()}
    return {"status":"NARWANI_METADATA_AND_COVERAGE_ONLY","source_md5":MD5,"source_sha256":SHA,
      "n_states":len(r),"n_columns":len(names),"design_columns":["Block","Rep","Dreisse","Macrophytes","Nutrients","TreatNR","Treatment"],
      "designs":{k:list(v) for k,v in designs.items()},"columns":cols,"week_numeric_coverage":cover}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run();a.out.write_text(json.dumps(z,indent=2,default=str,allow_nan=False)+"\n")
    print("NARWANI_DESIGN",json.dumps({k:v for k,v in z.items() if k!="columns"},default=str))
    for c in z["columns"]:print("NARWANI_COLUMN",json.dumps(c,default=str))
if __name__=="__main__":main()
