"""Bounded original-source screening for independent experimental role forecasts."""
from __future__ import annotations
import argparse,csv,hashlib,io,json,re,urllib.request,urllib.error
from pathlib import Path
from collections import Counter
from urllib.parse import urlparse
import openpyxl

def fetch(url,cap):
    if urlparse(url).scheme!="https":raise ValueError("public HTTPS only")
    req=urllib.request.Request(url,headers={"User-Agent":"AHIF academic source reproducibility","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=40) as r:
        if urlparse(r.geturl()).scheme!="https":raise ValueError("non HTTPS redirect")
        b=r.read(cap+1)
    if len(b)>cap:raise ValueError("file cap")
    return b

def workbook(b):
    w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    z=[]
    for s in w.worksheets:
        # Header and field cardinalities for treatment/time/id, not response values.
        first=list(s.iter_rows(min_row=1,max_row=min(6,s.max_row),values_only=True))
        heads=[[str(v)[:100] if v is not None else "" for v in row[:60]] for row in first[:2]]
        labels=[str(v or "") for v in first[0]] if first else []
        fields={i:Counter() for i,l in enumerate(labels) if re.search(r"pond|tank|mesocosm|trt|treat|date|day|week|time|species|tax|genus|group|rep|block|nutrient|mussel|plant|pulse",l,re.I)}
        for row in s.iter_rows(min_row=2,values_only=True):
            for i,c in fields.items():c[str(row[i] if i<len(row) else None)]+=1
        z.append({"sheet":s.title,"rows":s.max_row,"columns":s.max_column,"first_header_row":heads[0] if heads else [],"id_fields":{labels[i]:{"n_levels":len(c),"values":dict(c) if len(c)<=40 else None} for i,c in fields.items()}})
    w.close()
    return z

def figshare():
    a=json.loads(fetch("https://api.figshare.com/v2/collections/3300440/articles?page_size=100",2_000_000))
    if not isinstance(a,list):raise ValueError("collection schema")
    items=[]
    for item in a:
        ident=int(item["id"]);m=json.loads(fetch("https://api.figshare.com/v2/articles/"+str(ident),2_000_000))
        files=[]
        for f in m.get("files",[]):
            name=f["name"];q={"id":f["id"],"name":name,"size":f["size"],"md5":f.get("computed_md5") or f.get("supplied_md5")}
            # Only candidate structured tables/text, or small ZIP manifest deferred.
            if name.lower().endswith((".csv",".tsv",".txt",".xlsx",".xls",".zip")) and int(f["size"])<=5_000_000:
                u=f.get("downloadUrl") or f.get("download_url")
                if not u:raise ValueError("missing download URL")
                b=fetch(u,5_000_000)
                if len(b)!=int(f["size"]) or hashlib.md5(b).hexdigest()!=q["md5"]:raise ValueError("Figshare byte mismatch")
                q["sha256"]=hashlib.sha256(b).hexdigest()
                if name.lower().endswith(".xlsx"):q["sheets"]=workbook(b)
                elif name.lower().endswith(".zip"):
                    import zipfile
                    with zipfile.ZipFile(io.BytesIO(b)) as z:q["members"]=[{"name":i.filename,"bytes":i.file_size} for i in z.infolist()[:100]]
                else:
                    txt=b.decode("utf-8-sig",errors="replace")
                    q["n_lines"]=len(txt.splitlines());q["first_header"]=txt.splitlines()[0][:400] if txt else ""
            files.append(q)
        items.append({"id":ident,"title":m.get("title"),"doi":m.get("doi"),"description":re.sub("<[^>]+>"," ",str(m.get("description","")))[:3500],"files":files})
    return {"status":"PUBLIC_COLLECTION_AUDITED_NOT_SCORED","n_items":len(items),"items":items}

def narwani():
    b=fetch("https://zenodo.org/records/4941856/files/Narwani_et_al_2019_Proceedings_B_data.xlsx?download=1",2_000_000)
    if hashlib.md5(b).hexdigest()!="c5f52c568a1119429093d3db1603e59c":raise ValueError("Narwani source MD5")
    return {"status":"SOURCE_BYTES_VERIFIED_SCHEMA_ONLY","source_record":4941856,"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b),"sheets":workbook(b)}

def zenmeta(rec):
    z=json.loads(fetch(f"https://zenodo.org/api/records/{rec}",4_000_000))
    if int(z["id"])!=rec:raise ValueError("Zenodo record identity")
    m=z["metadata"]
    return {"status":"PUBLIC_MANIFEST_ONLY","id":rec,"title":m.get("title"),"doi":m.get("doi"),"description":re.sub("<[^>]+>"," ",str(m.get("description","")))[:4500],
      "files":[{"name":f["key"],"size":f["size"],"checksum":f["checksum"]} for f in z.get("files",[])]}

def run():
    result={}
    for k,f in (("jiang_pat2008",figshare),("narwani2019",narwani),("hayashi2026",lambda:zenmeta(19396153)),("wojcik",lambda:zenmeta(14773023))):
        try:result[k]=f()
        except urllib.error.HTTPError as e:result[k]={"status":"STOP_HTTP","code":e.code}
        except (ValueError,OSError,KeyError,TypeError) as e:result[k]={"status":"STOP_GATE","reason":str(e)[:400]}
    return {"status":"INDEPENDENT_SOURCE_SCREEN_NOT_ECOLOGICAL_SCORING","sources":result}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run();a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n");print("INDEPENDENT_ROLE_SOURCE_AUDIT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
