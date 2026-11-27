"""Audit publisher's experimental data workbook bytes and structure without ecological scoring."""
import argparse,hashlib,io,json,urllib.request,zipfile,xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
BASE="https://media.springernature.com/original/springer-static/esm/art%3A10.15252%2Fmsb.20178157/MediaObjects/"
FILES={"EV1":"44320_2018_BFMSB178157_MOESM2_ESM.xlsx","EV3":"44320_2018_BFMSB178157_MOESM4_ESM.xlsx"}
NS={"m":"http://schemas.openxmlformats.org/spreadsheetml/2006/main","r":"http://schemas.openxmlformats.org/officeDocument/2006/relationships","p":"http://schemas.openxmlformats.org/package/2006/relationships"}
def get(name):
    req=urllib.request.Request(BASE+name,headers={"User-Agent":"AHIF academic source audit","Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,application/octet-stream"})
    with urllib.request.urlopen(req,timeout=70) as r:b=r.read(40_000_001)
    assert len(b)<=40_000_000 and b.startswith(b"PK"),"Not a bounded XLSX"
    return b
def schema(b):
    with zipfile.ZipFile(io.BytesIO(b)) as z:
        assert len(z.infolist())<500 and sum(v.file_size for v in z.infolist())<200_000_000
        x=ET.fromstring(z.read("xl/workbook.xml"))
        rr=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rels={v.get("Id"):v.get("Target") for v in rr}
        ss=[]
        if "xl/sharedStrings.xml" in z.namelist():
            for v in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si",NS):ss.append("".join(t.text or "" for t in v.findall(".//m:t",NS)))
        result=[]
        for s in x.findall("m:sheets/m:sheet",NS):
            target=rels[s.get("{"+NS["r"]+"}id")]
            path=target.lstrip("/") if target.startswith("/") else "xl/"+target
            q=ET.fromstring(z.read(path));rows=q.findall("m:sheetData/m:row",NS);typ=Counter();formulas=0;cells=0;preview=[]
            for i,r in enumerate(rows):
                vals=[]
                for c in r.findall("m:c",NS):
                    v=c.find("m:v",NS);t=c.get("t","n");f=c.find("m:f",NS);typ[t]+=1;cells+=1;formulas+=f is not None
                    val=None if v is None else v.text
                    if t=="s" and val is not None:val=ss[int(val)]
                    elif t=="inlineStr":val="".join(t.text or "" for t in c.findall(".//m:t",NS))
                    if i<14 and val is not None and len(vals)<30:vals.append([c.get("r"),str(val)[:160] if t in ("s","str","inlineStr") else ("numeric" if t=="n" else val)])
                if i<14:preview.append(vals)
            result.append({"sheet":s.get("name"),"rows":len(rows),"cells":cells,"cell_types":dict(typ),"formulas":formulas,"header_preview":preview})
        return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    result={"status":"SOURCE_ONLY_NO_SCORES","publication":"10.15252/msb.20178157","files":{}}
    for key,name in FILES.items():
        try:
            b=get(name);result["files"][key]={"status":"VERIFIED_XLSX","url":BASE+name,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),"sheets":schema(b)}
        except Exception as e:result["files"][key]={"status":"STOP_SOURCE_OR_SCHEMA","url":BASE+name,"error":type(e).__name__+": "+str(e)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+"\n")
    print("VENTURELLI_SOURCE_RESULT",json.dumps(result),flush=True)
if __name__=="__main__":main()
