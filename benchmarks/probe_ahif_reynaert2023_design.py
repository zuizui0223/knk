"""Inspect exact plant workbook metadata and design only; no biomass effect scoring."""
import argparse,json,hashlib,io,collections
from pathlib import Path
import openpyxl
from benchmarks.probe_ahif_context_cross_taxa import get
URL="https://zenodo.org/records/10167759/files/Reynaertetal_Oikos_2023.xlsx?download=1"
SHA="7f51f1c3cfbd0400d7d60097df00fbf91dd60e3d98dcc5b49fcea5ce2bbd5d54"
def load():
    b,_=get(URL,118245)
    assert len(b)==118145 and hashlib.sha256(b).hexdigest()==SHA
    w=openpyxl.load_workbook(io.BytesIO(b),read_only=True,data_only=True)
    s=w["BM_per_species"];rows=list(s.values);h=rows[0]
    assert len(h)==10 and len(rows)==2307
    nonempty=[(i,x) for i,x in enumerate(rows[1:],2) if any(v is not None for v in x)]
    incomplete=[{"excel_row":i,"metadata":{str(k):str(v) for k,v in zip(h[:8],x[:8])},"nonempty_outcome_columns":sum(v is not None for v in x[8:])} for i,x in nonempty if x[2] is None or x[5] is None]
    print("REYNAERT_ROW_GRAIN",json.dumps({"nonempty_rows":len(nonempty),"blank_rows":len(rows)-1-len(nonempty),"incomplete_keys":incomplete[:25],"n_incomplete":len(incomplete)}),flush=True)
    assert not incomplete,"STOP incomplete nonempty sample keys"
    r=[dict(zip(h,x)) for i,x in nonempty]
    return w,r
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    w,r=load()
    readme=[[str(x) if x is not None else None for x in row] for row in w["README_Reynaertetal_2023_Oikos"].values]
    groups=collections.defaultdict(list)
    for x in r:groups[(x["Fati"],x["Pot"])].append(x)
    composition=collections.Counter();treatments=collections.Counter();dups=[];bad=[]
    for (f,pot),g in groups.items():
        sp=[x["Species"] for x in g];com="+".join(sorted(set(sp)))
        if len(set(sp))!=len(sp):dups.append([f,pot,sp])
        if len(set((x["Treatment"],x["Days"]) for x in g))!=1:bad.append([f,pot])
        composition[(str(f),str(g[0]["Treatment"]),str(g[0]["Days"]),com)]+=1
        treatments[(str(f),str(g[0]["Treatment"]),str(g[0]["Days"]))]+=1
    result={"status":"SOURCE_DESIGN_ONLY","source_sha256":SHA,"n_rows":len(r),"n_fati_pots":len(groups),"readme":readme,
      "metadata_values":{k:sorted(set(str(x[k]) for x in r)) for k in ("Fati","Treatment","Fctgrp","Species","Letter","Days")},
      "species_per_pot":dict(collections.Counter(len(g) for g in groups.values())),
      "dup_species_pots":dups[:20],"inconsistent_treatment_pots":bad[:20],
      "treatments":[list(k)+[v] for k,v in sorted(treatments.items())],
      "composition":[list(k)+[v] for k,v in sorted(composition.items())],
      "sample_metadata":[{k:str(x[k]) for k in ("Date","Fati","Pot","Treatment","Fctgrp","Species","Letter","Days")} for x in r[:8]],
      "missing_corrected_biomass_rows":sum(x["Dry weight (g) (corrected for # of individuals)"] is None for x in r)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2)+"\n")
    # Print metadata and membership, not biomass.
    brief={**result,"composition":result["composition"][:100],"n_composition_cells":len(composition)}
    print("REYNAERT_DESIGN_RESULT",json.dumps(brief),flush=True)
if __name__=="__main__":main()
