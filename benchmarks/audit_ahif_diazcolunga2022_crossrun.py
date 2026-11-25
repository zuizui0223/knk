"""Audit cross-run taxon support. No biological inference from disjoint raw ESV keys."""
import argparse,csv,io,json,hashlib,re
from pathlib import Path
from collections import Counter,defaultdict
import numpy as np
from benchmarks import score_ahif_diazcolunga2022_pool as base
TAX="data/taxonomy_silva.csv"
base.FILES[TAX]="b174b773a2887ee5fe2341bbfc01a5386ab57285"
def summarize():
    data,audit=base.load()
    rows=csv.reader(io.StringIO(base.get("data/OTU_table.csv").decode("utf-8-sig")));next(rows)
    seq=[r[0] for r in rows]
    rawtax=base.get(TAX).decode("utf-8-sig");tr=csv.reader(io.StringIO(rawtax));h=next(tr);t=list(tr)
    assert len(seq)==len(t) and all(r[0]==s for r,s in zip(t,seq)), "Taxonomy order must exactly match OTU"
    meta=list(csv.DictReader(io.StringIO(base.get("data/metaData.csv").decode("utf-8-sig"))))
    ret={"status":"AUDIT_ONLY","taxonomy_header":h,"sequence_id_lengths":dict(Counter(map(len,seq))),"sequence_keys_acgt_fraction":sum(bool(re.fullmatch("[ACGTNacgtn]+",s)) for s in seq)/len(seq),"taxonomy_first_record":t[0][1:],"medium":{}}
    for m,(p,c,mix) in data.items():
        parent=p.sum(0)>0;endpoint=(c.sum(0)+np.stack(list(mix.values())).sum(0))>0
        pp=[seq[i].upper() for i in np.flatnonzero(parent)]
        ee=[seq[i].upper() for i in np.flatnonzero(endpoint)]
        prefix,exp=base.EXP[m]
        ms=[r for r in meta if r["Experiment"] in {f"{prefix}-T12-X{x}" for x in base.MAP[m]} or (r["Experiment"]==exp and r["Carbon"]==m)]
        fields=["Sequencing_Run","LinkerPrimerSequence","ReversePrimer"]
        runs={}
        for role,z in (("pre",[r for r in ms if r["Transfer"]=="12"]),("endpoint",[r for r in ms if r["Transfer"]=="7"])):
            runs[role]={k:dict(Counter(r.get(k,"") for r in z)) for k in fields}
        exact=[any(a in b or b in a for a in pp) for b in ee] if all(re.fullmatch("[ACGTN]+",s) for s in pp+ee) else []
        r={"parent_taxa":int(parent.sum()),"endpoint_taxa":int(endpoint.sum()),"shared_raw_keys":int((parent&endpoint).sum()),"parent_length":dict(Counter(map(len,pp))),"endpoint_length":dict(Counter(map(len,ee))),"run_metadata":runs,"endpoint_exact_substring_matches":sum(exact),"endpoint_mass_substring_supported":None}
        if exact:
            ix=np.flatnonzero(endpoint);masks=np.zeros(len(seq),bool);masks[ix]=exact
            r["endpoint_mass_substring_supported"]=float(np.mean([v[masks].sum() for v in mix.values()]))
        r["taxonomic_rank_support"]={}
        for k in range(1,len(h)):
            keys=[tuple(v[1:k+1]) if all(v[j].strip() not in ("","NA","NaN","None","Unassigned","unclassified") for j in range(1,k+1)) else None for v in t]
            roster={v for v,x in zip(keys,parent) if x and v is not None}
            endkeys={v for v,x in zip(keys,endpoint) if x and v is not None}
            mask=np.array([v in roster for v in keys]);classified=np.array([v is not None for v in keys])
            r["taxonomic_rank_support"][h[k]]={"parent_keys":len(roster),"endpoint_keys":len(endkeys),"shared_keys":len(roster&endkeys),"endpoint_mass_supported":float(np.mean([v[mask].sum() for v in mix.values()])),"endpoint_mass_classified":float(np.mean([v[classified].sum() for v in mix.values()])),"parent_mass_classified":float(np.mean(p[:,classified].sum(1)))}
        ret["medium"][m]=r
    return ret
def main():
    ap=argparse.ArgumentParser();ap.add_argument("--out",type=Path,required=True);a=ap.parse_args();r=summarize()
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print("DIAZ_MARKER_AUDIT",json.dumps(r,allow_nan=False),flush=True)
if __name__=="__main__":main()
