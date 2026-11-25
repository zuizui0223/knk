"""Technical post-audit harmonization; preserve original score and report both resolution sensitivities."""
import argparse,csv,io,json,re
from pathlib import Path
from collections import Counter
import numpy as np
from benchmarks import score_ahif_diazcolunga2022_pool as base
TAX="data/taxonomy_silva.csv"
base.FILES[TAX]="b174b773a2887ee5fe2341bbfc01a5386ab57285"

def aggregate(data,keys):
    unique={k:i for i,k in enumerate(dict.fromkeys(keys))}
    indices=np.array([unique[k] for k in keys])
    def collapse(v):
        return np.bincount(indices,weights=v,minlength=len(unique))
    out={}
    for m,(p,c,mix) in data.items():
        q=np.stack([collapse(v) for v in p])
        r=np.stack([collapse(v) for v in c])
        mm={k:collapse(v) for k,v in mix.items()}
        assert np.allclose(p.sum(1),q.sum(1)) and np.allclose(c.sum(1),r.sum(1))
        assert all(np.isclose(mix[k].sum(),v.sum()) for k,v in mm.items())
        out[m]=(q,r,mm)
    return out,len(unique)

def brief(results):
    return {k:{**{j:v for j,v in r.items() if j not in ("pairs","source_disjoint_calibration")},"source_disjoint_calibration":{m:{x:y for x,y in z.items() if x!="by_pair"} for m,z in r["source_disjoint_calibration"].items()}} for k,r in results.items()}

def run():
    data,audit=base.load()
    rd=csv.reader(io.StringIO(base.get("data/OTU_table.csv").decode("utf-8-sig")));next(rd)
    seq=[r[0].upper() for r in rd]
    assert len(seq)==audit["n_taxa"] and len(set(seq))==len(seq)
    tr=csv.reader(io.StringIO(base.get(TAX).decode("utf-8-sig")));header=next(tr);tax=list(tr)
    assert len(tax)==len(seq) and all(r[0].upper()==s for r,s in zip(tax,seq))
    assert all(re.fullmatch("[ACGTN]+",s) for s in seq)
    offsets={};support={};parent_idx=set();endpoint_idx=set()
    for m,(p,c,mix) in data.items():
        a=np.flatnonzero(p.sum(0)>0);b=np.flatnonzero(c.sum(0)+np.stack(list(mix.values())).sum(0)>0)
        assert all(len(seq[i])==238 for i in a) and all(len(seq[i]) in (232,233) for i in b)
        parent_idx.update(map(int,a));endpoint_idx.update(map(int,b))
        cts=Counter();matched=set()
        for j in b:
            s=seq[j]
            if len(s)!=233:continue
            for i in a:
                for off in range(6):
                    if seq[i][off:off+233]==s:
                        cts[off]+=1;matched.add(int(j))
        offsets[m]=dict(cts)
        support[m]={"matched_233_endpoint_keys":len(matched),"endpoint_keys":len(b),"mean_mixed_mass_matched":float(np.mean([sum(v[i] for i in matched) for v in mix.values()]))}
    modes={};summaries={}
    valid=all(len(v)==1 and sum(v.values())>0 for v in offsets.values()) and len({next(iter(v)) for v in offsets.values() if v})==1
    result={"status":"POST_AUDIT_REPAIR_NOT_INDEPENDENT_CONFIRMATION","source_commit":base.REF,"offset_matches":offsets,"support":support,"common_segment_eligible":valid}
    if valid:
        off=next(iter(next(iter(offsets.values()))))
        keys=[("raw",i) for i in range(len(seq))]
        for i in parent_idx:keys[i]=("common233",seq[i][off:off+233])
        for i in endpoint_idx:
            if len(seq[i])==233:keys[i]=("common233",seq[i])
        result["common_offset_zero_based"]=off
        modes["common233"]=keys
    else:
        result["common_segment_stop"]="INCONSISTENT_OR_UNIDENTIFIED_OFFSETS"
    for rank in ("Family","Genus"):
        k=header.index(rank)
        invalid={"","na","nan","none","unassigned","unclassified","unknown"}
        keys=[("lineage",*r[1:k+1]) if all(x.strip().lower() not in invalid for x in r[1:k+1]) else ("unresolved",i) for i,r in enumerate(tax)]
        modes[rank.lower()]=keys
    for label,keys in modes.items():
        d,n=aggregate(data,keys)
        summaries[label]={"n_aggregated_keys":n,"results":base.run(d)}
    result["modes"]=summaries
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    r=run();a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(r,indent=2,allow_nan=False)+"\n")
    print("DIAZ_REPAIRED_RESULT",json.dumps({**r,"modes":{k:{"n_aggregated_keys":v["n_aggregated_keys"],"results":brief(v["results"])} for k,v in r["modes"].items()}},allow_nan=False),flush=True)
if __name__=="__main__":main()
