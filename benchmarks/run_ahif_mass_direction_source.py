"""Source-isolated runner for the already frozen mass/direction factorial estimand.

Acquisition failures in one source must not invalidate other sources.
Only source loading is adapted from the original PR14 loader; the scorer,
physical-unit partition seed, targets and archived per-source replay are unchanged.
"""
import argparse,itertools,json,hashlib
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_cross_source_concentration as src
from benchmarks import score_ahif_replication_scale as rs
from benchmarks import score_ahif_mass_direction_transfer as md

VAL=Path("validation/ahif_v1")
DB=Path("database/ahif_v1")

def receipt(name,run):
    z=json.loads((VAL/name).read_text())
    assert z["source_run"]==f"https://github.com/zuizui0223/knk/actions/runs/{run}"
    return z

def register(*args,**kwargs):
    return rs.previous(*args,**kwargs)

def load(st,a):
    if st=="pond":
        p=receipt("narwani2019_forecast_summary_v1.json",45557051796)
        x,y,arms,names,G,H,audit=src.pond_grid()
        for c,b in src.CONTRASTS:
            key=b+"-"+c
            register("pond",key,src.WEEKS,x[c],x[b],names,G,["cyanobacteria","other_named_algae"],False,
                     [v["contrast"] for v in p["contrasts"][key]["persistence"]],True)
    elif st=="kansas":
        k=receipt("kansas2020_published_sown_role_v1.json",45543133793)
        mapping=DB/"foster2009_curated_name_bridge.csv"
        assert src.gitblob(mapping)==src.MAP_BLOB
        a.root.mkdir(parents=True,exist_ok=True)
        man=src.acquire(a.root);states,roster,dup,drift=src.read(a.root,man)
        assert not dup[src.EXPERIMENT] and not drift[src.EXPERIMENT]
        z,blocks=src.kansas_grid(states[src.EXPERIMENT]);names=roster[src.EXPERIMENT]
        s,o,u,pre,ex=src.read_mapping(mapping,names);G=np.column_stack((s,o,u)).astype(float)
        register(st,"seeded-control",list(range(1,16)),z[:,:,0],z[:,:,1],names,G,["sown_list","other_named"],True,
                 [[v["pooled"][q] for q in ("sown_list","other_named","unresolved")] for v in k["by_age"]],True)
    elif st=="poem":
        po=receipt("poem_field_absolute_compensation_v1.json",45536266086)
        assert hashlib.md5(a.poem.read_bytes()).hexdigest()=="3faa5f7d1619a3beae5941ba763b7a26"
        pp,pa=src.poem_panels(a.poem);gr,names=src.poem_matrices(pp)
        assert names==po["taxon_names"]
        lookup={(v["cohort"],v["age"],v["early"]):v["species_delta"] for v in po["contrasts"]};assert len(lookup)==18
        for c,g in itertools.product(src.COHORTS,src.EARLY):
            G=src.role_matrix(names,g)[:,:3]
            A=np.array([gr[c,t,"S"] for t in src.AGES]);B=np.array([gr[c,t,g] for t in src.AGES])
            register(st,c+"_"+g+"_S",list(src.AGES),A,B,names,G,["early_sown","other_sown","nonsown"],False,
                     [lookup[c,t,g] for t in src.AGES])
    elif st=="field":
        fi=receipt("vansteijn2026_field_opposition_v1.json",45537486831)
        for n,md5 in src.FIELD_MD5.items():assert hashlib.md5((a.field/n).read_bytes()).hexdigest()==md5,n
        fp,fa=src.field_panels(a.field);gr,names,fast,slow,_=src.field_unpack(fp)
        assert names==fi["taxon_names"]
        lookup={(v["site"],v["year"],v["reference"],v["history"]):v["species_delta"] for v in fi["taxon_contrasts"]}
        assert len(lookup)==32
        for s,h,r in itertools.product(src.SITES,src.HISTORIES,src.REFS):
            A=np.array([gr[s,t,r] for t in src.YEARS]);B=np.array([gr[s,t,h] for t in src.YEARS]);G=src.rolemap(names,fast,slow,h)[:,:3]
            register(st,s+"_"+h+"_"+r,list(src.YEARS),A,B,names,G,["early_sown","other_sown","nonsown"],True,
                     [lookup[s,t,r,h] for t in src.YEARS])
    elif st=="mancuso":
        g,names,named,audit=rs.m.source(a.mancuso);mats,rosters,coverage,pairs=rs.m.panels(g,names,named)
        assert set(pairs)==rs.m.EXPECTED
        for r,f in pairs:
            A=mats[r,600];B=mats[r,f]
            register(st,f"{r}:{f}",rs.m.DAYS,A,B,names,np.ones((len(names),1)),["all"],False,B.mean(1)-A.mean(1))
    elif st=="hayashi":
        x,ids,taxa=rs.hy.read(a.hayashi,"otu99");x2,ids2,_=rs.hy.read(a.hayashi,"asv")
        assert ids==ids2 and len(taxa)==331 and sum(map(len,ids.values()))==760
        old=json.loads((VAL/"hayashi2026_direction_summary_20261125.json").read_text())
        assert {f"{s}:{d}":len(v) for (s,d),v in ids.items()}==old["coverage"]
        for s in rs.hy.SOURCES:
            for d in rs.hy.LEVELS[1:]:
                ka,kb=f"hayashi/{s}/1.0",f"hayashi/{s}/{d}"
                rs.register(st,f"{s}:{d}",rs.hy.DAYS,x[s,1.],x[s,d],ka,kb,(ka,kb))
    else:raise ValueError(st)
    assert set(rs.panels)=={st} and len(rs.panels[st])==rs.EXPECTED[st]
    if st!="hayashi":
        ref=json.loads((VAL/"one_lag_fungal_summary_20261125.json").read_text())
        v=rs.ol.alloc.nested([r for r in rs.ol.newrows if r["primary"]],"gain")
        assert set(v)=={st} and abs(v[st]-ref["previous_sources"]["primary"]["by_metric"]["gain"][st])<1e-10

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,choices=md.SOURCES)
    for k in ("root","poem","field","mancuso","hayashi","out"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args();md.selftest();load(a.source,a)
    prior=json.loads((VAL/"replication_scale_summary_20261125.json").read_text())
    r=md.source(a.source,prior)
    z={"status":"RETROSPECTIVE_FACTORIAL_TRANSFER_NOT_META","source":a.source,"panels":len(rs.panels[a.source]),"result":r}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    brief={**r,"comparisons":{k:{j:v for j,v in q.items() if j!="by_panel"} for k,q in r["comparisons"].items()}}
    print("MASS_DIRECTION_TRANSFER_RESULT",a.source,json.dumps(brief,allow_nan=False),flush=True)
if __name__=="__main__":main()
