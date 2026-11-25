"""Descriptive matched-stratum opposition in curated seed-addition synthesis."""
from __future__ import annotations
import argparse,hashlib,json,math
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from benchmarks.audit_ahif_seed_synthesis import acquire,load,number,audit

HASHES={"plot":("715c3ed27e6c11d41185a02a2ab07b90","0455c64afd914af5c7f85273ae950104ece31bf94fe9737f1a172aa6d9f437a3"),
        "species":("420281357484766c25df281a4cd2ed45","d9f09ea90b3bac6e9913000d93acf24e1c8b0293bebac563477c9e90a9fd590e")}
SPLITS=256
SEED=20261125
MODES=("direct","indirect")

def num(v):
    x=number(v)
    if x is None or x<0:raise ValueError("invalid nonnegative numeric value")
    return x

def metrics(v):
    v=np.asarray(v,float);g=float(abs(v).sum());n=float(abs(v.sum()));c=max(0.,g-n)
    return {"gross":g,"net":n,"cancellation":c,"fraction":c/g if g>1e-12 else None,
            "signed_total":float(v.sum())}

def read(root,manifest):
    for k,(md5,sha) in HASHES.items():
        b=(root/(k+".csv")).read_bytes()
        if hashlib.md5(b).hexdigest()!=md5 or hashlib.sha256(b).hexdigest()!=sha:raise ValueError("source byte drift")
    _,p,pe=load(root/"plot.csv");_,s,se=load(root/"species.csv")
    if len(p)!=3393 or len(s)!=58426 or (pe,se)!=("cp1252","cp1252"):raise ValueError("source count or encoding")
    idx={r["unique.id_"]:r for r in p}
    if len(idx)!=3393:raise ValueError("duplicate plot state")
    roster=defaultdict(set);groups=defaultdict(list);unique=set()
    for x in s:
        key=x["unique.id_"];q=idx[key]
        for field in ("Experiment_","site","block","plot","yr.trt"):
            if x[field]!=q[field]:raise ValueError("source join field mismatch "+field)
        if not math.isclose(num(x["seed.rich"]),num(q["seed.rich"])) or not math.isclose(num(x["biomass.plot"]),num(q["biomass.plot"]),rel_tol=1e-9,abs_tol=1e-9):
            raise ValueError("source numeric join")
        if (x["trt"].rstrip("s")!=q["trt"].rstrip("s")):raise ValueError("source treatment join")
        a=(key,x["species"])
        if a in unique or not x["species"].strip():raise ValueError("duplicate/missing taxon state")
        unique.add(a);roster[q["Experiment_"]].add(x["species"]);groups[key].append(x)
    if set(idx)!=set(groups) or len(roster)!=12:raise ValueError("source join coverage")
    roster={e:sorted(v) for e,v in roster.items()};states=defaultdict(list)
    pkey=Counter();assign=defaultdict(set)
    for id,q in idx.items():
        e=q["Experiment_"];a=num(q["yr.trt"]);dose=num(q["seed.rich"]);tr=q["trt"]
        if int(a)!=a or int(dose)!=dose or (tr=="Control")!=(dose==0) or tr not in ("Control","Seeds"):
            raise ValueError("source age/arm")
        block=(q["site"],q["block"]);physical=(e,*block,q["plot"])
        pkey[physical,int(a)]+=1;assign[physical].add((tr,dose))
        rows=groups[id];modes={v["biomass.measure"] for v in rows}
        if not modes<=set(MODES):raise ValueError("unknown biomass.measure")
        mode=next(iter(modes)) if len(modes)==1 else "mixed"
        x=np.zeros(len(roster[e]));col={s:i for i,s in enumerate(roster[e])}
        for v in rows:x[col[v["species"]]]=num(v["biomass.sp"])
        states[e].append({"age":int(a),"dose":int(dose),"tr":tr,"block":block,
           "mode":mode,"x":x,"total":num(q["biomass.plot"]),"rich":num(q["rich.plot"])})
    # Some author plot identifiers may repeat inside an age because code denotes
    # a grouping rather than a physical plot. Don't silently count as independent.
    duplicates={e:sum(n-1 for (p,a),n in pkey.items() if p[0]==e and n>1) for e in states}
    drift={e:sum(len(v)>1 for p,v in assign.items() if p[0]==e) for e in states}
    return states,roster,duplicates,drift

def summarize_strata(rows):
    delta=np.array([v["delta"] for v in rows]);d=delta.mean(0);m=metrics(d)
    within=[metrics(v) for v in delta]
    a={k:float(np.mean([v[k] for v in within])) for k in ("gross","net","cancellation")}
    a["fraction"]=a["cancellation"]/a["gross"] if a["gross"]>1e-12 else None
    bysite=defaultdict(list)
    for v in rows:bysite[v["block"][0]].append(v["delta"])
    se=metrics(np.mean([np.mean(v,axis=0) for v in bysite.values()],axis=0))
    return m,a,se

def crossfit(rows):
    x=np.array([r["delta"] for r in rows]);n=len(x)
    if n<4:return {"status":"STOP_FEWER_THAN_FOUR_BLOCK_GROUPS","n_blocks":n}
    # Roster-based deterministic ordering and seed; same roster across years
    # receives same partitions. This is a descriptive score only.
    ids=[r["block"] for r in rows]
    rng=np.random.default_rng(SEED+int(hashlib.sha256(json.dumps(ids).encode()).hexdigest()[:8],16))
    g=[];net=[]
    for _ in range(SPLITS):
        p=rng.permutation(n)
        for a,b in ((p[:n//2],p[n//2:]),(p[n//2:],p[:n//2])):
            d=x[a].mean(0);v=x[b].mean(0)
            g.append(float(np.sign(d)@v));net.append(float(np.sign(d.sum())*v.sum()))
    return {"status":"DESCRIPTIVE_HELD_BLOCK_SIGN_SCORE","n_blocks":n,"n_splits":SPLITS,
            "signed_gross":float(np.mean(g)),"signed_net":float(np.mean(net)),
            "gap":float(np.mean(g)-np.mean(net))}

def score(states,roster,duplicates,drift):
    out=[]
    for e,ss in sorted(states.items()):
        counts=Counter(v["mode"] for v in ss);rows=[]
        audit={"experiment":e,"n_states":len(ss),"n_taxa":len(roster[e]),
               "measurement_mode_counts":dict(counts),"duplicate_physical_key_age_rows":duplicates[e],
               "physical_keys_with_treatment_or_dose_drift":drift[e]}
        if duplicates[e]>0:
            out.append({**audit,"status":"STOP_AMBIGUOUS_PHYSICAL_KEY_AGE","comparisons":[]})
            continue
        for age in sorted({v["age"] for v in ss}):
            for mode in MODES:
                here=[v for v in ss if v["age"]==age and v["mode"]==mode]
                if not here:continue
                by=defaultdict(list)
                for v in here:by[v["block"]].append(v)
                rs=[];excluded=0;dose=set()
                for block,vs in sorted(by.items()):
                    c=[v for v in vs if v["tr"]=="Control"];t=[v for v in vs if v["tr"]=="Seeds"]
                    if not c or not t:excluded+=len(vs);continue
                    dose.update(v["dose"] for v in t)
                    delta=np.mean([v["x"] for v in t],axis=0)-np.mean([v["x"] for v in c],axis=0)
                    total=float(np.mean([v["total"] for v in t])-np.mean([v["total"] for v in c]))
                    rich=float(np.mean([v["rich"] for v in t])-np.mean([v["rich"] for v in c]))
                    rs.append({"block":block,"delta":delta,"reported_total_delta":total,
                               "reported_richness_delta":rich,"n_control":len(c),"n_seeded":len(t)})
                if not rs:
                    rows.append({"age":age,"mode":mode,"status":"STOP_NO_MATCHED_ARM_STRATUM",
                                 "states":len(here),"excluded_unmatched_states":excluded})
                    continue
                m,w,se=summarize_strata(rs)
                nt=sum(v["n_seeded"] for v in rs);nc=sum(v["n_control"] for v in rs)
                total=float(np.mean([v["reported_total_delta"] for v in rs]))
                rich=float(np.mean([v["reported_richness_delta"] for v in rs]))
                cf=({"status":"STOP_TREATMENT_DRIFT_PHYSICAL_KEYS","n_blocks":len(rs)}
                    if drift[e]>0 else crossfit(rs))
                rows.append({"age":age,"mode":mode,"status":"SCORED_DESCRIPTIVE",
                    "n_matched_strata":len(rs),"n_total_strata":len(by),
                    "n_control_states":nc,"n_seeded_states":nt,"excluded_unmatched_states":excluded,
                    "seeded_richness_levels":sorted(dose),
                    "pooled_species_metrics":m,"mean_within_stratum_metrics":w,
                    "equal_site_weighted_species_metrics":se,
                    "reported_total_delta":total,"reported_richness_delta":rich,
                    "species_sum_minus_reported_total_delta":m["signed_total"]-total,
                    "crossfit":cf,"n_matched_sites":len({r["block"][0] for r in rs})})
        for mode in MODES:
            good=[v for v in rows if v["mode"]==mode and v["status"]=="SCORED_DESCRIPTIVE"]
            if good:
                last=max(v["age"] for v in good)
                for v in good:v["final_matched_age"]=v["age"]==last
        out.append({**audit,"status":"AUDITED_DESCRIPTIVE","comparisons":rows})
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True)
    manifest=acquire(a.root)
    states,roster,dup,drift=read(a.root,manifest)
    ex=score(states,roster,dup,drift)
    z={"status":"CURATED_SYNTHESIS_DESCRIPTIVE_NOT_RAW_META_OR_RECRUITMENT_IDENTITY","source_manifest":manifest,
       "n_experiment_labels":len(ex),"n_plot_states":sum(x["n_states"] for x in ex),"experiments":ex,
       "limits":["Not independent raw source audits; possible source family overlap.",
                 "Indirect biomass allocated from cover, can mechanically produce opposition.",
                 "No seed membership, cannot identify recruited taxa.",
                 "Matched-stratum state differences descriptive; no causal or p-value claims.",
                 "Held-block crossfit is descriptive, not independent replication."]}
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SEED_OPPOSITION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
