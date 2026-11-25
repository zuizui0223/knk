"""Source-pinned fixed-role response contributor diagnostics; no intervention inference."""
from __future__ import annotations
import argparse,itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_narwani_forecast import grid as pond_grid,CONTRASTS,WEEKS,SPLITS
from benchmarks.score_ahif_kansas_forecast import grid as kansas_grid,EXPERIMENT
from benchmarks.audit_ahif_seed_synthesis import acquire
from benchmarks.score_ahif_seed_synthesis_opposition import read
from benchmarks.score_ahif_kansas_sown_roles import read_mapping
from benchmarks.score_ahif_kansas_group_forecast import gitblob,MAP_BLOB

FIRST=5
RULES=("early","previous","expanding","concurrent")
KS=(1,3)
EPS=1e-12

def top(q,names,k):
    ix=sorted((i for i,v in enumerate(q) if v>EPS),key=lambda i:(-q[i],names[i]))
    return ix[:k]
def overlap(a,b):
    a=np.maximum(a,0);b=np.maximum(b,0)
    if a.sum()<=EPS or b.sum()<=EPS:return None
    return float(np.minimum(a/a.sum(),b/b.sum()).sum())
def mean(vals):
    v=[x for x in vals if x is not None]
    return float(np.mean(v)) if v else None

def group_metrics(d,names):
    q=np.maximum(d,0);gross=float(q.sum());ix=top(q,names,3)
    return {"net":float(d.sum()),"gross_support":gross,"gross_opposing":float(np.maximum(-d,0).sum()),
       "top1":names[ix[0]] if ix else None,"top3":[names[i] for i in ix],
       "top1_share":float(q[ix[0]]/gross) if ix else None,
       "top3_share":float(q[ix].sum()/gross) if ix else None,
       "n_positive_contributors":int(np.sum(q>EPS))}
def delete(test,train,names,k):
    ix=top(train,names,k);gross=float(np.maximum(test,0).sum())
    rest=float(test.sum()-test[ix].sum())
    return {"selected":[names[i] for i in ix],"n_selected":len(ix),
       "capture_fraction":float(np.maximum(test[ix],0).sum()/gross) if gross>EPS else None,
       "random_expected_capture_fraction":len(ix)/len(names) if gross>EPS else None,
       "removed_signed_net":float(test[ix].sum()),"residual_net":rest}

def evaluate(name,unit,dates,source_names,G,full,split_pairs,prior=None,window=()):
    masks=[G[:,i].astype(bool) for i in (0,1)]
    gnames=[[n for n,b in zip(source_names,m) if b] for m in masks];sign=[1.,-1.]
    if len(full)!=15 or len(dates)!=15 or np.any(G.sum(1)!=1):raise ValueError("grid or partition")
    pool=[]
    for t,date in enumerate(dates):
        groups=[]
        for i,m in enumerate(masks):
            d=sign[i]*full[t,m];row=group_metrics(d,gnames[i])
            row["adjacent_support_overlap"]=overlap(sign[i]*full[t-1,m],d) if t else None
            groups.append(row)
        group_delta=(full[t]@G).tolist()
        if prior is not None and not np.allclose(group_delta,prior[t],atol=1e-9,rtol=1e-9):raise ValueError("archived pooled delta drift")
        pool.append({"time":date,"group_delta":group_delta,"opposition":bool(groups[0]["net"]>0 and groups[1]["net"]>0),"groups":groups})
    rows=[]
    for fold,train,held in split_pairs:
        if train.shape!=full.shape or held.shape!=full.shape:raise ValueError("split shape")
        for t in range(FIRST,15):
            group=[]
            for i,m in enumerate(masks):
                a=sign[i]*train[:,m];v=sign[i]*held[t,m];n=gnames[i]
                scores={"early":np.maximum(a[:FIRST],0).mean(0),"previous":np.maximum(a[t-1],0),
                        "expanding":np.maximum(a[:t],0).mean(0),"concurrent":np.maximum(a[t],0)}
                g=group_metrics(v,n)
                g["lag_overlap"]=overlap(a[t-1],v);g["concurrent_overlap"]=overlap(a[t],v)
                g["deletions"]={rule:{str(k):delete(v,scores[rule],n,k) for k in KS} for rule in RULES}
                group.append(g)
            rows.append({"time":dates[t],"fold":fold,"groups":group})
    def summarize(sel):
        if not sel:return None
        z={"n_cells":len(sel),"n_full_opposition":sum(all(q["net"]>0 for q in r["groups"]) for r in sel),"groups":[],"joint":{}}
        for i in (0,1):
            v=[r["groups"][i] for r in sel];common=[a for a in v if a["lag_overlap"] is not None and a["concurrent_overlap"] is not None]
            d={"mean_net":mean([a["net"] for a in v]),"mean_gross_support":mean([a["gross_support"] for a in v]),
               "n_positive_net":sum(a["net"]>0 for a in v),
               "mean_top1_share":mean([a["top1_share"] for a in v]),"mean_top3_share":mean([a["top3_share"] for a in v]),
               "top1_frequency":dict(Counter(a["top1"] for a in v if a["top1"] is not None)),
               "overlap":{"n_common":len(common),"lag":mean([a["lag_overlap"] for a in common]),
                  "concurrent":mean([a["concurrent_overlap"] for a in common]),
                  "temporal_gap":mean([a["concurrent_overlap"]-a["lag_overlap"] for a in common])},"deletions":{}}
            for rule in RULES:
                d["deletions"][rule]={}
                for k in KS:
                    v2=[a["deletions"][rule][str(k)] for a in v];positive=[b for a,b in zip(v,v2) if a["net"]>0]
                    capture=[b for b in v2 if b["capture_fraction"] is not None]
                    d["deletions"][rule][str(k)]={"mean_residual_net":mean([b["residual_net"] for b in v2]),
                       "n_residual_positive":sum(b["residual_net"]>0 for b in v2),
                       "n_full_positive":len(positive),"n_retained_among_full_positive":sum(b["residual_net"]>0 for b in positive),
                       "mean_support_capture":mean([b["capture_fraction"] for b in capture]),
                       "mean_random_capture":mean([b["random_expected_capture_fraction"] for b in capture]),
                       "n_capture_defined":len(capture),"selected_frequency":dict(Counter(n for b in v2 for n in b["selected"]))}
            z["groups"].append(d)
        for rule in RULES:
            z["joint"][rule]={}
            for k in KS:
                key=str(k);both=lambda r:all(g["deletions"][rule][key]["residual_net"]>0 for g in r["groups"])
                fullpos=[r for r in sel if all(q["net"]>0 for q in r["groups"])]
                z["joint"][rule][key]={"n_residual_opposition":sum(both(r) for r in sel),
                  "n_retained_among_full_opposition":sum(both(r) for r in fullpos),
                  "n_full_opposition":len(fullpos)}
        return z
    summary={"name":name,"unit":unit,"n_taxa":len(source_names),"group_sizes":G.sum(0).astype(int).tolist(),
        "dates":dates,"n_folds":len(split_pairs),"pooled":pool,"all_forecast_dates":summarize(rows),
        "known_persistence_window":list(window),"known_window_descriptive":summarize([r for r in rows if r["time"] in window]),
        "by_time":[{"time":t,**summarize([r for r in rows if r["time"]==t])} for t in dates[FIRST:]]}
    if len(rows)!=len(split_pairs)*10:raise ValueError("cell count")
    return summary,rows

def run(root,mapping,prior_pond,prior_kansas):
    pp=json.loads(prior_pond.read_text());pk=json.loads(prior_kansas.read_text())
    if pp["source_run"]!="https://github.com/zuizui0223/knk/actions/runs/45557051796" or pk["source_run"]!="https://github.com/zuizui0223/knk/actions/runs/45543133793":raise ValueError("prior lineage")
    x,y,arms,names,G,H,audit=pond_grid()
    if dict(zip(names,["cyanobacteria" if v[0] else "other_named_algae" if v[1] else "unresolved_other" for v in G]))!=pp["mapping"]:raise ValueError("pond grouping drift")
    summaries={};cells={}
    for a,b in CONTRASTS:
        name=b+"-"+a;folds=[]
        for ia,ib in itertools.product(SPLITS,repeat=2):
            ha=[i for i in range(4) if i not in ia];hb=[i for i in range(4) if i not in ib]
            tr=x[b][:,list(ib)].mean(1)-x[a][:,list(ia)].mean(1)
            te=x[b][:,hb].mean(1)-x[a][:,ha].mean(1)
            fold={"train_control":[arms[a][i] for i in ia],"train_treat":[arms[b][i] for i in ib],"held_control":[arms[a][i] for i in ha],"held_treat":[arms[b][i] for i in hb]}
            folds.append((fold,tr,te))
        prev=pp["contrasts"][name]["persistence"]
        key="pond_"+name
        summaries[key],cells[key]=evaluate(key,"mm^3/L",WEEKS,names,G,x[b].mean(1)-x[a].mean(1),folds,[v["contrast"] for v in prev],(12,13,14,15,16,19) if name=="MD-C" else ())
        if [v["opposition"] for v in summaries[key]["pooled"]]!=[v["opposite"] for v in prev]:raise ValueError("pond sign reproduction")
    if gitblob(mapping)!=MAP_BLOB:raise ValueError("Kansas map blob")
    manifest=acquire(root);states,roster,dup,drift=read(root,manifest)
    if dup[EXPERIMENT] or drift[EXPERIMENT]:raise ValueError("Kansas physical key")
    z,blocks=kansas_grid(states[EXPERIMENT]);names=roster[EXPERIMENT];s,o,u,pre,ex=read_mapping(mapping,names)
    G=np.column_stack((s,o,u)).astype(float);d=z[:,:,1]-z[:,:,0];folds=[]
    for ia in SPLITS:
        ha=[i for i in range(4) if i not in ia]
        folds.append(({"train_blocks":[blocks[i] for i in ia],"held_blocks":[blocks[i] for i in ha]},d[:,list(ia)].mean(1),d[:,ha].mean(1)))
    pr=[[q["pooled"][k] for k in ("sown_list","other_named","unresolved")] for q in pk["by_age"]]
    summaries["kansas"],cells["kansas"]=evaluate("kansas","g/m^2",list(range(1,16)),names,G,d.mean(1),folds,pr,tuple(range(7,16)))
    if [v["opposition"] for v in summaries["kansas"]["pooled"]]!=[v["pooled"]["opposite_named_signs"] for v in pk["by_age"]]:raise ValueError("Kansas sign reproduction")
    return {"status":"SOURCE_PINNED_FIXED_ROLE_CONTRIBUTOR_DIAGNOSTIC_NOT_CAUSAL","summaries":summaries,"cells":cells,
       "source_runs":[pp["source_run"],pk["source_run"]],"limits":["Oriented treatment contrast contributions, not realized abundance turnover or causal removals.","Selection from training half; contemporaneous comparator not a forecast.","All 10 forecast dates primary, known windows outcome-informed descriptive only.","Overlapping halves and repeated dates dependent, no p-values.","Kansas curated conditional roster; Narwani not an order experiment."]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--mapping",type=Path,required=True)
    p.add_argument("--prior-pond",type=Path,required=True);p.add_argument("--prior-kansas",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True);z=run(a.root,a.mapping,a.prior_pond,a.prior_kansas)
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    # Full per-time and split diagnostics in artifact; bounded log for archival.
    q={k:v for k,v in z.items() if k!="cells"}
    for v in q["summaries"].values():v.pop("by_time")
    print("ROLE_CONTRIBUTOR_RESULT",json.dumps(q,allow_nan=False))
if __name__=="__main__":main()
