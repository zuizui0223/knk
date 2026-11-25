"""Cross-source forecast of gross signed-response concentration."""
from __future__ import annotations
import argparse,itertools,json,hashlib
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_narwani_forecast import grid as pond_grid,CONTRASTS,WEEKS
from benchmarks.score_ahif_kansas_forecast import grid as kansas_grid,EXPERIMENT
from benchmarks.score_ahif_kansas_group_forecast import gitblob,MAP_BLOB
from benchmarks.audit_ahif_seed_synthesis import acquire
from benchmarks.score_ahif_seed_synthesis_opposition import read
from benchmarks.score_ahif_kansas_sown_roles import read_mapping
from benchmarks.score_poem_field_storage import source_panels as poem_panels
from benchmarks.score_ahif_poem_compensation import matrices as poem_matrices,role_matrix,COHORTS,AGES,EARLY
from benchmarks.score_ahif_field_opposition import unpack as field_unpack,rolemap,SITES,YEARS,HISTORIES,REFS
from benchmarks.score_ahif_vansteijn2026_field import source_panels as field_panels
from benchmarks.probe_ahif_vansteijn2026_field import MD5 as FIELD_MD5
EPS=1e-12
MODELS={"intercept":(), "richness":(0,1,2,3), "dominance":(0,1,2,3,4,5),
        "history":(0,1,2,3,6,7,8), "history_dominance":tuple(range(9))}
SOURCES=("kansas","pond","poem","field")
LAM=.1

def profile(v):
    n=int((v>EPS).sum());s=float(v.sum())
    if s<=EPS:return [0,1,0.,0.]
    p=v/s
    return [n,0,float(p.max()),float(p@p)]
def support(d):
    v=np.maximum(d,0);v=np.where(v>EPS,v,0.);s=float(v.sum())
    if s<=EPS:return None,0,None
    p=v/s
    return float(p.max()),int((v>0).sum()),float(1/(p@p))
def feats(a,b,t,sign,mode):
    # a/b: train-only time x taxon. t's data must not be read.
    u=a[t-1] if mode=="control" else (a[t-1]+b[t-1])/2
    n,absent,top,sim=profile(u)
    h,m,_=support(sign*(b[t-1]-a[t-1]))
    return [float(np.log(len(u))),float(np.log1p(n)),float(absent),float(sign),
            top,sim,0. if h is None else h,float(np.log1p(m)),float(h is None)]
def splits(n,paired):
    ii=list(itertools.combinations(range(n),2))
    jj=[(a,a) for a in ii] if paired else list(itertools.product(ii,repeat=2))
    for a,b in jj:
        c=tuple(i for i in range(n) if i not in a);d=tuple(i for i in range(n) if i not in b)
        yield a,b,c,d
        if n!=4:yield c,d,a,b

def panel(st,key,dates,a,b,names,G,roles,paired,prior,prior_is_group=False):
    if a.shape!=b.shape or a.ndim!=3 or a.shape[0]!=len(dates) or a.shape[2]!=len(names):raise ValueError("source panel shape")
    if np.any(a<0) or np.any(b<0) or not np.isfinite(a).all() or not np.isfinite(b).all() or not np.all(G.sum(1)==1):raise ValueError("source panel values")
    p=b.mean(1)-a.mean(1);q=p@G if prior_is_group else p
    if not np.allclose(q,prior,rtol=1e-9,atol=1e-9):raise ValueError("prior source contrast drift "+key)
    rows=[];n=a.shape[1]
    for fold,(ia,ib,ha,hb) in enumerate(splits(n,paired)):
        if set(ia)&set(ha) or set(ib)&set(hb):raise ValueError("fold leakage")
        tr0=a[:,list(ia)].mean(1);tr1=b[:,list(ib)].mean(1)
        te=b[:,list(hb)].mean(1)-a[:,list(ha)].mean(1)
        for gi,role in enumerate(roles):
            mask=G[:,gi].astype(bool);nn=[v for v,k in zip(names,mask) if k];size=len(nn)
            if size==0:raise ValueError("empty role")
            aa=tr0[:,mask];bb=tr1[:,mask]
            for t,date in enumerate(dates[1:],1):
                for sign in (1,-1):
                    v=sign*te[t,mask];y,m,eff=support(v)
                    lag=sign*(bb[t-1]-aa[t-1])
                    idx=sorted((i for i,x in enumerate(lag) if x>EPS),key=lambda i:(-lag[i],nn[i]))[:int(np.ceil(.1*size))]
                    net=float(v.sum());res=float(net-v[idx].sum())
                    rows.append({"source":st,"panel":key,"time":date,"role":role,"direction":sign,"fold":fold,
                      "roster":size,"target_top1":y,"support_n":m,"support_effective":eff,
                      "normalized_excess":(y-1/m)/(1-1/m) if m>=2 else None,
                      "control_features":feats(aa,bb,t,sign,"control"),"both_features":feats(aa,bb,t,sign,"both"),
                      "net_positive":net>EPS,"residual_positive":res>EPS,"n_excluded":len(idx)})
    return rows

def get_sources(root,poem_zip,field_dir,db,validation):
    rows=[];p=json.loads((validation/"narwani2019_forecast_summary_v1.json").read_text())
    k=json.loads((validation/"kansas2020_published_sown_role_v1.json").read_text())
    po=json.loads((validation/"poem_field_absolute_compensation_v1.json").read_text())
    fi=json.loads((validation/"vansteijn2026_field_opposition_v1.json").read_text())
    expect=[(p,45557051796),(k,45543133793),(po,45536266086),(fi,45537486831)]
    if any(z["source_run"]!=f"https://github.com/zuizui0223/knk/actions/runs/{n}" for z,n in expect):raise ValueError("prior lineage")
    x,y,arms,names,G,H,audit=pond_grid()
    for a,b in CONTRASTS:
        key=b+"-"+a
        rows+=panel("pond",key,WEEKS,x[a],x[b],names,G,["cyanobacteria","other_named_algae"],False,[v["contrast"] for v in p["contrasts"][key]["persistence"]],True)
    mapping=db/"foster2009_curated_name_bridge.csv"
    if gitblob(mapping)!=MAP_BLOB:raise ValueError("Kansas mapping")
    man=acquire(root);states,roster,dup,drift=read(root,man)
    if dup[EXPERIMENT] or drift[EXPERIMENT]:raise ValueError("Kansas ambiguous plot")
    z,blocks=kansas_grid(states[EXPERIMENT]);names=roster[EXPERIMENT];s,o,u,pre,ex=read_mapping(mapping,names)
    G=np.column_stack((s,o,u)).astype(float)
    rows+=panel("kansas","seeded-control",list(range(1,16)),z[:,:,0],z[:,:,1],names,G,["sown_list","other_named"],True,
                [[v["pooled"][q] for q in ("sown_list","other_named","unresolved")] for v in k["by_age"]],True)
    pp,pa=poem_panels(poem_zip);gr,names=poem_matrices(pp)
    if names!=po["taxon_names"]:raise ValueError("POEM roster")
    lookup={(v["cohort"],v["age"],v["early"]):v["species_delta"] for v in po["contrasts"]}
    if len(lookup)!=18:raise ValueError("POEM prior grain")
    for c,g in itertools.product(COHORTS,EARLY):
        G=role_matrix(names,g)[:,:3]
        a=np.array([gr[c,t,"S"] for t in AGES]);b=np.array([gr[c,t,g] for t in AGES])
        rows+=panel("poem",c+"_"+g+"_S",list(AGES),a,b,names,G,["early_sown","other_sown","nonsown"],False,[lookup[c,t,g] for t in AGES])
    for n,md5 in FIELD_MD5.items():
        if hashlib.md5((field_dir/n).read_bytes()).hexdigest()!=md5:raise ValueError("field source bytes")
    fp,fa=field_panels(field_dir);gr,names,fast,slow,_=field_unpack(fp)
    if names!=fi["taxon_names"]:raise ValueError("field roster")
    lookup={(v["site"],v["year"],v["reference"],v["history"]):v["species_delta"] for v in fi["taxon_contrasts"]}
    if len(lookup)!=32:raise ValueError("field prior grain")
    for s,h,r in itertools.product(SITES,HISTORIES,REFS):
        a=np.array([gr[s,t,r] for t in YEARS]);b=np.array([gr[s,t,h] for t in YEARS]);G=rolemap(names,fast,slow,h)[:,:3]
        rows+=panel("field",s+"_"+h+"_"+r,list(YEARS),a,b,names,G,["early_sown","other_sown","nonsown"],True,[lookup[s,t,r,h] for t in YEARS])
    return rows

def weights(rows):
    keys=[(r["source"],r["panel"],r["time"],r["role"],r["direction"]) for r in rows]
    counts=Counter(keys);ns=Counter(k[0] for k in counts)
    return np.array([1/(len(ns)*ns[k[0]]*counts[k]) for k in keys])
def fit(train,test,cols,mode):
    w=weights(train);y=np.array([r["target_top1"] for r in train])
    X=np.array([r[mode+"_features"] for r in train])[:,cols];Z=np.array([r[mode+"_features"] for r in test])[:,cols]
    mu=w@X;sd=np.sqrt(w@((X-mu)**2));sd=np.where(sd>EPS,sd,1.)
    a=np.column_stack((np.ones(len(X)),(X-mu)/sd));b=np.column_stack((np.ones(len(Z)),(Z-mu)/sd))
    pen=np.eye(a.shape[1])*LAM;pen[0,0]=0.
    coef=np.linalg.solve(a.T@(w[:,None]*a)+pen,a.T@(w*y))
    return np.clip(b@coef,0,1),coef.tolist()
def score(rows,mode,scheme,multi):
    r=[v for v in rows if v["target_top1"] is not None and (not multi or v["support_n"]>=2)]
    clusters=[(s,) for s in SOURCES] if scheme=="source" else [("kansas",),("pond",),("poem","field")]
    pred=[]
    for cl in clusters:
        tr=[v for v in r if v["source"] not in cl];te=[v for v in r if v["source"] in cl]
        if not tr or not te:raise ValueError("empty outer fold")
        pp={};coef={}
        for m,cols in MODELS.items():
            pp[m],coef[m]=fit(tr,te,cols,mode)
        for i,v in enumerate(te):
            pred.append({"source":v["source"],"panel":v["panel"],"time":v["time"],"role":v["role"],"direction":v["direction"],"fold":v["fold"],"target":v["target_top1"],"pred":{m:float(p[i]) for m,p in pp.items()}})
    if len(pred)!=len(r):raise ValueError("held coverage")
    def loss(p):
        w=weights(p);y=np.array([v["target"] for v in p])
        return {m:float(w@((np.array([v["pred"][m] for v in p])-y)**2)) for m in MODELS}
    by={s:loss([v for v in pred if v["source"]==s]) for s in SOURCES};agg=loss(pred)
    def gains(a):
        return {"dominance_over_richness":1-a["dominance"]/a["richness"],"dominance_over_history":1-a["history_dominance"]/a["history"],
                "richness_over_intercept":1-a["richness"]/a["intercept"],"history_over_richness":1-a["history"]/a["richness"]}
    return {"n_eligible":len(r),"n_by_source":dict(Counter(v["source"] for v in pred)),"mse":agg,"gains":gains(agg),
            "by_source":{s:{"mse":v,"gains":gains(v)} for s,v in by.items()}},pred

def descriptive(rows):
    out=[]
    for st in SOURCES:
        r=[v for v in rows if v["source"]==st]
        for role,sign in sorted({(v["role"],v["direction"]) for v in r}):
            z=[v for v in r if v["role"]==role and v["direction"]==sign];q=[v for v in z if v["target_top1"] is not None]
            pos=[v for v in z if v["net_positive"]]
            avg=lambda a,k:float(np.mean([v[k] for v in a])) if a else None
            bins={}
            for n in (1,2,3,4,5):
                w=[v for v in q if (v["support_n"]==n if n<5 else v["support_n"]>=5)]
                bins[str(n) if n<5 else "5+"]= {"n":len(w),"mean_top1":avg(w,"target_top1")}
            out.append({"source":st,"role":role,"direction":sign,"n_cells":len(z),"n_defined":len(q),
                        "mean_roster":avg(z,"roster"),"mean_top1":avg(q,"target_top1"),
                        "mean_support_count":avg(q,"support_n"),"mean_effective_support":avg(q,"support_effective"),
                        "mean_normalized_excess":avg([v for v in q if v["support_n"]>=2],"normalized_excess"),
                        "n_single_contributor":len([v for v in q if v["support_n"]==1]),"by_support_count":bins,
                        "n_positive_net":len(pos),"n_retained_10pct_roster":sum(v["residual_positive"] for v in pos),
                        "mean_removed_taxa":avg(pos,"n_excluded")})
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--poem",type=Path,required=True);p.add_argument("--field",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True);rows=get_sources(a.root,a.poem,a.field,Path("database/ahif_v1"),Path("validation/ahif_v1"))
    if len({(v["source"],v["panel"]) for v in rows})!=19:raise ValueError("contrast coverage")
    summary={};pred={}
    for mode,scheme,multi in (("control","source",False),("both","source",False),("control","cluster",False),("control","source",True)):
        key=f"{mode}_{scheme}_{'multi' if multi else 'all'}";summary[key],pred[key]=score(rows,mode,scheme,multi)
    z={"status":"SOURCE_LOCKED_CROSS_SOURCE_RESPONSE_CONCENTRATION_DESCRIPTIVE",
       "n_source_experiments":4,"n_contrast_panels":19,"n_all_cells":len(rows),
       "n_undefined_support":sum(v["target_top1"] is None for v in rows),
       "cells_by_source":dict(Counter(v["source"] for v in rows)),
       "descriptive":descriptive(rows),"forecast":summary,"lambda":LAM,
       "limits":["Conditional support concentration; not support occurrence prediction.","Lag control not necessarily pretreatment.",
                 "Only four source experiments, overlapping repeated units; descriptive, not formal meta-analysis.",
                 "Kansas original source subset not verified. No harmonized functional traits or causal claims."]}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps({**z,"feature_cells":rows,"prediction_cells":pred},allow_nan=False)+"\n")
    print("CROSS_SOURCE_CONCENTRATION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
