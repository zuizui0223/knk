"""Independent Narwani pond forecast; frozen source, mapping and held-time/pond tests."""
from __future__ import annotations
import argparse,itertools,json
from collections import Counter
from pathlib import Path
import numpy as np
from benchmarks.audit_ahif_narwani_metadata import load_source,isnum,MD5,SHA
from benchmarks.score_ahif_kansas_forecast import forecasts,signscore,MODELS as SPEC_MODELS
from benchmarks.score_ahif_kansas_group_forecast import model_predictions,summarize
GROUPS=("Cyanobacteria","Golden_algae","Diatoms","Dinoflagellates","Cryptophytes","Green_algae")
CYAN=set("Aphanocapsa Chroococus Gomphosphaeria Phormidium Synechococcales Unknown_1 Unknown_3 Unknown_4".split())
DIAT=set("Asterionella Cyclotella_1 Cyclotella_2 Diatoma Nitzschia Stephanodiscus_1 Stephanodiscus_2 Stephanodiscus_3 Synedra_1 Synedra_2 Synedra_3 Unknown_5".split())
GREEN=set("Carteria Coelastrum_1 Coelastrum_2 Coelastrum_3 Cosmarium_1 Cosmarium_2 Elakatothrix Kirchneriella Lagerheimia_1 Lagerheimia_2 Monoraphidium Nephrocitium Oocystis_1 Oocystis_2 Pediastrum Phacotus Rayssiella Scenedesmus_1 Scenedesmus_2 Scenedesmus_3 Scenedesmus_4 Tetraedron Willea Unknown_8 Unknown_9 Unknown_10 Unknown_12".split())
OTHER=set("Gymnodinium Unknown_6 Unknown_7 Unkonwn_11".split())
WEEKS=[1,3,5,8,9,11,12,13,14,15,16,19,24,28,32]
CONTRASTS=[("C","MD"),("C","M"),("C","D"),("O","C")]
FIRST=5
SPLITS=list(itertools.combinations(range(4),2))
def mapping(names):
    groups=[CYAN,DIAT|GREEN|OTHER]
    if len(set(names))!=66 or any(not s<=set(names) for s in groups) or groups[0]&groups[1]:raise ValueError("mapping roster")
    g=np.array([[n in groups[0],n in groups[1],n not in groups[0]|groups[1]] for n in names],float)
    if not np.all(g.sum(1)==1) or tuple(g.sum(0))!=(8,43,15):raise ValueError("partition coverage")
    return g,np.column_stack((g[:,0],g[:,1]+g[:,2]))
def grid():
    names,r,m=load_source();columns=[v for v in names if v.startswith("v_") and v[2:] not in GROUPS]
    taxa=[v[2:] for v in columns];src=["v_"+v for v in GROUPS];G,H=mapping(taxa)
    by={(v["Week"],v["Pond"]):v for v in r}
    if len(by)!=660:raise ValueError("duplicate pond-week")
    ponds=sorted({v["Pond"] for v in r});design={p:[v for v in r if v["Pond"]==p] for p in ponds}
    arms={a:[] for a in ("O","C","D","M","MD")}
    expect={"O":(0,0,0),"C":(0,0,1),"D":(1,0,1),"M":(0,1,1),"MD":(1,1,1)}
    for p,v in design.items():
        z={(q["Treatment"],q["Dreisse"],q["Macrophytes"],q["Nutrients"]) for q in v}
        if len(z)!=1:raise ValueError("allocation drift")
        a,d,m,n=next(iter(z))
        if (d,m,n)!=expect[a]:raise ValueError("treatment code")
        arms[a].append(p)
    if any(len(v)!=4 for v in arms.values()):raise ValueError("arm allocation")
    complete=[w for w in sorted({v["Week"] for v in r}) if all(isnum(by[w,p][c]) for p in ponds for c in columns+src)]
    if complete!=WEEKS:raise ValueError("whole-week complete-case gate")
    x={a:np.array([[ [by[w,p][c] for c in columns] for p in ps] for w in WEEKS],float)*1e-6 for a,ps in arms.items()}
    y={a:np.array([[ [by[w,p][c] for c in src] for p in ps] for w in WEEKS],float)*1e-6 for a,ps in arms.items()}
    if any(not np.isfinite(v).all() or np.any(v<0) for v in list(x.values())+list(y.values())):raise ValueError("invalid observed biovolume")
    if any(np.any(v.sum(2)<=0) for v in x.values()):raise ValueError("zero total")
    allx=np.concatenate(list(x.values()),axis=1);ally=np.concatenate(list(y.values()),axis=1)
    cy=allx@G;source=np.stack((ally[:,:,0],ally[:,:,1:].sum(2)),axis=2)
    sums=(cy[:,:,0]-source[:,:,0]).ravel();tot=(allx.sum(2)-ally.sum(2)).ravel()
    audit={"cyan_abs_error_mean":float(np.abs(sums).mean()),"cyan_abs_error_max":float(np.abs(sums).max()),
           "total_abs_error_mean":float(np.abs(tot).mean()),"total_abs_error_max":float(np.abs(tot).max()),
           "cyan_exact_rows":int(np.isclose(sums,0,atol=1e-8).sum()),"total_exact_rows":int(np.isclose(tot,0,atol=1e-8).sum())}
    return x,y,arms,taxa,G,H,audit
def opposite(v):
    return bool(v[0]>0 and v[1]<0)
def reverse(v):
    return bool(v[0]<0 and v[1]>0)
def persistence(x0,x1,G):
    d=(x1.mean(1)-x0.mean(1));z=d@G;rows=[]
    for k,w in enumerate(WEEKS):
        unresolved=G[:,2].astype(bool);u=d[k,unresolved]
        lower=z[k,0]+np.minimum(u,0).sum()
        upper=z[k,1]+np.maximum(u,0).sum()
        half=[];smin=[]
        for a,b in itertools.product(SPLITS,repeat=2):
            h=(x1[k,list(b)].mean(0)-x0[k,list(a)].mean(0))@G
            half.append(opposite(h));smin.append(min(h[0],-h[1]))
        rows.append({"week":w,"contrast":z[k].tolist(),"opposite":opposite(z[k]),"reverse":reverse(z[k]),
           "robust_unresolved":bool(lower>0 and upper<0),"held_half_opposition_fraction":float(np.mean(half)),
           "held_half_signed_min":float(np.mean(smin))})
    return rows
def simple_summary(rows):
    obs=np.array([v["observed"] for v in rows]);out={}
    for m in rows[0]["predictions"]:
        pred=np.array([v["predictions"][m] for v in rows])
        out[m]={"sse":float(((pred-obs)**2).sum()),"total_sse":float(((pred.sum(1)-obs.sum(1))**2).sum())}
    return out
def evaluate(x,y,arms,G,H,control,treat):
    x0,x1=x[control],x[treat];y0,y1=y[control],y[treat]
    labels={"conservative3":["cyanobacteria","other_named_algae","unresolved_other"],"cyan_other2":["cyanobacteria","all_other_morphotypes"]}
    groups={"conservative3":G,"cyan_other2":H};rows={n:[] for n in groups};sr=[];ss=[]
    src=lambda a:np.stack((a[...,0],a[...,1:].sum(-1)),axis=-1)
    for a,b in itertools.product(SPLITS,repeat=2):
        c=[i for i in range(4) if i not in a];d=[i for i in range(4) if i not in b]
        train=np.stack((x0[:,list(a)].mean(1),x1[:,list(b)].mean(1)),axis=1)
        test=x1[:,d].mean(1)-x0[:,c].mean(1)
        delta=train[:,1]-train[:,0]
        tsource=src(y1[:,list(b)].mean(1)-y0[:,list(a)].mean(1));hsource=src(y1[:,d].mean(1)-y0[:,c].mean(1))
        for k in range(FIRST,len(WEEKS)):
            sp=forecasts(train[:,0],delta,np.any(train[:k]>0,axis=(0,1)),k)
            key={"week":WEEKS[k],"train_control":[arms[control][i] for i in a],"train_treat":[arms[treat][i] for i in b],
                 "held_control":[arms[control][i] for i in c],"held_treat":[arms[treat][i] for i in d]}
            ss.append({**key,"scores":{m:float(((sp[m]-test[k])**2).sum()) for m in SPEC_MODELS},"concurrent":signscore(delta[k],test[k])})
            for n,g in groups.items():
                pred=model_predictions(train,k,g,sp);obs=test[k]@g
                rows[n].append({**key,"observed":obs.tolist(),"predictions":{m:v.tolist() for m,v in pred.items()}})
            sr.append({**key,"observed":hsource[k].tolist(),"predictions":{"zero":[0.,0.],"previous":tsource[k-1].tolist(),"cumulative":tsource[:k].mean(0).tolist()}})
    if any(len(v)!=360 for v in [*rows.values(),sr,ss]):raise ValueError("split count")
    summary={}
    for n,v in rows.items():
        summary[n]={"labels":labels[n],"group_sizes":groups[n].sum(0).astype(int).tolist(),
           "aggregate":summarize(v,labels[n],True),
           "by_week":[{"week":w,**summarize([r for r in v if r["week"]==w],labels[n],True)} for w in WEEKS[FIRST:]]}
    summary["morphotype_sse"]={m:float(sum(r["scores"][m] for r in ss)) for m in SPEC_MODELS}
    summary["source_group_totals"]={"aggregate":simple_summary(sr),"by_week":[{"week":w,"contrast":src(y1.mean(1)-y0.mean(1))[k].tolist(),"opposite":opposite(src(y1.mean(1)-y0.mean(1))[k])} for k,w in enumerate(WEEKS)]}
    summary["persistence"]=persistence(x0,x1,G)
    return summary,{"representations":rows,"source_groups":sr,"morphotypes":ss}
def run():
    x,y,arms,taxa,G,H,audit=grid();out={};cells={}
    for a,b in CONTRASTS:out[b+"-"+a],cells[b+"-"+a]=evaluate(x,y,arms,G,H,a,b)
    return {"status":"SOURCE_VERIFIED_NARWANI_POND_TIME_FORECAST_DESCRIPTIVE",
      "source_md5":MD5,"source_sha256":SHA,"unit":"mm^3/L","n_ponds":20,"n_morphotypes":66,"n_complete_pond_weeks":300,
      "weeks":WEEKS,"excluded_missing_weeks":[7,10],"warmup_observations":FIRST,"n_cells_per_contrast":360,
      "arm_ponds":arms,"mapping":dict(zip(taxa,["cyanobacteria" if v[0] else "other_named_algae" if v[1] else "unresolved_other" for v in G])),
      "source_group_reconciliation":audit,"primary":"MD-C","contrasts":out,
      "limits":["Independent ponds; overlapping train/test splits and repeated weeks dependent.","Not a history manipulation or mechanism test.","Morphotype biovolumes estimated; sums may include non-algae.","No p-values, model tuning or held-out pond/week leakage."],"cells":cells}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run();a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    # Full derived cells and per-week summaries in artifact. Keep log bounded.
    s={k:v for k,v in z.items() if k!="cells"}
    for v in s["contrasts"].values():
        for n in ("conservative3","cyan_other2"):v[n].pop("by_week")
    print("NARWANI_FORECAST_RESULT",json.dumps(s,allow_nan=False))
if __name__=="__main__":main()
