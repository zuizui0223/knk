"""Goldman endpoint source/dose forecasts, source-disjoint outcome holdouts.
Only observed pre-pseudopoint profiles. Separate stock proxy and p3 control lanes.
"""
import argparse,csv,io,json,zipfile,hashlib,re
from pathlib import Path
from collections import defaultdict,Counter
import numpy as np
from sklearn.linear_model import LogisticRegression
from benchmarks import probe_ahif_goldman2025 as source

COMBOS=("XBA-XBB","XCA-XCB","XDA-XDB","XFA-XFB","XBA-XCA","XCA-XDA","XDA-XFA","XFA-XBA")
PARENTS=sorted(set("-".join(COMBOS).split("-")))
RATIOS=("1000:1","100:1","10:1","1:1","1:10","1:100","1:1000")
METRICS=("fixed_pool","fixed_identity","fixed_dose","fit_pool","fit_identity","fit_dose","fit_dose_abundance")
MISSING={"","NA","NaN","nan","null","None"}
EPS=1e-6
def normkey(s):
    s=s.upper()
    assert re.fullmatch("[ACGTN]+",s), "Unexpected sequence key"
    return s
def rows(b):
    return csv.DictReader(io.StringIO(b.decode("utf-8-sig")),delimiter="\t")
def positive(s):
    if s in MISSING:return None
    x=float(s);assert np.isfinite(x) and x>=0
    return x if x>0 else None
def valid_profile(v):
    return bool(v) and abs(sum(v.values())-1)<=EPS
def add(d,key,seq,value):
    if value is None:return
    if seq in d[key]:raise ValueError("Duplicate observed taxon in profile "+str(key))
    d[key][seq]=value
def load():
    names=list(source.FILES);data={k:source.get(k) for k in names}
    receipt={k:{"blob":source.FILES[k],"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b)} for k,b in data.items()}
    with zipfile.ZipFile(io.BytesIO(data[names[0]])) as f:
        assert f.namelist()==["mixtureDataframe.txt"] and f.infolist()[0].file_size<=150_000_000
        actual=defaultdict(dict);control=defaultdict(dict)
        for r in rows(f.read(f.namelist()[0])):
            if r["combo"] not in COMBOS:raise ValueError("Unexpected combo")
            p=int(r["passage"]);ratio=r["ratio"];line=r["line"]
            if p not in (3,5):raise ValueError("Unexpected passage")
            a=positive(r["rel_abundance_old"])
            if a is None:continue
            seq=normkey(r["OTU"])
            if line.startswith("actual-") and ratio in RATIOS:
                assert r["mixtureType"]=="actual"
                add(actual,(r["combo"],p,ratio,int(line[-1])),seq,a)
            elif line in ("theoretical-1","theoretical-2") and ratio in ("1:0","0:1"):
                assert r["mixtureType"]=="theoretical"
                add(control,(r["combo"],p,ratio,int(line[-1])),seq,a)
    # Original metadata verifies physical stock IDs, ratios (orientation corrected by author code)
    meta=list(rows(data[names[1]]))
    observed=set()
    for r in meta:
        if r["group"]!="e0017" or r["extractionReplicate"]!="1" or r["biosample2"]=="NA":continue
        left,right=r["biosample1"].split("-")[0],r["biosample2"].split("-")[0]
        combo=left+"-"+right
        rev=False
        if combo not in COMBOS:combo=right+"-"+left;rev=True
        assert combo in COMBOS
        rat=r["ratio"]
        if rev:rat=":".join(rat.split(":")[::-1])
        assert rat in RATIOS
        key=(combo,int(r["passage"]),rat,int(r["inoculationReplicate"]))
        assert key not in observed
        observed.add(key)
    assert len(observed)==336
    assert set(actual)<=observed
    p5={k:v for k,v in actual.items() if k[1]==5}
    missing=sorted(k for k in observed if k[1]==5 and k not in p5)
    assert ("XFA-XBA",5,"100:1",1) in missing,"Author-excluded sample found"
    bad={str(k):sum(v.values()) for k,v in p5.items() if not valid_profile(v)}
    covered={c:{r:sum(k[0]==c and k[2]==r for k in p5) for r in RATIOS} for c in COMBOS}
    audit={"source_files":receipt,"planned_p5":168,"observed_p5":len(p5),"missing_p5":[list(k) for k in missing],"p5_nonunit_profiles":bad,"p5_ratio_counts":covered,"p5_sum_range":[min(sum(v.values()) for v in p5.values()),max(sum(v.values()) for v in p5.values())]}
    # Deduplicate actual measured parent controls reused by the author across pairs.
    ctrl={};ctrl_conflicts=[];ctrl_bad={}
    for (c,p,r,rep),v in control.items():
        if p!=3:continue
        parent=c.split("-")[0 if r=="1:0" else 1]
        key=(parent,rep)
        if not valid_profile(v):ctrl_bad[str((c,p,r,rep))]=sum(v.values());continue
        if key in ctrl:
            if set(ctrl[key])!=set(v) or max(abs(ctrl[key][s]-v[s]) for s in v)>1e-10:ctrl_conflicts.append(str((c,p,r,rep)))
        else:ctrl[key]=v
    # Prior separate stock experiment; verify exact sample crosswalk including passage and medium.
    pm={r["sample"]:r for r in rows(data[names[3]])}
    stocks=defaultdict(dict);stock_meta={};stock_conflicts=[]
    for r in rows(data[names[2]]):
        if r["passage"]!="15" or r["media"]!="mBHImucin" or r["PCRreplicate"]!="1":continue
        k=r["sample"];q=pm.get(k)
        if q is None or any(r[col]!=q[col] for col in ("biosample","media","passage","inoculationreplicate","PCRreplicate","filename")):
            stock_conflicts.append(k);continue
        par=r["biosample"].split("-")[0];assert par in PARENTS
        a=positive(r["rel_abundance"])
        if a is not None:add(stocks,k,normkey(r["OTU"]),a)
        stock_meta[k]=(par,r["inoculationreplicate"])
    audit["p3_control_profiles"]={p:sorted(rep for par,rep in ctrl if par==p) for p in PARENTS}
    audit["p3_control_conflicts"]=ctrl_conflicts;audit["p3_control_nonunit"]=ctrl_bad
    audit["stock_metadata_conflicts"]=sorted(set(stock_conflicts))
    audit["p15_stock_profiles"]={p:[k for k,v in stocks.items() if stock_meta[k][0]==p] for p in PARENTS}
    audit["p15_stock_nonunit"]={k:sum(v.values()) for k,v in stocks.items() if not valid_profile(v)}
    print("GOLDMAN_DOSE_AUDIT",json.dumps(audit,allow_nan=False),flush=True)
    assert not bad and all(n>=1 for x in covered.values() for n in x.values()),"STOP incomplete or nonunit p5 observed profile"
    lanes={}
    def mean_profiles(x):
        result={}
        for p in PARENTS:
            v=x[p];assert v
            m=defaultdict(float)
            for z in v:
                for s,w in z.items():m[s]+=w/len(v)
            assert valid_profile(m)
            result[p]=dict(m)
        return result
    if not ctrl_conflicts and not ctrl_bad and all(audit["p3_control_profiles"].values()):
        lanes["p3_control_proxy"]=mean_profiles({p:[v for (par,r),v in ctrl.items() if par==p] for p in PARENTS})
    if not stock_conflicts and not audit["p15_stock_nonunit"] and all(audit["p15_stock_profiles"].values()):
        lanes["p15_stock_proxy"]=mean_profiles({p:[v for k,v in stocks.items() if stock_meta[k][0]==p] for p in PARENTS})
    return p5,lanes,audit

def data_rows(mats,pre,threshold):
    ks=sorted(set().union(*(set(v) for v in pre.values())))
    K=len(ks)
    pp={p:np.array([pre[p].get(k,0.) for k in ks]) for p in PARENTS}
    det={p:(v>0 if threshold==0 else v>=threshold) for p,v in pp.items()}
    regional=np.any(np.stack(list(det.values())),axis=0)
    result={}
    for c in COMBOS:
        x=[];y=[];w=[];fixed=[];dose_groups=[];a,b=c.split("-")
        for r in RATIOS:
            pairs=sorted((k,v) for k,v in mats.items() if k[0]==c and k[2]==r)
            assert pairs
            a0,b0=map(float,r.split(":"))
            for k,z in pairs:
                val=np.array([z.get(s,0.) for s in ks]);target=(val>0 if threshold==0 else val>=threshold)
                for focal,donor,q in ((a,b,b0/(a0+b0)),(b,a,a0/(a0+b0))):
                    risk=regional&~det[focal]
                    assert risk.any()
                    P=np.mean([det[p] for p in PARENTS if p not in (focal,donor)],axis=0)[risk]
                    D=det[donor][risk].astype(float)
                    L=np.full(len(P),np.log10(q/(1-q))/3)*D
                    A=np.log10(np.maximum(q*pp[donor][risk],1e-8))/8
                    x.append(np.column_stack((P,D,L,A)));y.append(target[risk].astype(int))
                    fixed.append(np.column_stack((P,(P+D)/2,(P+q*D)/2)))
                    w.append(np.full(len(P),1/(7*len(pairs)*2*len(P))))
                    dose_groups.append(np.full(len(P),r+"|"+focal))
        X=np.concatenate(x);Y=np.concatenate(y);W=np.concatenate(w);F=np.concatenate(fixed)
        assert abs(W.sum()-1)<1e-10
        result[c]={"X":X,"Y":Y,"W":W,"F":F}
    return result,K,int(regional.sum())
def brier(y,p,w):
    return float(np.sum(w*(y-p)**2)/sum(w))
def evaluate(table,mode):
    out={}
    for held in COMBOS:
        a,b=held.split("-")
        def keep(c):
            u=c.split("-")
            if mode=="source_disjoint":return not set((a,b))&set(u)
            return not {a[:2],b[:2]}&{p[:2] for p in u}
        train=[c for c in COMBOS if keep(c)]
        h=table[held];y=h["Y"];w=h["W"]
        base={k:brier(y,h["F"][:,j],w) for j,k in enumerate(METRICS[:3])}
        r={"n_training_pairs":len(train),"training_pairs":train,"n_risk_cells":len(y),"target_prevalence":float(np.sum(w*y)),"fixed":base}
        if len(train)<2:r["status"]="STOP_INSUFFICIENT_TRAINING_PAIRS";out[held]=r;continue
        tx=np.concatenate([table[c]["X"] for c in train]);ty=np.concatenate([table[c]["Y"] for c in train]);tw=np.concatenate([table[c]["W"]*100 for c in train])
        if len(np.unique(ty))<2:r["status"]="STOP_ONE_TARGET_CLASS";out[held]=r;continue
        r["status"]="SCORED";r["fit"]={}
        for n,k in enumerate(METRICS[3:],1):
            model=LogisticRegression(C=1,solver="lbfgs",max_iter=1000,tol=1e-9)
            model.fit(tx[:,:n],ty,sample_weight=tw)
            assert model.n_iter_.max()<1000
            pred=model.predict_proba(h["X"][:,:n])[:,1]
            r["fit"][k]=brier(y,pred,w)
        out[held]=r
    valid=[r for r in out.values() if r["status"]=="SCORED"]
    means={}
    if valid:
        for k in METRICS:
            means[k]=float(np.mean([r["fixed" if k.startswith("fixed") else "fit"][k] for r in valid]))
        means["identity_gain"]=means["fit_pool"]-means["fit_identity"]
        means["dose_gain"]=means["fit_identity"]-means["fit_dose"]
        means["dose_abundance_gain"]=means["fit_dose"]-means["fit_dose_abundance"]
        means["fixed_identity_gain"]=means["fixed_pool"]-means["fixed_identity"]
        means["fixed_dose_increment"]=means["fixed_identity"]-means["fixed_dose"]
    return {"n_scored_pairs":len(valid),"n_total_pairs":len(out),"means":means,"n_positive_dose":sum(r["fit"]["fit_identity"]>r["fit"]["fit_dose"] for r in valid),"n_positive_identity":sum(r["fit"]["fit_pool"]>r["fit"]["fit_identity"] for r in valid),"by_pair":out}

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    mats,lanes,audit=load();result={"status":"EXPLORATORY_ENDPOINT_DOSE_SOURCE_PROXY_NOT_TEMPORAL_META","source_commit":source.REF,"audit":audit,"lanes":{}}
    for label in ("p15_stock_proxy","p3_control_proxy"):
        if label not in lanes:result["lanes"][label]={"status":"STOP_SOURCE_PROXY_CROSSWALK_OR_PROFILE"};continue
        pre=lanes[label];roster=set().union(*(set(v) for v in pre.values()))
        cov={c:float(np.mean([sum(w for s,w in v.items() if s in roster) for k,v in mats.items() if k[0]==c])) for c in COMBOS}
        z={"source_roster":len(roster),"parent_lengths":dict(Counter(map(len,roster))),"endpoint_lengths":dict(Counter(map(len,set().union(*(set(v) for v in mats.values()))))),"coverage_by_pair":cov}
        if min(cov.values())<.9:
            z["status"]="STOP_CROSSRUN_OR_ROSTER_SUPPORT_LT_0.90";result["lanes"][label]=z;continue
        z["status"]="SCORED_CONDITIONAL_SOURCE_PROXY";z["results"]={}
        for mode,threshold in (("threshold_0.001",.001),("raw_detection",0)):
            table,k,n=data_rows(mats,pre,threshold)
            z["results"][mode]={"source_roster":n,"source_disjoint":evaluate(table,"source_disjoint"),"household_disjoint":evaluate(table,"household_disjoint")}
        result["lanes"][label]=z
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False)+"\n")
    print("GOLDMAN_DOSE_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
