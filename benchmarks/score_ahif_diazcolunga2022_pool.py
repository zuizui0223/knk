"""Independent crossed-parent forward detection audit (endpoint only)."""
import argparse,csv,hashlib,io,itertools,json,re,urllib.request
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
REF="c21bfb9bfe8cbe7dbcff52e76c4d85e3ee7e7f6d"
BASE="https://raw.githubusercontent.com/jdiazc9/coalescence/"+REF+"/"
FILES={"data/metaData.csv":"4039c5ed6f1aa87519d846b1ddc952eb07fb5327","data/OTU_table.csv":"00e527c890bc324c2e1b3d3f8a842ef484b22b0b","data/analyzeData.R":"88f31dc57777f0c249840d67f6cd1fc9e8f4ef2f"}
MAP={"Glutamine":[1,3,5,7,8,10,11,12],"Citrate":[1,2,3,7,8,9,10,11]}
EXP={"Glutamine":("Gln","T7-C1-2"),"Citrate":("Citr","T7-C4-2")}
PAIRS=list(itertools.combinations(range(8),2))
def get(p):
    with urllib.request.urlopen(urllib.request.Request(BASE+p,headers={"User-Agent":"AHIF source audit"}),timeout=90) as r:b=r.read(20_000_001)
    assert len(b)<=20_000_000 and hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==FILES[p],p
    return b
def load():
    data={p:get(p) for p in FILES}
    receipt={p:{"blob":FILES[p],"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest()} for p,b in data.items()}
    script=data["data/analyzeData.R"].decode("utf-8-sig")
    for m,ids in MAP.items():
        pat=r"original_communities\[\['"+m+r"'\]\]\s*<-\s*setNames\(paste\('Community',c\(([^)]*)\)"
        q=re.search(pat,script)
        assert q and [int(v) for v in q.group(1).split(",")]==ids,(m,"parent mapping")
    meta=list(csv.DictReader(io.StringIO(data["data/metaData.csv"].decode("utf-8-sig"))))
    assert len(meta)==1550 and len({r["Sample"] for r in meta})==1550
    samples={};needed=[]
    for m,(prefix,exp) in EXP.items():
        def unique(rows):
            assert len(rows)==1,(m,"nonunique metadata",len(rows))
            return rows[0]["Sample"]
        parent=[unique([r for r in meta if r["Experiment"]==f"{prefix}-T12-X{x}" and r["Comm1"]==f"Community{x}" and r["Carbon"]==m and r["Transfer"]=="12"]) for x in MAP[m]]
        block=[r for r in meta if r["Experiment"]==exp and r["Carbon"]==m and r["Transfer"]=="7" and r["StabilizingCarb1"]==m]
        assert len(block)==72
        control=[unique([r for r in block if r["Comm1"]==f"Community{i+1}" and r["Comm2"]=="None"]) for i in range(8)]
        mix={(i,j):unique([r for r in block if r["Comm1"]==f"Community{i+1}" and r["Comm2"]==f"Community{j+1}" and r["StabilizingCarb2"]==m]) for i in range(8) for j in range(8) if i!=j}
        assert len(mix)==56
        samples[m]=(parent,control,mix)
        needed+=parent+control+list(mix.values())
    assert len(needed)==144 and len(set(needed))==144
    rd=csv.reader(io.StringIO(data["data/OTU_table.csv"].decode("utf-8-sig")))
    h=next(rd);assert len(h)==len(set(h)), "duplicate header"
    assert set(needed)<=set(h),("missing columns",sorted(set(needed)-set(h)))
    idx=[h.index(s) for s in needed];seen=set();mat=[]
    for row in rd:
        assert len(row)==len(h) and row[0] not in seen
        seen.add(row[0]);v=np.array([float(row[i]) for i in idx])
        assert np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all()
        mat.append(v)
    x=np.array(mat).T
    depth=x.sum(1);audit={"n_taxa":len(seen),"n_profiles":len(needed),"n_parent_profiles":16,"n_control_profiles":16,"n_mixed_profiles":112,"min_reads":int(depth.min()),"max_reads":int(depth.max()),"below_1000":{s:int(d) for s,d in zip(needed,depth) if d<1000},"source_files":receipt,"parent_mapping":MAP}
    print("DIAZ_SOURCE_AUDIT",json.dumps(audit),flush=True)
    if audit["below_1000"]:raise ValueError("STOP selected profiles below prespecified read gate")
    p=x/depth[:,None];d=dict(zip(needed,p))
    out={}
    for m,(par,c,mix) in samples.items():
        out[m]=(np.stack([d[s] for s in par]),np.stack([d[s] for s in c]),{k:d[s] for k,s in mix.items()})
    return out,audit
def detection(x,threshold):
    return x>0 if threshold==0 else x>=threshold
def losses(P,D,y):
    H=(P+D)/2
    a=float(np.mean((P-y)**2));b=float(np.mean((D-y)**2));c=float(np.mean((H-y)**2))
    return np.array([a,b,c,a-c,a-b])
METRICS=("pool_loss","donor_loss","hybrid_loss","hybrid_gain","donor_gain")
def metric_dict(a):
    return dict(zip(METRICS,map(float,a)))
def safe_corr(x,y):
    return float(spearmanr(x,y).statistic) if len(x)>=3 and np.std(x)>0 and np.std(y)>0 else None
def avg(x):
    return np.mean(x,axis=0) if x else None
def predict_medium(p,c,mix,threshold):
    d=detection(p,threshold);ctl=detection(c,threshold);mx={k:detection(v,threshold) for k,v in mix.items()}
    pool=d.any(0);rows=[];subgroups={s:[] for s in ("shared_pool","partner_only")}
    sub_n={s:0 for s in subgroups};outside=[]
    for i,j in PAIRS:
        gains=[];controls=[];wrong=[];subs={s:[] for s in subgroups}
        for well in ((i,j),(j,i)):
            outside.append(float(mix[well][~pool].sum()))
            for focal,partner in ((i,j),(j,i)):
                risk=pool&~d[focal]
                assert risk.any()
                P=np.mean(d[[k for k in range(8) if k not in (focal,partner)]],axis=0)[risk]
                D=d[partner,risk].astype(float);Y=mx[well][risk].astype(float);C=ctl[focal,risk].astype(float)
                a=losses(P,D,Y);b=losses(P,D,C)
                gains.append(a);controls.append(b)
                alt=[]
                for k in range(8):
                    if k in (focal,partner):continue
                    Q=np.mean(d[[l for l in range(8) if l not in (focal,k)]],axis=0)[risk]
                    alt.append(losses(Q,d[k,risk].astype(float),Y))
                assert len(alt)==6
                wrong.append(float(np.mean(alt,axis=0)[2]-a[2]))
                for name,mask in (("shared_pool",P>0),("partner_only",(P==0)&(D==1))):
                    if mask.any():
                        subs[name].append(losses(P[mask],D[mask],Y[mask]))
                        sub_n[name]+=1
        x=avg(gains);y=avg(controls)
        row={"i":i,"j":j,"mismatch_TV":float(np.abs(p[i]-p[j]).sum()/2),"mixed":metric_dict(x),"unmixed_control":metric_dict(y),"mixed_minus_control_gain":float(x[3]-y[3]),"actual_vs_alternative_hybrid_gain":float(np.mean(wrong))}
        for name in subs:
            z=avg(subs[name])
            if z is not None:subgroups[name].append(z)
        rows.append(row)
    assert len(rows)==28 and len(outside)==56
    summary={"pairs":len(rows),"mixed_wells":len(mx),"source_roster_taxa":int(pool.sum()),"mixed":{k:float(np.mean([r["mixed"][k] for r in rows])) for k in METRICS},"unmixed_control":{k:float(np.mean([r["unmixed_control"][k] for r in rows])) for k in METRICS},"mixed_minus_control_gain":float(np.mean([r["mixed_minus_control_gain"] for r in rows])),"actual_vs_alternative_hybrid_gain":float(np.mean([r["actual_vs_alternative_hybrid_gain"] for r in rows])),"mismatch_spearman":safe_corr([r["mismatch_TV"] for r in rows],[r["mixed"]["hybrid_gain"] for r in rows]),"n_positive_hybrid_pairs":sum(r["mixed"]["hybrid_gain"]>0 for r in rows),"mean_endpoint_mass_outside_roster":float(np.mean(outside)),"subgroups":{k:{"n_defined_pair_orient_wells":sub_n[k],"n_defined_pairs":len(v),"metrics":metric_dict(avg(v)) if v else None} for k,v in subgroups.items()}}
    return summary,rows
def calibrate(p,c,mix,threshold):
    # Train only on pairs disjoint in both parent communities. Aggregate each
    # training pair once per taxon, separately for donor-detection states.
    d=detection(p,threshold);mx={k:detection(v,threshold) for k,v in mix.items()}
    pool=d.any(0);K=d.shape[1];rows=[]
    pair_stats={}
    for i,j in PAIRS:
        n=np.zeros((3,K));s=np.zeros((3,K))
        for well in ((i,j),(j,i)):
            for f,o in ((i,j),(j,i)):
                risk=pool&~d[f]
                for state,mask in ((0,risk),(1,risk&~d[o]),(2,risk&d[o])):
                    n[state]+=mask;s[state]+=mask*mx[well]
        valid=n>0;rate=np.divide(s,n,out=np.zeros_like(s),where=valid)
        pair_stats[i,j]=(valid.astype(float),rate)
    for i,j in PAIRS:
        train=[k for k in PAIRS if i not in k and j not in k]
        assert len(train)==15
        n=sum(pair_stats[k][0] for k in train);s=sum(pair_stats[k][1] for k in train)
        base=(s[0]+1)/(n[0]+2)
        cond=np.divide(s[1:]+1,n[1:]+2)
        cond=np.where(n[1:]>0,cond,base[None,:])
        z=[];fb=0;total=0
        for well in ((i,j),(j,i)):
            for f,o in ((i,j),(j,i)):
                risk=pool&~d[f];state=d[o].astype(int);q=cond[state,np.arange(K)]
                y=mx[well][risk].astype(float);a=float(np.mean((base[risk]-y)**2));b=float(np.mean((q[risk]-y)**2))
                z.append([a,b,a-b]);fb+=int((n[1:][state,np.arange(K)][risk]==0).sum());total+=int(risk.sum())
        a,b,g=avg(z);rows.append({"i":i,"j":j,"base_loss":float(a),"donor_conditional_loss":float(b),"gain":float(g),"fallback_cells":fb,"risk_cells":total})
    return {"pairs":len(rows),"train_pairs_per_held":15,"base_loss":float(np.mean([r["base_loss"] for r in rows])),"donor_conditional_loss":float(np.mean([r["donor_conditional_loss"] for r in rows])),"gain":float(np.mean([r["gain"] for r in rows])),"n_positive_pairs":sum(r["gain"]>0 for r in rows),"fallback_cells":sum(r["fallback_cells"] for r in rows),"risk_cells":sum(r["risk_cells"] for r in rows),"by_pair":rows}
def run(data):
    out={}
    for label,threshold in (("raw_detected",0),("relative_0.001",.001)):
        summary={};pair_rows={};cal={}
        for m,(p,c,mix) in data.items():
            summary[m],pair_rows[m]=predict_medium(p,c,mix,threshold)
            cal[m]=calibrate(p,c,mix,threshold)
        pooled={k:float(np.mean([v["mixed"][k] for v in summary.values()])) for k in METRICS}
        centered_x=[];centered_y=[]
        for rows in pair_rows.values():
            x=np.array([r["mismatch_TV"] for r in rows]);y=np.array([r["mixed"]["hybrid_gain"] for r in rows])
            centered_x.extend(x-x.mean());centered_y.extend(y-y.mean())
        centered=float(np.corrcoef(centered_x,centered_y)[0,1])
        deleted={}
        for origin in sorted(set(sum(MAP.values(),[]))):
            means=[]
            for m,rows in pair_rows.items():
                a=[r["mixed"]["hybrid_gain"] for r in rows if origin not in (MAP[m][r["i"]],MAP[m][r["j"]])]
                assert a;means.append(float(np.mean(a)))
            deleted[str(origin)]=float(np.mean(means))
        out[label]={"by_medium":summary,"equal_medium_mixed":pooled,"medium_centered_mismatch_pearson":centered,"delete_original_id_gain":deleted,"equal_medium_control_gain":float(np.mean([v["unmixed_control"]["hybrid_gain"] for v in summary.values()])),"equal_medium_mixed_minus_control_gain":float(np.mean([v["mixed_minus_control_gain"] for v in summary.values()])),"equal_medium_actual_vs_alternative_hybrid_gain":float(np.mean([v["actual_vs_alternative_hybrid_gain"] for v in summary.values()])),"source_disjoint_calibration":cal,"pairs":pair_rows}
    return out
def main():
    a=argparse.ArgumentParser();a.add_argument("--out",type=Path,required=True);a=a.parse_args()
    data,audit=load();out={"status":"INDEPENDENT_ENDPOINT_POOL_DIAGNOSTIC_NOT_TEMPORAL_META","source_commit":REF,"publication":"10.1073/pnas.2111261119","audit":audit,"results":run(data)}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(out,allow_nan=False)+"\n")
    brief={**out,"results":{k:{**{j:v for j,v in r.items() if j not in ("pairs","source_disjoint_calibration")},"source_disjoint_calibration":{m:{x:y for x,y in z.items() if x!="by_pair"} for m,z in r["source_disjoint_calibration"].items()}} for k,r in out["results"].items()}}
    print("DIAZ_POOL_RESULT",json.dumps(brief,allow_nan=False),flush=True)
if __name__=="__main__":main()
