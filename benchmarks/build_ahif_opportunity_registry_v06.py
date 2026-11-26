"""Build AHIF opportunity database: derived risk cells, exposure assignment and measured coculture endpoints."""
import argparse,csv,hashlib,io,itertools,json,sqlite3,subprocess,sys
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_villarreal2025_external as vv
from benchmarks import score_ahif_diazcolunga2022_harmonized as dh
from benchmarks import score_ahif_diazcolunga2022_pool as dz
from benchmarks import score_ahif_goldman2025_dose as gg
DB=Path("database/ahif_mechanism_v0")
VAL=Path("validation/ahif_v1")
TOL=1e-10
SOURCES=[
 ("villarreal2025","10.1002/ecy.70256","10.6084/m9.figshare.28651985.v1","one-way repeated transfer","20 linked donor-recipient tank pairs, two source lakes","fixed-volume count sample","Preceding donor detection, assigned transfer, subsequent detection; not ancestry or actual transmitted dose","no"),
 ("diazcolunga2022","10.1073/pnas.2111261119",dz.REF,"crossed parental coalescence and isolate coculture","eight parent labels per medium, shared controls and parent origins","relative read composition at endpoint","Common-marker detection and controlled coculture endpoint; no signed growth effect; analytical orientation","no"),
 ("goldman2025","10.1073/pnas.2322440122",gg.source.REF,"crossed parent coalescence at seven nominal ratios","eight parent labels, four household origins, shared parental profiles","relative read composition at p5 endpoint","Stock or earlier unmixed proxy, nominal mixing fraction; not actual transferred taxon counts or temporal order","no")]
def read(p):return json.loads((VAL/p).read_text())
def check(a,b,label):
    if not abs(float(a)-float(b))<TOL:raise AssertionError((label,a,b))
def local_key(s):
    return hashlib.sha256(str(s).encode("utf-8")).hexdigest()
def detect(x,t):return x>0 if t==0 else x>=t
def rows_csv(con,table,out):
    cur=con.execute('SELECT * FROM "'+table+'"')
    with (out/(table+".csv")).open("w",newline="",encoding="utf-8") as f:
        w=csv.writer(f);w.writerow([r[0] for r in cur.description]);w.writerows(cur)
class Registry:
    def __init__(self,con):
        self.c=con;self.taxa=set()
    def lane(self,l,s,res,t,proxy,measurement,risk,role):
        self.c.execute("INSERT INTO opportunity_lanes VALUES (?,?,?,?,?,?,?,?)",(l,s,res,t,proxy,measurement,risk,role))
    def event(self,l,s,e,context,dep,target,focal,partner,pre,post,exposure,orientation,q,weight,keys,labels,res,pn,pool_n,d,a,y,c=None):
        n=len(keys)
        assert n>0 and len(set(keys))==n and len(labels)==n
        assert all(len(z)==n for z in (pn,d,a,y)) and (c is None or len(c)==n)
        assert np.isfinite(a).all() and np.all((a>=0)&(a<=1+1e-12))
        self.c.execute("INSERT INTO opportunity_events VALUES ("+",".join("?" for _ in range(14))+")",
                       (l,s,e,context,dep,target,focal,partner,pre,post,exposure,orientation,q,weight))
        trows=[];rrows=[]
        for i,k in enumerate(keys):
            tax=(s,k)
            if tax not in self.taxa:
                trows.append((s,k,labels[i],res));self.taxa.add(tax)
            rrows.append((l,s,e,k,int(pn[i]),pool_n,int(d[i]),float(a[i]),int(y[i]),None if c is None else int(c[i])))
        self.c.executemany("INSERT INTO opportunity_taxa VALUES (?,?,?,?)",trows)
        self.c.executemany("INSERT INTO opportunity_risk VALUES (?,?,?,?,?,?,?,?,?,?)",rrows)
    def event_scores(self,l,context=None,dep=None):
        cols="n_risk,event_weight,pool_loss,donor_loss,hybrid_loss,dose_loss,control_pool_loss,control_hybrid_loss"
        query="SELECT "+cols+" FROM opportunity_event_brier WHERE lane=?";args=[l]
        if context is not None:query+=" AND context=?";args.append(context)
        if dep is not None:query+=" AND dependency_group=?";args.append(dep)
        return self.c.execute(query,args).fetchall()
def villarreal(R,root):
    s="villarreal2025";l="villarreal_raw";res="author_taxon_bin"
    refs=vv.acquire(root);mats,audit=vv.load(root,refs)
    prior=read("villarreal2025_donor_pool_summary_20261126.json")
    external=read("villarreal2025_external_summary_20261126.json")
    # JSON roundtrip converts source parser tuple pairs to the archived JSON lists.
    audit=json.loads(json.dumps(audit,allow_nan=False))
    assert audit==external["source_audit"], [k for k in audit if audit[k]!=external["source_audit"].get(k)]
    R.lane(l,s,res,0,"observed_preceding_sample","fixed_volume_count_fraction","all 63 author bins undetected in actual recipient", "primary")
    labels=audit["taxon_codes"];keys=[local_key(k) for k in labels];assert len(keys)==63
    for donor,rec in (("A0","Aa"),("B0","Ab")):
        src=mats[donor];dst=mats[rec]
        for i in range(10):
            a=vv.EXPECTED[donor][i];b=vv.EXPECTED[rec][i]
            for t in range(4):
                pre,post=vv.DATES[t:t+2];interval=pre+"-"+post;context=rec+"|"+interval
                risk=dst[t,i]==0;idx=np.flatnonzero(risk)
                P=(np.delete(src[t],i,axis=0)[:,risk]>0).sum(0);D=src[t,i,risk]>0
                A=src[t,i,risk]/src[t,i].sum();Y=dst[t+1,i,risk]>0
                R.event(l,s,f"{a}:{b}:{interval}",context,f"pair:{a}:{b}",f"tank:{b}:{post}",f"tank:{b}:{pre}",f"tank:{a}:{pre}",
                        pre,post,"no_transfer_placebo" if t==0 else "assigned_transfer","actual_recipient",None,.1,
                        [keys[k] for k in idx],[labels[k] for k in idx],res,P,9,D,A,Y)
        for t in range(4):
            interval=vv.DATES[t]+"-"+vv.DATES[t+1];r=R.event_scores(l,rec+"|"+interval)
            assert len(r)==10
            m=np.array([x[:5] for x in r],float).mean(0);q=prior["results"][rec]["intervals"][interval]
            for i,k in ((2,"pool_loss"),(3,"binary_loss"),(4,"hybrid_loss")):check(m[i],q[k],(s,rec,interval,k))
    return {"source_receipt":"villarreal2025_external_summary_20261126.json","replayed_run":45716249512,"verified_profiles":200,"events":80,"source_files":audit["source_files"]}

def diaz(R):
    s="diazcolunga2022";res="common233"
    prior=read("diazcolunga2022_harmonized_endpoint_20261126.json")
    q=prior["output"];assert q["common_segment_eligible"] and q["common_offset_zero_based"]==0
    data,audit=dz.load()
    raw=dz.get("data/OTU_table.csv")
    rd=csv.reader(io.StringIO(raw.decode("utf-8-sig")));h=next(rd);allrows=list(rd)
    seq=[r[0].upper() for r in allrows];assert len(seq)==audit["n_taxa"]==4328
    assert len(set(seq))==len(seq)
    keys=[("raw",i) for i in range(len(seq))];pa=set();ep=set()
    for m,(p,c,mix) in data.items():
        a=np.flatnonzero(p.sum(0)>0);b=np.flatnonzero(c.sum(0)+np.stack(list(mix.values())).sum(0)>0)
        assert all(len(seq[i])==238 for i in a) and all(len(seq[i]) in (232,233) for i in b)
        pa.update(map(int,a));ep.update(map(int,b))
        offsets=Counter()
        for j in b:
            if len(seq[j])!=233:continue
            for i in a:
                for off in range(6):
                    if seq[i][off:off+233]==seq[j]:offsets[off]+=1
        assert {str(k):v for k,v in offsets.items()}==q["offset_matches"][m]
        assert set(offsets)=={0}
    for i in pa:keys[i]=("common233",seq[i][:233])
    for i in ep:
        if len(seq[i])==233:keys[i]=("common233",seq[i])
    d,n=dh.aggregate(data,keys)
    unique=list(dict.fromkeys(keys));assert n==len(unique)==q["modes"]["common233"]["n_aggregated_keys"]
    marker_keys=[local_key(json.dumps(k,separators=(",",":"))) for k in unique]
    # Obtain actual sample IDs and same-medium physical well linkage from exact metadata.
    meta=list(csv.DictReader(io.StringIO(dz.get("data/metaData.csv").decode("utf-8-sig"))))
    samples={};allneeded=set()
    def one(rows):
        assert len(rows)==1,("sample mapping",len(rows))
        return rows[0]
    for m,(prefix,exp) in dz.EXP.items():
        block=[r for r in meta if r["Experiment"]==exp and r["Carbon"]==m and r["Transfer"]=="7" and r["StabilizingCarb1"]==m]
        assay_exp={"Glutamine":"T7-C3-2","Citrate":"T7-C6-2"}[m]
        ablock=[r for r in meta if r["Experiment"]==assay_exp and r["Carbon"]==m and r["Transfer"]=="7" and r["StabilizingCarb1"]==m]
        assert len(block)==len(ablock)==72
        par=[one([r for r in meta if r["Experiment"]==f"{prefix}-T12-X{x}" and r["Comm1"]==f"Community{x}" and r["Carbon"]==m and r["Transfer"]=="12"])["Sample"] for x in dz.MAP[m]]
        def controls(b,prefix):
            return [one([r for r in b if r["Comm1"]==f"{prefix}{i+1}" and r["Comm2"]=="None"]) for i in range(8)]
        def pairs(b,prefix):
            return {(i,j):one([r for r in b if r["Comm1"]==f"{prefix}{i+1}" and r["Comm2"]==f"{prefix}{j+1}" and r["StabilizingCarb2"]==m]) for i in range(8) for j in range(8) if i!=j}
        c=controls(block,"Community");mx=pairs(block,"Community")
        ac=controls(ablock,"Champion");am=pairs(ablock,"Champion")
        for k in mx:
            assert (mx[k]["Row"],mx[k]["Column"],mx[k]["Sequencing_Run"])==(am[k]["Row"],am[k]["Column"],am[k]["Sequencing_Run"])
        samples[m]=(par,c,mx,ac,am)
        allneeded.update(r["Sample"] for r in ac+list(am.values()))
    assert len(allneeded)==128
    # Assay measurements: same-run relative read signatures, no isolate-level strain deconvolution.
    idx={k:h.index(k) for k in sorted(allneeded)}
    mat={k:np.array([float(r[j]) for r in allrows]) for k,j in idx.items()}
    depth={k:int(v.sum()) for k,v in mat.items()}
    assert all(np.isfinite(v).all() and (v>=0).all() and (v==np.floor(v)).all() for v in mat.values())
    assert min(depth.values())>=1000
    mat={k:v/depth[k] for k,v in mat.items()}
    for m,(par,c,mx,ac,am) in samples.items():
        rawp,rawc,rawmix=data[m]
        mono=np.stack([mat[x["Sample"]] for x in ac])
        for (i,j),r in sorted(am.items()):
            a=mono[i];b=mono[j];v=mat[r["Sample"]]
            mi=(a>=.001)&(b<.001);mj=(b>=.001)&(a<.001)
            u=float(v[mi].sum());w=float(v[mj].sum());rem=1-u-w
            tv=float(np.abs(v-(a+b)/2).sum()/2)
            community_tv=float(np.abs(rawmix[i,j]-(rawc[i]+rawc[j])/2).sum()/2)
            aid=f"{m}:{r['Sample']}";pair_id=f"{m}:{min(i,j)+1}:{max(i,j)+1}"
            R.c.execute("INSERT INTO opportunity_assays VALUES ("+",".join("?" for _ in range(18))+")",
                (aid,s,m,f"Champion{i+1}",f"Champion{j+1}",pair_id,r["Sample"],r["Row"]+r["Column"],ac[i]["Sample"],ac[j]["Sample"],depth[r["Sample"]],int(mi.sum()),int(mj.sum()),u,w,rem,u/(u+w) if u+w>1e-12 else None,tv,"controlled_coculture_relative_read_endpoint"))
            R.c.execute("INSERT INTO opportunity_assay_links VALUES (?,?,?,?,?,?,?)",(aid,s,mx[i,j]["Sample"],c[i]["Sample"],c[j]["Sample"],community_tv,"same_medium_origin_label_and_well_not_taxon_identity"))
    for t,mode,role in ((0,"raw_detected","primary"),(.001,"relative_0.001","sensitivity")):
        l="diaz_common233_"+mode
        R.lane(l,s,res,t,"representative_original_parent_R0_harmonized","relative_read_fraction","regional pre-parent roster undetected in analytical focal parent",role)
        for m,(p,c,mix) in d.items():
            par,ctl,targets,_,_=samples[m]
            dp=detect(p,t);cp=detect(c,t);pool=dp.any(0)
            for i,j in dz.PAIRS:
                dep=f"{m}:{i+1}:{j+1}"
                for well in ((i,j),(j,i)):
                    target=targets[well]["Sample"];z=detect(mix[well],t)
                    for f,o in ((i,j),(j,i)):
                        risk=pool&~dp[f];ix=np.flatnonzero(risk)
                        P=dp[[k for k in range(8) if k not in (f,o)]][:,risk].sum(0)
                        R.event(l,s,f"{target}:focal{f+1}",m,dep,target,par[f],par[o],"T12_parent","T7_postmix","assigned_coalescence","analytical",None,.25,
                            [marker_keys[k] for k in ix],["study_local_marker" for k in ix],res,P,6,dp[o,risk],p[o,risk],z[risk],cp[f,risk])
            r=R.event_scores(l,m);assert len(r)==112
            a=np.asarray([[x[2],x[3],x[4],x[6],x[7]] for x in r],float).mean(0)
            old=q["modes"]["common233"]["results"][mode]["by_medium"][m]
            for k,v in zip(("pool_loss","donor_loss","hybrid_loss"),a[:3]):check(v,old["mixed"][k],(s,mode,m,k))
            check(a[3],old["unmixed_control"]["pool_loss"],(m,"control_pool"))
            check(a[4],old["unmixed_control"]["hybrid_loss"],(m,"control_hybrid"))
    return {"source_receipt":"diazcolunga2022_harmonized_endpoint_20261126.json","replayed_run":45717729776,"verified_detection_profiles":144,
            "additional_assay_profiles":128,"assay_monocultures":16,"assay_mixed_wells":112,"assay_reads_range":[min(depth.values()),max(depth.values())],"source_files":audit["source_files"]}
def goldman(R):
    s="goldman2025";res="author_ASV"
    mats,lanes,audit=gg.load()
    prior=read("goldman2025_dose_summary_20261126.json")
    assert audit==prior["audit"] and len(mats)==167 and len(lanes)==2
    for lane,pre in sorted(lanes.items()):
        for mode,t in (("threshold_0.001",.001),("raw_detection",0)):
            l="goldman_"+lane+"_"+mode
            R.lane(l,s,res,t,lane,"relative_read_fraction","regional source-proxy roster undetected in analytical focal proxy","primary" if lane=="p15_stock_proxy" and t==.001 else "sensitivity")
            ks=sorted(set().union(*(set(v) for v in pre.values())))
            ids=[local_key(k) for k in ks];pp={p:np.array([pre[p].get(k,0.) for k in ks]) for p in gg.PARENTS}
            det={p:detect(v,t) for p,v in pp.items()};regional=np.any(np.stack(list(det.values())),axis=0)
            old=prior["lanes"][lane]["results"][mode]["source_disjoint"]
            assert old["n_scored_pairs"]==8 and int(regional.sum())==prior["lanes"][lane]["results"][mode]["source_roster"]
            for combo in gg.COMBOS:
                a,b=combo.split("-")
                for ratio in gg.RATIOS:
                    g=[(k,v) for k,v in sorted(mats.items()) if k[0]==combo and k[2]==ratio]
                    assert len(g) in (2,3)
                    ra,rb=map(float,ratio.split(":"))
                    for k,v in g:
                        target=f"{combo}:p5:{ratio}:rep{k[3]}"
                        y=detect(np.array([v.get(x,0.) for x in ks]),t)
                        for f,o,q in ((a,b,rb/(ra+rb)),(b,a,ra/(ra+rb))):
                            risk=regional&~det[f];idx=np.flatnonzero(risk)
                            P=np.stack([det[p] for p in gg.PARENTS if p not in (f,o)])[:,risk].sum(0)
                            R.event(l,s,target+":"+f,combo,combo,target,f"{lane}:{f}",f"{lane}:{o}",
                                "p15_prior_stock_experiment" if lane=="p15_stock_proxy" else "p3_separate_unmixed_culture","p5_mixed_endpoint","assigned_coalescence","analytical",q,1/(7*len(g)*2),
                                [ids[i] for i in idx],["study_local_ASV" for i in idx],res,P,6,det[o][risk],pp[o][risk],y[risk])
                rr=R.event_scores(l,dep=combo);assert len(rr) in (40,42)
                a=np.array([x[:6] for x in rr],float)
                check(a[:,1].sum(),1,(combo,"weight"))
                assert int(a[:,0].sum())==old["by_pair"][combo]["n_risk_cells"]
                for j,k in ((2,"fixed_pool"),(4,"fixed_identity"),(5,"fixed_dose")):
                    check(np.sum(a[:,1]*a[:,j]),old["by_pair"][combo]["fixed"][k],(s,l,combo,k))
    return {"source_receipt":"goldman2025_dose_summary_20261126.json","replayed_run":45718594570,"verified_endpoint_profiles":167,
      "stock_proxy_profiles":24,"earlier_control_profiles":14,"source_files":audit["source_files"]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);p.add_argument("--scratch",type=Path,required=True);a=p.parse_args()
    subprocess.run([sys.executable,"-m","benchmarks.build_ahif_mechanism_registry_v05","--out",str(a.out)],check=True)
    con=sqlite3.connect(a.out/"ahif_mechanism_v0.sqlite");con.execute("PRAGMA foreign_keys=ON")
    prior_tables={k:con.execute('SELECT COUNT(*) FROM "'+k+'"').fetchone()[0] for (k,) in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    con.executescript((DB/"opportunity_schema_v06.sql").read_text())
    con.executemany("INSERT INTO opportunity_sources VALUES (?,?,?,?,?,?,?,?)",SOURCES);reg=Registry(con)
    receipts={"villarreal2025":villarreal(reg,a.scratch/"villarreal"),"diazcolunga2022":diaz(reg),"goldman2025":goldman(reg)}
    con.commit()
    assert con.execute("PRAGMA integrity_check").fetchone()[0]=="ok" and con.execute("PRAGMA foreign_key_check").fetchall()==[]
    assert all(con.execute('SELECT COUNT(*) FROM "'+k+'"').fetchone()[0]==v for k,v in prior_tables.items())
    assert con.execute("SELECT COUNT(*) FROM opportunity_sources").fetchone()[0]==3
    assert con.execute("SELECT COUNT(*) FROM opportunity_lanes").fetchone()[0]==7
    assert con.execute("SELECT COUNT(*) FROM opportunity_assays").fetchone()[0]==112
    assert con.execute("SELECT COUNT(*) FROM opportunity_assay_links").fetchone()[0]==112
    assert con.execute("SELECT COUNT(*) FROM opportunity_signed_interaction_effects").fetchone()[0]==0
    coverage=con.execute("""SELECT e.lane,COUNT(*),COUNT(DISTINCT physical_target),SUM(b.n_risk)
        FROM opportunity_events e JOIN opportunity_event_brier b USING(lane,event_id,source) GROUP BY e.lane ORDER BY e.lane""").fetchall()
    for l,n,target,risk in coverage:
        assert n== (80 if l=="villarreal_raw" else 224 if l.startswith("diaz") else 334)
    assert sum(x[1] for x in coverage)==1864
    # Full data in derived artifact; brief ledger has counts and reproducible source receipts only.
    tables=[x[0] for x in con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name LIKE 'opportunity_%' ORDER BY name")]
    views=["opportunity_detection_summary","opportunity_pairwise_assay_bridge","opportunity_event_brier"]
    for t in tables+views:rows_csv(con,t,a.out)
    counts={t:con.execute('SELECT COUNT(*) FROM "'+t+'"').fetchone()[0] for t in tables}
    assay=con.execute("""SELECT medium,COUNT(*),AVG(monoculture_reference_tv),AVG(community_reference_tv),AVG(unresolved_mass),SUM(signature_share_i IS NOT NULL)
       FROM opportunity_assays JOIN opportunity_assay_links USING(assay_id,source) GROUP BY medium ORDER BY medium""").fetchall()
    primary=con.execute("""SELECT lane,COUNT(*) FROM opportunity_events JOIN opportunity_lanes USING(lane,source) WHERE role='primary' GROUP BY lane""").fetchall()
    counts["opportunity_primary_events"]=sum(x[1] for x in primary)
    con.close()
    z={"status":"VALIDATED_DERIVED_OPPORTUNITY_REGISTRY_NOT_SIGNED_NETWORK_NOT_META","version":"0.6.0","counts":counts,
      "coverage":[dict(zip(("lane","events","physical_targets_within_lane","dependent_risk_rows"),r)) for r in coverage],
      "assay_summary":[dict(zip(("medium","wells","mean_isolate_reference_tv","mean_community_reference_tv","mean_unresolved_signature_mass","defined_signature_shares"),r)) for r in assay],
      "source_receipts":receipts,"prior_table_counts_unchanged":True,"integrity":"ok","foreign_key_check":"ok","meta_eligible":0,
      "notes":["Three existing endpoint study families, no newly recruited publication.","Risk rows and analytical orientations dependent; thresholds/proxies duplicate physical targets.","Coculture relative-read endpoint and linked original community pair, not signed pairwise coefficient or ancestry."]}
    (a.out/"opportunity_audit.json").write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    base=json.loads((a.out/"registry_audit.json").read_text());base["version"]="0.6.0";base["counts"].update(counts);base["opportunity_audit"]="opportunity_audit.json"
    (a.out/"registry_audit.json").write_text(json.dumps(base,indent=2)+"\n")
    print("AHIF_OPPORTUNITY_V06_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
