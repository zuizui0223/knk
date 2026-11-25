"""Prescored Portal capture-response allocation pilot; single site, no causal claims."""
import argparse
from pathlib import Path
import hashlib, json
import itertools
import numpy as np
import pandas as pd
SOURCE = "a4dc266c0068d26370cd88bec5272a02623d5470"
FILES = {
    "captures": ("Rodents/Portal_rodent.csv","b2c6da8e3c18ca93a2d57083c09f2eb39ed28bb5"),
    "traps": ("Rodents/Portal_rodent_trapping.csv","2a4d18388a4296dd81e8bc16f8dadd92a63e0214"),
    "species": ("Rodents/Portal_rodent_species.csv","bcecdbeab3fb364a7712d17486fa12b5ceda835f"),
    "plots": ("SiteandMethods/Portal_plots.csv","20a4b46d97509787d510df60d2bd898a2fe0ad8a"),
    "moons": ("Rodents/moon_dates.csv","9f48e59d38514f135a97c06f200e33efb1ccba37"),
}
def norm(a):
    a=np.asarray(a,dtype=float); s=a.sum()
    return a/s if s>0 else None
def mix(a,b):
    if a is None: return None
    return a if b is None else (a+b)/2
def tv(a,b):
    return float(np.abs(a-b).sum()/2) if a is not None and b is not None else np.nan
def overlap(a,b):
    if a.sum()<=0 or b.sum()<=0: return np.nan
    aa=a==a.max(); bb=b==b.max()
    return float((aa&bb).sum()/(aa.sum()*bb.sum()))
def forecast(A,B,start,end,part,orientation):
    """A/B maps integer moon to plot-by-taxon counts, full balanced roster."""
    rows=[]
    for k in range(start,end-2,4):
        if not all(t in A for t in range(k,k+4)): continue
        older=[t for t in range(k-12,k) if t in A]
        pairs=[(t,t+1) for t in range(k-12,k,2) if t in A and t+1 in A]
        if len(older)<8 or len(pairs)<3: continue
        meanA={t:A[t].mean(axis=0) for t in older+list(range(k,k+4))}
        da=np.abs(meanA[k+1]-meanA[k])
        db=B[k+3].mean(axis=0)-B[k+2].mean(axis=0)
        dx=A[k+3].mean(axis=0)-A[k+2].mean(axis=0)
        gross=np.abs(db).sum()
        q=norm(np.abs(db))
        mu=np.mean([meanA[t] for t in older],axis=0)
        var=np.mean([A[t].var(axis=0,ddof=1) for t in older],axis=0)
        prior=np.array([np.abs(meanA[v]-meanA[u]) for u,v in pairs])
        H=norm(prior.mean(axis=0)); V=norm(np.sqrt(var)); L=norm(da)
        M1=norm(mu); M2=norm(np.sqrt(mu))
        pred={"M1":M1,"M2":M2,"V":V,"H":H,"L":L,
              "AH":mix(M2,L),"VH":mix(V,L),"HH":mix(H,L)}
        scores={name:tv(p,q) for name,p in pred.items()}
        alternatives=[tv(mix(H,norm(p)),q) for p in prior]
        alt=float(np.mean(alternatives)) if all(np.isfinite(alternatives)) else np.nan
        sg=(np.dot(np.sign(dx),db)+np.dot(np.sign(db),dx))/2
        sn=(np.sign(dx.sum())*db.sum()+np.sign(db.sum())*dx.sum())/2
        rows.append(dict(partition=part,orientation=orientation,lag_start=k,
          target_end=k+3,older_periods=len(older),older_pairs=len(pairs),gross=float(gross),
          lag_defined=L is not None,dominant_overlap=overlap(meanA[k+1],np.abs(db)),
          lag_overlap=overlap(da,np.abs(db)),cancel=float(1-abs(db.sum())/gross) if gross>0 else np.nan,
          signed_gap=float(sg-sn),H_oldpair_mix=alt,**scores))
    return rows
def summarize(rows):
    cols=["M1","M2","V","H","L","AH","VH","HH","H_oldpair_mix"]
    rec=[]
    for (p,o),g in rows.groupby(["partition","orientation"],sort=True):
        good=g[np.isfinite(g[["H","HH","V","VH","M2","AH"]]).all(axis=1)]
        d={"partition":p,"orientation":o,"windows":len(g),"paired_scored":len(good),
           "lag_defined":int(g.lag_defined.sum())}
        for c in cols+["dominant_overlap","lag_overlap","cancel","signed_gap"]:
            d[c]=float(good[c].mean()) if len(good) else np.nan
        for key,a,b in [("history_gain","H","HH"),("variance_gain","V","VH"),
                        ("sqrt_gain","M2","AH"),("locality_gain","H_oldpair_mix","HH")]:
            dd=(good[a]-good[b]).dropna()
            d[key]=float(dd.mean());d[key+"_wins"]=int((dd>0).sum());d[key+"_n"]=len(dd)
        rec.append(d)
    return pd.DataFrame(rec)
def selftest():
    assert tv(norm([1,2]),norm([2,1]))==1/3 or abs(tv(norm([1,2]),norm([2,1]))-1/3)<1e-12
    a={t:np.array([[5+t%3,4],[6+t%3,3]]) for t in range(1,30)}
    b={t:np.array([[3+t%3,6],[4+t%3,5]]) for t in range(1,30)}
    r=forecast(a,b,13,28,0,0)
    assert len(r)==4
    b2={t:x.copy() for t,x in b.items()};b2[16][:]=[200,0]
    assert forecast(a,b2,13,28,0,0)[0]["gross"]!=r[0]["gross"]
    assert all(0<=r[0][k]<=1 for k in ["H","V","HH","VH"])
    print("Synthetic checks PASS")
def load(root):
    ds={}; hashes={}
    for name,(file,sha) in FILES.items():
        p=root/file; data=p.read_bytes()
        assert hashlib.sha1(b"blob "+str(len(data)).encode()+b"\0"+data).hexdigest()==sha,file
        hashes[name]=hashlib.sha256(data).hexdigest()
        d=pd.read_csv(p,low_memory=False);d.columns=d.columns.str.lower()
        ds[name]=d
        print(name,len(d),list(d.columns))
    return ds,hashes
def main():
    p=argparse.ArgumentParser();p.add_argument("root",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True);selftest()
    ds,hashes=load(a.root)
    c,tr,sp,pl,mo=(ds[k] for k in ["captures","traps","species","plots","moons"])
    required=[(c,["period","plot","species","id","note1"]),(tr,["period","plot","sampled","effort","qcflag"]),
        (sp,["speciescode","censustarget","unidentified","rodent"]),(pl,["year","month","plot","treatment","resourcetreatment","anttreatment"]),
        (mo,["newmoonnumber","newmoondate","period"])]
    for d,cols in required: assert set(cols)<=set(d.columns),str(set(cols)-set(d.columns))
    pl=pl[pl.year.between(2016,2024)].copy()
    assert not pl.duplicated(["year","month","plot"]).any()
    cond=(pl.treatment.str.lower()=="control")&(pl.resourcetreatment.str.lower()=="none")&(pl.anttreatment.str.lower()=="none")
    pl["eligible"]=cond
    pc=pl.groupby("plot").agg(n=("month","size"),ok=("eligible","all"))
    roster=sorted(pc.index[(pc.n==108)&pc.ok].astype(int).tolist())
    assert len(roster)>=4, f"too few constant controls {roster}; labels {pl.treatment.unique()}"
    roster=roster[:len(roster)//2*2]; print("Fixed control plots",roster)
    spp=sorted(sp.loc[(sp.rodent==1)&(sp.censustarget==1)&(sp.unidentified==0),"speciescode"].unique())
    mo=mo.copy();mo["date"]=pd.to_datetime(mo.newmoondate)
    mo=mo[mo.date.between("2016-01-01","2024-12-31")].copy()
    assert not mo.newmoonnumber.duplicated().any()
    mp=mo[mo.period.notna()&(mo.period>0)][["period","newmoonnumber"]].copy()
    assert not mp.period.duplicated().any()
    tr=tr[tr.plot.isin(roster)&tr.period.isin(mp.period)].copy()
    # stop rather than resolving duplicated trap effort records
    assert not tr.duplicated(["plot","period"]).any()
    tr=tr[(tr.sampled==1)&(tr.effort==49)&(tr.qcflag==1)]
    complete=tr.groupby("period").plot.nunique()
    periods=complete.index[complete==len(roster)]
    tr=tr[tr.period.isin(periods)]
    raw=c[c.period.isin(periods)&c.plot.isin(roster)].copy()
    raw=raw[raw.species.isin(spp)&~pd.to_numeric(raw.note1,errors="coerce").isin([2,4])]
    if "qcflag" in raw: raw=raw[raw.qcflag==1]
    ids=raw.id.astype("string").str.strip()
    missing=ids.isna()|ids.isin(["","NA","nan","0"])
    keyed=raw.loc[~missing].copy()
    ambiguous=keyed.groupby(["period","plot","id"]).species.nunique()
    bad=ambiguous[ambiguous>1].index
    if len(bad):
        keyed=keyed[~pd.MultiIndex.from_frame(keyed[["period","plot","id"]]).isin(bad)]
    n_keyed=len(keyed);keyed=keyed.drop_duplicates(["period","plot","id"])
    counted=pd.concat([keyed,raw[missing]],ignore_index=True)
    ct=counted.groupby(["period","plot","species"]).size()
    grid=pd.MultiIndex.from_product([periods,roster,spp],names=["period","plot","species"])
    matrix=ct.reindex(grid,fill_value=0).to_numpy().reshape(len(periods),len(roster),len(spp))
    moon_map=dict(zip(mp.period,mp.newmoonnumber))
    mats={int(moon_map[t]):matrix[i] for i,t in enumerate(periods)}
    start,end=int(mo.newmoonnumber.min()),int(mo.newmoonnumber.max())
    base=np.arange(len(roster)); n=len(base)//2
    splits=[(base[::2],base[1::2])]
    rng=np.random.default_rng(20261125)
    for j in range(3):
        order=rng.permutation(base);splits.append((order[:n],order[n:]))
    rows=[]
    for part,(left,right) in enumerate(splits):
        for orientation,(i,j) in enumerate([(left,right),(right,left)]):
            rows+=forecast({t:v[i] for t,v in mats.items()},{t:v[j] for t,v in mats.items()},start,end,part,orientation)
    audit=dict(source_commit=SOURCE,sha256=hashes,raw_capture_rows=len(c),source_trap_rows=len(ds["traps"]),
        control_plots=roster,species_roster=spp,moons_in_scope=len(mo),complete_periods=len(periods),
        eligible_capture_rows=len(raw),missing_id_rows=int(missing.sum()),ambiguous_ids=len(bad),
        removed_duplicate_id_rows=n_keyed-len(keyed),counted_rows=len(counted),scope=[start,end],
        partition_plots=[[[roster[i] for i in x],[roster[i] for i in y]] for x,y in splits])
    (a.out/"audit.json").write_text(json.dumps(audit,indent=2))
    print(json.dumps(audit,indent=2))
    if not rows: raise SystemExit("STOP: no eligible complete disjoint windows")
    out=pd.DataFrame(rows);out.to_csv(a.out/"window_diagnostics.csv",index=False)
    summary=summarize(out);summary.to_csv(a.out/"summary.csv",index=False)
    print(summary.to_string(index=False))
if __name__=="__main__": main()
