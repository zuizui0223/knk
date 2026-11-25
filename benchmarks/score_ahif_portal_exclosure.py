"""Source-locked Portal held-plot response-direction pilot; no causal treatment inference."""
import argparse,hashlib,json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from benchmarks.score_ahif_replication_scale import scores,METRICS
SOURCE="a4dc266c0068d26370cd88bec5272a02623d5470"
FILES={
 "captures":("Rodents/Portal_rodent.csv","b2c6da8e3c18ca93a2d57083c09f2eb39ed28bb5"),
 "traps":("Rodents/Portal_rodent_trapping.csv","2a4d18388a4296dd81e8bc16f8dadd92a63e0214"),
 "species":("Rodents/Portal_rodent_species.csv","bcecdbeab3fb364a7712d17486fa12b5ceda835f"),
 "plots":("SiteandMethods/Portal_plots.csv","20a4b46d97509787d510df60d2bd898a2fe0ad8a"),
 "moons":("Rodents/moon_dates.csv","9f48e59d38514f135a97c06f200e33efb1ccba37"),
 "history":("SiteandMethods/Portal_plot_treatments.csv","deb0c0615291074d0c3f8c16f5b24d354017f6bd"),
 "methods":("SiteandMethods/Methods.md","c196c58d5a9a71685fe739ba9af127bed455197d"),
 "flags":("Rodents/Portal_rodent_datanotes.csv","47868b824fbeebd54b90f6718d3224b972e11b0a"),
}
C=[4,5,6,7,11,13,14,17,18,24];E=[2,3,8,15,19,20,21,22];PLOTS=C+E
SEED=2026112511;N=256;SIZES=(1,2,3,4)
BAD_PERIOD={4,5,6,10,14,18};BAD_ROW={2,4,5,6,7,8,9,10,11,13,14,15,18}
def load(root):
    ds={};h={}
    for k,(name,sha) in FILES.items():
        p=root/name;b=p.read_bytes()
        assert hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==sha,(name,"git hash")
        h[name]=hashlib.sha256(b).hexdigest()
        if p.suffix==".csv":
            d=pd.read_csv(p,low_memory=False);d.columns=d.columns.str.lower()
            if k in ("species","captures"):
                col="speciescode" if k=="species" else "species"
                d[col]=pd.read_csv(p,usecols=[col],keep_default_na=False,dtype=str)[col]
            ds[k]=d
    return ds,h
def prepare(root):
    ds,h=load(root);c,tr,sp,pl,mo,hi,fl=(ds[k] for k in ("captures","traps","species","plots","moons","history","flags"))
    assert set(fl.note1.astype(int))==set(range(1,19))
    hi=hi.set_index("plot")
    col="apr2015_present"
    assert sorted(hi.index[hi[col]=="Unmanipulated control"].tolist())==C
    assert sorted(hi.index[hi[col]=="Kangaroo rats removed"].tolist())==E
    pl=pl[pl.year.between(2016,2024)&pl["plot"].isin(PLOTS)].copy()
    assert len(pl)==18*108 and not pl.duplicated(["year","month","plot"]).any()
    assert pl.resourcetreatment.isna().all() and pl.anttreatment.isna().all()
    assert set(pl[pl["plot"].isin(C)].treatment)=={"control"} and set(pl[pl["plot"].isin(E)].treatment)=={"exclosure"}
    spp=sorted(sp.loc[(sp.rodent==1)&(sp.censustarget==1)&(sp.unidentified==0),"speciescode"].unique())
    removed=sorted(sp.loc[(sp.rodent==1)&(sp.censustarget==1)&(sp.unidentified==0)&sp.scientificname.str.startswith("Dipodomys "),"speciescode"].unique())
    assert removed==["DM","DO","DS"] and "NA" in spp
    mo=mo.copy();mo["date"]=pd.to_datetime(mo.newmoondate);mo=mo[mo.date.between("2016-01-01","2024-12-31")].copy()
    mp=mo[mo.period.notna()&(mo.period>0)][["period","newmoonnumber"]]
    assert not mp.period.duplicated().any() and not mp.newmoonnumber.duplicated().any()
    tr=tr[tr["plot"].isin(PLOTS)&tr.period.isin(mp.period)].copy()
    assert not tr.duplicated(["plot","period"]).any()
    tr=tr[(tr.sampled==1)&(tr.effort==49)&(tr.qcflag==1)]
    full=tr.groupby("period")["plot"].nunique()
    full=set(full.index[full==18])
    r=c[c.period.isin(full)&c["plot"].isin(PLOTS)].copy()
    notes=pd.to_numeric(r.note1,errors="coerce")
    bad=set(r.loc[notes.isin(BAD_PERIOD),"period"])
    eligible=sorted(full-bad)
    assert eligible
    raw=r[r.period.isin(eligible)&r.species.isin(spp)&~notes.isin(BAD_ROW)].copy()
    if "qcflag" in raw:raw=raw[raw.qcflag==1]
    sid=raw.id.astype("string").str.strip()
    miss=sid.isna()|sid.isin(["","NA","nan","0"])
    keyed=raw[~miss].copy();keyed["sid"]=sid[~miss]
    ambigu=keyed.groupby(["period","plot","sid"]).species.nunique()
    ambiguous=ambigu[ambigu>1].index
    if len(ambiguous):keyed=keyed[~pd.MultiIndex.from_frame(keyed[["period","plot","sid"]]).isin(ambiguous)]
    before=len(keyed);keyed=keyed.drop_duplicates(["period","plot","sid"])
    cross=keyed.groupby(["period","sid"])["plot"].nunique();multi=cross[cross>1].index
    cross_mask=pd.MultiIndex.from_frame(keyed[["period","sid"]]).isin(multi)
    cnt=pd.concat([keyed,raw[miss]],ignore_index=True)
    cnt2=pd.concat([keyed[~cross_mask],raw[miss]],ignore_index=True)
    period_to_moon={int(p):int(t) for p,t in mp.itertuples(index=False,name=None)}
    moons=[period_to_moon[p] for p in eligible];order=np.argsort(moons)
    moonlist=[moons[i] for i in order]
    loc={t:i for i,t in enumerate(moonlist)}
    pairs=[(loc[t-1],loc[t]) for t in moonlist if t-1 in loc and t>=int(mo.newmoonnumber.min())+2]
    if not pairs:raise ValueError("STOP no consecutive full-effort lunar pairs")
    def mat(rows):
        v=rows.groupby(["period","plot","species"]).size()
        grid=pd.MultiIndex.from_product([eligible,PLOTS,spp],names=["period","plot","species"])
        return v.reindex(grid,fill_value=0).to_numpy().reshape(len(eligible),18,len(spp))[order].astype(float)
    audit=dict(source_commit=SOURCE,hashes=h,raw_capture_rows=len(c),control=C,exclosure=E,treatment_history={str(k):{c:(None if pd.isna(v) else str(v)) for c,v in hi.loc[k].items()} for k in PLOTS},
               n_moons_scope=len(mo),n_full_effort_periods=len(full),bad_flag_periods=sorted(map(int,bad)),
               n_eligible_periods=len(eligible),n_consecutive_pairs=len(pairs),eligible_moons=moonlist,
               pairs=[[moonlist[i],moonlist[j]] for i,j in pairs],species=spp,excluded_target_species=removed,
               n_eligible_capture_rows=len(raw),n_missing_id_rows=int(miss.sum()),n_ambiguous_plot_period_ids=len(ambiguous),
               n_duplicate_rows_removed=before-len(keyed),n_counted=len(cnt),n_crossplot_individual_periods=len(multi),n_crossplot_rows=int(cross_mask.sum()),n_counted_no_crossplot=len(cnt2))
    return mat(cnt),mat(cnt2),spp,removed,pairs,audit
def score_lane(x,pairs,oc,oe):
    vals={};oks={}
    lag=[i for i,j in pairs];target=[j for i,j in pairs]
    for n in SIZES:
        ds=[];hs=[]
        for o in (0,1):
            ic=oc[:,:n] if o==0 else oc[:,-n:];jc=oc[:,-n:] if o==0 else oc[:,:n]
            ie=oe[:,:n] if o==0 else oe[:,-n:];je=oe[:,-n:] if o==0 else oe[:,:n]
            assert not np.any((ic[:,:,None]==jc[:,None,:])) and not np.any((ie[:,:,None]==je[:,None,:]))
            # time x iteration x sample x species -> iteration x time x species
            d=(x[:,len(C)+ie,:].mean(2)-x[:,ic,:].mean(2)).transpose(1,0,2)
            h=(x[:,len(C)+je,:].mean(2)-x[:,jc,:].mean(2)).transpose(1,0,2)
            ds.append(d);hs.append(h)
        d=np.stack(ds,axis=1);h=np.stack(hs,axis=1)
        v,ok=scores(d[:,:,lag],d[:,:,target],h[:,:,target])
        vals[n]=v;oks[n]=ok
    common=np.logical_and.reduce([oks[n] for n in SIZES])
    def summ(v,mask):
        dates=[i for i in range(v.shape[2]) if mask[:,:,i].any()]
        if not dates:return None
        row=dict(zip(METRICS,np.mean([v[:,:,i][mask[:,:,i]].mean(0) for i in dates],axis=0).tolist()))
        o=1-row["allocation_tv"]
        row.update(O=o,Astar=.5+row["H"]/o if o>1e-12 else None,J=row["H"]-row["G"],n_dates=len(dates),n_defined=int(mask.sum()))
        if row["Astar"] is not None:assert -1e-10<=row["Astar"]<=1+1e-10
        assert row["J"]>=-1e-10
        return row
    return {str(n):dict(n_total=int(oks[n].size),common=summ(vals[n],common),all_defined=summ(vals[n],oks[n])) for n in SIZES}
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    x,x2,spp,rem,pairs,audit=prepare(a.root)
    rng=np.random.default_rng(SEED);oc=np.array([rng.permutation(len(C)) for _ in range(N)]);oe=np.array([rng.permutation(len(E)) for _ in range(N)])
    mask=[i for i,s in enumerate(spp) if s not in rem]
    result=dict(status="PORTAL_SINGLE_SITE_NON_TARGET_DESCRIPTIVE_NOT_META",seed=SEED,n_splits=N,audit=audit,
      primary_nontarget=score_lane(x[:,:,mask],pairs,oc,oe),all_species=score_lane(x,pairs,oc,oe),
      nontarget_no_crossplot_ID=score_lane(x2[:,:,mask],pairs,oc,oe))
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,allow_nan=False)+"\n")
    print("PORTAL_EXCLOSURE_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
