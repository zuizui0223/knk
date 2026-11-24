"""Source-locked independent replay of bacterial community terminal states.
No introduced-species intervention; prediction only, not treatment moderation.
"""
from __future__ import annotations
import argparse,csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
import numpy as np
from benchmarks.score_ahif_saine2025_succession import holm,SEED

METADATA="metadata.tsv"
MATRIX="seqtable.csv"
BLOBS={METADATA:"7394c6545e0a138ba1ba54a259580d8c63c7b505",
        MATRIX:"eaa43f867b0df78d124c742ca0b23f971c099ea4"}
EXPERIMENTS=["0D","7D_rep1","7D_rep2","7D_rep3","7D_rep4"]
LOCATIONS={"AshridgeEstate":11,"BurhamBeeches":32,"BeechWood":6,"LongWood":11,
 "OldWindsorWood":14,"PullingshillWood":21,"SilwoodPark":1,"StAlbansCity":14,
 "CombeHill":21,"WindsorForest":6,"GreenbroomCovert":50,"WindsorGreatPark":1,
 "WaterersPark":14,"WythamWoods":73}

def gitblob(b):
    return hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()

def get_meta(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        m=list(csv.DictReader(f,delimiter="\t"))
    if len(m)!=1375 or len({r["sampleid"] for r in m})!=1375:
        raise ValueError("matched metadata count")
    parents=defaultdict(dict)
    for r in m:
        key=r["parent"];e=r["Experiment"]
        if e not in EXPERIMENTS or e in parents[key]:
            raise ValueError("unmatched or duplicate trajectory")
        parents[key][e]=r
    ids=sorted(parents)
    if len(ids)!=275:raise ValueError("source parent count")
    locations=[];sampleids=[]
    for i in ids:
        g=parents[i]
        if set(g)!=set(EXPERIMENTS) or len({r["Location"] for r in g.values()})!=1:
            raise ValueError("incomplete or location-inconsistent parent")
        locations.append(g["0D"]["Location"])
        sampleids.extend(g[e]["sampleid"] for e in EXPERIMENTS)
    if Counter(locations)!=LOCATIONS:raise ValueError("source collection location counts")
    return ids,np.array(locations),sampleids

def read_matrix(path,want):
    # Author file may be tab-separated despite .csv suffix. Infer only
    # delimiter/orientation from exact matched sample ID set, not outcomes.
    with path.open(encoding="utf-8-sig",newline="") as f:
        head=f.readline()
        sep=max(("\t",",",";"),key=head.count)
        f.seek(0)
        rows=list(csv.reader(f,delimiter=sep))
    if not rows or len(rows[0])<2:raise ValueError("empty matrix")
    if len({len(z) for z in rows})!=1:raise ValueError("nonrectangular matrix")
    columns=rows[0][1:];rowids=[z[0] for z in rows[1:]]
    if len(set(columns))!=len(columns) or len(set(rowids))!=len(rowids):
        raise ValueError("duplicate source matrix identifiers")
    w=set(want);ascol=w.issubset(columns);asrow=w.issubset(rowids)
    if ascol==asrow:raise ValueError("exact sample orientation unresolved")
    if asrow:
        by={name:i for i,name in enumerate(rowids)}
        x=np.asarray([rows[by[s]+1][1:] for s in want],float)
        names=columns;orientation="samples_by_asv"
    else:
        by={name:i for i,name in enumerate(columns)}
        i=[by[s]+1 for s in want]
        x=np.asarray([[row[j] for j in i] for row in rows[1:]],float).T
        names=rowids;orientation="asv_by_samples"
    if not np.isfinite(x).all() or np.any(x<0) or np.any(x!=np.floor(x)) or np.any(x.sum(1)<=0):
        raise ValueError("noninteger/invalid source counts or empty matched state")
    return x.reshape(275,5,-1),{"orientation":orientation,"original_rows":len(rowids),
      "original_columns":len(columns),"n_asvs":len(names),"separator":repr(sep)}

def features_and_outcomes(x):
    if x.ndim!=3 or x.shape[1]!=5 or np.any(x.sum(2)<=0):
        raise ValueError("trajectory shape or zero row")
    b=x[:,0];n=b.sum(1);p=b/n[:,None]
    f=np.column_stack((np.log1p(n),np.log1p((b>0).sum(1)),np.sum(p*p,1),p.max(1)))
    end=x[:,1:]/x[:,1:].sum(2,keepdims=True)
    named=end.mean(1);rank=np.sort(end,axis=2)[:,:,::-1].mean(1)
    # Zero trailing columns do not change the geometry.
    k=int(np.max(np.count_nonzero(x[:,1:],axis=2)))
    return f,{"named":named,"rank":rank[:,:k]}

def maps(f,locations):
    n=len(f);m0=np.zeros((n,n));m1=np.zeros((n,n))
    for site in np.unique(locations):
        tr=np.flatnonzero(locations!=site);te=np.flatnonzero(locations==site)
        mu=f[tr].mean(0);sd=f[tr].std(0);z=(f-mu)/np.where(sd>1e-12,sd,1.)
        base=np.column_stack((np.ones(n),z[:,0],z[:,0]**2))
        full=np.column_stack((base,z[:,1:]))
        for a,H in ((base,m0),(full,m1)):
            penalty=np.ones(a.shape[1]);penalty[0]=0.
            H[np.ix_(te,tr)]=a[te]@np.linalg.solve(a[tr].T@a[tr]+np.diag(penalty),a[tr].T)
    return m0,m1

def restricted(f,locations,y):
    d=(f[:,0]-f[:,0].mean())/max(f[:,0].std(),1e-12)
    a=np.column_stack(((locations[:,None]==np.unique(locations)).astype(float),d,d*d))
    if np.linalg.matrix_rank(a)!=a.shape[1]:raise ValueError("restricted design singular")
    pinv=np.linalg.pinv(a);pred=a@(pinv@y);h=np.sum(a*pinv.T,axis=1)
    resid=y-pred
    sat=h>=1-1e-10
    if np.any(np.abs(resid[sat])>1e-8):raise ValueError("nonzero saturated residual")
    hc2=np.divide(resid,np.sqrt(np.maximum(1-h,0))[:,None],out=np.zeros_like(resid),where=~sat[:,None])
    return pred,hc2,int(sat.sum())

def gain(a,b):
    return float((a-b)/a) if a>1e-15 else 0.

def score(f,locations,outcomes,nboot=999):
    H0,H1=maps(f,locations);n=len(f)
    E0=np.eye(n)-H0;E1=np.eye(n)-H1
    Q0=E0.T@E0;Q1=E1.T@E1
    rng=np.random.default_rng(SEED+2301)
    signs=rng.choice([-1.,1.],size=(nboot,n))
    result=[]
    for key,y in outcomes.items():
        err0=np.sum((E0@y)**2,axis=1);err1=np.sum((E1@y)**2,axis=1)
        a,b=float(err0.sum()),float(err1.sum());g=gain(a,b)
        fits,res,sat=restricted(f,locations,y)
        ff=fits@fits.T;fr=fits@res.T;rr=res@res.T
        null=[]
        for s in signs:
            cross=fr*s[None,:]
            gram=ff+cross+cross.T+rr*s[:,None]*s[None,:]
            null.append(gain(np.sum(Q0*gram),np.sum(Q1*gram)))
        null=np.asarray(null)
        folds=[]
        for site in np.unique(locations):
            mask=locations==site;v0=float(err0[mask].sum());v1=float(err1[mask].sum())
            folds.append({"location":str(site),"n_parents":int(mask.sum()),"depth_sse":v0,
                "biology_sse":v1,"gain":gain(v0,v1),"gain_without_location":gain(a-v0,b-v1)})
        keep=np.array([Counter(locations)[site]>1 for site in locations])
        result.append({"endpoint":key,"n_features":y.shape[1],"depth_sse":a,"biology_sse":b,
            "held_location_gain":g,"p_approx":float((1+np.count_nonzero(null>=g-1e-12))/(nboot+1)),
            "null_gain_95":np.quantile(null,[.025,.975]).tolist(),"n_saturated_restricted_parents":sat,
            "gain_excluding_singletons_from_scoring":gain(err0[keep].sum(),err1[keep].sum()),"folds":folds})
    for row,p in zip(result,holm([r["p_approx"] for r in result])):row["p_holm2_approx"]=p
    return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True);p.add_argument("--bootstrap",type=int,default=999)
    a=p.parse_args()
    for name,sha in BLOBS.items():
        if gitblob((a.source/name).read_bytes())!=sha:raise ValueError("original Git blob mismatch "+name)
    ids,loc,want=get_meta(a.source/METADATA)
    x,audit=read_matrix(a.source/MATRIX,want)
    f,ys=features_and_outcomes(x)
    z={"status":"SOURCE_LOCKED_SPONTANEOUS_SUCCESSION_PREDICTION_NOT_INTERVENTION",
       "n_parents":275,"n_final_replicates":1100,"n_matched_states":1375,"n_locations":14,
       "source_asv_matrix":audit,"n_bootstrap":a.bootstrap,
       "source_blob_sha":BLOBS,"original_location_counts":dict(Counter(loc)),
       "baseline_read_quantiles":np.quantile(x[:,0].sum(1),[0,.25,.5,.75,1]).tolist(),
       "results":score(f,loc,ys,a.bootstrap),
       "limits":["No introduced-species intervention. Not replication of Saine treatment moderation.",
                 "Original authors already show initial-state dependence.",
                 "Fixed linear predictive model, processed read profiles not biomass.",
                 "Approximate restricted-model wild reference; parent errors assumed independent conditional on locations.",
                 "All four final replicates averaged within original parent; no independent replication inflation."]}
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("REPLAY_STATE_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
