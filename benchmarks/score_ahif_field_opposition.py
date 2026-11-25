"""Unclosed plant-hit role opposition with matched-date block contrasts."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.stats import t
from benchmarks.score_ahif_vansteijn2026_field import source_panels,SITES,YEARS
from benchmarks.probe_ahif_vansteijn2026_field import MD5

HISTORIES=("FtS","StF")
REFS=("SFlate","SFearly")
ROLES=("early_sown_group","other_sown","nonsown","total_hits")
SEED=20261125
NSPLIT=256

def metrics(d):
    d=np.asarray(d,float);g=float(np.abs(d).sum());n=float(abs(d.sum()));c=max(0.,g-n)
    return {"gross":g,"net":n,"cancellation":c,"fraction":c/g if g>1e-12 else None,"signed_total":float(d.sum())}

def gap(a,b):
    return float(np.sign(a)@b-np.sign(a.sum())*b.sum())

def rolemap(names,fast,slow,h):
    early=set(fast if h=="FtS" else slow);other=set(slow if h=="FtS" else fast)
    if early&other or len(early)!=10 or len(other)!=10 or not (early|other)<=set(names):
        raise ValueError("role roster")
    return np.array([[float(s in early),float(s in other),float(s not in early|other),1.] for s in names])

def unpack(panels):
    first=panels[SITES[0],YEARS[0]]
    names=sorted(set(first["names"])|set(first["sown_target"]))
    if len(names)!=78 or len(first["names"])!=76 or len(first["sown_target"])!=20:
        raise ValueError("roster changed")
    fast=first["fast"];slow=first["slow"];groups={};coverage=[]
    for site in SITES:
        ids=panels[site,YEARS[0]]["plots"]
        for year in YEARS:
            p=panels[site,year]
            if p["plots"]!=ids or p["fast"]!=fast or p["slow"]!=slow or p["names"]!=first["names"]:
                raise ValueError("persistent plot or roster drift")
            x=np.zeros((20,78));x[:,[names.index(s) for s in p["names"]]]=p["X"]
            if np.any(x<0) or not np.isfinite(x).all() or np.any(x.sum(1)<=0):raise ValueError("hit counts")
            d=x.sum(1)
            coverage.append({"site":site,"year":year,"min_total":float(d.min()),"max_total":float(d.max()),
              "mean_total":float(d.mean()),"n_distinct_totals":len(np.unique(d)),"constant_total":bool(np.ptp(d)<1e-12)})
            for h in (*HISTORIES,*REFS):
                v=[]
                for b in range(5):
                    ix=np.flatnonzero(p["history"][4*b:4*b+4]==h)
                    if len(ix)!=1:raise ValueError("block allocation")
                    v.append(x[4*b+ix[0]])
                groups[site,year,h]=np.asarray(v)
    return groups,names,fast,slow,coverage

def contrasts(groups,names,fast,slow):
    tax={};roles={}
    for site in SITES:
        for year in YEARS:
            for ref in REFS:
                rr=[]
                for h in HISTORIES:
                    d=groups[site,year,h]-groups[site,year,ref]
                    tax[site,year,ref,h]=d
                    rr.append(d@rolemap(names,fast,slow,h))
                roles[site,year,ref]=np.mean(rr,axis=0)
    return tax,roles

def welch(xs,nfamily):
    # fixed sites equal-weighted; five independent block-level differences each.
    x=np.asarray(xs,float)
    if x.shape!=(2,5,4):raise ValueError("block grain")
    est=x.mean((0,1));v=x.var(axis=1,ddof=1)/5/4
    vs=v.sum(0);den=(v*v/4).sum(0)
    rows=[]
    for j,k in enumerate(ROLES):
        if vs[j]<=1e-15 or den[j]<=0:
            rows.append({"role":k,"estimate":float(est[j]),"ci95":None,"status":"ZERO_VARIANCE"})
            continue
        df=vs[j]**2/den[j];se=np.sqrt(vs[j]);half=t.ppf(1-.05/(2*nfamily),df)*se
        rows.append({"role":k,"estimate":float(est[j]),"se":float(se),"df":float(df),
                     "ci95":[float(est[j]-half),float(est[j]+half)],"status":"APPROX_WELCH"})
    return rows

def crossfit(tax,roles,nsplit):
    rng=np.random.default_rng(SEED)
    values={k:[] for k in tax};rv={(y,r):[] for y in YEARS for r in REFS}
    for _ in range(nsplit):
        ix={s:rng.permutation(5) for s in SITES}
        for tr,te in ((slice(0,2),slice(2,5)),(slice(2,5),slice(0,2))):
            for k,d in tax.items():
                p=ix[k[0]];a=d[p[tr]].mean(0);b=d[p[te]].mean(0)
                values[k].append(gap(a,b))
            for y,r in rv:
                a=np.mean([roles[s,y,r][ix[s][tr],:3].mean(0) for s in SITES],axis=0)
                b=np.mean([roles[s,y,r][ix[s][te],:3].mean(0) for s in SITES],axis=0)
                rv[y,r].append(gap(a,b))
    return ({k:float(np.mean(v)) for k,v in values.items()},{k:float(np.mean(v)) for k,v in rv.items()})

def support(rows):
    d={x["role"]:x["ci95"] for x in rows}
    neg=lambda s:d[s] is not None and d[s][1]<0
    pos=lambda s:d[s] is not None and d[s][0]>0
    return {"other_sown_negative_nonsown_positive":bool(neg("other_sown") and pos("nonsown")),
            "strict_early_positive_other_or_nonsown_negative":bool(pos("early_sown_group") and (neg("other_sown") or neg("nonsown")))}

def run(root,nsplit=NSPLIT):
    for name,h in MD5.items():
        if hashlib.md5((root/name).read_bytes()).hexdigest()!=h:raise ValueError("source hash")
    panels,audit=source_panels(root)
    groups,names,fast,slow,coverage=unpack(panels)
    tax,roles=contrasts(groups,names,fast,slow)
    cf,rcf=crossfit(tax,roles,nsplit)
    rows=[]
    for (s,y,r,h),d in sorted(tax.items()):
        m=d.mean(0)
        rows.append({"site":s,"year":y,"reference":r,"history":h,"species_delta":m.tolist(),
                     "taxon_metrics":metrics(m),"crossfit_taxon_signed_gap":cf[s,y,r,h],
                     "role_delta":(m@rolemap(names,fast,slow,h)).tolist()})
    result=[]
    for y in YEARS:
        for r in REFS:
            x=np.array([roles[s,y,r] for s in SITES])
            primary=y=="2022" and r=="SFlate"
            intervals=welch(x,4 if primary else 1)
            v=x.mean((0,1));rr=[a for a in rows if a["year"]==y and a["reference"]==r]
            if len(rr)!=4 or not np.isclose(v[:3].sum(),v[3]):raise ValueError("role aggregation")
            gm={k:float(np.mean([a["taxon_metrics"][k] for a in rr])) for k in ("gross","net","cancellation")}
            gm["fraction"]=gm["cancellation"]/gm["gross"] if gm["gross"]>1e-12 else None
            gm["mean_crossfit_taxon_signed_gap"]=float(np.mean([a["crossfit_taxon_signed_gap"] for a in rr]))
            result.append({"year":y,"reference":r,"primary":primary,"interval_family_size":4 if primary else 1,
                           "role_estimates":intervals,"role_metrics":metrics(v[:3]),
                           "by_site_role_means":{s:x[i].mean(0).tolist() for i,s in enumerate(SITES)},
                           "pooled_role_crossfit_gap":rcf[y,r],"mean_within_contrast_taxon_metrics":gm,
                           "primary_directional_gate":support(intervals) if primary else None})
    return {"status":"SOURCE_LOCKED_FIELD_HIT_OPPOSITION_NOT_BIOMASS_OR_MEDIATION",
            "audit":audit,"source_md5":MD5,"n_physical_plots":40,"n_blocks":10,"n_splits":nsplit,
            "taxon_names":names,"n_sown":20,"n_nonsown":58,"role_names":list(ROLES),
            "total_hit_coverage":coverage,"comparisons":result,"taxon_contrasts":rows,
            "limits":["Same source experiment as earlier field scores; not new source count.",
                "Pin hits not biomass, density, or recruitment rate; nonsown includes possible preexisting taxa.",
                "Different plot-level totals do not remove sampling constraints or detection effects.",
                "Two fixed sites, five blocks/site; approximate independent-block simultaneous intervals.",
                "Shared synchronous reference captured within block; repeated years not independent.",
                "No recruitment intervention, no causal compensation or equivalence."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    z=run(a.source);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("FIELD_OPPOSITION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
