"""Fixed-site routing contrasts and deterministic richness/concentration rank residuals."""
from __future__ import annotations
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t
from benchmarks.score_ahif_vansteijn2026_field import (
    source_panels, blocked_assignment_stream, score_matrix, SourceStop, SITES, SEED
)
from benchmarks.score_ahif_nontarget_temporal import roster_data, changes, holm_values, pair_differences
from benchmarks.probe_ahif_vansteijn2026_field import MD5

TREAT=("FtS","StF","SFearly","SFlate")
# Fast early vs late, slow early vs late, same vs staggered.
COEF=np.asarray([[1,-1,1,-1],[-1,1,1,-1],[-1,-1,1,1]],float)/2

def fitted_curve(p,kind):
    p=np.asarray(p,float)
    if p.ndim!=1 or np.any(p<0) or not np.isclose(p.sum(),1.) or np.any(p[:-1]<p[1:]-1e-12):
        raise SourceStop("Invalid sorted positive-normalized vector")
    k=int(np.count_nonzero(p>0)); c=float(p@p)
    out=np.zeros_like(p)
    if k==1: out[0]=1.; return out
    if not 1/k-1e-12<=c<=1+1e-12: raise SourceStop("Invalid concentration")
    if kind=="one_dominant":
        a=(1+np.sqrt(max(0.,(k-1)*(k*c-1))))/k
        out[:k]=(1-a)/(k-1);out[0]=a
    elif kind=="geometric":
        j=np.arange(k,dtype=float)
        lo=0.;hi=1.
        for _ in range(70):
            q=(lo+hi)/2
            z=q**j;z/=z.sum()
            if z@z>c: lo=q
            else: hi=q
        z=((lo+hi)/2)**j;out[:k]=z/z.sum()
    else: raise ValueError("unknown curve")
    if not np.isclose(out@out,c,atol=1e-10) or not np.isclose(out.sum(),1.):
        raise SourceStop("Curve concentration mismatch")
    return out

def shape_residual(p,kind):
    return p-np.asarray([fitted_curve(x,kind) for x in p])

def block_contrasts(x,labels):
    if len(labels)!=20 or len(x)!=20: raise SourceStop("Expected 20 complete plots per site")
    x=np.asarray(x,float)
    out=[]
    for i in range(0,20,4):
        ls=list(labels[i:i+4])
        if set(ls)!=set(TREAT): raise SourceStop("Incomplete source block")
        z=np.asarray([x[i+ls.index(t)] for t in TREAT])
        out.append((COEF@z).ravel())
    return np.asarray(out)

def contrast_stat(a,b):
    d=a.mean(0)-b.mean(0)
    var=float(a.var(0,ddof=1).sum()/len(a)+b.var(0,ddof=1).sum()/len(b))
    q=float(d@d)
    return (q/var if var>1e-15 else None),q,var

def fixed_site_wild(a,b):
    """Approximate block-error wild bootstrap; not randomization inference."""
    if a.shape!=b.shape or a.shape[0]!=5 or a.ndim!=2:
        raise SourceStop("Five source blocks per site with matching response width")
    obs,q,var=contrast_stat(a,b)
    out={"difference_l2":float(np.sqrt(q)),"difference_squared_noise_corrected":q-var,
         "variance_trace":var,"n_blocks_per_site":5,"status":"APPROXIMATE_WILD_BOOTSTRAP_NOT_EXACT"}
    if obs is None:
        return {**out,"status":"STOP_ZERO_VARIANCE","statistic":None,"p_approx":None}
    ar=(a-a.mean(0))*np.sqrt(5/4); br=(b-b.mean(0))*np.sqrt(5/4)
    null=[]
    for signs in itertools.product((-1.,1.),repeat=10):
        s=np.asarray(signs)
        z,_,_=contrast_stat(ar*s[:5,None],br*s[5:,None])
        if z is None: return {**out,"status":"STOP_ZERO_NULL_VARIANCE","statistic":obs,"p_approx":None}
        null.append(z)
    return {**out,"statistic":obs,"n_sign_patterns":1024,
            "p_approx":float((1+np.count_nonzero(np.asarray(null)>=obs-1e-12))/1025)}

def welch_scalar(a,b):
    a=np.asarray(a,float).ravel(); b=np.asarray(b,float).ravel()
    if a.shape!=(5,) or b.shape!=(5,):raise SourceStop("Expected scalar block contrast")
    d=float(a.mean()-b.mean()); va=float(a.var(ddof=1)/5);vb=float(b.var(ddof=1)/5)
    if va+vb<=1e-15:return {"status":"STOP_ZERO_VARIANCE","difference_MAL_minus_UME":d,"p_approx":None}
    se=np.sqrt(va+vb); df=(va+vb)**2/(va**2/4+vb**2/4)
    crit=student_t.ppf(.975,df)
    return {"status":"WELCH_APPROXIMATE_FIXED_SITE_CONTRAST","difference_MAL_minus_UME":d,
            "MAL_mean":float(a.mean()),"UME_mean":float(b.mean()),"se":float(se),
            "df":float(df),"ci95_unadjusted":[float(d-crit*se),float(d+crit*se)],
            "p_approx":float(2*student_t.sf(abs(d/se),df))}

def adjust(rows,key,pout):
    # Keep original family size even if a stop occurs; stopped tests count as p=1
    ps=[1. if x[key] is None else x[key] for x in rows]
    for row,p in zip(rows,holm_values(ps)):
        row[pout]=None if row[key] is None else p

def analyse(panel,n_perm=9999):
    rank=[]; fit=[]; geo=[]; alt=[]; delta={}
    for i,s in enumerate(SITES):
        a=panel[s,"2021"];b=panel[s,"2024"]
        if a["plots"]!=b["plots"] or a["names"]!=b["names"] or not np.array_equal(a["history"],b["history"]):
            raise SourceStop("Repeated plot mismatch")
        x=roster_data(a);y=roster_data(b);labels=a["history"]
        d,l=changes(panel,s)
        delta[s]=(d,l)
        perm=blocked_assignment_stream(labels,n_perm=n_perm,seed=SEED+501+i)
        rank.append({"site":s,**score_matrix(d["nontarget_rank"],labels,perm,metric="euclidean")})
        for kind,dest in (("geometric",geo),("one_dominant",alt)):
            ra=shape_residual(x["nontarget_rank"],kind)
            rb=shape_residual(y["nontarget_rank"],kind)
            dest.append({"site":s,"curve":kind,**score_matrix(rb-ra,labels,perm,metric="euclidean")})
            for year,z,res in (("2021",x,ra),("2024",y,rb)):
                fit.append({"site":s,"year":year,"curve":kind,
                    "mean_squared_residual_l2":float(np.mean(np.sum(res*res,axis=1))),
                    "richness_range":[float(z["nontarget_richness"].min()),float(z["nontarget_richness"].max())],
                    "simpson_range":[float(z["nontarget_simpson_concentration"].min()),float(z["nontarget_simpson_concentration"].max())]})
    for rows in (geo,alt):adjust(rows,"p_upper","p_holm_family")
    wild=[]; scalar=[]
    for endpoint in ("nontarget_rank","nontarget_richness","nontarget_simpson_concentration"):
        a=block_contrasts(delta["MAL"][0][endpoint],delta["MAL"][1])
        b=block_contrasts(delta["UME"][0][endpoint],delta["UME"][1])
        z={"endpoint":endpoint,**fixed_site_wild(a,b)}
        if endpoint!="nontarget_rank":
            z["scalar_contrast_means"]={"MAL":a.mean(0).tolist(),"UME":b.mean(0).tolist(),
                       "contrast_order":["fast_early","slow_early","same_vs_staggered"]}
            aa=pair_differences(delta["MAL"][0][endpoint],delta["MAL"][1])
            bb=pair_differences(delta["UME"][0][endpoint],delta["UME"][1])
            scalar.append({"endpoint":endpoint,"contrast":"FtS_minus_StF",**welch_scalar(aa,bb)})
        wild.append(z)
    adjust(wild,"p_approx","p_holm_approx_family");adjust(scalar,"p_approx","p_holm_approx_family")
    return {"status":"EXPLORATORY_FIXED_SITE_AND_SHAPE_DIAGNOSTICS_NOT_MEDIATION",
            "n_physical_plots":40,"n_publications":1,"raw_rank_change_reference":rank,
            "primary_geometric_residual":geo,"secondary_one_dominant_residual":alt,
            "curve_fit_descriptives":fit,"approximate_fixed_site_contrasts":wild,
            "approximate_welch_order_sensitivity":scalar,
            "limits":["Same source plots; prior outcome motivated endpoints.",
                      "Shape residual is derived post-treatment endpoint, not controlled direct effect.",
                      "Wild bootstrap and Welch with five blocks/site are approximate, not exact randomization or site causation."]}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--original-dir",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--permutations",type=int,default=9999)
    a=p.parse_args()
    for name,md5 in MD5.items():
        if hashlib.md5((a.original_dir/name).read_bytes()).hexdigest()!=md5:
            raise SourceStop("Original hash mismatch: "+name)
    panel,_=source_panels(a.original_dir)
    z=analyse(panel,a.permutations)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("ROUTING_SHAPE_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
