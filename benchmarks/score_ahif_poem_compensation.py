"""POEM field absolute biomass compensation, preserving shared-control and repeated-plot grain."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
from scipy.stats import t
from benchmarks.score_poem_field_storage import source_panels,matrix_for,FOUR,SOURCE_SOWN

COHORTS=("2020","2021")
AGES=(1,2,3)
EARLY=("F","G","L")
GUILDS={"G":("Agrostis capillaris","Anthoxanthum odoratum","Bromus hordeaceus","Festuca ovina"),
        "L":("Lathyrus pratensis","Lotus corniculatus","Trifolium arvense","Trifolium campestre"),
        "F":("Dianthus deltoides","Jasione montana","Pilosella officinarum","Pimpinella saxifraga","Potentilla argentea","Silene vulgaris")}
ROLES=("early_sown_group","other_sown","nonsown","total_identified_shoot")
SEED=20261125
NSPLIT=256

def metrics(d):
    d=np.asarray(d,float)
    g=float(np.abs(d).sum());n=float(abs(d.sum()));c=max(0.,g-n)
    return {"gross":g,"net":n,"cancellation":c,"fraction":c/g if g>1e-12 else None,
            "signed_total":float(d.sum())}

def sign_gap(a,b):
    return float(np.sign(a)@b-np.sign(a.sum())*b.sum())

def role_matrix(names,g):
    first=set(GUILDS[g]);sown=set(SOURCE_SOWN)
    return np.asarray([[float(x in first),float(x in sown-first),float(x not in sown),1.] for x in names])

def matrices(panels):
    groups={};names0=None;ids0={}
    if set().union(*map(set,GUILDS.values()))!=set(SOURCE_SOWN) or sum(map(len,GUILDS.values()))!=14:
        raise ValueError("guild roster mismatch")
    for c in COHORTS:
        for a in AGES:
            x,l,ids,names=matrix_for(panels[c,a],FOUR)
            if names0 is None:names0=names
            if names!=names0 or len(names)!=49 or not set(SOURCE_SOWN)<=set(names):
                raise ValueError("roster mismatch")
            if x.shape!=(20,49) or not np.isfinite(x).all() or (x<0).any():
                raise ValueError("mass matrix mismatch")
            for g in (*EARLY,"S"):
                ix=np.flatnonzero(l==g)
                k=(c,g);v=[ids[i] for i in ix]
                if len(v)!=5:raise ValueError("plot count")
                if k in ids0 and v!=ids0[k]:raise ValueError("plot order changed")
                ids0[k]=v
                groups[c,a,g]=x[ix]
    return groups,names0

def role_contributions(groups,names,age):
    # Independent treatment means. S common control is role-averaged per plot.
    out=[]
    for c in COHORTS:
        for g in EARLY:
            out.append((groups[c,age,g]@role_matrix(names,g),1/6))
        control=sum(groups[c,age,"S"]@role_matrix(names,g) for g in EARLY)/3
        out.append((control,-1/2))
    return out

def combine(contrib):
    return sum(w*x.mean(0) for x,w in contrib)

def welch(contrib,nfamily):
    obs=combine(contrib)
    vs=np.array([w*w*x.var(0,ddof=1)/len(x) for x,w in contrib])
    dfs=np.array([len(x)-1 for x,w in contrib])
    var=vs.sum(0)
    den=(vs*vs/dfs[:,None]).sum(0)
    rows=[]
    for j,k in enumerate(ROLES):
        if var[j]<=1e-15 or den[j]<=0:
            rows.append({"role":k,"estimate":float(obs[j]),"ci95":None,"status":"ZERO_VARIANCE"})
            continue
        df=var[j]**2/den[j];se=np.sqrt(var[j]);half=t.ppf(1-.05/(2*nfamily),df)*se
        rows.append({"role":k,"estimate":float(obs[j]),"se":float(se),"df":float(df),
                     "ci95":[float(obs[j]-half),float(obs[j]+half)],"status":"APPROX_WELCH"})
    return rows

def by_age(groups,names):
    rows=[]
    for a in AGES:
        for c in COHORTS:
            for g in EARLY:
                d=groups[c,a,g].mean(0)-groups[c,a,"S"].mean(0)
                rr=d@role_matrix(names,g)
                rows.append({"cohort":c,"age":a,"early":g,"species_delta":d.tolist(),
                             "taxon_metrics":metrics(d),"role_delta":rr.tolist(),
                             "role_metrics":metrics(rr[:3])})
    return rows

def crossfit(groups,names,nsplit=NSPLIT):
    rng=np.random.default_rng(SEED)
    keys=[(c,a,g) for c in COHORTS for a in AGES for g in EARLY]
    values={k:[] for k in keys};roles={a:[] for a in AGES}
    for _ in range(nsplit):
        # Shared plot allocation across years, and shared S for all comparisons.
        ix={(c,g):rng.permutation(5) for c in COHORTS for g in (*EARLY,"S")}
        half=[{(c,a,g):groups[c,a,g][p[sl]] for (c,g),p in ix.items() for a in AGES}
              for sl in (slice(0,2),slice(2,5))]
        for train,test in (half,half[::-1]):
            for c,a,g in keys:
                d=train[c,a,g].mean(0)-train[c,a,"S"].mean(0)
                v=test[c,a,g].mean(0)-test[c,a,"S"].mean(0)
                values[c,a,g].append(sign_gap(d,v))
            for a in AGES:
                d=combine(role_contributions(train,names,a))[:3]
                v=combine(role_contributions(test,names,a))[:3]
                roles[a].append(sign_gap(d,v))
    return ({k:float(np.mean(v)) for k,v in values.items()},
            {a:float(np.mean(v)) for a,v in roles.items()})

def run(path,nsplit=NSPLIT):
    panels,audit=source_panels(path)
    if audit["total_source_rows"]!=5210 or audit["n_unique_independent_plots_across_cohorts"]!=50:
        raise ValueError("source grain")
    groups,names=matrices(panels)
    rows=by_age(groups,names)
    cf,rcf=crossfit(groups,names,nsplit)
    for r in rows:r["crossfit_taxon_signed_gap"]=cf[r["cohort"],r["age"],r["early"]]
    ages=[]
    for a in AGES:
        cs=role_contributions(groups,names,a)
        role=combine(cs)
        rr=welch(cs,4 if a==2 else 1)
        r=[v for v in rows if v["age"]==a]
        gm={k:float(np.mean([v["taxon_metrics"][k] for v in r])) for k in ("gross","net","cancellation")}
        gm["fraction"]=gm["cancellation"]/gm["gross"] if gm["gross"]>1e-12 else None
        gm["mean_crossfit_taxon_signed_gap"]=float(np.mean([v["crossfit_taxon_signed_gap"] for v in r]))
        early=rr[0]["ci95"];other=[v["ci95"] for v in rr[1:3]]
        support=early is not None and early[0]>0 and any(x is not None and x[1]<0 for x in other)
        ages.append({"age":a,"primary":a==2,"interval_family_size":4 if a==2 else 1,
            "role_estimates":rr,"role_metrics":metrics(role[:3]),
            "pooled_role_crossfit_gap":rcf[a],"mean_within_contrast_taxon_metrics":gm,
            "opposing_role_signs_supported":bool(support) if a==2 else None})
    return {"status":"SOURCE_LOCKED_POEM_ABSOLUTE_BIOMASS_COMPENSATION_NOT_MEDIATION",
        "audit":audit,"n_primary_plots":40,"n_source_cohorts":2,"n_split_draws":nsplit,
        "taxon_names":names,"role_names":list(ROLES),"ages":ages,"contrasts":rows,
        "limits":["Same POEM source experiment, no new independent study.",
            "Primary age2 chosen from previous results; exploratory, not independent confirmation.",
            "Shared synchronous controls combined at plot level; approximate independent-group intervals.",
            "Identified shoot grams per equal-area sampled quadrat, not whole plot production.",
            "Arithmetic cancellation and cross-fitted signed gap do not identify causal interactions."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--zip",type=Path,required=True);p.add_argument("--out",type=Path,required=True);p.add_argument("--splits",type=int,default=NSPLIT);a=p.parse_args()
    z=run(a.zip,a.splits);a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("POEM_COMPENSATION_RESULT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
