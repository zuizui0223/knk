#!/usr/bin/env python3
"""Source-audited two-species identity-vs-rank and interaction data audit.

Original: Zou & Rudolf 2026, Ecology DOI 10.1002/ecy.70454.
Public author source: hengxingzou/Zou2026Ecology pinned
commit 5f28e5da78853bb31de1e1fa86f20f4892d7f2bb.
This is a POST HOC representation reanalysis, NOT source-model refitting or
independently validated mediation.

Pre-score protocol: docs/AHIF_ZOU2026_SOURCE_AND_POSTHOC_REPRESENTATION_AUDIT_20261124.md
"""
from __future__ import annotations
import argparse, csv, hashlib, itertools, json, math, random
from collections import Counter, defaultdict
from pathlib import Path

COUNT_BLOB="8fef4ec5818319452e1febb4c51042530b39a676"
COEFF_BLOB="5071d3d4bfdfb831326e5eb49d3960c0902c29fc"
LAMBDA_BLOB="3380ce79795a8ca2419524480cf6a70e0acdfbd3"
SOURCE_COMMIT="5f28e5da78853bb31de1e1fa86f20f4892d7f2bb"
DAYS=(-18,-12,-6,0,6,12,18)
GEN=(1,2,3)
DENSITY=(15,30)
REP=(1,2,3)
BOOTSTRAPS=20000
SEED=20261124

def pinned_csv(path:Path, git_sha:str):
    raw=path.read_bytes()
    sha=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    if sha!=git_sha:
        raise ValueError(f"Original author source hash mismatch: {path.name} got {sha}")
    with path.open(encoding="utf-8-sig",newline="") as f:
        rows=list(csv.DictReader(f))
    if not rows:raise ValueError("Empty original CSV")
    return rows,hashlib.sha256(raw).hexdigest()

def source_data(countpath:Path,coeffpath:Path,lambdapath:Path):
    raw,sha=pinned_csv(countpath,COUNT_BLOB)
    c,sc=pinned_csv(coeffpath,COEFF_BLOB)
    lm,sl=pinned_csv(lambdapath,LAMBDA_BLOB)
    if len(raw)!=414 or len(c)!=42 or len(lm)!=6:
        raise ValueError("Expected 414 whole vials, 42 posterior parameter summaries and 6 intrinsic growth summaries")
    expected_cols=("Vial_Name","Num_Gen","Arriv_Time","Start_T","Start_F",
                   "Replication","First_Addition","Second_Addition","Check_Date",
                   "Adults_T","Adults_F","Relative_Arriv_Time")
    if tuple(raw[0])!=expected_cols:raise ValueError("Original vial CSV schema changed")
    rows=[];keys=set()
    for x in raw:
        name=x["Vial_Name"]
        if name in keys:raise ValueError("Repeated independent vial "+name)
        keys.add(name)
        g=int(x["Num_Gen"]);day=int(x["Relative_Arriv_Time"])
        t=int(x["Start_T"]);f=int(x["Start_F"])
        r=int(x["Replication"]);adultT=int(x["Adults_T"]);adultF=int(x["Adults_F"])
        if g not in GEN or day not in DAYS or r not in REP or min(t,f,adultT,adultF)<0:
            raise ValueError("Invalid original vial factors or counts")
        if t>0 and f>0 and adultT+adultF==0:
            raise ValueError("Undefined mixed-vial species proportions; STOP")
        rows.append({"vial":name,"gen":g,"day":day,"startT":t,"startF":f,
                     "rep":r,"adultT":adultT,"adultF":adultF})
    configs={ (r["startT"],r["startF"]) for r in rows}
    expected_configs={(0,30),(0,60),(10,20),(15,15),(20,10),(20,40),(30,0),(30,30),(40,20),(60,0)}
    if configs!=expected_configs:raise ValueError("Unexpected original initial-density configurations")
    mixed=[r for r in rows if r["startT"]>0 and r["startF"]>0]
    controls=[r for r in rows if r["startT"]==0 or r["startF"]==0]
    if len(mixed)!=378 or len(controls)!=36:
        raise ValueError("Wrong mixed and solo vial count")
    exact=Counter((r["gen"],r["day"],r["startT"],r["startF"]) for r in rows)
    if any(v!=3 for v in exact.values()) or len(exact)!=138:
        raise ValueError("Expected each generation×arrival×initial-density source cell 3 vials")
    for g in GEN:
        for day in DAYS:
            expected_n=30 if day==0 else 18
            if sum(r["gen"]==g and r["day"]==day for r in rows)!=expected_n:
                raise ValueError("Missing generation×arrival vial block")
    params={}
    for x in c:
        key=(int(x["Num_Gen"]),int(x["Relative_Arriv_Time"]),x["Coeff"])
        if key in params or key[0] not in GEN or key[1] not in DAYS or key[2] not in ("alpha12","alpha21"):
            raise ValueError("Missing/extra coefficient summary factor")
        lo,hi,med=(float(x[k]) for k in ("Q10","Q90","Median"))
        if not all(math.isfinite(v) for v in (lo,hi,med)) or not 0<lo<=med<=hi:
            raise ValueError("Invalid author posterior parameter quantile")
        params[key]={"q10":lo,"q90":hi,"median":med}
    if len(params)!=42:raise ValueError("Incomplete author fitted alpha table")
    if sorted((x["Species"],int(x["Num_Gen"])) for x in lm)!=[(s,g) for s in ("F","T") for g in GEN]:
        raise ValueError("Incomplete source intrinsic growth table")
    return rows,params,{
        "study":"Zou and Rudolf 2026, Ecology, 10.1002/ecy.70454",
        "source_commit":SOURCE_COMMIT,
        "count_source_git_blob":COUNT_BLOB,"coeff_source_git_blob":COEFF_BLOB,
        "lambda_source_git_blob":LAMBDA_BLOB,
        "count_source_sha256":sha,"coeff_source_sha256":sc,
        "lambda_source_sha256":sl,
        "n_original_vials":len(rows),"n_mixed_vials":len(mixed),
        "n_single_species_controls":len(controls),
        "n_equal_initial_density_mixed_vials":sum(r["startT"]==r["startF"] and r["startT"]>0 for r in rows),
        "n_posterior_parameter_summary_rows":len(params),
        "replication_unit":"independent experimental vial, NOT adults and not model posterior samples"
    }

def adult_share(r):
    return r["adultT"]/(r["adultT"]+r["adultF"])

def two_outcomes(row):
    p=adult_share(row)
    return p,max(p,1-p)

def group_contrasts(rows:list[dict],boot=BOOTSTRAPS,seed=SEED):
    rng=random.Random(seed)
    out=[];means=[]
    for g in GEN:
        for day in DAYS:
            sample=[r for r in rows if r["gen"]==g and r["day"]==day
                    and r["startT"]==r["startF"] and r["startT"] in DENSITY]
            if len(sample)!=6 or Counter(r["startT"] for r in sample)!={15:3,30:3}:
                raise ValueError("One of 21 balanced-initial-frequency conditions is unbalanced")
            y=[two_outcomes(r) for r in sample]
            means.append({"n_gen":g,"relative_arrival_days":day,"n_vials":6,
                          "mean_T_share":sum(v[0] for v in y)/6,
                          "mean_top1_share_unlabeled":sum(v[1] for v in y)/6,
                          "n_T_final_mass_winner":sum(v[0]>.5 for v in y),
                          "n_F_final_mass_winner":sum(v[0]<.5 for v in y),
                          "n_ties":sum(v[0]==.5 for v in y)})
        for absday in (6,12,18):
            a=[r for r in rows if r["gen"]==g and r["day"]==-absday and
                r["startT"]==r["startF"] and r["startT"] in DENSITY]
            b=[r for r in rows if r["gen"]==g and r["day"]==absday and
                r["startT"]==r["startF"] and r["startT"] in DENSITY]
            if len(a)!=6 or len(b)!=6:
                raise ValueError("One first-arrival direction missing")
            observed=[]
            for metric in (0,1):
                va=[two_outcomes(r)[metric] for r in a]
                vb=[two_outcomes(r)[metric] for r in b]
                observed.append(sum(va)/6-sum(vb)/6)
            blocks={}
            for level in DENSITY:
                aa=[r for r in a if r["startT"]==level]
                bb=[r for r in b if r["startT"]==level]
                if len(aa)!=3 or len(bb)!=3:raise ValueError("Broken density block")
                blocks[level]=[two_outcomes(r) for r in aa+bb]
            indices=list(itertools.combinations(range(6),3))
            if len(indices)!=20:raise ValueError("Expected 20 block allocations")
            distributions=[[],[]]
            for ia in indices:
                for ib in indices:
                    diffs=[0.,0.]
                    for level,idx in zip(DENSITY,(ia,ib)):
                        block=blocks[level]
                        selected=set(idx)
                        for metric in (0,1):
                            mean_p=sum(v[metric] for i,v in enumerate(block) if i in selected)
                            mean_q=sum(v[metric] for i,v in enumerate(block) if i not in selected)
                            diffs[metric]+=(mean_p-mean_q)/6
                    for metric in (0,1):distributions[metric].append(diffs[metric])
            if any(len(z)!=400 for z in distributions):raise ValueError("Exact allocation count !=400")
            p=[sum(abs(t)>=abs(obs)-1e-12 for t in d)/400
               for obs,d in zip(observed,distributions)]
            boot_diffs=[[],[]]
            for _ in range(boot):
                vals=[[],[]]
                for level in DENSITY:
                    for direction, members in (("T",a),("F",b)):
                        v=[two_outcomes(r) for r in members if r["startT"]==level]
                        if len(v)!=3:raise ValueError("Bootstrap source wrong")
                        draws=[v[rng.randrange(3)] for _ in range(3)]
                        weight=+1 if direction=="T" else -1
                        for metric in (0,1):
                            vals[metric].append(weight*sum(x[metric] for x in draws)/6)
                for metric in (0,1):boot_diffs[metric].append(sum(vals[metric]))
            ci=[]
            for d in boot_diffs:
                d.sort()
                ci.append([d[max(0,int(.025*boot)-1)],d[min(boot-1,int(.975*boot))]])
            out.append({"n_generations":g,"lead_days":absday,
                "n_T_first_vials":len(a),"n_F_first_vials":len(b),
                "n_vials_per_initial_density_and_direction":3,
                "delta_T_identity_share_Tfirst_minus_Ffirst":observed[0],
                "delta_unlabeled_top1_share_Tfirst_minus_Ffirst":observed[1],
                "delta_T_share_bootstrap95":ci[0],
                "delta_unlabeled_top1_bootstrap95":ci[1],
                "exact_400_within_density_permutation_p_two_sided_T_share":p[0],
                "exact_400_within_density_permutation_p_two_sided_top1":p[1],
                "n_exact_permutations":400,
                "note":"Within-source posthoc difference, not multi-species ranked abundance or causal mediation"})
    return means,out

def pearson(a,b):
    if len(a)!=len(b) or len(a)<2:raise ValueError("Correlation needs matching arrays")
    ma=sum(a)/len(a);mb=sum(b)/len(b)
    den1=sum((x-ma)**2 for x in a)
    den2=sum((y-mb)**2 for y in b)
    if den1<=0 or den2<=0:raise ValueError("Degenerate coefficient")
    return sum((x-ma)*(y-mb) for x,y in zip(a,b))/math.sqrt(den1*den2)

def coefficient_context(means:list[dict],params:dict):
    result=[]
    for row in means:
        gen=row["n_gen"];day=row["relative_arrival_days"]
        a21=params[gen,day,"alpha21"]["median"]
        a12=params[gen,day,"alpha12"]["median"]
        result.append({**row,
            "source_fitted_alpha12_median_F_impact_on_T":a12,
            "source_fitted_alpha21_median_T_impact_on_F":a21,
            "log_alpha21_to_alpha12_asymmetry":math.log(a21/a12),
            "coeff_source_is_same_count_dataset":True})
    return {"n_cells":len(result),
        "source_fitted_coefficient_asymmetry_vs_mean_T_share_pearson_in_sample":pearson(
            [r["log_alpha21_to_alpha12_asymmetry"] for r in result],
            [r["mean_T_share"] for r in result]),
        "correlation_by_generation_IN_SAMPLE_no_p":{
            str(g):pearson(
                [r["log_alpha21_to_alpha12_asymmetry"] for r in result if r["n_gen"]==g],
                [r["mean_T_share"] for r in result if r["n_gen"]==g])
            for g in GEN},
        "condition_level":result,
        "no_independent_prediction":"The fitted parameters are derived from the SAME adult population counts; correlation has no independent validation or mediation interpretation."}

def analyze(countpath:Path,coeffpath:Path,lambdapath:Path,bootstrap=BOOTSTRAPS):
    rows,params,src=source_data(countpath,coeffpath,lambdapath)
    means,contrasts=group_contrasts(rows,boot=bootstrap)
    if len(means)!=21 or len(contrasts)!=9:raise ValueError("Unexpected score count")
    coefficient=coefficient_context(means,params)
    return {
        "status":"ORIGINAL_SOURCE_AUTHENTIC_POSTHOC_TWO_SPECIES_STRUCTURAL_CONTRAST",
        "source_audit":src,
        "protocol":"docs/AHIF_ZOU2026_SOURCE_AND_POSTHOC_REPRESENTATION_AUDIT_20261124.md",
        "analysis":{"n_balanced_initial_vials":126,
            "relative_starting_frequency":"T:F = 1:1 at 15:15 or 30:30",
            "n_gen":3,"arrival_time_levels":list(DAYS),
            "symmetric_time_contrast_days":[6,12,18],
            "n_symmetric_contrasts":9,
            "bootstrap_whole_vials_within_generation_density_and_arrival":bootstrap,
            "exact_conditionally_exchangeable_label_allocations":400,
            "seed":SEED},
        "condition_means":means,
        "symmetric_arrival_reversal_contrasts":contrasts,
        "fitted_interaction_parameter_context":coefficient,
        "limits":[
            "Only two competing species; sorted ranks reduce to max adult share and cannot test >=11 taxa law.",
            "Coefficients are POSTERIOR FITTED from same counts and cannot independently predict the source group means.",
            "Original authors already showed history-dependent competitive coefficients, non-novel to AHIF.",
            "Nine contrasts are post hoc exploratory, not independent source study effects.",
            "Exact conditional randomization p values assume exchangeability within original initial-density vials; no claim of randomized cohort block without original methods audit.",
            "The different relative abundance and max-share response transformations are not information-theoretic bits or direct mediation.",
            "Original source solo vials excluded from abundance contrast but available for coefficient model source.",
            "No changes to frozen KNK v28-v38 or Ecology Letters submission."]
    }

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts",type=Path,required=True)
    p.add_argument("--coefficients",type=Path,required=True)
    p.add_argument("--lambdas",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--bootstrap",type=int,default=BOOTSTRAPS)
    a=p.parse_args()
    res=analyze(a.counts,a.coefficients,a.lambdas,bootstrap=a.bootstrap)
    a.out.parent.mkdir(exist_ok=True,parents=True)
    a.out.write_text(json.dumps(res,ensure_ascii=False,indent=2,allow_nan=False)+"\n")
    print("ZOU2026_SCORE_START")
    print("SOURCE",json.dumps(res["source_audit"]))
    for r in res["symmetric_arrival_reversal_contrasts"]:
        print("SYM_CONTRAST",json.dumps(r))
    corr=res["fitted_interaction_parameter_context"]
    print("SOURCE_COEFF_CORR",json.dumps({
        k:v for k,v in corr.items() if k!="condition_level"}))
    print("ZOU2026_SCORE_END",a.out)

if __name__=="__main__":main()
