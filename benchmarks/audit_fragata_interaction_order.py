#!/usr/bin/env python3
"""KNK source-locked mechanism anchor: 2022 spider mite order experiment.

This is a source-replication audit, NOT a new discovery of coexistence:
Fragata et al. (2022) Ecology Letters doi:10.1111/ele.14021
Source: https://github.com/irfragata/order_arrival_niche_preemption

Independent units are whole experimental BOXES. Source is read only in ephemeral
CI storage; outputs contain only derived group-level metrics and protocol.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

SOURCE_COMMIT = "f1568dfb4a0da7bc9811985e416daec86447d0f3"
RAW_BLOB = "dd734646364eb75464ab16b0169afb4b06b6e6d4"
AGE_MAP_BLOB = "354dd41ece632b8e5f8bd0266233f5a88f1517e5"
ORDER_LABELS = {
    "10Te10Tu": "simultaneous",
    "10Te10Tu48h": "Tu_first",
    "10Te48h10Tu": "Te_first",
}
N_PERMUTATIONS_EXACT = 252 * 252
ORDER = ("10Te10Tu", "10Te10Tu48h", "10Te48h10Tu")
SEED = 20261123


class SourceStop(ValueError):
    pass


def check_git_blob(path: Path, expected: str) -> tuple[bytes,str]:
    data=path.read_bytes()
    hash_=hashlib.sha1(f"blob {len(data)}\0".encode()+data).hexdigest()
    if hash_!=expected:
        raise SourceStop(f"Source Git blob mismatch: {path.name}: {hash_}")
    return data,hashlib.sha256(data).hexdigest()


def read_source(source:Path, age_map:Path):
    raw,sha=check_git_blob(source,RAW_BLOB)
    raw_age,age_sha=check_git_blob(age_map,AGE_MAP_BLOB)
    maps={}
    for r in csv.DictReader(raw_age.decode("utf-8-sig").splitlines(),delimiter="\t"):
        key=(r["Treatment"].strip(),int(r["Box"]),int(r["Leaf"]))
        if key in maps or r["Age"] not in ("old","new"):
            raise SourceStop("Invalid or repeated leaf age mapping")
        maps[key]=r["Age"]
    if len(maps)!=360:
        raise SourceStop(f"Expected 360 leaf age maps, got {len(maps)}")
    rows=list(csv.DictReader(raw.decode("utf-8-sig").splitlines()))
    if len(rows)!=3214:
        raise SourceStop(f"Source row-count differs from original release: {len(rows)}")
    expected_cols={"Box","Treatment","Block","Leaf","Leaflet","Direction","Te","Tu"}
    if not expected_cols.issubset(rows[0]):
        raise SourceStop("Source header mismatch")
    boxes={}
    for row in rows:
        treatment=row["Treatment"].strip()
        block=int(row["Block"])
        box=int(row["Box"])
        leaf=int(row["Leaf"])
        if block not in (1,2) or leaf not in (2,3,4,5):
            raise SourceStop("Invalid block/leaf ID")
        if box not in range(1,11):
            raise SourceStop("Unexpected replicate box")
        leaf_age=maps.get((treatment,box,leaf))
        if leaf_age is None:
            raise SourceStop("Source leaf_age join unmatched")
        value={}
        for sp in ("Te","Tu"):
            s=row[sp].strip()
            if s in ("NA",""):
                if not ((treatment=="20Te" and sp=="Tu") or
                        (treatment=="20Tu" and sp=="Te")):
                    raise SourceStop(f"Unexpected NA in mixed treatment {treatment}/{sp}")
                s="0"
            v=int(s)
            if v<0:
                raise SourceStop("Negative adult female counts")
            value[sp]=v
        key=(treatment,block,box)
        if key not in boxes:
            boxes[key]={"treatment":treatment,"block":block,"box":box,
                        "Te":0,"Tu":0,"leaves":{},
                        "new_Te":0,"new_Tu":0,
                        "n_source_subrecords":0}
        b=boxes[key]
        b["n_source_subrecords"]+=1
        for sp in ("Te","Tu"):
            b[sp]+=value[sp]
        ll=b["leaves"].setdefault(leaf,{"Te":0,"Tu":0})
        for sp in ("Te","Tu"):
            ll[sp]+=value[sp]
            if leaf_age=="new":
                b[f"new_{sp}"]+=value[sp]
    if len(boxes)!=90:
        raise SourceStop(f"Expected 90 independent boxes, got {len(boxes)}")
    treatments=Counter(k[0] for k in boxes)
    if len(treatments)!=9 or any(v!=10 for v in treatments.values()):
        raise SourceStop(f"Unbalanced experiment treatments: {treatments}")
    grouped=Counter((t,b) for t,b,rep in boxes)
    if any(grouped[(t,b)]!=5 for t in treatments for b in (1,2)):
        raise SourceStop("Expected each treatment to have 5 boxes per block")
    result=[]
    for b in boxes.values():
        if len(b["leaves"])!=4:
            raise SourceStop("Box did not have all 4 leaf strata")
        te,tu=b["Te"],b["Tu"]
        if b["treatment"] not in ("20Te","20Tu") and (te<=0 or tu<=0) \
            and b["treatment"] in ORDER:
            raise SourceStop("Zero count in balanced 10:10 treatment")
        denom=te+tu
        ll=b["leaves"]
        leaf_ov=sum(min((v["Te"]/te if te else 0),(v["Tu"]/tu if tu else 0))
                    for v in ll.values()) if te>0 and tu>0 else None
        result.append({
            "treatment":b["treatment"],"block":b["block"],"box":b["box"],
            "Te":te,"Tu":tu,
            "fraction_Te":te/denom if denom>0 else None,
            "fraction_Te_leaf4":ll[4]["Te"]/te if te>0 else None,
            "fraction_Te_leaf3":ll[3]["Te"]/te if te>0 else None,
            "fraction_Te_new_leaf":b["new_Te"]/te if te>0 else None,
            "fraction_Tu_new_leaf":b["new_Tu"]/tu if tu>0 else None,
            "spatial_leaf_overlap":leaf_ov,
        })
    return result,{"coexistence_sha256":sha,"leafmap_sha256":age_sha,
                   "raw_count_rows":len(rows),"n_independent_boxes":len(boxes),
                   "n_treatments":len(treatments),"boxes_per_treatment":dict(treatments),
                   "n_blocks":2,"status":"SOURCE_AUDIT_PASS"}


def exact_blocked_pairwise(data:list[dict], a:str,b:str,metric:str)->dict:
    one=[r for r in data if r["treatment"]==a]
    other=[r for r in data if r["treatment"]==b]
    if len(one)!=10 or len(other)!=10:
        raise SourceStop("Unbalanced head-to-head exact contrast")
    observed=np.mean([r[metric] for r in one])-np.mean([r[metric] for r in other])
    differences=[]
    for block in (1,2):
        pool=[r[metric] for r in data if r["block"]==block and
              r["treatment"] in (a,b)]
        if len(pool)!=10 or not all(math.isfinite(x) for x in pool):
            raise SourceStop("Nonfinite response or unexpected block sizes")
        vals=np.asarray(pool,dtype=float)
        total=float(np.sum(vals))
        differences.append(np.array([(2*float(np.sum(vals[list(ix)]))-total)/10
                                     for ix in itertools.combinations(range(10),5)]))
    dist=differences[0][:,None]+differences[1][None,:]
    p=float(np.mean(np.abs(dist)>=abs(observed)-1e-12))
    return {"comparison":f"{ORDER_LABELS[a]} minus {ORDER_LABELS[b]}",
            "raw_treatments":[a,b],
            "observed_mean_difference":float(observed),
            "exact_two_sided_p":p,
            "exact_assignments":int(dist.size),
            "n_boxes_each_group":10,
            "blocks":2}


def bootstrap_blocked(data:list[dict], a:str,b:str,metric:str,
                      n:int=20000,seed:int=SEED)->dict:
    rng=np.random.default_rng(seed)
    ma=np.zeros(n)
    mb=np.zeros(n)
    for treatment, target in ((a,ma),(b,mb)):
        for block in (1,2):
            x=np.array([r[metric] for r in data if r["treatment"]==treatment and
                        r["block"]==block],dtype=float)
            if len(x)!=5:
                raise SourceStop("Cannot bootstrap nonbalanced cells")
            target[:]+=rng.choice(x,size=(n,5),replace=True).sum(axis=1)/10
    delta=ma-mb
    return {"resample_unit":"whole box within treatment and block",
            "bootstrap_n":n,
            "bootstrap_seed":seed,
            "difference_percentile95":np.quantile(delta,[.025,.975]).tolist()}


def holm(p:list[float])->list[float]:
    p=np.asarray(p,dtype=float)
    order=np.argsort(p)
    sorted_q=np.maximum.accumulate([min(1,p[ii]*(len(p)-i))
                                    for i,ii in enumerate(order)])
    results=np.empty(len(p))
    results[order]=sorted_q
    return results.tolist()


def score(data:list[dict],audit:dict)->dict:
    primary=[r for r in data if r["treatment"] in ORDER]
    means={}
    keys=("fraction_Te","Te","Tu","fraction_Te_leaf4",
          "fraction_Te_leaf3","fraction_Te_new_leaf",
          "fraction_Tu_new_leaf","spatial_leaf_overlap")
    for treatment in ORDER:
        rows=[r for r in primary if r["treatment"]==treatment]
        means[ORDER_LABELS[treatment]]={
            "source_treatment":treatment,
            "n_independent_boxes":len(rows),
            "block_counts":{str(bl):sum(r["block"]==bl for r in rows)
                            for bl in (1,2)},
            "means":{k:float(np.mean([r[k] for r in rows])) for k in keys}
        }
    pairs=[(ORDER[1],ORDER[0]),(ORDER[2],ORDER[0]),(ORDER[2],ORDER[1])]
    metric_tests={}
    for key in ("fraction_Te","fraction_Te_leaf4","fraction_Te_leaf3",
                "Te","Tu"):
        tests=[exact_blocked_pairwise(primary,a,b,key) for a,b in pairs]
        adjusted=holm([x["exact_two_sided_p"] for x in tests])
        for i,q in enumerate(adjusted):
            tests[i]["holm_three_pairwise_p"]=q
        metric_tests[key]=tests
    primary_contrast=exact_blocked_pairwise(
        primary,ORDER[2],ORDER[1],"fraction_Te")
    primary_contrast["bootstrap"]=bootstrap_blocked(
        primary,ORDER[2],ORDER[1],"fraction_Te")
    return {
        "status":"REANALYSIS_OF_PUBLISHED_INTERACTION_EXPERIMENT_NOT_NEW_DISCOVERY",
        "original_source_paper":"Fragata et al. (2022), Ecology Letters, DOI 10.1111/ele.14021",
        "original_repository":"https://github.com/irfragata/order_arrival_niche_preemption",
        "pinned_source_commit":SOURCE_COMMIT,
        "source_audit":audit,
        "balanced_treatments":means,
        "primary_10_10_arrival_contrast":primary_contrast,
        "posthoc_exact_tests":metric_tests,
        "interpretation_limitations":[
            "90 experimental boxes total, but the primary contrast is 10+10 boxes; leaf-level counts are subsamples",
            "Original source already established history dependence of abundance and predicted coexistence",
            "Cross-sectional abundance does NOT directly identify competition coefficients",
            "Leaf 4 secondary targeted contrasts may fail even if published full-model spatial niche pre-emption is supported",
            "Source measured other initial-frequency levels and plant gene expression; this pilot does not refit published coexistence model",
            "Do not pool two-species mite effect with multi-taxon KNK rank-architecture effect or overwrite v28-v38",
        ]
    }


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--coexistence",type=Path,required=True)
    ap.add_argument("--leaf-map",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    data,audit=read_source(args.coexistence,args.leaf_map)
    result=score(data,audit)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print("INTERACTION_FORMATION_AUDIT_START")
    print("SOURCE_AUDIT",json.dumps(audit))
    for k,v in result["balanced_treatments"].items():
        print("TREATMENT_MEAN",k,json.dumps(v["means"]))
    print("PRIMARY_CONTRAST",json.dumps(result["primary_10_10_arrival_contrast"]))
    print("SECONDARY_LEAF4_TESTS",json.dumps(result["posthoc_exact_tests"]["fraction_Te_leaf4"]))
    print("INTERACTION_FORMATION_AUDIT_DONE")

if __name__=="__main__":
    main()
