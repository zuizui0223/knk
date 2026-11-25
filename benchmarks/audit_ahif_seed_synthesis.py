"""Official Figshare v1 source access and schema/lineage gate; no ecological scores."""
from __future__ import annotations
import argparse,csv,hashlib,json,math,re,urllib.request,urllib.error
from collections import Counter,defaultdict
from pathlib import Path
from urllib.parse import urlparse

ARTICLE=12319682
URL=f"https://api.figshare.com/v2/articles/{ARTICLE}/versions/1"
WANTED={"seedadd_plot_level.csv":"plot","seedadd_sp_level.csv":"species"}
MAX=50_000_000
NA={"","NA","NaN","nan","NULL","null"}

def request(url,limit):
    if urlparse(url).scheme!="https":raise ValueError("non-HTTPS public source")
    req=urllib.request.Request(url,headers={"User-Agent":"AHIF academic reproducibility","Accept":"*/*"})
    with urllib.request.urlopen(req,timeout=90) as res:
        data=res.read(limit+1)
    if len(data)>limit:raise ValueError("public file exceeds cap")
    return data

def acquire(root):
    meta=json.loads(request(URL,2_000_000))
    if int(meta.get("id",-1))!=ARTICLE or int(meta.get("version",-1))!=1:
        raise ValueError("not pinned public Figshare version")
    doi=str(meta.get("doi","")).lower()
    if doi not in ("10.6084/m9.figshare.12319682","10.6084/m9.figshare.12319682.v1"):
        raise ValueError("Figshare article DOI mismatch")
    out={}
    for item in meta.get("files",[]):
        name=item["name"]
        lane=WANTED.get(name.lower())
        if lane is None:continue
        if lane in out:raise ValueError("duplicate expected file")
        declared=int(item["size"]);md5=item.get("computed_md5") or item.get("supplied_md5")
        if not md5 or not re.fullmatch("[a-f0-9]{32}",md5) or not 0<declared<=MAX:
            raise ValueError("missing/inconsistent declared size/checksum")
        u=item["downloadUrl"]
        if urlparse(u).hostname not in ("figshare.com","api.figshare.com","ndownloader.figshare.com","s3-eu-west-1.amazonaws.com","s3.amazonaws.com"):
            raise ValueError("unexpected official manifest file host")
        if item.get("is_link_only"):raise ValueError("external link not an archived file")
        b=request(u,MAX)
        if len(b)!=declared or hashlib.md5(b).hexdigest()!=md5:
            raise ValueError("original Figshare byte mismatch")
        p=root/(lane+".csv");p.write_bytes(b)
        out[lane]={"source_name":name,"file_id":item["id"],"size":len(b),"md5":md5,"sha256":hashlib.sha256(b).hexdigest()}
    if set(out)!={"plot","species"}:raise ValueError("required original files absent")
    return {"article_id":ARTICLE,"version":1,"doi":doi,"title":meta.get("title"),"license":meta.get("license"),
            "source_description":re.sub("<[^>]+>"," ",str(meta.get("description","")))[:12000],
            "files":out}

def load(p):
    with p.open(encoding="utf-8-sig",newline="") as f:
        rd=csv.DictReader(f)
        fields=rd.fieldnames or []
        if len(fields)!=len(set(fields)):raise ValueError("duplicate header")
        rows=list(rd)
    if not rows or any(None in r or any(v is None for v in r.values()) for r in rows):
        raise ValueError("nonrectangular CSV")
    return fields,rows

def number(x):
    if str(x).strip() in NA:return None
    try:v=float(x)
    except (ValueError,TypeError):return None
    return v if math.isfinite(v) else None

def field_profile(fields,rows):
    ans={}
    for f in fields:
        c=Counter(r[f] for r in rows);v={"distinct":len(c),"missing":sum(n for k,n in c.items() if k.strip() in NA)}
        if len(c)<=24 and not re.search(r"unique|plot.?id|species",f,re.I):
            v["values"]={str(k):n for k,n in sorted(c.items())}
        ans[f]=v
    return ans

def audit(root,manifest):
    pf,p=load(root/"plot.csv");sf,s=load(root/"species.csv")
    required={"unique.id_","Experiment_","site","block","yr.trt","seed.rich","biomass.plot","rich.plot"}
    if not required<=set(pf) or not (required|{"species","biomass.sp"})<=set(sf):
        return {"status":"STOP_SCHEMA_RECONCILIATION","manifest":manifest,"plot_fields":pf,"species_fields":sf,
                "n_plot_rows":len(p),"n_species_rows":len(s),
                "missing_plot_fields":sorted(required-set(pf)),"missing_species_fields":sorted((required|{"species","biomass.sp"})-set(sf))}
    groups=defaultdict(list)
    for r in s:groups[r["unique.id_"]].append(r)
    idx=defaultdict(list)
    for r in p:idx[r["unique.id_"]].append(r)
    common=set(idx)&set(groups)
    md=[];a=[];dups=0
    for key in sorted(common):
        pp=idx[key]
        if len(pp)!=1:continue
        q=pp[0];rr=groups[key]
        bad=[f for f in required-{"unique.id_"} if any(r[f]!=q[f] and not (number(r[f]) is not None and number(q[f]) is not None and math.isclose(number(r[f]),number(q[f]),rel_tol=1e-9,abs_tol=1e-9)) for r in rr)]
        if bad:md.append({"experiment":q["Experiment_"],"fields":bad})
        ns=len({r["species"] for r in rr});dups+=len(rr)-ns
        vals=[number(r["biomass.sp"]) for r in rr];total=number(q["biomass.plot"]);rich=number(q["rich.plot"])
        mass_ok=all(v is not None and v>=0 for v in vals)
        pos=sum(v is not None and v>0 for v in vals)
        ratio=sum(vals)/total if mass_ok and total is not None and total>0 else None
        a.append({"experiment":q["Experiment_"],"metadata_match":not bad,"species_count":ns,"positive_weights":pos,
                  "reported_richness":rich,"valid_mass":mass_ok,"reported_total":total,"weight_sum_ratio":ratio})
    ex=[]
    for exp in sorted({x["Experiment_"] for x in p}):
        pp=[r for r in p if r["Experiment_"]==exp];ss=[r for r in s if r["Experiment_"]==exp];aa=[r for r in a if r["experiment"]==exp]
        strata=defaultdict(list)
        for row in pp:strata[row["site"],row["block"],row["yr.trt"]].append(row)
        both=0
        for v in strata.values():
            x=[number(r["seed.rich"]) for r in v]
            if any(z==0 for z in x) and any(z is not None and z>0 for z in x):both+=1
        ratio=sorted(x["weight_sum_ratio"] for x in aa if x["weight_sum_ratio"] is not None)
        ex.append({"experiment":exp,"plot_rows":len(pp),"unique_plot_state_ids":len({x["unique.id_"] for x in pp}),
            "species_rows":len(ss),"distinct_species":len({x["species"] for x in ss}),
            "seed_richness_levels":dict(sorted(Counter(x["seed.rich"] for x in pp).items())),
            "treatment_labels":dict(Counter(x.get("trt","MISSING_COLUMN") for x in pp)),
            "site_count":len({x["site"] for x in pp}),"site_block_yeartrt_strata":len(strata),"strata_with_seed0_and_positive":both,
            "n_joined_unique_states":len(aa),"metadata_mismatch_states":sum(not x["metadata_match"] for x in aa),
            "valid_species_mass_states":sum(x["valid_mass"] for x in aa),
            "species_count_matches_reported_richness":sum(x["reported_richness"]==x["species_count"] for x in aa),
            "positive_weights_match_reported_richness":sum(x["reported_richness"]==x["positive_weights"] for x in aa),
            "species_weight_sum_ratio_to_plot_mass_quantiles":([ratio[0],ratio[len(ratio)//2],ratio[-1]] if ratio else []),
            "states_mass_sum_within_one_percent":sum(abs(v-1)<=.01 for v in ratio)})
    return {"status":"SOURCE_BYTES_VERIFIED_SCHEMA_AUDITED_NOT_ECOLOGICAL_SCORING","manifest":manifest,
            "plot_rows":len(p),"species_rows":len(s),"plot_unique_state_ids":len(idx),
            "duplicate_plot_state_rows":sum(len(v)-1 for v in idx.values()),
            "species_unique_state_ids":len(groups),"duplicate_state_species_rows":dups,
            "plot_only_state_ids":len(set(idx)-set(groups)),"species_only_state_ids":len(set(groups)-set(idx)),
            "metadata_mismatch_fields":dict(Counter(f for r in md for f in r["fields"])),
            "n_experiments":len(ex),"experiments":ex,"plot_field_profile":field_profile(pf,p),
            "species_field_profile":field_profile(sf,s),
            "potential_original_physical_id_fields":[f for f in pf if re.search("plot|rep|subplot|year",f,re.I)],
            "potential_seed_membership_fields":[f for f in sf if re.search("sown|seed|colonist|resident|invad",f,re.I)],
            "limits":["Secondary synthesis, not original raw-source audit.",
               "Species masses may derive from cover and measured total.",
               "State IDs may not identify independent plots across years.",
               "No treatment inference or recruitment identity assignment."]}

def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    a.root.mkdir(parents=True,exist_ok=True)
    try:
        m=acquire(a.root);z=audit(a.root,m)
    except urllib.error.HTTPError as e:
        z={"status":"STOP_PUBLIC_SOURCE_HTTP_ERROR","http_status":e.code,"source":e.url,
           "limits":["No access workaround or ecological scores."]}
    except (ValueError,KeyError,OSError) as e:
        z={"status":"STOP_SOURCE_GATE","reason":str(e)[:600],"limits":["No ecological scores."]}
    a.out.write_text(json.dumps(z,indent=2,allow_nan=False)+"\n")
    print("SEED_SYNTHESIS_AUDIT",json.dumps(z,allow_nan=False))
if __name__=="__main__":main()
