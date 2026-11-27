"""Finite, observed density response contrasts from author-processed five-species CFU ratios."""
import argparse,csv,hashlib,io,itertools,json,math,statistics as st,urllib.request
from pathlib import Path
from collections import Counter
REF="14c9ac7e3b8ab05447fbe19519b1e3c8f94afea3"
BASE="https://raw.githubusercontent.com/padpadpadpad/five_species_community_paper/"+REF+"/"
FILES={"data/interaction_data.csv":"3d3c347465cb5694e30f5c4cef814416bed3814a","scripts/phenotypic_assays.R":"aa7b38b96990ec8ef437d929491dff3e39e9c49e"}
SP=tuple("aopsv")
def subsets(xs):
    return [tuple(v) for n in range(len(xs)+1) for v in itertools.combinations(xs,n)]
def get(path):
    with urllib.request.urlopen(urllib.request.Request(BASE+path,headers={"User-Agent":"AHIF academic source audit"}),timeout=60) as r:b=r.read(1_000_001)
    assert len(b)<=1_000_000 and hashlib.sha1(b"blob "+str(len(b)).encode()+b"\0"+b).hexdigest()==FILES[path],path
    return b
def load():
    data={p:get(p) for p in FILES}
    script=data["scripts/phenotypic_assays.R"].decode("utf-8-sig")
    assert "rel_fit: measured relative fitness of focal species compared to growth in monoculture" in script
    r=list(csv.DictReader(io.StringIO(data["data/interaction_data.csv"].decode("utf-8-sig"))))
    w={};rows=[]
    assert len(r)==75
    for a in r:
        i=a["focal_spp"];S=tuple(sorted(a["other_spp"]))
        assert i in SP and i not in S and all(j in SP for j in S) and len(set(S))==len(S)
        key=(i,S);assert key not in w
        x=float(a["rel_fit"]);assert math.isfinite(x) and x>=0
        w[key]=x
        if len(S)>1:
            e=float(a["expected_fitness_tot"]);h=float(a["tot_hoi"])
            assert math.isfinite(e) and math.isfinite(h) and abs(x-e-h)<1e-10
        else:
            assert a["expected_fitness_tot"] in ("NA","") and a["tot_hoi"] in ("NA","")
        rows.append({"focal":i,"partners":"".join(S),"n_partners":len(S),"relative_density":x,
                     "finite_response_vs_mono":x-1,"measurement":"author_processed_CFU_relative_density",
                     "meta_eligible":"no"})
    for i in SP:
        expected=subsets([j for j in SP if j!=i])
        assert all((i,S) in w for S in expected if S)
        w[i,()]=1.
    receipt={p:{"git_blob":FILES[p],"sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b)} for p,b in data.items()}
    return w,rows,receipt
def sign(x,eps=0):
    return 1 if x>eps else -1 if x< -eps else 0
def compare_sign(rows,eps):
    counts=Counter()
    for r in rows:
        a=sign(r["pair_marginal"],eps);b=sign(r["context_marginal"],eps)
        if not a or not b:counts["unresolved"]+=1
        else:counts["same" if a==b else "opposite"]+=1
    n=counts["same"]+counts["opposite"]
    return {"n_rows":len(rows),"same":counts["same"],"opposite":counts["opposite"],"unresolved":counts["unresolved"],
            "agreement_defined":counts["same"]/n if n else None}
def error_summary(rows,models,group):
    if not rows:return None
    groups={}
    for r in rows:groups.setdefault(tuple(r[k] for k in group),[]).append(r)
    means={m:st.mean(st.mean(abs(r["target"]-r[m]) for r in g) for g in groups.values()) for m in models}
    return {"n_rows":len(rows),"n_groups":len(groups),"MAE":means}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    w,endpoint,receipt=load();context=[];full=[]
    for i in SP:
        other=[j for j in SP if j!=i]
        g=math.prod(w[i,(j,)] for j in other)**.25
        for S in subsets(other):
            if len(S)<2:continue
            z=w[i,S]
            full.append({"focal":i,"partners":"".join(S),"n_partners":len(S),"target":z,
                         "no_response":1.,"actual_product":math.prod(w[i,(j,)] for j in S),"anonymous_product":g**len(S)})
        for j in other:
            pair=w[i,(j,)];base_delta=pair-1
            for S in subsets([k for k in other if k!=j]):
                t=tuple(sorted(S+(j,)));b=w[i,S];y=w[i,t]
                cand=[k for k in other if k not in S]
                context.append({"focal":i,"added":j,"background":"".join(S),"background_size":len(S),
                                "pair_marginal":base_delta,"background_density":b,"target":y,
                                "context_marginal":y-b,"zero":b,"additive":max(0.,b+base_delta),
                                "multiplicative":b*pair,"anonymous":b*st.mean(w[i,(k,)] for k in cand)})
    assert len(context)==160 and len(full)==55
    bg=[r for r in context if r["background_size"]>0]
    assert len(bg)==140
    idrows=[r for r in bg if r["background_size"]<3];assert len(idrows)==120
    fullid=[r for r in full if r["n_partners"]<4];assert len(fullid)==50
    sm={}
    for eps in (0,.05):
        sm[str(eps)]={"overall":compare_sign(bg,eps),"by_background_size":{str(n):compare_sign([r for r in bg if r["background_size"]==n],eps) for n in (1,2,3)}}
    edge=[]
    for i in SP:
        for j in SP:
            if i==j:continue
            r=[x for x in bg if x["focal"]==i and x["added"]==j];assert len(r)==7
            c=compare_sign(r,0)
            edge.append({"focal":i,"added":j,**c,"sign_reversal_in_any_background":c["opposite"]>0,
                         "threshold_0.05":compare_sign(r,.05)})
    models=("zero","additive","multiplicative","anonymous")
    fm=("no_response","actual_product","anonymous_product")
    result={"status":"AUTHOR_PROCESSED_FINITE_DENSITY_RESPONSE_NOT_GLV_NOT_META",
      "source_commit":REF,"publication":"10.1099/mic.0.001489","source_files":receipt,
      "coverage":{"focal_species":5,"author_relative_density_rows":75,"pair_alone_rows":20,
                  "background_present_rows":140,"primary_identity_rows":120,"higher_richness_endpoints":55},
      "zero_endpoints":sum(r["relative_density"]==0 for r in endpoint),
      "sign":sm,"edges_with_any_reversal":sum(r["sign_reversal_in_any_background"] for r in edge),
      "edges_with_any_reversal_0.05":sum(r["threshold_0.05"]["opposite"]>0 for r in edge),
      "pair_alone_signs":dict(Counter(str(sign(r["pair_marginal"])) for r in context if r["background_size"]==0)),
      "marginal_forecasts":{"all_backgrounds":error_summary(bg,models,("focal","added")),
          "primary_identity":error_summary(idrows,models,("focal","added")),
          "by_background_size":{str(n):error_summary([r for r in bg if r["background_size"]==n],models,("focal","added")) for n in (1,2,3)},
          "by_focal":{i:error_summary([r for r in idrows if r["focal"]==i],models,("focal","added")) for i in SP}},
      "full_endpoint_forecasts":{"all":error_summary(full,fm,("focal",)),
          "identity_comparable":error_summary(fullid,fm,("focal",)),
          "by_partner_count":{str(n):error_summary([r for r in full if r["n_partners"]==n],fm,("focal",)) for n in (2,3,4)}},
      "edge_summaries":edge}
    a.out.mkdir(parents=True,exist_ok=True)
    for name,r in (("finite_density_endpoints",endpoint),("context_marginal_responses",context),("full_endpoint_forecasts",full)):
        with (a.out/(name+".csv")).open("w",newline="") as f:
            wr=csv.DictWriter(f,fieldnames=list(r[0]));wr.writeheader();wr.writerows(r)
    (a.out/"summary.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("CASTLEDINE_CONTEXT_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
