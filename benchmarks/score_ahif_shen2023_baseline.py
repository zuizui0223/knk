"""Post-score baseline-state sensitivity, same Shen cultures and fixed dates."""
import argparse,json
from pathlib import Path
from collections import Counter
import numpy as np
from benchmarks import score_ahif_shen2023_context as s
VAL=Path("validation/ahif_v1/shen2023_context_summary_20261127.json")
def eq(a,b,path=""):
    if isinstance(a,dict):
        assert a.keys()==b.keys(),path
        for k in a:eq(a[k],b[k],path+"/"+str(k))
    elif isinstance(a,list):
        assert len(a)==len(b),path
        for i,(x,y) in enumerate(zip(a,b)):eq(x,y,path+"/"+str(i))
    elif isinstance(a,(float,int)) and not isinstance(a,bool):assert abs(a-b)<1e-10,(path,a,b)
    else:assert a==b,(path,a,b)
def contrasts(means,scale):
    result=[]
    for i,f in enumerate(s.FOCAL):
        for comp,b,d in (("P",1,4),("S",2,5)):
            a,B,c,D=means[i,[0,b,3,d]]/scale[i]
            for typ,u,v in (("competitor",B-a,D-c),("predator",c-a,D-B)):
                result.append({"focal":f,"competitor":comp,"added_agent":typ,"base_marginal":float(u),"background_marginal":float(v),"difference":float(v-u),"direction":s.relation(u,v),"direction_threshold_005":s.relation(u,v,.05)})
    return {"rows":result,"counts":dict(Counter(r["direction"] for r in result)),"counts_threshold_005":dict(Counter(r["direction_threshold_005"] for r in result))}
def main():
    p=argparse.ArgumentParser();p.add_argument("--csv",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    mat,common,strict,audit=s.load(a.csv)
    old=json.loads(VAL.read_text());eq(audit,old["audit"])
    ref=mat(common)[:,0].mean((1,2))
    for lane,ds in (("primary_common",common),("after_second_predator_addition",[d for d in common if d>=9])):
        oldlane=old["results"][lane];eq(ds,oldlane["dates"]);eq(s.rows_for(mat(ds),ds,ref),oldlane["marginals"]);eq(s.forecast(mat(ds)),oldlane["forecast"])
    ds=[d for d in common if d>=9];assert ds==[9,11,13,15,17]
    x=mat(ds).mean(3);ini=mat([1])[:,:,:,0];day3=mat([3])[:,:,:,0];assert (ini>0).all()
    means=x.mean(2);b=ini.mean(2);b3=day3.mean(2)
    out={}
    for name,v in (("raw_post",means),("day1_pre_assignment_diagnostic",b),("day3_pre_predator_diagnostic",b3),("change_from_day1",(x-ini).mean(2)),("change_from_day3",(x-day3).mean(2))):
        out[name]=contrasts(v,ref)
    out["ratio_of_culture_post_to_day1"]=contrasts((x/ini).mean(2),np.ones(2))
    out["ratio_of_arm_means_post_to_day1"]=contrasts(means/b,np.ones(2))
    z=ini-b[:,:,None];y=x-means[:,:,None]
    den=(z*z).sum((1,2));assert (den>0).all()
    slope=(z*y).sum((1,2))/den
    g=ini.mean((1,2));adj=means-slope[:,None]*(b-g[:,None])
    out["within_arm_common_slope_ANCOVA"]=contrasts(adj,ref)
    support=[]
    for i,f in enumerate(s.FOCAL):
        for comp,j,k in (("P",1,4),("S",2,5)):
            ids=[0,j,3,k];r=ini[i,ids,:]
            lo=float(r.min(1).max());hi=float(r.max(1).min())
            support.append({"focal":f,"competitor":comp,"max_of_min":lo,"min_of_max":hi,"common_support":lo<=hi,"grand_mean_day1":float(g[i]),"grand_mean_inside_all_four_arms":bool(lo<=g[i]<=hi)})
    result={"status":"POST_SCORE_BASELINE_DIAGNOSTIC_NOT_CAUSAL_NOT_META","dates":ds,"n_cultures":24,"source_replay":True,"normalization_CD":dict(zip(s.FOCAL,ref.tolist())),"day1_range":[float(ini.min()),float(ini.max())],"within_arm_slopes":dict(zip(s.FOCAL,slope.tolist())),"ANCOVA_support":support,"results":out}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print("SHEN_BASELINE_RESULT",json.dumps(result,allow_nan=False),flush=True)
if __name__=="__main__":main()
