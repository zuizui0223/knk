"""Source-locked independent zooplankton mass/sign validation, linked donor clusters."""
import argparse,csv,hashlib,io,json,re
from pathlib import Path
import numpy as np
from benchmarks import score_ahif_mass_direction_transfer as md
from benchmarks import score_ahif_replication_scale as rs
from benchmarks.probe_ahif_villarreal2025 import get,ARTICLE,DOI
N=256;SEED=2026112613
RAW="Dadosbrutosfinais1.1.txt";CODE="1. Analyses-v7-06Oct2025.R"
EXPECTED={
 "A0":[1,5,6,7,8,18,20,28,30,38],
 "Aa":[16,15,4,26,22,2,23,25,21,13],
 "B0":[3,11,12,17,19,24,29,32,36,40],
 "Ab":[31,35,14,34,9,27,33,39,10,37]
}
DATES=("T1","T2","T3","T4","T5")
PANELS=(("A0","Aa"),("A0","Ab"),("Aa","Ab"))
def acquire(root):
    receipt=json.loads(Path("validation/ahif_v1/villarreal2025_source_receipt_20261126.json").read_text())
    assert receipt["n_checksum_verified"]==6 and receipt["doi"]==DOI
    refs={f["name"]:f for f in receipt["files"]}
    m=json.loads(get(f"https://api.figshare.com/v2/articles/{ARTICLE}/versions/1",2_000_000))
    assert int(m["id"])==ARTICLE and int(m["version"])==1 and m["doi"]==DOI
    fs={f["name"]:f for f in m["files"]};assert set(fs)==set(refs)
    root.mkdir(parents=True,exist_ok=True)
    for name in (RAW,CODE):
        f=fs[name];r=refs[name]
        assert f["id"]==r["file_id"] and int(f["size"])==r["size"]
        assert (f.get("computed_md5") or f.get("supplied_md5"))==r["md5"]
        b=get(f.get("downloadUrl") or f.get("download_url"),r["size"]+100)
        assert len(b)==r["size"] and hashlib.md5(b).hexdigest()==r["md5"] and hashlib.sha256(b).hexdigest()==r["sha256"]
        (root/name).write_bytes(b)
    return refs
def vector(script,name):
    # Source code is inspected, never executed.
    pat=r"(?m)^\s*"+name+r"\s*(?:<-|=)\s*c\s*\(([^)]*)\)"
    v=re.search(pat,script)
    if not v:
        raise ValueError("missing explicit author "+name+" vector; source lines 20–32: "+repr(script.splitlines()[19:32]))
    raw=re.sub(r"#.*","",v.group(1))
    assert re.fullmatch(r"[\d,\s]+",raw),("unexpected pair expression",name,raw)
    return [int(x) for x in raw.split(",") if x.strip()]
def load(root,refs):
    script=(root/CODE).read_text(encoding="utf-8-sig",errors="strict")
    donor=vector(script,"donor");recipient=vector(script,"recipient")
    assert len(donor)==len(recipient)==20 and len(set(donor+recipient))==40
    pairs=list(zip(donor,recipient))
    a=[(d,r) for d,r in pairs if d in EXPECTED["A0"]]
    b=[(d,r) for d,r in pairs if d in EXPECTED["B0"]]
    # Author's tank lists specify pair ordering; explicit parsed vectors must agree.
    assert set(a)==set(zip(EXPECTED["A0"],EXPECTED["Aa"])),a
    assert set(b)==set(zip(EXPECTED["B0"],EXPECTED["Ab"])),b
    assert len(a)==len(b)==10
    raw=(root/RAW).read_bytes()
    # The original file has a non-UTF8 Portuguese header. Use lossless decoding; all numeric columns remain unchanged.
    try:
        decoded=raw.decode("utf-8-sig");encoding="utf-8-sig"
    except UnicodeDecodeError:
        decoded=raw.decode("latin-1");encoding="latin-1";assert decoded.encode(encoding)==raw
    r=csv.reader(io.StringIO(decoded,newline=""),delimiter="\t");header=next(r);rows=list(r)
    assert len(header)==len(set(header))==70 and header[:3]==["Tanques","Data","Tratamento"]
    assert len(rows)==200 and all(len(x)==70 for x in rows)
    removed={27,28,57,58}
    keep=[i for i in range(3,70) if i not in removed];taxa=[header[i] for i in keep]
    assert len(taxa)==63
    expected={int(t):arm for arm,ts in EXPECTED.items() for t in ts}
    dat={}
    for r in rows:
        tank=int(r[0]);day=r[1];arm={"AA":"Aa","AB":"Ab"}.get(r[2],r[2])
        assert day in DATES and expected[tank]==arm and (tank,day) not in dat
        x=np.asarray([float(r[i]) for i in keep])
        assert np.isfinite(x).all() and (x>=0).all() and (x==np.floor(x)).all()
        dat[tank,day]=x
    assert set(dat)=={(t,d) for t in range(1,41) for d in DATES}
    allx=np.stack(list(dat.values()))
    assert allx.sum()==88774 and allx.sum(1).min()==126 and allx.sum(1).max()==1225
    assert (allx.sum(0)>0).all()
    mats={arm:np.array([[dat[t,d] for t in ids] for d in DATES]) for arm,ids in EXPECTED.items()}
    audit={"source_doi":DOI,"source_files":{k:{"md5":refs[k]["md5"],"sha256":refs[k]["sha256"]} for k in (RAW,CODE)},
       "tank_pairs":pairs,"n_tanks":40,"n_linked_pairs":20,"n_profiles":200,"n_post_profiles":120,
       "n_taxa":len(taxa),"source_encoding":encoding,"taxon_codes":taxa,"excluded_columns":[header[i] for i in sorted(removed)],
       "individuals":int(allx.sum()),"sample_total_range":[int(allx.sum(1).min()),int(allx.sum(1).max())],
       "post_dates_days":[65,86,101],"arm_sample_counts":{k:int(v.shape[0]*v.shape[1]) for k,v in mats.items()}}
    return mats,audit
def split(mats,orders,n,panel):
    ca,cb=panel;g=lambda arm: "A" if arm in ("A0","Aa") else "B"
    a,b=mats[ca],mats[cb];oa,ob=orders[g(ca)],orders[g(cb)]
    assert n+n<=10
    ds=[];hs=[]
    for o in (0,1):
        ia,ha=(oa[:,:n],oa[:,-n:]) if o==0 else (oa[:,-n:],oa[:,:n])
        ib,hb=(ob[:,:n],ob[:,-n:]) if o==0 else (ob[:,-n:],ob[:,:n])
        assert not np.any(ia[:,:,None]==ha[:,None,:]) and not np.any(ib[:,:,None]==hb[:,None,:])
        if g(ca)==g(cb):assert np.array_equal(ia,ib) and np.array_equal(ha,hb)
        ds.append((b[:,ib,:].mean(2)-a[:,ia,:].mean(2)).transpose(1,0,2))
        hs.append((b[:,hb,:].mean(2)-a[:,ha,:].mean(2)).transpose(1,0,2))
    return np.stack(ds,axis=1),np.stack(hs,axis=1)
def score(mats):
    rng=np.random.default_rng(SEED)
    orders={k:np.array([rng.permutation(10) for _ in range(N)]) for k in ("A","B")}
    out={}
    for panel in PANELS:
        vv={};oo={};ll={}
        for n in (2,4):
            d,h=split(mats,orders,n,panel)
            v,ok=md.factorial(d[:,:,2:3],d[:,:,3:4],h[:,:,4:5])
            w,q=rs.scores(d[:,:,3:4],d[:,:,4:5],h[:,:,4:5])
            vv[n]=v;ll[n]=w;oo[n]=ok&q
        common=oo[2]&oo[4];assert common.any()
        pp={}
        for n in (2,4):
            def sum_(mask):
                x=dict(zip(md.KEYS,vv[n][mask].mean(0).tolist()))
                y=dict(zip(rs.METRICS,ll[n][mask].mean(0).tolist()))
                o=1-y["allocation_tv"]
                y.update(O=o,Astar=.5+y["H"]/o if o>1e-12 else None,J=y["H"]-y["G"])
                assert abs(x["total_gain"]-x["mass_symmetric"]-x["sign_symmetric"])<1e-10
                assert abs(x["RR"]-y["lag_signed_tv"])<1e-10
                return {"n":int(mask.sum()),"factorial":x,"one_lag":y}
            pp[str(n)]={"n_total":int(common.size),"common":sum_(common),"all_defined":sum_(oo[n])}
        out[panel[1]+"-"+panel[0]]=pp
    primary={}
    for n in ("2","4"):
        p=[out["Aa-A0"][n],out["Ab-A0"][n]]
        a=[r["common"] for r in p]
        primary[n]={"n_common":sum(r["n"] for r in a),"n_total":sum(r["n_total"] for r in p),
           "factorial":{k:float(np.mean([r["factorial"][k] for r in a])) for k in md.KEYS},
           "one_lag":{k:float(np.mean([r["one_lag"][k] for r in a])) for k in rs.METRICS}}
        y=primary[n]["one_lag"];o=1-y["allocation_tv"];y.update(O=o,Astar=.5+y["H"]/o if o>1e-12 else None,J=y["H"]-y["G"])
    return primary,out
def main():
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    md.selftest();rs.selftest()
    refs=acquire(a.root);mats,audit=load(a.root,refs)
    print("VILLARREAL_COUNT_GRAIN_PASS",json.dumps({k:v for k,v in audit.items() if k!="taxon_codes"}),flush=True)
    p,q=score(mats)
    z={"status":"INDEPENDENT_ZOOPLANKTON_EXTERNAL_VALIDATION_NOT_META","source_audit":audit,"seed":SEED,"n_splits":N,"primary":p,"by_panel":q,
      "notes":["Source-recipient linked cluster split; current A0 donor is reference, not unmanipulated water control","Ab-Aa is secondary; no p-values from dependent folds"]}
    a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(z,allow_nan=False)+"\n")
    print("VILLARREAL_EXTERNAL_RESULT",json.dumps(z,allow_nan=False),flush=True)
if __name__=="__main__":main()
