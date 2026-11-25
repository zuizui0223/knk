"""Illustrative synthetic controls: no change in latent ecological abundance.
Parameters fixed a priori, not calibrated to any empirical source.
"""
import json
import numpy as np
SEED=2026112512
N=6000;K=40;T=6
rng=np.random.default_rng(SEED)
def norm(x):return x/x.sum(-1,keepdims=True)
def signmass(m,s):
    return np.concatenate([m*np.where(s>0,1,np.where(s<0,0,.5)),m*np.where(s<0,1,np.where(s>0,0,.5))],axis=-1)
def tv(a,b):return np.abs(a-b).sum(-1)/2
def compare(old,rec,target):
    ok=(np.abs(old).sum(-1)>0)&(np.abs(rec).sum(-1)>0)&(np.abs(target).sum(-1)>0)
    o,r,y=old[ok],rec[ok],target[ok];a,b,q=norm(np.abs(o)),norm(np.abs(r)),signmass(norm(np.abs(y)),y)
    oo=tv(signmass(a,o),q);ro=tv(signmass(b,o),q);orr=tv(signmass(a,r),q);rr=tv(signmass(b,r),q)
    metrics={'OO':oo,'RO':ro,'OR':orr,'RR':rr,'recent_mass_only':tv(a,norm(np.abs(y)))-tv(b,norm(np.abs(y))), 'recent_mass_old_sign':oo-ro,'recent_sign_recent_mass':ro-rr,'hybrid_over_recent':rr-ro,'total_recent':oo-rr}
    return {'n':int(ok.sum()),'metrics':{k:{'mean':float(v.mean()),'mc_se':float(v.std(ddof=1)/np.sqrt(len(v)))} for k,v in metrics.items()}}
res={}
for scenario in ('fixed_observation','temporally_correlated_detection'):
    base=rng.lognormal(np.log(12),.65,(N,K));sign=rng.choice([-1,1],(N,K));treatment=base*np.exp(.35*sign)
    if scenario=='fixed_observation':
        w=np.ones((N,T,K))
    else:
        x=np.empty((N,T,K));x[:,0]=rng.normal(0,.55,(N,K));phi=.85
        for t in range(1,T):x[:,t]=phi*x[:,t-1]+np.sqrt(1-phi**2)*rng.normal(0,.55,(N,K))
        w=np.exp(x-.55**2/2)
    # Training and held observations are independent conditional on the shared detection factor.
    d=(rng.poisson(2*treatment[:,None]*w)-rng.poisson(2*base[:,None]*w))/2
    h=(rng.poisson(2*treatment*w[:,-1])-rng.poisson(2*base*w[:,-1]))/2
    recent=d[:,4];old=d[:,3];cumulative=d[:,:4].mean(1)
    res[scenario]={'equal_one':compare(old,recent,h),'cumulative_vs_recent':compare(cumulative,recent,h)}
print(json.dumps({'status':'ILLUSTRATIVE_OBSERVATION_NULL_NOT_EMPIRICAL','seed':SEED,'n':N,'taxa':K,'times':T,'train_per_arm':2,'held_per_arm':2,'results':res},indent=2))
