"""Synthetic negative control, not BioTIME or empirical evidence.
Independent event counts conditional on fixed species-level means/dispersion.
Two 2-event halves and disjoint 4-year windows, matching M2/M4 prescore.
"""
import numpy as np
import json
rng=np.random.default_rng(20261125)
N,S=4000,30
mu=np.exp(rng.normal(1.5,1.4,(N,S)))
def norm(x):
    total=x.sum(-1,keepdims=True)
    return np.divide(x,total,out=np.zeros_like(x,dtype=float),where=total>0)
def sim(case):
    m=mu[:,None,None,None,:]
    if case=='poisson':
        counts=rng.poisson(np.broadcast_to(m,(N,4,2,2,S)))
        variance=mu
    else:
        alpha=np.full_like(mu,.6) if case=='nb_common' else np.exp(rng.normal(-.5,1.5,(N,S)))
        scale=m*alpha[:,None,None,None,:]
        shape=1/alpha[:,None,None,None,:]
        counts=rng.poisson(rng.gamma(shape,scale,size=(N,4,2,2,S)))
        variance=mu+alpha*mu**2
    half=counts.mean(3)
    A,B=half[:,:,0,:],half[:,:,1,:]
    h=np.abs(A[:,1]-A[:,0]); d=np.abs(B[:,3]-B[:,2]); avg=A[:,:3].mean(1)
    p2=norm(np.sqrt(avg)); p1=norm(avg); p=norm(h); ok=h.sum(1)>0
    p4=np.where(ok[:,None],.5*(p2+p),p2); p3=np.where(ok[:,None],.5*(p1+p),p1)
    q=norm(d); valid=(d.sum(1)>0)&(avg.sum(1)>0)
    def tv(pred): return (.5*np.abs(q-pred).sum(1))[valid]
    scores={k:float(tv(v).mean()) for k,v in [('M1',p1),('M2',p2),('M3',p3),('M4',p4),('oracle_sqrt_variance',norm(np.sqrt(variance)))]}
    gain=tv(p2)-tv(p4)
    se=gain.std(ddof=1)/np.sqrt(len(gain))
    return {'n_scored':int(valid.sum()),'mean_TV':scores,'M4_minus_M2_improvement_pp':float(100*gain.mean()),'MC95_pp':[float(100*(gain.mean()-1.96*se)),float(100*(gain.mean()+1.96*se))], 'fraction_improved':float((gain>0).mean())}
result={case:sim(case) for case in ['poisson','nb_common','nb_species']}
print(json.dumps(result,indent=2))
