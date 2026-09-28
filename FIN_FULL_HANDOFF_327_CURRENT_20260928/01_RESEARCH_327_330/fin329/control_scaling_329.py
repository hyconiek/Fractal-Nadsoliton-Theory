import numpy as np, math, json
from scipy.optimize import brentq

def stats(N,k):
 z=np.load(f'/mnt/data/defect_research_inputs/base_N{N}.npz',allow_pickle=True);st=z['states'];lab=z['labels'];pi=z['pi'];j0=np.where(lab==0)[0];p0=st[j0,0]/N;D=N-st[j0,0];b=pi[j0];b=b/b.sum();a=k*p0;M=a.max();logZ=M+math.log(float(np.sum(b*np.exp(a-M))));w=b*np.exp(a-logZ);m=float(w@p0);md=float(w@D);kl=k*m-logZ;return {'mean_p0':m,'mean_D':md,'KL_nats':kl,'KL_bits':kl/math.log(2)}
anchorN=6;anchorK=12.;anchor=stats(anchorN,anchorK);KLt=anchor['KL_nats'];rows=[]
for N in range(3,11):
 fixedK=stats(N,12.0);theta=anchorK/anchorN;fixedTheta=stats(N,theta*N)
 def f(k):return stats(N,k)['KL_nats']-KLt
 hi=12.
 while f(hi)<0 and hi<1e5:hi*=2
 kkl=brentq(f,0,hi); fixedKL=stats(N,kkl)
 rows.append({'N':N,'fixed_kappa':{'kappa':12.0,'theta':12/N,**fixedK},'fixed_theta':{'kappa':theta*N,'theta':theta,**fixedTheta},'fixed_KL_anchor':{'kappa':kkl,'theta':kkl/N,**fixedKL}})
out={'anchor':{'N':anchorN,'kappa':anchorK,'theta':anchorK/anchorN,**anchor},'identity':'mu_{N,kappa}/mu_{N,0} is proportional to exp[-(kappa/N) D]; therefore theta=kappa/N is the direct defect-tilt control parameter.','rows':rows}
open('/mnt/data/defect_research_inputs/CONTROL_SCALING_329.json','w').write(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
