import numpy as np, math
from scipy.optimize import minimize
Z=np.load('/mnt/data/FIN_NEXT_RESEARCH_20260925/quartic_tensor_certificate.npz')
alphas=np.array(Z['alphas'],int); g=3.7183448981203875
c=Z['pred'][0]+g*Z['pred'][1]
D=7

def value_grad(z):
    z=np.asarray(z,float)
    val=0.; gr=np.zeros(D)
    for coef,a in zip(c,alphas):
        if abs(coef)<1e-18: continue
        # monomial
        base=coef
        for i,p in enumerate(a):
            if p: base*=z[i]**p
        val+=base
        for i,p in enumerate(a):
            if p:
                # robust derivative without division
                t=coef*p
                for j,q in enumerate(a):
                    qq=q-(1 if j==i else 0)
                    if qq:t*=z[j]**qq
                gr[i]+=t
    return val,gr

def fun(z,sign): return sign*value_grad(z)[0]
def jac(z,sign): return sign*value_grad(z)[1]
cons={'type':'eq','fun':lambda z:np.dot(z,z)-1,'jac':lambda z:2*z}
rng=np.random.default_rng(123456)
for mode,sign in [('max',-1),('min',1)]:
    sols=[]
    seeds=[]
    # axes, pair mixtures, and random starts
    for i in range(D):
        x=np.zeros(D);x[i]=1;seeds.append(x);seeds.append(-x)
    for i in range(D):
      for j in range(i+1,D):
        x=np.zeros(D);x[i]=x[j]=1/np.sqrt(2);seeds.append(x)
    for _ in range(1500):
        x=rng.normal(size=D);x/=np.linalg.norm(x);seeds.append(x)
    for x in seeds:
        r=minimize(fun,x,args=(sign,),jac=jac,constraints=[cons],method='SLSQP',options={'ftol':1e-13,'maxiter':500,'disp':False})
        z=r.x/np.linalg.norm(r.x);v=value_grad(z)[0]
        sols.append((v,z,r.success,r.nit))
    sols.sort(key=lambda t:t[0],reverse=(mode=='max'))
    best=sols[0]
    vals=np.array([s[0] for s in sols])
    print(mode,'best',repr(float(best[0])),'success',best[2],'nit',best[3],'z',best[1].tolist())
    print(mode,'quantiles',np.quantile(vals,[0,.01,.1,.5,.9,.99,1]).tolist())
    # sector radii and hidden split using direct construction
    z=best[1]; print(mode,'sector_norms',[float(np.hypot(z[0],z[1])),float(np.hypot(z[2],z[3])),float(np.hypot(z[4],z[5])),abs(float(z[6]))])
