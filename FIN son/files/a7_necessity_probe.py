"""How many independent numbers does the strict kernel (omega, phi, eta, beta) feed into A7?
Maps kernel parameters -> (lam4/lam3, lam5/lam3, lam6/lam3)  and reports Jacobian rank + robustness of ordering."""
import numpy as np
from fin_core import strict_L
p0=np.array([0.18575,0.1625,1.8,1.0])
def ratios(p):
    L=strict_L(omega=p[0],phi=p[1],eta=p[2],beta=p[3])
    return np.array([L[4]/L[3],L[5]/L[3],L[6]/L[3]]), L
r0,L0=ratios(p0)
print("kernel params (omega,phi,eta,beta) =",p0)
print("ratios lam4/lam3, lam5/lam3, lam6/lam3 =",r0)
print("L_k k=1..6:",np.round(L0[1:7],5), " retained=largest four (k=3..6); dropped=k1,k2")
J=np.zeros((3,4))
for i in range(4):
    h=1e-6*max(1,abs(p0[i])); pp=p0.copy(); pm=p0.copy(); pp[i]+=h; pm[i]-=h
    J[:,i]=(ratios(pp)[0]-ratios(pm)[0])/(2*h)
print("Jacobian d(ratios)/d(omega,phi,eta,beta):\n",np.round(J,5))
u,s,vt=np.linalg.svd(J); print("singular values:",s)
print("null direction in (omega,phi,eta,beta) space (ratios invariant to 1st order):",np.round(vt[-1],4))
# ordering / monotonicity robustness over a broad random box
rng=np.random.default_rng(0); ok=0; tot=0; spans=[]
for _ in range(4000):
    p=np.array([rng.uniform(0.0,0.6),rng.uniform(-0.5,0.5),rng.uniform(0.8,3.0),rng.uniform(0.3,3.0)])
    L=strict_L(omega=p[0],phi=p[1],eta=p[2],beta=p[3])
    tot+=1
    if L[3]<L[4]<L[5]<L[6] and L[1]<L[2]<L[3]: ok+=1; spans.append(L[6]/L[3])
print(f"random kernel box: strictly increasing L_1<...<L_6 in {ok}/{tot} = {100*ok/tot:.1f}% ; lam6/lam3 range [{min(spans):.3f},{max(spans):.3f}]")
# how much do ratios move over the box vs. the 3 measured ratios' natural scale
R=[];
for _ in range(2000):
    p=np.array([rng.uniform(0.0,0.6),rng.uniform(-0.5,0.5),rng.uniform(0.8,3.0),rng.uniform(0.3,3.0)])
    R.append(ratios(p)[0])
R=np.array(R); print("ratio ranges over box (min,max):",np.round(R.min(0),3),np.round(R.max(0),3))
