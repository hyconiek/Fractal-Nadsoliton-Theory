"""Is the harmonic cutoff k>=3 (rank 7) reproduced by a parameter-free spectral threshold rule on L_k of the strict kernel?"""
import numpy as np
from fin_core import strict_L
L=strict_L(); Lk=L[1:7]
print("L_1..L_6 =",np.round(Lk,4)); 
mult=np.array([2,2,2,2,2,1])
rules={"L_k >= mean(L_1..L_6)":Lk.mean(),
       "L_k >= multiplicity-weighted mean (11 modes)":float((mult*Lk).sum()/mult.sum()),
       "L_k >= median of 11 modes (=L_3, boundary case)":float(np.median(np.repeat(Lk,mult))),
       "L_k >= 0.5*max":0.5*Lk.max(), "L_k >= 0.8*max":0.8*Lk.max()}
for n,t in rules.items():
    keep=[k for k in range(1,7) if Lk[k-1]>=t-1e-12]; print(f"  {n:52s} thr={t:.4f} keeps modes {keep}")
rng=np.random.default_rng(3); p0=np.array([0.18575,0.1625,1.8,1.0])
def frac(box,n=4000):
    ok=0;tot=0
    for _ in range(n):
        p=box(rng); Lq=strict_L(omega=p[0],phi=p[1],eta=p[2],beta=p[3])[1:7]
        tot+=1; ok+= [k for k in range(1,7) if Lq[k-1]>=Lq.mean()]==[3,4,5,6]
    return ok/tot
print("rule 'L_k>=mean' returns exactly {3,4,5,6}:")
print("   wide kernel box (omega 0-0.6, phi +-0.5, eta 0.8-3, beta 0.3-3): %.1f%%"%(100*frac(lambda r:np.array([r.uniform(0,.6),r.uniform(-.5,.5),r.uniform(.8,3),r.uniform(.3,3)]))))
print("   +-20%% around FIN parameters:                                  %.1f%%"%(100*frac(lambda r:p0*(1+0.2*r.uniform(-1,1,4)))))
print("   +-50%% around FIN parameters:                                  %.1f%%"%(100*frac(lambda r:p0*(1+0.5*r.uniform(-1,1,4)))))
