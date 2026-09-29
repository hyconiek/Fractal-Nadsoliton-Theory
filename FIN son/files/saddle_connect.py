"""Where do the index-1 critical points at odd d and at d=4 lead?  Integrate the N->inf mean-field drift
   dp/dt = softmax(g A p) - p   from  p* +/- eps*v  (v = unstable tangent eigenvector)."""
import numpy as np
from scipy.integrate import solve_ivp
from fin_core import *
from shell_ladder import sym_crit
A7,L=A7_default(); G=G_FROZEN
def drift(t,p): 
    h=G*(A7@p); h-=h.max(); q=np.exp(h); q/=q.sum(); return q-p
def V(p):
    pc=np.clip(p,1e-300,None); return float(np.sum(pc*np.log(Q*pc))-0.5*G*p@A7@p)
U=np.linalg.svd(np.ones((1,Q)))[2][1:].T
def describe(p):
    j=int(np.argmax(p)); return f"peak label {j:2d} (p={p[j]:.3f}) |p-u|1={np.abs(p-1/Q).sum():.3f} V={V(p):+.5f}"
vloc=V_localized(G,A7)[0]
print("V_loc =",round(vloc,5))
for d in (1,3,4,5):
    v0,p=sym_crit(G,A7,d)
    H=np.diag(1/np.clip(p,1e-300,None))-G*A7
    w,vec=np.linalg.eigh(U.T@H@U); vneg=U@vec[:,0]; vneg/=np.linalg.norm(vneg)
    print(f"\n-- d={d}: critical V={v0:+.5f}  B_d={v0-vloc:.4f}  lowest tangent Hessian eigenvalue={w[0]:.4f} (2nd={w[1]:.3f})")
    for s in (+1,-1):
        eps=1e-4; p0=p+s*eps*vneg
        if p0.min()<=0: p0=np.clip(p0,1e-12,None); p0/=p0.sum()
        sol=solve_ivp(drift,[0,4000],p0,method="LSODA",rtol=1e-10,atol=1e-13)
        pe=sol.y[:,-1]; print(f"   side {s:+d}: -> {describe(pe)}")
# how many distinct local minima exist at G?  multistart on the mean-field flow
rng=np.random.default_rng(7); found={}
for _ in range(300):
    p0=rng.dirichlet(np.full(Q,0.4))
    sol=solve_ivp(drift,[0,4000],p0,method="LSODA",rtol=1e-9,atol=1e-12); pe=sol.y[:,-1]
    key=(round(V(pe),4), int(np.argmax(pe)) if np.abs(pe-1/Q).sum()>0.05 else -1)
    found.setdefault(key[0],set()).add(key[1])
print("\ndistinct attractor energies of the mean-field flow at G (V -> set of peak labels):")
for v in sorted(found): print(f"   V={v:+.4f}  labels={sorted(found[v])}")
