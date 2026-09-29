"""Barrier ladder B_d between localized states at label distance d (d=1..6).
Critical point = minimum of V_g on the reflection-symmetric subspace p_j = p_{(d-j) mod 12}; index (# negative Hessian
directions on the simplex tangent space) is reported: index 1 => genuine communication saddle."""
import numpy as np, math
from scipy.optimize import minimize
from fin_core import *
def sym_crit(g,A,d,starts=6):
    orbits=[];seen=set()
    for a in range(Q):
        if a in seen: continue
        orb=sorted({a,(d-a)%Q}); orbits.append(orb); seen|=set(orb)
    def expand(eta):
        s=np.exp(eta-eta.max()); p=np.zeros(Q)
        for si,orb in zip(s,orbits):
            for a in orb: p[a]=si
        return p/p.sum()
    def f(eta):
        p=expand(eta); pc=np.clip(p,1e-300,None)
        return float(np.sum(pc*np.log(Q*pc))-0.5*g*p@A@p)
    best=None
    rng=np.random.default_rng(d)
    for s in range(starts):
        eta0=np.full(len(orbits),-6.0)
        for i,orb in enumerate(orbits):
            if 0 in orb: eta0[i]=0.0
        if s: eta0+=0.5*rng.standard_normal(len(orbits))
        r=minimize(f,eta0,method="BFGS",options={"gtol":1e-11})
        p=expand(r.x)
        if best is None or r.fun<best[0]: best=(r.fun,p)
    return best
def hess_index(p,g,A):
    # Hessian of V wrt p on tangent space sum dp=0:  diag(1/p) - g A ; restrict to tangent
    H=np.diag(1/np.clip(p,1e-300,None))-g*A
    U=np.linalg.svd(np.ones((1,Q)))[2][1:].T      # basis of sum=0 subspace
    w=np.linalg.eigvalsh(U.T@H@U)
    return int((w<-1e-8).sum())
def ladder(g,A,label):
    vl=V_localized(g,A)[0]; out={}
    for d in range(1,7):
        v,p=sym_crit(g,A,d); out[d]=(v-vl,hess_index(p,g,A),float(np.abs(p-1/Q).sum()))
    print(label,f"(V_loc={vl:.5f})")
    print("   d : "+"  ".join(f"{d:6d}" for d in out))
    print("   B_d: "+"  ".join(f"{out[d][0]:6.3f}" for d in out))
    print("   idx: "+"  ".join(f"{out[d][1]:6d}" for d in out))
    print("   |p-u|1: "+"  ".join(f"{out[d][2]:6.3f}" for d in out),flush=True)
    return out
if __name__=="__main__":
    A7,L=A7_default(); Ap,lbar=projector_countermodel(); G=G_FROZEN
    ge=3.7183448981203777; gep=g_eq(Ap,lo=2.0,hi=8.0,n=40); gp=G*gep/ge
    ladder(G,A7,"FIN A7 @G")
    ladder(gp,Ap,"projector k>=3 @ same G/g_eq")
    # standard Potts reference: all jumps equivalent
    Ppotts=A_from_lams({k:1.0 for k in range(1,7)}); gpo=5.275369600142116*G/ge
    ladder(gpo,Ppotts,"standard Potts @ same G/g_eq (should be d-independent)")
    # robustness of the argmin shell under eigenvalue perturbation
    rng=np.random.default_rng(1); win={d:0 for d in range(1,7)}; n=0
    for _ in range(14):
        lam=np.array([L[3],L[4],L[5],L[6]])*(1+0.15*rng.uniform(-1,1,4))
        A=A_from_lams({3:lam[0],4:lam[1],5:lam[2],6:lam[3]})
        try:
            ge_=g_eq(A,lo=2.0,hi=10.0,n=40)
            if ge_ is None: continue
            g_=1.3837*ge_; vl=V_localized(g_,A)[0]
            B={d:sym_crit(g_,A,d,starts=4)[0]-vl for d in range(1,7)}
            win[min(B,key=B.get)]+=1; n+=1
        except Exception as e: pass
    print("argmin_d B_d over",n,"random +-15% eigenvalue perturbations (at G=1.3837 g_eq):",win)
