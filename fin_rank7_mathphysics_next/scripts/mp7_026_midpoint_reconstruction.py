#!/usr/bin/env python3
import json, math
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares, root

ROOT=Path(__file__).resolve().parents[1]
INP=ROOT/'results/MP7-015_g37_local_interval_roots.json'
OUT=ROOT/'results/MP7-026_midpoint_reconstruction.json'

d=json.loads(INP.read_text())
g37=3.7
roots=[np.array([float(x) for x in c['center']],dtype=float) for c in d['certificates']]
j=np.arange(12,dtype=float)
T=np.stack([
    np.cos(2*np.pi*3*j/12),
    np.cos(2*np.pi*4*j/12),
    np.cos(2*np.pi*5*j/12),
    (-1.0)**j,
],axis=1)
norm=np.array([6.,6.,6.,12.])

def stationary_residual(loglam):
    lam=np.exp(loglam); scale=np.sqrt(lam/norm); X=T*scale
    rr=[]
    for s in roots:
        h=X@s; h-=h.max(); p=np.exp(h); p/=p.sum()
        mu=X.T@p
        rr.extend(g37*mu-s)
    return np.array(rr)

ls=least_squares(stationary_residual,np.zeros(4),xtol=1e-14,ftol=1e-14,gtol=1e-14,max_nfev=10000)
lam=np.exp(ls.x); scale=np.sqrt(lam/norm); X=T*scale

def stats(s,g):
    h=X@s; m=h.max(); p=np.exp(h-m); p/=p.sum()
    mu=X.T@p; Y=X-mu
    M=(Y*p[:,None]).T@Y
    H=np.eye(4)/g-M
    return p,mu,H

def F(s,g): return s/g-stats(s,g)[1]

def aug(x):
    s=x[:4]; g=x[4]; v=x[5:]
    p,mu,H=stats(s,g)
    return np.r_[s/g-mu,H@v,np.dot(v,v)-1]

v0=np.array([0.50736753967725,0.528685655313811,0.553760694068750,0.395498105244206])
v0/=np.linalg.norm(v0)
s0=(roots[0]+roots[1])/2
x0=np.r_[s0,3.5156447168395917,v0]
rr=root(aug,x0,method='hybr',tol=1e-12)
if not rr.success: raise RuntimeError(rr.message)
sF=rr.x[:4]; gF=float(rr.x[4]); v=rr.x[5:]; v/=np.linalg.norm(v)
pF,muF,HF=stats(sF,gF)
z=X@v; zc=z-(pF@z)
a_mid=-float(v@sF)/(gF*gF)
b_mid=-float(pF@(zc**3))
eigs=np.linalg.eigvalsh(HF)

# Use accepted midpoint a,b for navigation constants.
a=a_mid; b=b_mid
xi_pref=math.sqrt(-2*a/b)
barrier_pref=(2**2.5/3)*(abs(a)**1.5)/math.sqrt(b)
soft_pref=math.sqrt(2*abs(a)*b)

def Phi(s,g):
    h=X@s; m=h.max(); logmean=m+math.log(np.exp(h-m).mean())
    return float(s@s/(2*g)-logmean)

def solve_branch(eps,sign):
    g=gF+eps
    init=sF+sign*xi_pref*math.sqrt(eps)*v
    sol=root(lambda s:F(s,g),init,method='lm',tol=1e-13)
    if np.max(np.abs(sol.fun))>1e-9:
        raise RuntimeError(f'branch solve failed eps={eps} sign={sign}: {sol.message}')
    return sol.x

branch=[]
for eps in np.geomspace(1e-7,1e-2,26):
    sm=solve_branch(float(eps),-1); sp=solve_branch(float(eps),+1)
    xim=float(v@(sm-sF)); xip=float(v@(sp-sF))
    if xim>xip:
        sm,sp=sp,sm; xim,xip=xip,xim
    em=np.linalg.eigvalsh(stats(sm,gF+eps)[2]); ep=np.linalg.eigvalsh(stats(sp,gF+eps)[2])
    pm=Phi(sm,gF+eps); pp=Phi(sp,gF+eps)
    branch.append({
        'epsilon':float(eps),
        'xi_minus_over_sqrt_epsilon':xim/math.sqrt(eps),
        'xi_plus_over_sqrt_epsilon':xip/math.sqrt(eps),
        'symmetric_separation_ratio':(xip-xim)/(2*xi_pref*math.sqrt(eps)),
        'local_energy_difference_ratio':abs(pp-pm)/(barrier_pref*eps**1.5),
        'soft_abs_minus_ratio':abs(float(em[0]))/(soft_pref*math.sqrt(eps)),
        'soft_abs_plus_ratio':abs(float(ep[0]))/(soft_pref*math.sqrt(eps)),
        'stationarity_residual_minus':float(np.max(np.abs(F(sm,gF+eps)))),
        'stationarity_residual_plus':float(np.max(np.abs(F(sp,gF+eps)))),
    })

out={
  'task':'MP7-026 midpoint reconstruction / navigation only',
  'input':'MP7-015 certified root centers plus Fourier model contract',
  'recovered_lambda_midpoints':{str(k):float(x) for k,x in zip([3,4,5,6],lam)},
  'least_squares_max_stationarity_residual_at_two_g37_roots':float(np.max(np.abs(ls.fun))),
  'fold_midpoint':{'s':[float(x) for x in sF],'g':gF,'v':[float(x) for x in v]},
  'fold_augmented_residual_max':float(np.max(np.abs(aug(np.r_[sF,gF,v])))),
  'fold_H4_eigenvalues_midpoint':[float(x) for x in eigs],
  'a_midpoint_reconstructed':a_mid,
  'b_midpoint_reconstructed':b_mid,
  'leading_prefactors_midpoint':{'xi':xi_pref,'local_barrier':barrier_pref,'soft_curvature':soft_pref},
  'branch_navigation':branch,
  'scientific_state':'NUMERICAL_EVIDENCE_ONLY_NOT_INTERVAL_PROOF',
  'proof_blocker':'The checkpoint omits the original outward spectral intervals and the full R7P-031 fold box/source checker required to enclose higher derivatives uniformly. Midpoint reconstruction cannot replace those proof inputs.'
}
OUT.write_text(json.dumps(out,indent=2))
print(json.dumps(out,indent=2))
