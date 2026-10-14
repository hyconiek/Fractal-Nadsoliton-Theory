#!/usr/bin/env python3
"""Bounded numerical checks for UF finite model. No trajectory simulation or LDP proof."""
import json, math, sys
from pathlib import Path
import numpy as np
from scipy.linalg import null_space, solve_continuous_lyapunov
BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE/'source_context'))
import theory_core as tc
Q=12
SIGMA=(-1.)**np.arange(Q)
U=np.ones(Q)/Q
PFAST=np.eye(Q)-np.ones((Q,Q))/Q-np.outer(SIGMA,SIGMA)/Q

def soft(z):
 z=np.asarray(z); v=np.exp(z-z.max()); return v/v.sum()

def rate_matrix(p,g,rule,loo=False,N=None):
 out=np.zeros((Q,Q))
 for i in range(Q):
  h=g*tc.A@(p-(np.eye(Q)[i]/N if loo else 0))
  if rule=='HB': r=soft(h)
  else:
   delta=h-h[i]
   if rule=='MET':r=np.exp(np.minimum(0,delta))/Q
   elif rule=='BAR':r=(2/Q)/(1+np.exp(-delta))
   else:raise ValueError(rule)
  out[i]=r
  out[i,i]=0
 return out

def drift(p,g,rule):
 R=rate_matrix(p,g,rule)
 flux=p[:,None]*R
 return flux.sum(axis=0)-flux.sum(axis=1)

def parity_drift(p,g,rule): return float(SIGMA@drift(p,g,rule))

def parity_bracket_coeff(p,g,rule):
 R=rate_matrix(p,g,rule)
 return float(np.sum((p[:,None]*R)*((SIGMA[None,:]-SIGMA[:,None])**2)))

def finite_generator_characteristics(n,g,rule):
 N=int(sum(n));p=n/N;R=rate_matrix(p,g,rule,True,N)
 flux=p[:,None]*R
 inc=SIGMA[None,:]-SIGMA[:,None]
 return {'N':N,'m':float(SIGMA@p), 'drift':float(np.sum(flux*inc)), 'bracket_rate_m':float(np.sum(flux*inc**2)/N), 'third_abs_m':float(np.sum(flux*np.abs(inc)**3)/(N*N)), 'total_real_offdiagonal_rate':float(N*np.sum(flux))}

def entropy_free_energy(p,g):
 p=np.asarray(p); inds=p>0
 return float(np.sum(p[inds]*np.log(12*p[inds]))-g/2*p@tc.A@p)

def main():
 lap=tc.lap;L6=float(lap[6]); chi=L6/12; gc=1/chi
 vals=np.linalg.eigvalsh(PFAST@(gc/12*tc.A-np.eye(Q))@PFAST)
 # zero eigenvalues in projected 2-d complement; evaluate fast via Fourier eigenvalues
 fast_eig=[float(gc*lap[k]/12-1) for k in [1,2,3,3,4,4,5,5]]+[-1,-1]
 # parent spectrum for rank7 features: active k3,4,5,6; hidden k1,2 (not in A)
 # Correct rank: hidden k1,2 4 dimensions, k3..5 6 dimensions = ten
 fast_eig=[-1]*4+[float(gc*lap[k]/12-1) for k in [3,4,5] for _ in (0,1)]
 out={'arithmetic_status':'binary64 proposal, NOT outward enclosure', 'lambda6':L6,'chi':chi,'critical_gain_proposal':gc,
 'active_spectrum':{str(k):float(lap[k]) for k in [3,4,5,6]},
 'uniform_fast_eigenvalues':sorted(fast_eig), 'fast_gap_proposal':-max(fast_eig),
 'at_gc_full_gibbs_competitor':{'F_uniform':entropy_free_energy(U,gc),'F_delta0':entropy_free_energy(np.eye(Q)[0],gc),'F_interior_99pct_delta':entropy_free_energy(.99*np.eye(Q)[0]+.01*U,gc)},
 'rules':{}}
 v=np.linalg.svd(PFAST)[0][:,:10]
 # use a true null-space basis of (1,sigma) instead
 T=null_space(np.vstack([np.ones(Q), SIGMA]))
 for rule in ['HB','MET','BAR']:
  RU=rate_matrix(U,gc,rule); c=U[:,None]*RU
  K=np.zeros((Q,Q));D=np.zeros((Q,Q))
  for i in range(Q):
   for j in range(i+1,Q):
    e=np.eye(Q)[j]-np.eye(Q)[i];K+=c[i,j]*np.outer(e,e)
    D+=(c[i,j]+c[j,i])*np.outer(e,e)
  matrix=gc/12*tc.A-np.eye(Q)
  B=T.T@matrix@T
  P=solve_continuous_lyapunov(B.T,-np.eye(10))
  mvals=[]
  for m in [1e-3,-1e-3,2e-3,-2e-3]:
   p=U+(m/12)*SIGMA
   b=parity_drift(p,gc,rule)
   mvals.append({'m':m,'drift':b,'cubic_coefficient':b/m**3})
  n=np.array([120]*12)
  n[0]+=5;n[2]-=5
  fin=finite_generator_characteristics(n,gc,rule)
  # exact parity derivative in orthogonal fast directions at small nonzero m
  eps=1e-6
  pm=U+1e-3*SIGMA/12
  fd=max(abs((parity_drift(pm+eps*T[:,j],gc,rule)-parity_drift(pm-eps*T[:,j],gc,rule))/(2*eps)) for j in range(10))
  # Nonlinearity coupling bound: at m=0, p=u+eps w (tiny perturbation); diagnostic only
  w=T[:,0]
  nonlinear_b=parity_drift(U+.001*w,gc,rule)
  out['rules'][rule]={'uniform_parity_drift':parity_drift(U,gc,rule),'uniform_m_bracket_coeff':parity_bracket_coeff(U,gc,rule),
  'max_g0_rate_diff':float(np.max(np.abs(rate_matrix(U,0,rule)-rate_matrix(U,0,'HB')))),
  'K_minus_Sigma_u_max':float(np.max(np.abs(K-(np.diag(U)-np.outer(U,U))))),
  'D_minus_2K':float(np.max(np.abs(D-2*K))),
  'fast_B_gap':float(-np.max(np.real(np.linalg.eigvals(B)))),'lyapunov_residual':float(np.linalg.norm(B.T@P+P@B+np.eye(10))),
  'lyapunov_min_eig':float(np.linalg.eigvalsh(P)[0]),'cubic_symmetric_numerics':mvals,
  'directional_fast_linear_response_at_m0001_nominal':fd,
  'fast_perturbation_drft_at_m0_nominal':nonlinear_b,
  'finite_loo_diagnostic':fin}
 (BASE/'DIAGNOSTICS.json').write_text(json.dumps(out,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
 print(json.dumps({'gc':gc,'gap':out['fast_gap_proposal'],'F_delta':out['at_gc_full_gibbs_competitor']['F_delta0'],'models':{k:{'cov':v['uniform_m_bracket_coeff'],'cubic':v['cubic_symmetric_numerics'][0]['cubic_coefficient']} for k,v in out['rules'].items()}},indent=2))
if __name__=='__main__': main()
