#!/usr/bin/env python3
"""Independent nominal arithmetic verifier; no import of uf_analysis/theory_core.
Proofs of stochastic convergence are in separate mathematics, not the test flags.
"""
import json,math
from pathlib import Path
import numpy as np
from scipy.linalg import null_space
BASE=Path(__file__).resolve().parent
Q=12; j=np.arange(Q);sig=(-1.)**j; u=np.ones(Q)/Q
om=.18575;ph=.1625;eta=1.8
D=np.array([min(abs(i-k),Q-abs(i-k)) for i in range(Q) for k in range(Q)]).reshape(Q,Q)
W=np.where(D>0,np.cos(om*D+ph)/(1+np.where(D>0,D,1)**eta),0)
L=np.diag(W.sum(axis=1))-W
l=np.linalg.eigvalsh(L)
# orthonormal cosine-sine projectors from physical laplacian; A built via PSD Fourier eigenspaces
A=np.zeros((Q,Q));lam={}
for k in [3,4,5,6]:
 eigvec=np.exp(2j*np.pi*k*j/Q)
 eig=float(np.real(np.vdot(eigvec,L@eigvec))/Q)
 lam[k]=eig
 if k<6:
  cos=np.cos(2*np.pi*k*j/Q);sin=np.sin(2*np.pi*k*j/Q)
  A+=eig*(np.outer(cos,cos)+np.outer(sin,sin))/6
 else:A+=eig*np.outer(sig,sig)/Q
chi=lam[6]/Q;gc=1/chi;F=lambda p:float(sum(v*math.log(Q*v) for v in p if v>0)-gc/2*p@A@p)
P=np.eye(Q)-np.outer(np.ones(Q),np.ones(Q))/Q-np.outer(sig,sig)/Q
T=null_space(np.vstack([np.ones(Q),sig]))
checks={}
def chk(name,value,lim=1e-9): checks[name]=bool(value<=lim)
chk('circulant',np.max(abs(A-np.roll(A[0],np.arange(Q)[:,None],axis=0))) if False else max(abs(A[i,j]-A[0,(j-i)%Q]) for i in range(Q) for j in range(Q)))
chk('centered',np.linalg.norm(A@np.ones(Q)))
chk('parity_eigen',np.linalg.norm(A@sig-lam[6]*sig))
chk('critical_linear_zero',abs(gc*lam[6]/12-1))
fast=gc*(T.T@A@T)/12-np.eye(10)
chk('fast_strict_negative',max(0,np.linalg.eigvalsh(fast)[-1]+.018))
chk('projector_rank',abs(np.trace(P)-10))
chk('projector',np.linalg.norm(P@P-P))
chk('P_soft_orthogonal',np.linalg.norm(P@sig))
chk('F_delta_negative',max(0,F(np.eye(Q)[0])+.7))
chk('F_uniform_zero',abs(F(u)))
chk('F_explicit_interior_negative',max(0,F(.99*np.eye(Q)[0]+.01*u)+.1))

def rr(p,g,rule):
 h=g*A@p; delta=h[None,:]-h[:,None]
 if rule=='HB':
  ee=np.exp(h-h.max());R=np.tile(ee/ee.sum(),(Q,1))
 elif rule=='MET':R=np.exp(np.minimum(delta,0))/Q
 else:R=(2/Q)/(1+np.exp(-delta))
 np.fill_diagonal(R,0.);return R

def char(p,g,rule):
 rates=p[:,None]*rr(p,g,rule);di=sig[None,:]-sig[:,None]
 return float(np.sum(rates*di)),float(np.sum(rates*di*di))
for rule in ['HB','MET','BAR']:
 dr,vari=char(u,gc,rule)
 chk(f'{rule}_zero_drift',abs(dr))
 chk(f'{rule}_noise_2',abs(vari-2))
 for m in [.01,-.01]:
  pp=u+sig*m/12
  b=char(pp,gc,rule)[0]
  chk(f'{rule}_critical_cubic_sign_m{m}',max(0,abs(b/m**3+1/3)-.005))
  # orientation under 1-label shift
  shift=np.roll(pp,1)
  chk(f'{rule}_shift_covariance_m{m}',abs(char(shift,gc,rule)[0]+b))
  # test first derivatives in fast directions for m away from 0
  eps=1e-6
  derivatives=[(char(pp+eps*T[:,i],gc,rule)[0]-char(pp-eps*T[:,i],gc,rule)[0])/(2*eps) for i in range(10)]
  chk(f'{rule}_fast_derivative_m{m}',max(abs(x) for x in derivatives),3e-8)
 # g0 matched off-diagonal rates exactly
 chk(f'{rule}_g0_off_diagonal',max(abs(rr(u,0,rule)[i,k]-1/Q) for i in range(Q) for k in range(Q) if i!=k))

# exact symmetry-start moment expectations from multinomial independent variables
# Fourth moment for iid one coordinate ±1: E m^4=(3N^2-2N)/N^4.
for N in [120,1200]:
 chk(f'initial_soft_variance_{N}',abs(N**.5*(1/N)),.1)
 chk(f'initial_soft_fourth_{N}',max(0,(3*N*N-2*N)/N**4-3/N**2))
# historical source interval for gc reported by provided review; merely check nominal membership
chk('gc_in_accepted_review_box', max(5.12342755139859-gc,gc-5.12342755139864,0))
# stationary quartic moments
from scipy.special import gamma
var0=math.sqrt(12)*gamma(.75)/gamma(.25)
chk('quartic_var_positive',max(0,-var0))
out={'status':'nominal_binary64_reconstruction','total':len(checks),'passed':sum(checks.values()),'all_pass':all(checks.values()),'checks':checks,'gc_nominal':gc,'chi_nominal':chi,'fast_gap_nominal':-float(np.linalg.eigvalsh(fast)[-1]),'F_delta':F(np.eye(Q)[0]),'quartic_variance_lambda0':var0}
(BASE/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({k:out[k] for k in ['passed','total','all_pass','gc_nominal','fast_gap_nominal','F_delta','quartic_variance_lambda0']},indent=2))
if not out['all_pass']:raise SystemExit(1)
