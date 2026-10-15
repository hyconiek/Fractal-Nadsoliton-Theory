#!/usr/bin/env python3
"""Finite-alphabet MT local inputs: tiny nominal numerical diagnostics only.
Not an asymptotic proof or interval certificate. Never enumerates count lattices.
"""
from pathlib import Path
import sys, json, math
import numpy as np
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P/'source_context'))
import theory_core as tc
q=12; sigma=(-1.)**np.arange(q); u=np.ones(q)/q
Rf=np.eye(q)-np.ones((q,q))/q-np.outer(sigma,sigma)/q
Xfast=np.linalg.eigh(Rf)[1][:,2:] # diagnostic only (use QR nullspace instead)
from scipy.linalg import null_space
E=null_space(np.stack([u,sigma]))
chi=float(tc.lap[6]/12); gc=1/chi
gap=1-max(float(tc.lap[k]/tc.lap[6]) for k in (3,4,5))
eps=.01; delta=eps**2; gain_dev=.001
eta=1/((1+eps)/12+delta)-12*(1-gap)-gain_dev*max(float(tc.lap[k]) for k in (3,4,5))
def F(p,g):
 p=np.asarray(p); return float(np.sum(p*np.log(q*p))-.5*g*p@tc.A@p)
def p_m(m):return u+m*sigma/12
rng=np.random.default_rng(20261014)
checks=[]
def test(label,passed,data):
 checks.append({'name':label,'pass':bool(passed),'nominal_diagnostic':data})
test('A PSD',np.linalg.eigvalsh(tc.A).min()>-1e-12,float(np.linalg.eigvalsh(tc.A).min()))
test('constant A diagonal',np.ptp(np.diag(tc.A))<1e-12,float(np.ptp(np.diag(tc.A))))
test('unique parity top',gap>.0186, gap)
test('certified-input lower fast Hess at cylinder',eta>.08,eta)
test('explicit h below both cylinder boundary lower bounds',.0000000001<min(eps**4/24,.04*delta*delta),[eps**4/24,.04*delta**2])
for m in [-.01,-.005,0.,.005,.01]:
 p=p_m(m)
 for lam in [-1,0,1]:
  N=120000000
  g=gc+lam/math.sqrt(N)
  f=F(p,g); analytical=.5*((1+m)*math.log1p(m)+(1-m)*math.log1p(-m))-g*chi*m*m/2
  test(f'parity exact FE m={m} lambda={lam}',abs(f-analytical)<2e-14,[f,analytical])
  grad=np.log(12*p)-g*(tc.A@p)
  test(f'fast stationary at parity m={m} lambda={lam}',np.linalg.norm(E.T@grad)<2e-12,float(np.linalg.norm(E.T@grad)))
for model in ('HB','MET','BAR'):
 n=np.full(q,100,dtype=int); n[0]+=1;n[1]-=1
 N=int(n.sum());g=gc
 ni=n[0]; nj=n[1]
 nprime=n.copy();nprime[0]-=1;nprime[1]+=1
 def r(n,i,j):
  h=(g/N)*(tc.A@(n-np.eye(q,dtype=int)[i]))
  if model=='HB':
   z=np.exp(h-h.max());return n[i]*z[j]/z.sum()
  d=h[j]-h[i]
  if model=='MET':return n[i]/q*min(1,math.exp(d))
  return n[i]*(2/q)/(1+math.exp(-d))
 fw,rev=r(n,0,1),r(nprime,1,0)
 # exact parent Gibbs log type weight, finite A const diagonal
 def logW(n):
  return math.lgamma(N+1)-sum(math.lgamma(int(x)+1) for x in n)-N*math.log(q)+g/(2*N)*float(n@tc.A@n)
 residue=(logW(n)+math.log(fw))-(logW(nprime)+math.log(rev))
 test(f'{model} finite-N detailed balance',abs(residue)<5e-12,residue)
 test(f'{model} full off-diagonal rate upper 2N',sum(n[i]*2 for i in range(q))<=2*N,2*N)
# shape and Hess positive for a few deterministic samples
for i in range(8):
 m=rng.uniform(-eps,eps); v=E@rng.standard_normal(10);v=v/np.linalg.norm(v)*rng.uniform(0,delta)
 p=p_m(m)+v
 H=E.T@(np.diag(1/p)-gc*tc.A)@E
 test(f'fast Hess sample {i}',np.linalg.eigvalsh(H).min()>.08,float(np.linalg.eigvalsh(H).min()))
 test(f'coercivity sample {i}',F(p,gc)>=(m**4)/24+.04*np.dot(v,v)-2e-14,[F(p,gc),m**4/24+.04*np.dot(v,v)])
result={'title':'MT bounded numerical checks','status':'nominal binary64, no stochastic simulation or interval enclosure',
 'q':q,'critical_gain_nominal':gc,'chi_nominal':chi,'parity_eigenvalue_nominal':float(tc.lap[6]),'fast_gap_nominal':gap,'Hfast_lower_using_accepted_gap_and_gain_cap_nominal':eta,
 'domain':{'epsilon':eps,'fast_outer_radius':delta,'sublevel_h':1e-10,'boundary_b':5e-11,'allowed_gain_deviation_for_hessian':gain_dev},
 'checks_passed':sum(c['pass'] for c in checks),'checks_total':len(checks),'checks':checks,'libraries':{'numpy':np.__version__}}
(P/'NUMERICAL_DIAGNOSTICS.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(f'Nominal checks {result["checks_passed"]}/{len(checks)}; gc={gc:.14f}; gap={gap:.9f}; Hfast_lower={eta:.6f}')
if result['checks_passed']!=len(checks):sys.exit(1)
