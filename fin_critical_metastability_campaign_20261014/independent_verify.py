#!/usr/bin/env python3
"""Independent small-checker: no theory_core or mt_analysis imported.
Rational domain checks are exact. Spectra/sample generator probes are binary64.
No code-based proof of global gradient/canonical path/stochastic claims.
"""
import math, json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
BASE=Path(__file__).resolve().parent
checks=[]
def ck(name,good,details=None,scope='binary64_diagnostic'):
 checks.append({'name':name,'pass':bool(good),'scope':scope,'details':details})
EPS=Q(1,100); DT=EPS**2; H=Q(1,10**10)
GAMMA=Q(93,5000)
LOW=1/((1+EPS)/12+DT) - 12*(1-GAMMA) - Q(1,1000)*Q(12,5)
ck('exact rational Hessian floor >0.08',LOW>Q(2,25),str(LOW),'exact_fraction_given_spectral_premise')
ck('exact parity boundary > h', EPS**4/24>H, str(EPS**4/24),'exact_fraction')
ck('exact fast boundary > h', Q(1,25)*DT**2>H, str(Q(1,25)*DT**2),'exact_fraction')
ck('positive b=h/2',H/2>0,str(H/2),'exact_fraction')
q=12;j=np.arange(q);omega=.18575;phi=.1625;eta=1.8
# independently constructed circulant Laplacian and Fourier projector
row=np.array([0. if d==0 else math.cos(omega*min(d,q-d)+phi)/(1+min(d,q-d)**eta) for d in range(q)])
lap0=np.sum(row)
lam=np.array([lap0-np.sum(row*np.cos(2*math.pi*k*j/q)) for k in range(7)])
A=sum((lam[k]/6)*(np.outer(np.cos(2*math.pi*k*j/q),np.cos(2*math.pi*k*j/q))+
  np.outer(np.sin(2*math.pi*k*j/q),np.sin(2*math.pi*k*j/q))) for k in (3,4,5))+(lam[6]/12)*np.outer((-1.)**j,(-1.)**j)
chi=lam[6]/q;gc=1/chi;u=np.ones(q)/q;sig=(-1.)**j
E=np.linalg.svd(np.column_stack([np.ones(q),sig]).T)[2][2:].T
ck('A positive semidefinite',np.linalg.eigvalsh(A).min()>-2e-14,float(np.linalg.eigvalsh(A).min()))
ck('A constant diagonal',np.ptp(np.diag(A))<1e-14,float(np.ptp(np.diag(A))))
ck('uniform critical top parity', all(lam[6]>lam[k] for k in (3,4,5)),[lam[k] for k in (3,4,5,6)])
ck('critical gain inherited interval',5.12342755139859<gc<5.12342755139864,gc)
ck('fast gap inherited lower',1-max(lam[k]/lam[6] for k in (3,4,5))>.0186,float(1-max(lam[k]/lam[6] for k in (3,4,5))))
ck('uniform fast Hess positive',np.linalg.eigvalsh(E.T@(12*np.eye(q)-gc*A)@E).min()>.22,float(np.linalg.eigvalsh(E.T@(12*np.eye(q)-gc*A)@E).min()))
def FE(p,g):return float(sum(x*math.log(q*x) for x in p if x>0)-g/2*p@A@p)
def lgmu(n):
 n=np.array(n,int);N=int(n.sum());return math.lgamma(N+1)-sum(math.lgamma(k+1) for k in n)-N*math.log(q)
def rates(n,rule,g):
 n=np.asarray(n,int);N=int(n.sum());R=np.zeros((q,q))
 for i in range(q):
  z=g*(A@(n-np.eye(q,dtype=int)[i]))/N
  if rule=='HB':
   t=np.exp(z-max(z));t/=sum(t)
  else:
   d=z-z[i]
   if rule=='MET':t=np.minimum(1.,np.exp(np.minimum(d,0.)))/q
   elif rule=='BAR':t=(2/q)/(1+np.exp(-d))
  R[i]=n[i]*t;R[i,i]=0
 return R
for rule in ('HB','MET','BAR'):
 for move in [(0,1),(5,8),(3,6),(11,0)]:
  n=np.full(q,800,dtype=int);n[move[0]]+=2;n[move[1]]-=2;N=sum(n)
  nn=n.copy();nn[move[0]]-=1;nn[move[1]]+=1
  fw=rates(n,rule,gc)[move[0],move[1]]
  bw=rates(nn,rule,gc)[move[1],move[0]]
  ww=lambda a:lgmu(a)+(gc/(2*sum(a)))*float(a@A@a)
  residue=ww(n)+math.log(fw)-ww(nn)-math.log(bw)
  ck(f'{rule} detailed balance edge {move}',abs(residue)<1e-10,residue)
  p=n/N
  D=sum(v*math.log(v*q) for v in p)
  ck(f'{rule} count type method bound {move}', lgmu(n)<=-N*D+1e-8 and lgmu(n)>=-N*D-12*math.log(N+1)-1e-8,[lgmu(n),-N*D])
for m in [-.009,-.006,0,.006,.009]:
 p=u+m*sig/q
 g=gc
 f=FE(p,g)
 I=(1+m)/2*math.log1p(m)+(1-m)/2*math.log1p(-m)
 ck(f'parity exact F m={m}',abs(f-(I-g*chi*m*m/2))<1e-13,[f,I-g*chi*m*m/2])
 ck(f'parity gradient orthogonal to fast m={m}',np.linalg.norm(E.T@(np.log(12*p)-g*A@p))<2e-13,None)
rng=np.random.default_rng(7)
for i in range(12):
 m=rng.uniform(-.009,.009);v=E@rng.standard_normal(10);v=v/np.linalg.norm(v)*rng.uniform(0,1e-4)
 p=u+m*sig/q+v
 bound=(m**4)/24 +.04*(v@v)
 ck(f'coercivity local probe {i}',FE(p,gc)+2e-14>=bound,[FE(p,gc),bound])
# analytic stationary density quartic has positive variance, gamma formula
variance=math.sqrt(12)*math.gamma(.75)/math.gamma(.25)
ck('quartic limiting variance 1.17',1.17<variance<1.172,variance)
# no finite-N parent Gibbs implication from numerical check
ck('global Gibbs competitor at gc',FE(np.eye(q)[0],gc)<-.75, FE(np.eye(q)[0],gc))
output={'tool':'independent_verify.py','status':'bounded_checks_only_not_a_theorem_prover', 'num_checks':len(checks),'passed':sum(t['pass'] for t in checks),'checks':checks,
  'critical_input_nominal':{'gc':gc,'lambda6':float(lam[6])}, 'rational_fast_hessian_floor':str(LOW)}
(BASE/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
print(f'Independent checks {output["passed"]}/{len(checks)}; rational Hess lower {float(LOW):.7f}')
if output['passed']!=len(checks):raise SystemExit(2)
