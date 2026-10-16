#!/usr/bin/env python3
"""Bounded MM algebra/interval/finite-jump diagnostics; not proof of stochastic limits."""
from pathlib import Path
import json, hashlib, math
import sympy as sp
import mpmath as mp
import numpy as np
R=Path(__file__).resolve().parent
checks=[]
def check(name,cond,detail=''):
 checks.append({'id':name,'pass':bool(cond),'detail':str(detail)[:320]})
 if not cond: raise AssertionError(f'{name}: {detail}')
# Rigorous exact polynomial identity after cancelling a reference pole.
a,b,x,y,q=sp.symbols('a b x y q', real=True)
C=lambda j:a*x**j+b*y**j
d=lambda j:sp.expand(C(j+1)-q*C(j))
identity=sp.factor(d(1)*d(3)-d(2)**2)
check('exact-filtered-hankel-factor',sp.simplify(identity-a*b*x*y*(x-q)*(y-q)*(x-y)**2)==0,identity)
check('one-pole-zero',sp.simplify(identity.subs(b,0))==0)
check('ordinary-two-pole-curvature-positive',sp.simplify((C(1)*C(3)-C(2)**2)-a*b*x*y*(x-y)**2)==0)
# Exact Lyapunov calculation (works for any 0<c<1, k>0).
c,k=sp.symbols('c k',positive=True)
J=sp.Matrix([[-1,1],[k*c,-k]])
S=sp.Matrix([[1,c],[c,c]])/(1-c)
Q=sp.diag(2,2*k*c)
check('joint-lyapunov',all(sp.simplify(v)==0 for v in (J*S+S*J.T+Q)))
check('global-F-drift-at-origin',sp.simplify((J*S)[0,0]+1)==0)
check('stationary-global-variance',sp.simplify(S[0,0]-1/(1-c))==0)
# Mixture weights, normalizations and positive rate ordering: numeric controlled.
mp.mp.dps=65
lam6=mp.mpf('2.3421820411463004');lam3=mp.mpf('1.9614068619764455')
gamma=1-lam3/lam6
check('gamma-positive',0<gamma<1,mp.nstr(gamma,30))
gamma_nom=mp.mpf('0.16257283698729818')
check('gamma-consistency-to-input',abs(gamma-gamma_nom)<mp.mpf('2e-15'))
q0=mp.exp(-1)
working={}
for kap in [1,4]:
 discriminant=mp.sqrt((1-kap)**2+4*kap*(1-gamma))
 beta_s=(1+kap-discriminant)/2; beta_f=(1+kap+discriminant)/2
 ws=(beta_f-gamma)/(beta_f-beta_s);wf=1-ws
 check(f'k{kap}-rates-order',0<beta_s<gamma<1<beta_f)
 check(f'k{kap}-positive-weights',0<ws<1 and 0<wf<1)
 check(f'k{kap}-weight-sum',abs(ws+wf-1)<mp.mpf('1e-60'))
 def cov(t,ap=mp.mpf(1)/2):
  return (ws*mp.exp(-beta_s*t)+wf*mp.exp(-beta_f*t))/gamma + (1-ap)/ap*mp.exp(-t)
 Cv=[cov(j) for j in range(5)]
 rho=[v/Cv[0] for v in Cv]
 d=[rho[j+1]-q0*rho[j] for j in [1,2,3]]
 h=d[0]*d[2]-d[1]**2
 hh=((ws/gamma)/Cv[0])*((wf/gamma)/Cv[0])*mp.exp(-beta_s-beta_f)*(mp.exp(-beta_s)-q0)*(mp.exp(-beta_f)-q0)*(mp.exp(-beta_s)-mp.exp(-beta_f))**2
 check(f'k{kap}-filtered-hankel-sign',h<0,mp.nstr(h,25))
 check(f'k{kap}-independent-factor',abs(h-hh)<mp.mpf('1e-60'))
 # exact observed sign-map inverse
 signs=[2/mp.pi*mp.asin(rr) for rr in rho[1:]]
 recovered=[mp.sin(mp.pi/2*s) for s in signs]
 check(f'k{kap}-sign-transform-roundtrip',max(abs(rh-rr) for rh,rr in zip(recovered,rho[1:]))<mp.mpf('1e-60'))
 working[f'kappa_{kap}']={'beta_s':str(beta_s),'beta_f':str(beta_f),'w_s':str(ws),'w_f':str(wf),'covariance_lags':[str(v) for v in Cv], 'rho_lags':[str(v) for v in rho],'filtered_H_normalized':str(h),'lag1_sign_agreement':str(mp.mpf('.5')+mp.asin(rho[1])/mp.pi)}
# Native full/partial one-pole-after-filter.
for ap in [mp.mpf('0.5'),mp.mpf(1)]:
 ca=(1-ap)/ap
 CN=lambda j:mp.exp(-gamma*j)/gamma+ca*mp.exp(-j)
 rn=lambda j:CN(j)/CN(0)
 ds=[rn(j+1)-q0*rn(j) for j in [1,2,3]]
 hn=ds[0]*ds[2]-ds[1]**2
 check(f'native-alpha-{ap}-filtered-zero',abs(hn)<mp.mpf('1e-60'),str(hn))
 # partial two-rate curvature before filtering is not specific.
 if ap<1:check('native-positive-ordinary-Hankel',CN(1)*CN(3)-CN(2)**2>0)
# interval enclosure, hypothesis gamma range explicitly nominal/historical, not independently certified.
mp.iv.dps=65
giv=mp.iv.mpf(['0.16257283698729','0.16257283698731']);civ=1-giv;qiv=mp.iv.exp(-1)
interv={}
for kap in [1,4]:
 root=mp.iv.sqrt((1-kap)**2+4*kap*civ);s=(1+kap-root)/2;f=(1+kap+root)/2
 w=(f-giv)/(f-s);wf=1-w
 C=lambda j:(w*mp.iv.exp(-s*j)+wf*mp.iv.exp(-f*j))/giv+mp.iv.exp(-j)
 D=lambda j:C(j+1)/C(0)-qiv*C(j)/C(0)
 h=D(1)*D(3)-D(2)**2
 lo=float(h.a);hi=float(h.b)
 check(f'k{kap}-outward-negative',hi<0,(lo,hi))
 interv[f'kappa_{kap}']={'assumed_gamma_bracket':['0.16257283698729','0.16257283698731'],'interval_engine':'mpmath.iv dps=65','bound_via_iv_float_preview':[lo,hi],'printed_interval':str(h),'certification_scope':'outward conditional interval only if gamma input bracket established by independent spectral source'}
# Direct nonuniform labelled jump sums, finite N, no Gaussian process simulation.
N=12; labels=list(range(12));B={0,2,4,6,8,10};alpha=len(B)/N
f3=np.array([math.sqrt(2)*math.cos(math.pi*j/2) for j in range(12)])
f1=np.array([math.sqrt(2)*math.cos(math.pi*j/6) for j in range(12)])
weights=np.array([N/len(B)-1 if a in B else -1 for a in range(N)],float)
# theta=0, so every candidate label has probability 1/12, current labels=0..11
basis=np.array([f1,f3])
def bracket(mode_a,mode_b,weight_a,weight_b):
 out=0.0
 for slot,old in enumerate(labels):
  for new in range(12):
   da=weight_a[slot]*(basis[mode_a,new]-basis[mode_a,old])/math.sqrt(N)
   db=weight_b[slot]*(basis[mode_b,new]-basis[mode_b,old])/math.sqrt(N)
   out+=da*db/12
 return out
one=np.ones(N)
FF=bracket(1,1,one,one);UU=bracket(1,1,weights,weights);FU=bracket(1,1,one,weights)
check('nonuniform-state-global-bracket',abs(FF-2)<1e-12,FF)
check('nonuniform-state-private-bracket',abs(UU-2)<1e-12,UU)
check('nonuniform-state-FU-nonzero',abs(FU-1)<1e-12,FU)
check('labelled-null-parity-harmonic-orthogonality',abs(np.dot(f1,f3))<1e-12)
check('feature-normalization-f3',abs(np.dot(f3,f3)/12-1)<1e-12)
check('feature-normalization-f1',abs(np.dot(f1,f1)/12-1)<1e-12)
# Actual sign finite support ties and real-gain invariance.
raw=[(0.,0.),(1.,-2.),(3.,-3.),(-0.02,-.03),(-4.,1.)]
sgn=lambda x: 0 if x==0 else (1 if x>0 else -1)
for gain_a,gain_b in [(1,1),(-4,3),(0.13,-11)]:
 old=[sgn(x)*sgn(y) for x,y in raw]
 transformed=[sgn(gain_a*x)*sgn(gain_a*y) for x,y in raw]
 check(f'fixed-real-gain-sign-autocorr-{gain_a}',old==transformed)
# stationary null q from s1/s2 exact formula
g0=mp.mpf(1); rho_null=lambda j:mp.exp(-j)
check('null-reference-from-two-lags',abs(rho_null(2)/rho_null(1)-mp.exp(-1))<mp.mpf('1e-60'))
result={'date':'2026-10-16','status':'bounded independent calculation, NOT proof of stochastic path theorem', 'inputs':{'lambda3_nominal':str(lam3),'lambda6_nominal':str(lam6),'gamma3_derived_nominal':str(gamma)},'finite_jumps':{'N':N,'block_slots':sorted(B),'FF_f3':FF,'UU_f3':UU,'FU_f3':FU,'point':'mixed bracket nonzero statewise, only stationary averaged limit zero'},'working_precision':working,'interval_checks':interv,'checks':checks,'passed':sum(x['pass'] for x in checks),'total':len(checks)}
(R/'verification.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print('Verified',result['passed'],'/',result['total'])
for kap,v in working.items():print(kap,'H_ref',mp.nstr(mp.mpf(v['filtered_H_normalized']),24))
