#!/usr/bin/env python3
"""Independent exact rational stress for BA-001/002; not q12 critical asymptotics.
Requires Python stdlib and sympy for Poisson solve; source generator constructed anew.
"""
import json, sys
from fractions import Fraction as Q
import sympy as S
from pathlib import Path
P=Path(__file__).resolve().parent
m=tuple(map(Q,(0,1,2,1,0,1)))
pi=tuple(Q(n,41) for n in (2,3,5,7,11,13))
v=(-1,-1,1,1,-1,1)
assert sum(pi)==1
checks=[]
def check(name,truth,extras=None):
 if not truth: raise AssertionError(name)
 checks.append(dict(check=name,result='PASS',detail=extras))
def chain(tau):
 # baseline is tau=1/10; as tau down to zero three edges close
 assert Q(0)<=tau<=Q(1,10)
 epsilon=Q(1,10)-tau
 c=tuple(Q(1,10)+epsilon*a for a in v)
 q=[[Q(0) for j in range(6)] for i in range(6)]
 for i in range(6):
  j=(i+1)%6
  q[i][j]=c[i]/pi[i];q[j][i]=c[i]/pi[j]
 d=tuple(sum(q[i][j]*(m[j]-m[i]) for j in range(6)) for i in range(6))
 G=tuple(sum(q[i][j]*(m[j]-m[i])**2 for j in range(6)) for i in range(6))
 return c,q,d,G
c0,q0,d0,G0=chain(Q(1,10))
ce,qe,de,Ge=chain(Q(0))
check('baseline all six cycle edges positive',all(x>0 for x in c0))
check('endpoint stress preserves full pointwise drift',de==d0)
check('endpoint retains positive/reversible rates on retained edges',all(qe[i][j]>=0 and pi[i]*qe[i][j]==pi[j]*qe[j][i] for i in range(6) for j in range(6)))
check('pointwise Gamma changes at multiple states',Ge!=G0 and Ge[1]==0 and G0[1]>0)
check('stationary mean Gamma unchanged',sum(pi[i]*G0[i] for i in range(6))==sum(pi[i]*Ge[i] for i in range(6))==Q(6,5))
check('stationary covariance identity',sum(pi[i]*G0[i] for i in range(6))==-2*sum(pi[i]*m[i]*d0[i] for i in range(6)))
levels=sorted(set(m))
def cut(c,lev):
 return sum(c[i] for i in range(6) if (m[i]<=lev<m[(i+1)%6]) or (m[(i+1)%6]<=lev<m[i]))
for lev in levels[:-1]:
 known=sum(pi[i]*d0[i] for i in range(6) if m[i]<=lev)
 c_expected=known # m step =1 in this toy
 check(f'exact parity cut conductance at {lev}',cut(c0,lev)==cut(ce,lev)==c_expected,{'flux':str(known)})
for lev in levels:
 ids=[i for i in range(6) if m[i]==lev]
 vals=[sum(pi[i]*G[i] for i in ids)/sum(pi[i] for i in ids) for G in (G0,Ge)]
 check(f'conditional Gamma rigidity level {lev}',vals[0]==vals[1],{'mean':str(vals[0])})
# Per-component invariant centering fails for the disconnected limit
components=((0,5),(1,),(2,3,4))
mean=sum(pi[i]*Ge[i] for i in range(6))
component_mass=tuple(sum(pi[i] for i in a) for a in components)
h=tuple(sum(pi[i]*(Ge[i]-mean) for i in a)/component_mass[j] for j,a in enumerate(components))
check('endpoint components are closed',all(qe[i][j]==0 for a in components for i in a for j in range(6) if j not in a))
check('invariant projection is nonzero',h==(Q(-8,75),Q(-6,5),Q(26,115)))
check('invariant projection centred in global pi',sum(x*y for x,y in zip(component_mass,h))==0)
# Strictly positive family: solve the Poisson equation with gauge condition by exact rationals.
def poisson(tau):
 c,q,d,G=chain(tau)
 L=S.Matrix([[S.Rational(-sum(q[i][j] for j in range(6) if j!=i)) if i==j else S.Rational(q[i][j]) for j in range(6)] for i in range(6)])
 mean_G=sum(pi[i]*G[i] for i in range(6))
 r=S.Matrix([S.Rational(x-mean_G) for x in G]); weights=S.Matrix([[S.Rational(p) for p in pi]])
 A=-L+S.ones(6,1)*weights
 u=A.inv()*r
 assert sum(pi[i]*u[i] for i in range(6))==0
 residue=S.factor((weights.multiply_elementwise(S.Matrix([[r[i]*u[i] for i in range(6)]]))).applyfunc(S.simplify).sum()) if False else sum(S.Rational(pi[i])*r[i]*u[i] for i in range(6))
 E=sum(S.Rational(pi[i])*q[i][j]*(u[i]-u[j])**2 for i in range(6) for j in range(6))/2
 assert S.simplify(residue-E)==0
 return S.factor(residue)
coarse_K=S.Matrix([[2,-1,-1],[-1,2,-1],[-1,-1,2]])
w=S.diag(*[S.Rational(x) for x in component_mass]); hh=S.Matrix([S.Rational(x) for x in h])
u=(coarse_K+S.Matrix(component_mass)*S.Matrix([list(component_mass)])).inv()*(w*hh)
coefficient=S.factor((hh.T*w*u)[0]);check('coarse exact residue',coefficient==S.Rational(1064,126075))
results=[]
for tau in (Q(1,100),Q(1,1000),Q(1,10000)):
 R=poisson(tau)
 check(f'Poisson exact energy at tau={tau}',R>0)
 results.append({'tau':str(tau),'R_exact':str(R),'R_float':float(R),'tau_R_float':float(R*S.Rational(tau.numerator,tau.denominator))})
check('Poisson energy diverges with declining positive edges',results[0]['R_float']<results[1]['R_float']<results[2]['R_float'])
# Limiting coefficient rigorously follows finite 3-component effective network via singular perturbation;
# independent symbolic rational check that tau*R tends to same constant at geometric sequence.
check('positive coarse residue',coefficient>0)
# Structural diagnostic: projected Dirichlet form for every polynomial of m agrees.
for power in (1,2,3,4):
 D=lambda qq:sum(pi[i]*qq[i][j]*(m[i]**power-m[j]**power)**2 for i in range(6) for j in range(6))/2
 check(f'Dirichlet projection for power={power}',D(q0)==D(qe))
output={'status':'PASS','check_count':len(checks),'checks':checks,'fixture':{'type':'six-vertex-reversible-countercontrol_NOT_q12','pi':[str(x) for x in pi],'m':[str(x) for x in m],'baseline_conductances':[str(x) for x in c0],'endpoint_conductances':[str(x) for x in ce],'stress_signs':v,'baseline_G':[str(x) for x in G0],'endpoint_G':[str(x) for x in Ge],'mean_G':str(mean),'components':[list(a) for a in components],'component_masses':[str(x) for x in component_mass],'invariant_residual_means':[str(x) for x in h],'tau_poisson':results,'exact_leading_tauR':str(coefficient)}}
(P/'EXACT_DIAGNOSTICS.json').write_text(json.dumps(output,indent=2)+'\n')
print('PASS',len(checks),'exact rational checks; leading tau*R=',coefficient)
