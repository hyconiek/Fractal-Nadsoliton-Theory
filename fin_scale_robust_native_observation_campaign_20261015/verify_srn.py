#!/usr/bin/env python3
"""SRN bounded independent diagnostics: exact projectors, discrete block generator,
small finite hypergeometric laws, exact common-gain invariance, angle and exponents.
Not a mechanized proof of the infinite-N native/MT theorems.
"""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import pi,cos,sqrt,atan,exp,comb,floor,isclose
import json
import mpmath as mp
import numpy as np
from sympy import integer_nthroot  # exact floors, never floating pow
OUT=Path(__file__).resolve().parent
checks=[]
def check(cond,name,details=None):
    checks.append({'name':name,'pass':bool(cond),'details':details})
    if not cond: raise AssertionError(name)

q=12
p3=[[F(1,6)*[F(1),F(0),F(-1),F(0)][(i-j)%4] for j in range(q)] for i in range(q)]
p4=[[F(1,6)*[F(1),F(-1,2),F(-1,2)][(i-j)%3] for j in range(q)] for i in range(q)]
p34=[[p3[i][j]-p4[i][j] for j in range(q)] for i in range(q)]
prod=lambda a,b:[[sum(a[i][m]*b[m][j] for m in range(q)) for j in range(q)] for i in range(q)]
zero=[[F(0) for j in range(q)] for i in range(q)]
check(prod(p3,p3)==p3,'P3 exact idempotent')
check(prod(p4,p4)==p4,'P4 exact idempotent')
check(prod(p3,p4)==zero,'P3P4 orthogonal')
check(prod(p4,p3)==zero,'P4P3 orthogonal')
check(all(p34[i][i]==0 for i in range(q)),'P3-P4 diagonal fixed')
check(all(sum(p34[i])==0 for i in range(q)),'P3-P4 rows centered')
check(all(sum(p34[i][j]*(-1)**j for j in range(q))==0 for i in range(q)),'parity eigenvector fixed')
check(all(p34[i][j]==p34[j][i] for i in range(q) for j in range(q)),'P3-P4 symmetric')
check(all(p34[i][j]==p34[(i+1)%q][(j+1)%q] for i in range(q) for j in range(q)),'P3-P4 circulant')

mp.mp.dps=60
omega=mp.mpf('0.18575');phi=mp.mpf('0.16250');eta=mp.mpf('1.8')
weight=lambda d:mp.cos(omega*d+phi)/(1+mp.mpf(d)**eta)
lam=lambda k:2*sum(weight(d)*(1-mp.cos(mp.pi*k*d/6)) for d in range(1,6))+weight(6)*(1-mp.cos(mp.pi*k))
L={k:lam(k) for k in [3,4,5,6]}
gamma=1-L[3]/L[6]; gamma_p=gamma-mp.mpf('0.01')/L[6]
check(L[6]>max(L[3],L[4],L[5]),'spectral top fixed under epsilon0')
check(0<gamma_p<gamma<1,'positive fast gaps and ordering')
check(L[6]>L[5]+mp.mpf('.01'),'top stays simple control')
check(L[3]+mp.mpf('.01')<L[5],'mode3 not top transverse')
check(abs(12/L[6]-mp.mpf('5.123427551399'))<mp.mpf('0.000000001'),'known critical gain match')

f1=np.array([sqrt(2)*cos(pi*j/6) for j in range(q)])
f3=np.array([sqrt(2)*cos(pi*j/2) for j in range(q)])
check(abs((f3[0]+f1[0])-(f3[6]+f1[6])-4*sqrt(2))<1e-12,'plus pair0-6 spacing')
check(abs((f3[8]-f1[8])-(f3[2]-f1[2])-3*sqrt(2))<1e-12,'minus pair2-8 spacing')
check(abs(sum(f1)/q)<1e-12 and abs(sum(f3)/q)<1e-12,'uniform features centered')
check(abs(sum(f1*f3)/q)<1e-12,'uniform cross orthogonal')
check(abs(sum(f1*f1)/q-1)<1e-12 and abs(sum(f3*f3)/q-1)<1e-12,'normalized uniform variance1')

H=lambda a,b: int(abs(b)>abs(a))
for w in [(1,2),(2,-1),(0,0),(3,-3),(0,1),(1,0)]:
    for aa in [F(2),F(-3,2),F(1,3)]:
        check(H(float(aa)*w[0],float(aa)*w[1])==H(*w),f'exact samplewise common gain {w},{aa}')
check(H(0,0)==0 and H(1,-1)==0,'ties retained')
check(H(2,1)!=H(2,3),'separate-gain invariance not assumed')
check(H(1,0)!=H(1,3),'offset invariance not assumed')

prob=lambda gg,a:2/mp.pi*mp.atan(mp.sqrt(1+a*(1/gg-1)))
d=(1/gamma_p-1/gamma)/(2*mp.pi)
check(d>0,'angular d positive')
for alpha in [mp.mpf('0.5'),mp.mpf('0.25'),mp.mpf('0.01'),mp.mpf('0.00001')]:
    a=prob(gamma,alpha);b=prob(gamma_p,alpha)
    check(b>a, f'angular gap positive alpha{alpha}')
    if alpha<=mp.mpf('0.01'):
        check(b-a>=d*alpha/mp.mpf('2'),f'local angular gap bound alpha{alpha}')
check(abs((prob(gamma_p,mp.mpf('1e-8'))-prob(gamma,mp.mpf('1e-8')))/mp.mpf('1e-8')-d)<mp.mpf('2e-7'),'d slope finite difference tiny alpha')
# exact conditional two-type hypergeom sum for a=b=8,H16,t7 (uniform fully enumerated 2-type draws)
A,B,t=8,8,7;Hpop=A+B
support=[x for x in range(t+1) if x<=A and t-x<=B]
probs=[F(comb(A,x)*comb(B,t-x),comb(Hpop,t)) for x in support]
check(sum(probs)==1,'two-type hypergeom exact mass 1')
check(all(x>=0 for x in probs),'two-type HG nonnegative')
check(max(probs)<F(1,2),'two-type HG maximum atom small')
# N12,k6 full hypergeom (sampling without replacement) conditional covariance cross
vals=[]
for ss in combinations(range(q),6):
    vals.append([sum(f1[list(ss)])/sqrt(6),sum(f3[list(ss)])/sqrt(6)])
vv=np.array(vals)
check(vv.shape==(924,2),'N12,k6 samples 924')
check(np.max(abs(vv.mean(0)))<1e-12,'finite HG sample means zero')
cov=vv.T@vv/924
check(np.max(abs(cov-(6/11)*np.eye(2)))<1e-12,'finite HG variance k^-1/2 block normalizer')
check(abs(cov[0,1])<1e-12,'finite HG cross feature covariance')

# independently reconstruct two A matrices via projector eigen-decomposition
P6=np.outer(np.array([(-1)**j for j in range(q)]),np.array([(-1)**j for j in range(q)]))/12
P5=np.array([[2/12*cos(2*pi*5*(i-j)/12) for j in range(q)] for i in range(q)])
AA=sum(float(L[r])*mat for r,mat in [(3,np.array(p3,dtype=float)),(4,np.array(p4,dtype=float)),(5,P5),(6,P6)])
Aplus=AA+0.01*np.array(p34,dtype=float)
check(np.max(abs(AA-AA.T))<1e-13,'A7 symmetric')
check(np.max(abs(AA@np.ones(q)))<1e-12,'A7 centered')
check(np.max(abs(np.diag(Aplus)-np.diag(AA)))<1e-14,'comparator diagonal preserved')
check(np.linalg.matrix_rank(AA,tol=1e-10)==7 and np.linalg.matrix_rank(Aplus,tol=1e-10)==7,'both rank seven')
# real jump generator and bracket checks N24, k8 and nonuniform block
B=np.array([2,1,1,1,1,0,0,0,1,1,0,0],dtype=int)
O=np.array([1,2,1,2,1,2,1,1,1,1,2,1],dtype=int)
N=int(sum(B+O));k=int(sum(B));alpha=k/N
check(N==24 and k==8,'native nonuniform test N24 k8')
count=B+O;FEATURES=np.stack([f1,f3],axis=1)
for model_id,MX in [('A7',AA),('control',Aplus)]:
  def chart(bb,oo):
      n=bb+oo
      Fr=n@FEATURES/sqrt(N)
      Rr=(bb-alpha*n)@FEATURES/sqrt(k)
      return np.array([Fr[0],Rr[0],Fr[1],Rr[1]])
  x=chart(B,O);drift=np.zeros(4);bracket=np.zeros((4,4));manual=np.zeros((4,4));legacyR=np.zeros(2)
  def smooth(n):
      v= (float(12/L[6])/N)*(MX@n)
      y=np.exp(v-max(v));return y/sum(y)
  for i in range(q):
      nminus=count.copy();nminus[i]-=1
      probs=smooth(nminus)
      for block,oldc in [(True,int(B[i])),(False,int(O[i]))]:
          if oldc==0:continue
          for jj in range(q):
              rate=oldc*probs[jj]
              b2=B.copy();o2=O.copy()
              if block:b2[i]-=1;b2[jj]+=1
              else:o2[i]-=1;o2[jj]+=1
              delta=chart(b2,o2)-x
              drift+=rate*delta; bracket+=rate*np.outer(delta,delta)
              df=FEATURES[jj]-FEATURES[i]
              factors=np.array([1/sqrt(N),(1-alpha if block else -alpha)/sqrt(k)])
              delta2=np.array([df[0]*factors[0],df[0]*factors[1],df[1]*factors[0],df[1]*factors[1]])
              manual+=rate*np.outer(delta2,delta2)
              check(np.max(abs(delta-delta2))<1e-12,f'native chart delta {model_id}-{block}-{i}-{jj}')
  check(np.max(abs(bracket-manual))<1e-12,'actual native 4D brackets '+model_id)
  qbar=smooth(count)
  expect=np.zeros(2)
  for i in range(q):
      nminus=count.copy();nminus[i]-=1
      expect+=((1-alpha)*B[i]-alpha*O[i])*(FEATURES.T@smooth(nminus)-FEATURES.T@qbar)
  boundcheck=drift[[1,3]]+x[[1,3]]-expect/sqrt(k)
  check(np.max(abs(boundcheck))<1e-12,'exact private drift identity '+model_id)

from fractions import Fraction
exponents={'martingale':F(2,32)-F(1,4)+F(2,64),'generator_N':F(3,32)-F(1,4)+F(1,64),'generator_k':F(3,32)-F(63,128)+F(1,64),'smoothing_eta':-F(1,32)+F(1,64),'smoothing_lattice':-F(63,128)+F(1,64)}
check(exponents['martingale']==F(-5,32),'martingale power -5/32')
check(exponents['generator_N']==F(-9,64),'global generator power -9/64')
check(exponents['generator_k']==F(-49,128),'private generator power -49/128')
check(exponents['smoothing_eta']==F(-1,64),'eta/alpha power -1/64')
check(exponents['smoothing_lattice']==F(-61,128),'lattice/alpha power -61/128')
for Nn in [12*100,12*100000,12*(2**61)]:
    kk=int(integer_nthroot(Nn**63,64)[0]);mm=int(integer_nthroot(Nn,4)[0]);aa=kk/Nn
    check(1<=kk<=Nn and mm>=1, f'exact integer-root floor fixture {Nn}')
    check(kk**64<=Nn**63<(kk+1)**64, f'exact floor(k) brackets {Nn}')
    check(mm**4<=Nn<(mm+1)**4,f'exact floor(M) brackets {Nn}')
    if Nn>=2**64: check(kk<=Nn//2, 'k<=N/2 guaranteed only eventually')

payload={'date':'2026-10-15','scientific_scope':'bounded reproducible diagnostics only; not a mechanized proof of asymptotic claims','checks_passed':sum(x['pass'] for x in checks),'checks_total':len(checks),'checks':checks,'spectral_nominal':{f'Lambda{k}':str(L[k]) for k in L},'gamma_A7':str(gamma),'gamma_control':str(gamma_p),'ideal_angular_gap_coefficient_nominal':str(d),'angular_probabilities':{str(a):{'A7':str(prob(gamma,mp.mpf(a))),'control':str(prob(gamma_p,mp.mpf(a)))} for a in ['0.5','0.25','0.01']},'power_ledger':{k:str(v) for k,v in exponents.items()},'finite_hypergeom_covariance':cov.tolist(),'claims_not_tested':['uniform k,N anti-concentration for arbitrary source counts','stochastic MT domain/mixing/trapping theorem','all-N four-dimensional generator estimates','finite-parent asymptotic error o(1)','physical calibration or instrument source']}
(OUT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(payload,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
print(f'PASS {len(checks)}/{len(checks)} bounded checks; nominal d_ang={mp.nstr(d,18)}')
