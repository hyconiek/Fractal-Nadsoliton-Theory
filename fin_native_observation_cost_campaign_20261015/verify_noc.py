#!/usr/bin/env python3
"""Reproduce bounded algebraic, generator and corrector checks for NOC.
Uses independent formulas and an explicit 12-mode spectral construction.
Does not constitute a mechanized proof of MT or large-N conclusions.
"""
import json, math
from pathlib import Path
from fractions import Fraction
import numpy as np
import mpmath as mp
ROOT=Path(__file__).resolve().parent
mp.mp.dps=70
q=12
wf=[mp.mpf(0)]+[mp.cos(mp.mpf('.18575')*d+mp.mpf('.16250'))/(1+mp.power(d,mp.mpf(9)/5)) for d in range(1,7)]
def lam(r):
 return 2*sum((wf[d]*(1-mp.cos(2*mp.pi*r*d/12)) for d in range(1,6)),mp.mpf(0))+wf[6]*(1-mp.cos(mp.pi*r))
L={r:lam(r) for r in [3,4,5,6]};gamma0=1-L[3]/L[6];gamma1=1-(L[3]+mp.mpf('.01'))/L[6]
cosf=np.array([math.sqrt(2)*math.cos(math.tau*3*j/12) for j in range(q)])
sig=np.array([(-1)**j for j in range(q)])
def A(eps):
 t=np.arange(q)
 P={r:np.array([[math.cos(math.tau*r*(i-j)/q)/6 for j in range(q)] for i in range(q)]) for r in (3,4,5)}
 return sum((float(L[r]+(mp.mpf(str(eps)) if r==3 else -mp.mpf(str(eps)) if r==4 else 0))*P[r] for r in (3,4,5)),np.zeros((q,q)))+float(L[6])*np.outer(sig,sig)/12

def softmax(v):
 w=np.exp(v-v.max());return w/w.sum()
def actual(n,b,A0):
 N=int(n.sum());k=int(b.sum());al=k/N;gc=12/float(L[6]);O=n-b
 F=float(np.dot(n,cosf)/math.sqrt(N));R=float(np.dot(b-al*n,cosf)/math.sqrt(k));
 drift=np.zeros(2);bracket=np.zeros((2,2));third=0
 for counts,inside in ((b,True),(O,False)):
  coeff=np.array([1/math.sqrt(N),(1-al if inside else -al)/math.sqrt(k)])
  for j in range(q):
   if counts[j]==0:continue
   nj=n.copy();nj[j]-=1;prob=softmax(gc/N*(A0@nj))
   deltas=cosf-cosf[j]
   drift+=counts[j]*(prob@deltas)*coeff
   bracket+=counts[j]*float(prob@(deltas*deltas))*np.outer(coeff,coeff)
   third+=counts[j]*float(prob@np.abs(deltas)**3)*float(np.linalg.norm(coeff)**3)
 return F,R,drift,bracket,third

checks=[]
def chk(name,cond):
 if not bool(cond):raise AssertionError(name)
 checks.append(name)
for eps in (0,.01):
 aa=A(eps)
 chk(f'operator symmetric {eps}',np.max(abs(aa-aa.T))<1e-13)
 chk(f'operator centered {eps}',np.max(abs(aa@np.ones(12)))<1e-12)
 chk(f'operator rank 7 {eps}',np.count_nonzero(np.linalg.eigvalsh(aa)>1e-7)==7)
 chk(f'operator parity eigen {eps}',np.max(abs(aa@sig-float(L[6])*sig))<1e-11)
chk('operator frozen diagonal',np.max(abs(np.diag(A(.01))-np.diag(A(0))))<1e-13)
chk('positive distinct gammas',0<gamma1<gamma0<1)
chk('feature centered',abs(np.mean(cosf))<1e-15)
chk('feature var1',abs(np.mean(cosf**2)-1)<1e-15)
# Exact rational projector checks P3-P4
p3=[[Fraction((1,0,-1,0)[(i-j)%4],6) for j in range(12)] for i in range(12)]
p4=[[Fraction((2,-1,-1)[(i-j)%3],12) for j in range(12)] for i in range(12)]
for p,name in ((p3,'P3'),(p4,'P4')):
 chk(f'{name} idempotent',all(sum(p[i][k]*p[k][j] for k in range(12))==p[i][j] for i in range(12) for j in range(12)))
chk('projectors orthogonal',all(sum(p3[i][k]*p4[k][j] for k in range(12))==0 for i in range(12) for j in range(12)))
chk('rational equal diag',all(p3[i][i]==p4[i][i]==Fraction(1,6) for i in range(12)))
# Generator exact (b,o) rate summation; balanced N120, various fixed blocks
for k in (12,60):
 n=np.full(12,10,dtype=int);b=np.zeros(12,dtype=int)
 for s in range(k):b[s%12]+=1
 f,r,d,B,third=actual(n,b,A(0))
 chk(f'balanced F drift {k}',abs(d[0])<2e-12)
 chk(f'balanced R drift {k}',abs(d[1]+r)<1e-10)
 chk(f'balanced FF approx2 {k}',abs(B[0,0]-2)<.2)
 chk(f'balanced RR approx2(1-alpha) {k}',abs(B[1,1]-2*(1-k/120))<.2)
 chk(f'balanced FR zero {k}',abs(B[0,1])<1e-12)
 chk(f'third finite positive {k}',math.isfinite(third) and third>0)
# nonuniform n,b tests global drift and private exact correction
n=np.array([13,9,12,7,11,10,9,10,9,11,10,9])
assert n.sum()==120
b=np.array([1,0,2,1,1,1,0,1,1,1,2,1]);assert b.sum()==12 and np.all(b<=n)
f,r,d,B,third=actual(n,b,A(0))
chk('nonuniform R controlled deviation',abs(d[1]+r)<.1)
chk('nonuniform generator PSD bracket',np.linalg.eigvalsh(B).min()>=-1e-14)
chk('nonuniform finite FF RR FR',all(math.isfinite(x) for x in B.flatten()))
# exact hypergeometric conditional moments from a population with n balanced
# E_b R=0, Var R=(N-k)/(N-1)Var_p(f)
for k in (12,60):
 N=120;alpha=k/N;n=np.full(12,10);p=n/N
 target=(N-k)/(N-1)*(p@cosf**2-(p@cosf)**2)
 chk(f'exact finite hypergeometric variance positive {k}',target>0)
 chk(f'R variance formula balanced {k}',abs(target-(N-k)/(N-1))<1e-12)
# stationary gaussian and corrector Poisson identity, using truncation and exact OU formula
for gamma in (float(gamma0),float(gamma1)):
 for alpha in (.01,.5):
  var=alpha/gamma+1-alpha;mu=math.exp(-var/2)
  def sem(t,ff,rr):
   a=math.exp(-gamma*t);b=math.exp(-t)
   vt=alpha/gamma*(1-a*a)+(1-alpha)*(1-b*b)
   return math.exp(-vt/2)*math.cos(math.sqrt(alpha)*a*ff+b*rr)
  psis=lambda ff,rr: sum(sem(t,ff,rr)-mu for t in range(240))
  f0,r0=1.2,-.35
  tail=abs(sem(240,f0,r0)-mu)/max(1e-16,1-math.exp(-gamma))
  chk(f'corrector series remainder {gamma:.4f} {alpha}',tail<5e-15)
  lhs=psis(f0,r0) - sum((sem(t+1,f0,r0)-mu for t in range(240)))
  chk(f'OU discrete Poisson telescope {gamma:.4f} {alpha}',abs(lhs-(math.cos(math.sqrt(alpha)*f0+r0)-mu))<3e-14)
  chk(f'corrector stationary mean {gamma:.4f} {alpha}',0<mu<1)
# model contrast and asymptotic cost
linear=mp.exp(-mp.mpf('.5'))*(1/gamma1-1/gamma0)/2
chk('strict positive slope of mean gap',linear>0)
for N in (120000,1200000,12000000):
 k=math.floor(N**(15/16));M=math.floor(N**.25);alpha=k/N
 chk(f'frozen block eventually <=N/2 N{N}',k<=N//2)
 mu0=math.exp(-(1-alpha+alpha/float(gamma0))/2)
 mu1=math.exp(-(1-alpha+alpha/float(gamma1))/2)
 chk(f'frozen mean control lower N{N}',mu0>mu1)
 chk(f'frozen window powered N{N}',M>1)
output={'status':'PASS','checks':checks,'n_checks':len(checks),'scope':'Exact rational P3,P4 and finite count-generator identities; floating spectral, Poisson and finite fixtures; not MT/UF/KO source proof or machine proof of NOC asymptotics','spectral':{'Lambda3':mp.nstr(L[3],50),'Lambda6':mp.nstr(L[6],50),'gamma_A7':mp.nstr(gamma0,50),'gamma_control':mp.nstr(gamma1,50),'mu_gap_linear_coefficient':mp.nstr(linear,50)},'example_N':120000,'example_k':math.floor(120000**(15/16)),'example_M':math.floor(120000**.25)}
(ROOT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(output,indent=2)+'\n')
print(f"PASS {len(checks)}/{len(checks)} checks; mu gap linear coefficient {float(linear):.14f}")
