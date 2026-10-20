#!/usr/bin/env python3
"""Independent small numerical checks. Proofs are in THEOREM.md."""
import json, math
import numpy as np
from scipy.linalg import expm
from scipy.integrate import quad
K=np.array([0, .46998567264502017, .19204355169010282, .09142861427792497, .04702916874565040, .02413122336363006, .011070817321442113])
n=12
W=np.zeros((n,n))
for i in range(n):
 for j in range(n):
  d=min((i-j)%n,(j-i)%n)
  if i!=j: W[i,j]=K[d]
s=W.sum(axis=1)
assert np.max(np.abs(s-s[0]))<1e-13
A=np.diag(s)-W
q=-A
assert np.linalg.norm(A-A.T)<1e-13
assert np.linalg.norm(A@np.ones(n))<1e-13
pi=np.ones(n)/n
alpha=.35
p0=(1-alpha)*pi
p0[0]+=alpha
T=1.
def prob(t):return expm(t*q)@p0
def inf(t):
 p=prob(t)
 v=q@p
 return float(np.sum(v*v/p))
S_heat=0.5*quad(inf,0,T,epsabs=5e-13,epsrel=5e-13,limit=200)[0]
pT=prob(T)
B=float(np.dot(np.sqrt(p0),np.sqrt(pT)))
theta=math.acos(max(-1,min(1,B)))
S_min=2*theta*theta/T
D=[]
for m in (8,16,32,64,128,256):
 dt=T/m
 H=expm(dt*q)
 p=p0.copy();val=0
 for step in range(m):
  next_p=H@p
  val+=float(np.sum(p*np.log(p/next_p)))/dt
  p=next_p
 D.append([m,val,abs(val-S_heat)])
assert D[-1][2]<D[0][2] and D[-1][2]<D[-2][2]
assert S_heat>S_min>0
assert abs(sum(p0)-1)<1e-14
assert np.min(pT)>0
# Nonconstant-speed analytic obstacle illustrated by Fisher speeds at t=0, .5, 1.
I_samples=[[t,inf(t)] for t in [0,.25,.5,.75,1,2,5]]
assert I_samples[0][1]>I_samples[-1][1]
# Static constraints alone do not select a unique 12-state clock: use a distance-dependent
# perturbation keeping row escape s and invariant uniform; NOT an adopted FIN law.
a=.1
kk=K.copy();kk[1]*=1+a
kk[2]-=K[1]/K[2]*a*K[2]  # simple algebra: w2 -> w2 - a*w1
assert kk[2]>0
WW=W.copy()
for i in range(n):
 for j in range(n):
  if i!=j:
   d=min((i-j)%n,(j-i)%n)
   WW[i,j]=kk[d]
assert np.max(abs(WW.sum(axis=1)-s))<1e-13
Qalt=WW-np.diag(s)
assert np.linalg.norm(Qalt-Qalt.T)<1e-13
assert np.linalg.norm(Qalt.T@pi)<1e-13
assert np.linalg.norm(Qalt-q)>0
# Exact path KL-rate at stationarity, as source Q relative to alternative Qalt.
rate=float(sum(pi[i]*sum(W[i,j]*math.log(W[i,j]/WW[i,j])-W[i,j]+WW[i,j]
 for j in range(n) if j!=i) for i in range(n)))
assert rate>0
# Fisher--Rao entropy gradient: all tangent eigenvalues equal to -1, while strict
# heat uses -lambda_k, which has nonconstant positive spectrum.
lams=np.linalg.eigvalsh(A)
assert abs(lams[0])<1e-13 and lams[-1]-lams[1]>1.5
u=pi
x=np.zeros(n);x[0]=1; x[1]=-1

def fisher_entropy_gradient(p):
 lograt=np.log(p/pi)
 return p*(lograt-np.dot(p,lograt))
for h in (1e-5,1e-6):
 fd=(fisher_entropy_gradient(pi+h*x)-fisher_entropy_gradient(pi-h*x))/(2*h)
 assert np.linalg.norm(fd-x)<1e-8
# Graph logarithmic-mean mobility reconstructs strict heat after supplying
# that precise mobility. This is representation, not independent origin.
p=prob(.3)
lograt=np.log(p/pi)
mob_drift=np.zeros(n)
for i in range(n):
 for j in range(n):
  if i==j: continue
  lmean=(p[i]-p[j])/(math.log(p[i])-math.log(p[j])) if abs(p[i]-p[j])>1e-14 else p[i]
  mob_drift[i]-=W[i,j]*lmean*(lograt[i]-lograt[j])
assert np.linalg.norm(mob_drift-q@p)<1e-12
# New generic finite clock-normalized nonuniqueness corroborates previous FIN source obstructions.
output={'strict_kernel_row_sum':float(s[0]),'n':n,'p0_alpha':alpha,'T':T,'S_heat':S_heat,'S_Fisher_geodesic_min':S_min,'strict_gap':float(lams[1]),'strict_laplacian_eigenvalues':[float(v) for v in lams],'Fisher_entropy_linearization_nonzero_spectrum':[1.]*11,'ratio_heat_to_min':S_heat/S_min,'Hellinger_B':B,'angle_radians':theta,'discrete_KL_checks':D,'fisher_speed_squared':I_samples,'comparison_activity_multiplier_d1':1+a,'comparison_activity_multiplier_d2':float(kk[2]/K[2]),'comparison_same_mean_escape_rate':float(s[0]),'comparison_path_KL_rate':rate,'graph_mobility_reconstruction_error':float(np.linalg.norm(mob_drift-q@p)),'assessment':'bounded float64 diagnostics, not formal proof or physical test'}
with open('/mnt/data/fin_intrinsic_kl_history_research_20261019/verification.json','w') as f: json.dump(output,f,indent=2,ensure_ascii=False)
print(json.dumps(output,indent=2))
