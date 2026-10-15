#!/usr/bin/env python3
"""Independent bounded checks for DKO: exact Fourier projectors, actual q12
LOO-HB jump characteristics, limiting Gaussian sign event and resource KL.
Not a machine proof of the large-N/MT stochastic process theorem.
"""
from __future__ import annotations
import json,math
from pathlib import Path
from fractions import Fraction as Q
import numpy as np
import mpmath as mp

ROOT=Path(__file__).resolve().parent
mp.mp.dps=70
W=[mp.mpf(0)]+[mp.cos(mp.mpf('0.18575')*d+mp.mpf('0.16250'))/(1+mp.power(d,mp.mpf(9)/5)) for d in range(1,7)]

def lam(k):
 return 2*sum((W[d]*(1-mp.cos(2*mp.pi*k*d/12)) for d in range(1,6)),mp.mpf(0))+W[6]*(1-mp.cos(mp.pi*k))
L={r:lam(r) for r in (3,4,5,6)}
Q12=12
u=np.ones(12)/12
features={r:np.array([math.sqrt(2)*math.cos(math.tau*r*j/12) for j in range(12)]) for r in (1,3,4)}
sig=np.array([(-1)**j for j in range(12)])

def proj(k):
 if k==3:cf=lambda d:[Q(1),Q(0),Q(-1),Q(0)][d%4]
 else:cf=lambda d:[Q(1),Q(-1,2),Q(-1,2)][d%3]
 return [[cf(i-j)/6 for j in range(12)] for i in range(12)]
P3,P4=proj(3),proj(4)
def matmul(M1,M2):return [[sum(M1[i][k]*M2[k][j] for k in range(12)) for j in range(12)] for i in range(12)]
zero=[[Q(0) for j in range(12)] for i in range(12)]

def operator(eps):
 t=np.arange(12);p={k:np.array([[math.cos(math.tau*k*(i-j)/12)/6 for j in range(12)] for i in range(12)]) for k in (3,4,5)}
 p6=np.outer(sig,sig)/12
 return float(L[3]+eps)*p[3]+float(L[4]-eps)*p[4]+float(L[5])*p[5]+float(L[6])*p6

def softmax(x):
 xx=np.exp(x-np.max(x));return xx/xx.sum()

def rates(n,A):
 N=n.sum();g=12/float(L[6]); Qs=[]
 for i in range(12):
  if n[i]==0:continue
  bg=n.copy();bg[i]-=1
  Qs.append((i,float(n[i]),softmax((g/N)*(A@bg))))
 return Qs

def actual_characteristics(n,block,A,rs=(1,3,4)):
 # compressed block count and outside count. Actual finite jump compensators,
 # no empirical simulation and no assumption of independent refreshes.
 N=int(n.sum());k=int(block.sum());a=k/N
 cB=(N/k-1); cO=-1.0
 Fdr=np.zeros(len(rs));Udr=np.zeros(len(rs));Fbr=np.zeros((len(rs),len(rs))); Ubr=Fbr.copy(); cross=Fbr.copy()
 for typ,typblock in [(block,cB),(n-block,cO)]:
  for i in range(12):
   if typ[i]==0:continue
   ng=n.copy();ng[i]-=1
   q=softmax((12/float(L[6])/N)*(A@ng))
   d=np.stack([[features[r][j]-features[r][i] for r in rs] for j in range(12)])
   md=q@d; bd=np.einsum('j,ji,jk->ik',q,d,d)
   count=float(typ[i]); Fdr+=count*md/math.sqrt(N);Udr+=count*typblock*md/math.sqrt(N)
   Fbr+=count*bd/N;Ubr+=count*typblock**2*bd/N;cross+=count*typblock*bd/N
 return {'Fdrift':Fdr,'Udrift':Udr,'FF':Fbr,'UU':Ubr,'FU':cross}

def corr(gamma,alpha,h):
 c=(1-alpha)/alpha;return (math.exp(-gamma*h)/gamma+c*math.exp(-h))/(1/gamma+c)

def agree(gamma,alpha=.5,h=1):
 return .5+math.asin(corr(gamma,alpha,h))/math.pi

def gaussian_K(gamma,alpha,M,h=1):
 ii=np.arange(M); d=np.abs(ii[:,None]-ii[None,:]);return (alpha/gamma)*np.exp(-gamma*h*d)+(1-alpha)*np.exp(-h*d)

def kl(A,B):
 c=np.linalg.cholesky(B); x=np.linalg.solve(c,A); c2=np.linalg.solve(c,x.T).T
 eig=np.linalg.eigvalsh((c2+c2.T)/2)
 return .5*sum(float(x-1-math.log(x)) for x in eig)

def main():
 tests=[]
 def check(name,predicate, data=None):
  if not predicate:raise AssertionError(name)
  tests.append(name)
 def close(a,b,tol=1e-10):return abs(a-b)<tol
 check('P3 rational idempotent',matmul(P3,P3)==P3)
 check('P4 rational idempotent',matmul(P4,P4)==P4)
 check('P3P4 rational orthogonality',matmul(P3,P4)==zero)
 check('P4P3 rational orthogonality',matmul(P4,P3)==zero)
 check('equal diagonal projector trace',all(P3[i][i]==P4[i][i]==Q(1,6) for i in range(12)))
 check('no projector parity response',all(sum(P3[i][j]*Q(int(sig[j])) for j in range(12))==0 and sum(P4[i][j]*Q(int(sig[j])) for j in range(12))==0 for i in range(12)))
 check('no projector uniform response',all(sum(P3[i])==0 and sum(P4[i])==0 for i in range(12)))
 A=operator(0);Ac=operator(.01)
 check('A epsilon parity eigen',np.max(np.abs((Ac-A)@sig))<1e-12)
 check('A epsilon rank7',sum(np.linalg.eigvalsh(Ac)>1e-6)==7)
 check('A epsilon same diagonal',np.max(np.abs(np.diag(Ac)-np.diag(A)))<1e-12)
 check('A epsilon centered',np.max(np.abs(Ac@np.ones(12)))<1e-12)
 check('f feature orthogonal uniform',max(abs(np.mean(features[r])) for r in features)<1e-15)
 check('f normalized covariance',np.max(np.abs(np.stack(list(features.values()))@np.stack(list(features.values())).T/12-np.eye(3)))<1e-14)
 gam={r:1-float(L[r])/float(L[6]) for r in (3,4,5)}
 gamma1=1.;g3=gam[3];gc3=1-float(L[3]+mp.mpf('.01'))/float(L[6]);g4=gam[4];gc4=1-float(L[4]-mp.mpf('.01'))/float(L[6]);
 check('strict transverse gaps',all(gam[r]>0 for r in gam) and gc3>0 and gc4>0)
 check('control decreases gamma3',gc3<g3)
 check('control increases gamma4',gc4>g4)
 # Exact finite-generator signatures for balanced N=120; initial full counts and fixed half-block labelled assignments.
 n=np.array([10]*12);bn=np.array([5]*12)
 d=actual_characteristics(n,bn,A)
 check('Fdr uniform symmetry',np.max(abs(d['Fdrift']))<1e-12)
 check('Udr uniform symmetry',np.max(abs(d['Udrift']))<1e-12)
 check('exact FU balanced block',np.max(abs(d['FU']))<1e-12)
 check('exact UU=FF for alpha half at balanced condition',np.max(abs(d['FF']-d['UU']))<1e-12)
 check('FF bracket near 2I',np.max(abs(d['FF']-2*np.eye(3)))<0.1)
 # inspect small perturbation to check actual drift formula for U, O(N^-1/2) finite correction.
 n2=n.copy();bn2=bn.copy();n2[1]+=3;n2[7]-=3;bn2[1]+=2;bn2[7]-=2
 d2=actual_characteristics(n2,bn2,A)
 # U has features opposite +/- at 1 and 7 for r3; identity residual bounded.
 Fu=np.array([sum(features[r]*n2/math.sqrt(120)) for r in (1,3,4)])
 Vu=np.array([sum(features[r]*bn2*math.sqrt(120)/60) for r in (1,3,4)])
 Uu=Vu-Fu
 check('private drift close to -U actual LOO generator',np.max(abs(d2['Udrift']+Uu))<0.05)
 check('private bracket near 2I nonuniform',np.max(abs(d2['UU']-2*np.eye(3)))<0.3)
 pA=agree(g3);pC=agree(gc3)
 check('gamma3 time sign model contrast',pC>pA)
 check('gamma4 response reversed',agree(gc4)<agree(g4))
 check('r1 clock mode',close(corr(1,.5,1),math.exp(-1)))
 check('positive bounded signature',0<pC-pA<.02)
 check('amplitude-units drop from sign formula',close(agree(g3),.5+math.asin(corr(g3,.5,1))/math.pi))
 # sanity checks resource KL and alpha^2 scaling inside Gaussian family, not finite FIN
 K0=gaussian_K(g3,.02,25);K1=gaussian_K(gc3,.02,25);delta=kl(K0,K1)
 K02=gaussian_K(g3,.01,25);K12=gaussian_K(gc3,.01,25);d2kl=kl(K02,K12)
 check('Gaussian KL strictly positive',delta>0 and d2kl>0)
 check('Gaussian KL approx alpha squared',2.5<delta/d2kl<5.5)
 check('Gaussian covariance posdef',np.linalg.eigvalsh(K0).min()>0)
 check('normalized null mode K(alpha)',np.max(abs(gaussian_K(1,.5,5)-np.exp(-np.abs(np.arange(5)[:,None]-np.arange(5)[None,:]))))<1e-12)
 # check a valid second threshold no physical significance
 sample_cost=[{'alpha':a,'M':m,'KL0_1':kl(gaussian_K(g3,a,m),gaussian_K(gc3,a,m))} for a,m in [(.02,12),(.02,24),(.01,12),(.01,24)]]
 data={'verification_type':'bounded analytic-algebraic and floating q12 generator checks','n_checks':len(tests),'tests':tests,
       'spectral':{f'Lambda{r}':mp.nstr(L[r],45) for r in L}, 'gamma':{'3_A7':g3,'3_control':gc3,'4_A7':g4,'4_control':gc4,'1':1},
       'frozen_alpha':.5,'h':1,'snapshot_event':{'p_A7':pA,'p_control':pC,'p_difference':pC-pA,'equal_prior_risk':(1-(pC-pA))/2},
       'sign_modes':{str(r):{'p_A7':agree(1 if r==1 else gam[r]),'p_control':agree(1 if r==1 else (gc3 if r==3 else gc4))} for r in (1,3,4)},
       'balanced_N120':{'FF':d['FF'].tolist(),'UU':d['UU'].tolist(),'FU':d['FU'].tolist()},
       'gaussian_resource_fixtures':sample_cost,'scope_warning':'Does not mechanically prove MT transfer, martingale FCLT, iid parent, or shrinking-alpha finite-N result.'}
 (ROOT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n')
 print('PASS',len(tests),'checks');print('sign p',pA,pC,'diff',pC-pA)
 print('KL resources',sample_cost)

if __name__=='__main__':main()
