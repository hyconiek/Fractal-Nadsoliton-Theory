#!/usr/bin/env python3
import math,itertools,numpy as np
from scipy.special import softmax, logsumexp
Q=12;g=3.7183448981203875;jj=np.arange(Q)
W=np.array([[0.0 if i==j else math.cos(.18575*min(abs(i-j),Q-abs(i-j))+.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)])
L=np.diag(W.sum(1))-W;lam=np.fft.fft(L[0]).real[:7]
cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
X=np.column_stack(cols);A=X@X.T

def comps(N,q,prefix=()):
    if q==1: yield prefix+(N,);return
    for a in range(N+1):yield from comps(N-a,q-1,prefix+(a,))

def check(N):
    states=list(comps(N,Q)); idx={s:i for i,s in enumerate(states)};M=len(states)
    Gen=np.zeros((M,M)); logpi=np.empty(M)
    for a,n in enumerate(states):
        nv=np.array(n,float);p=nv/N;field=g*A@p; q=softmax(field)
        logZ=logsumexp(field)
        logmult=math.lgamma(N+1)-sum(math.lgamma(v+1) for v in n)
        logpi[a]=logmult+N*g*.5*(p@A@p)+logZ
        for i in range(Q):
          if n[i]==0: continue
          for j in range(Q):
            if i==j:continue
            nn=list(n);nn[i]-=1;nn[j]+=1;b=idx[tuple(nn)]
            Gen[a,b]+=n[i]*q[j]
        Gen[a,a]=-Gen[a].sum()
    pi=np.exp(logpi-logsumexp(logpi))
    stat=np.max(np.abs(pi@Gen)); db=0.
    for a in range(M):
      for b in range(a+1,M):
        if Gen[a,b] or Gen[b,a]: db=max(db,abs(pi[a]*Gen[a,b]-pi[b]*Gen[b,a]))
    print('N',N,'states',M,'stationarity_max',stat,'detailed_balance_max',db,'rowsum',np.max(abs(Gen.sum(1))))
    assert stat<3e-14 and db<3e-14
for N in (1,2,3):check(N)

# exact algebraic single-edge log-ratio check at random interior p-like count states
N=24;n=np.ones(Q,dtype=int)*2
for i,j in [(0,1),(0,5),(3,8)]:
 p=n/N; field=g*A@p; q=softmax(field)
 nn=n.copy();nn[i]-=1;nn[j]+=1;pp=nn/N;fieldp=g*A@pp;qp=softmax(fieldp)
 lhs=math.log(n[i]*q[j]/(nn[j]*qp[i]))
 def logweight(nv):
   p=nv/N; f=g*A@p
   return math.lgamma(N+1)-sum(math.lgamma(int(v)+1) for v in nv)+N*g*.5*(p@A@p)+logsumexp(f)
 rhs=logweight(nn)-logweight(n)
 print('edge',i,j,'rate_log_ratio',lhs,'weight_log_ratio',rhs,'resid',lhs-rhs)
 assert abs(lhs-rhs)<2e-14
