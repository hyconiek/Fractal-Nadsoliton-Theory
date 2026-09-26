#!/usr/bin/env python3
import math,numpy as np
from scipy.special import softmax
Q=12;g=3.7183448981203875
j=np.arange(Q)
W=np.array([[0.0 if a==b else math.cos(.18575*min(abs(a-b),Q-abs(a-b))+.1625)/(1+min(abs(a-b),Q-abs(a-b))**1.8) for b in range(Q)] for a in range(Q)])
L=np.diag(W.sum(1))-W; lam=np.fft.fft(L[0]).real[:7]
cols=[]
for k in (3,4,5): cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols);A=X@X.T

def q_naive(p): return softmax(g*A@p)
def q_loo(p,i): return softmax(g*A@p-(g/len_scale)*A[:,i])

def tri_affinity(N,i=0,jj=1,k=2,loo=False):
    global len_scale;len_scale=N
    n=np.ones(Q,dtype=int)*(N//Q); p=n/N
    Astate=p.copy();B=Astate.copy();B[i]-=1/N;B[jj]+=1/N;C=Astate.copy();C[i]-=1/N;C[k]+=1/N
    def target(p,dep,targ):
      return (q_loo(p,dep) if loo else q_naive(p))[targ]
    f=target(Astate,i,jj)*target(B,jj,k)*target(C,k,i)
    r=target(Astate,i,k)*target(C,k,jj)*target(B,jj,i)
    return math.log(f/r)

for N in (12,24,48,96,192,384):
    an=tri_affinity(N,loo=False); al=tri_affinity(N,loo=True)
    exact=g/N*(A[0,1]-A[0,0])
    print(N,'naive',repr(an),'N*aff',repr(N*an),'formula',repr(exact),'loo',repr(al))
    assert abs(an-exact)<2e-14
    assert abs(al)<2e-14
print('A00',A[0,0],'A01',A[0,1],'g*(A01-A00)',g*(A[0,1]-A[0,0]))
