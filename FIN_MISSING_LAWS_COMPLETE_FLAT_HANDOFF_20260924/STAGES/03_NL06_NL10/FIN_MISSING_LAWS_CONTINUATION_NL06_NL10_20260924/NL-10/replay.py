#!/usr/bin/env python3
import numpy as np, math
from scipy.optimize import root
N=12
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/(1+min(abs(i-j),N-abs(i-j))**1.8) for j in range(N)] for i in range(N)])
A=np.diag(W.sum(1))-W;L=np.fft.fft(A[0]).real[:7];j=np.arange(12)
C=np.column_stack([math.sqrt(L[3]/6)*np.cos(2*np.pi*3*j/12),math.sqrt(L[4]/6)*np.cos(2*np.pi*4*j/12),math.sqrt(L[5]/6)*np.cos(2*np.pi*5*j/12),math.sqrt(L[6]/12)*((-1.)**j)])
g=3.7183448981203875
sl=np.array([1.8199035812800827,1.913989554668724,1.914569132546848,1.367203280195504])
def P(s):
 h=C@s;m=h.max();return float(s@s/(2*g)-(m+math.log(np.mean(np.exp(h-m)))))
def G(s):
 h=C@s;w=np.exp(h-h.max());p=w/w.sum();return s/g-p@C
def H(s):
 h=C@s;w=np.exp(h-h.max());p=w/w.sum();mu=p@C;Y=C-mu
 return np.eye(4)/g-Y.T@(p[:,None]*Y)
ts=np.linspace(0,1,10001);vv=np.array([P(t*sl) for t in ts]);tm=ts[vv.argmax()]
r=root(G,tm*sl,jac=H)
print("uniform/localized",P(np.zeros(4)),P(sl))
print("line_barrier",tm,vv.max())
print("saddle",r.x.tolist(),P(r.x),np.linalg.eigvalsh(H(r.x)).tolist(),np.linalg.norm(G(r.x)))
assert abs(P(sl))<1e-10
assert vv[1:-1].min()>0
assert sum(np.linalg.eigvalsh(H(r.x))<0)==1
