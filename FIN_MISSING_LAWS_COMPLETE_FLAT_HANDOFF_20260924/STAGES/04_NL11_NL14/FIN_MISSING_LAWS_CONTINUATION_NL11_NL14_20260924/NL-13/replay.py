#!/usr/bin/env python3
# See results.json for the L=15/20/25 campaign. This replay checks the
# coexistence endpoints and the central index-one barrier, while the full BVP
# code is recorded in REPORT/results provenance.
import numpy as np, math
from scipy.optimize import root
N=12
W=np.array([[0 if i==j else math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/(1+min(abs(i-j),N-abs(i-j))**1.8) for j in range(N)] for i in range(N)],float)
A=np.diag(W.sum(1))-W;L=np.fft.fft(A[0]).real[:7];j=np.arange(12)
C=np.column_stack([math.sqrt(L[3]/6)*np.cos(2*np.pi*3*j/12),math.sqrt(L[4]/6)*np.cos(2*np.pi*4*j/12),math.sqrt(L[5]/6)*np.cos(2*np.pi*5*j/12),math.sqrt(L[6]/12)*((-1.)**j)])
g=3.7183448981203875;s=np.array([1.8199035812800827,1.913989554668724,1.914569132546848,1.367203280195504])
def P(x):
 h=C@x;m=h.max();return float(x@x/(2*g)-(m+math.log(np.mean(np.exp(h-m)))))
def G(x):
 h=C@x;w=np.exp(h-h.max());p=w/w.sum();return x/g-p@C
print(P(np.zeros(4)),P(s),np.linalg.norm(G(s)))
assert abs(P(s))<1e-10 and np.linalg.norm(G(s))<1e-8
