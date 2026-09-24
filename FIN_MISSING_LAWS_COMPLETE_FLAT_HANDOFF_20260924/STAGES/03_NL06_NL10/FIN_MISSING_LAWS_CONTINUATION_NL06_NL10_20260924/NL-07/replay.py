#!/usr/bin/env python3
import numpy as np, math
from scipy.linalg import eigh
N=12
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/(1+min(abs(i-j),N-abs(i-j))**1.8) for j in range(N)] for i in range(N)])
A=np.diag(W.sum(1))-W; L=np.fft.fft(A[0]).real[:7]
c=np.array([math.sqrt(L[k]/6) for k in (3,4,5)])
s=np.array([3.5475952773959873,3.783140506071962,3.677956492814376]);a=c*s
M=65536;t=2*np.pi*np.arange(M)/M;k=np.array([3.,4.,5.]);B=t[:,None]*k
h=np.sum(a*np.cos(B),1);w=np.exp(h-h.max());p=w/w.sum()
D=-np.sin(B)*a;md=p@D;F=(D*p[:,None]).T@D-np.outer(md,md)
H=np.diag(a*(p@np.cos(B)))-F
J=np.array([[3.,0.,1.],[4.,0.,1.],[5.,1.,1.]])
Fy=J.T@F@J; Hy=J.T@H@J
Fq=Fy[1:,1:]-np.outer(Fy[0,1:],Fy[0,1:])/Fy[0,0]
Hs=Hy[1:,1:]
print("Fq",Fq.tolist())
print("Fq_eigs",np.linalg.eigvalsh(Fq).tolist())
print("omega2_if_Fisher",eigh(Hs,Fq,eigvals_only=True).tolist())
assert np.min(np.linalg.eigvalsh(Fq))>0
