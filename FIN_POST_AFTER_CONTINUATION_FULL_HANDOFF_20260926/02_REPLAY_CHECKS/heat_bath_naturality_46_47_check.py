#!/usr/bin/env python3
import math, numpy as np
from scipy.special import softmax

N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
L=np.fft.fft(A[0]).real[:7]
j=np.arange(N)

cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/N),
        np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/N)
    ]
cols += [np.sqrt(L[6]/12)*(-1.)**j]
X=np.column_stack(cols)

s=np.array([
  1.8199035812800826819822480178121342,
  1.9139895546687240691680856665430389,
  1.9145691325468480852004194017926659,
  1.3672032801955039010722866317201310
])
theta=np.array([s[0],0,s[1],0,s[2],0,s[3]])
p=softmax(X@theta)
P=np.diag(p)

mu=p@X
Xc=X-np.ones((N,1))*mu
F=Xc.T@P@Xc

PV=Xc@np.linalg.inv(F)@Xc.T@P
PC=np.eye(N)-np.ones((N,1))*p[None,:]
PH=PC-PV

assert np.linalg.norm(PV@PV-PV)<1e-12
assert np.linalg.norm(PH@PH-PH)<1e-12
assert np.linalg.norm(PV@PH)<1e-12

# Heat-bath/reset at c=1.
Q1=-PC
assert np.max(np.abs(Q1-(np.ones((N,1))*p[None,:]-np.eye(N))))<1e-12

# Exact positivity interval for Qc=-PH-cPV.
lo=-np.inf
hi=np.inf
for i in range(N):
    for jj in range(N):
        if i==jj: continue
        v=PV[i,jj]
        pj=p[jj]
        if v>1e-14:
            hi=min(hi,1+pj/v)
        elif v<-1e-14:
            lo=max(lo,1+pj/v)

assert lo < 1 < hi
assert hi-lo > 0.04

c=1.01
Q=-PH-c*PV
off=Q.copy()
np.fill_diagonal(off,np.nan)
assert np.nanmin(off)>0
assert np.linalg.norm(Q@np.ones(N))<1e-12
assert np.linalg.norm(P@Q-Q.T@P)<1e-12

# Spectral rates in p-self-adjoint similarity representation.
SqrtP=np.diag(np.sqrt(p))
InvSqrtP=np.diag(1/np.sqrt(p))
Qs=SqrtP@Q@InvSqrtP
e=np.linalg.eigvalsh(Qs)
# 7 retained rates -c, 4 hidden rates -1, one zero.
assert np.sum(np.isclose(e,-c,atol=1e-10))==7
assert np.sum(np.isclose(e,-1,atol=1e-10))==4
assert np.sum(np.isclose(e,0,atol=1e-10))==1

# D12 covariance for translated states/features.
def perm(a,eps=1):
    R=np.zeros((N,N))
    for x in range(N):
        R[x,(eps*x+a)%N]=1
    return R

maxcov=0.0
for a in range(N):
    R=perm(a,1)
    pp=R@p
    XX=R@X
    PP=np.diag(pp)
    mm=pp@XX
    XXc=XX-np.ones((N,1))*mm
    FF=XXc.T@PP@XXc
    PV2=XXc@np.linalg.inv(FF)@XXc.T@PP
    PC2=np.eye(N)-np.ones((N,1))*pp[None,:]
    PH2=PC2-PV2
    Q2=-PH2-c*PV2
    maxcov=max(maxcov,np.linalg.norm(Q2-R@Q@R.T))
assert maxcov<1e-11

print("PASS")
print("valid_c_interval",[float(lo),float(hi)])
print("test_c",c)
print("minimum_offdiag_rate",float(np.nanmin(off)))
print("reversibility_error",float(np.linalg.norm(P@Q-Q.T@P)))
print("D12_covariance_error",maxcov)
print("rates_hidden",-1.0)
print("rates_retained",-c)
print("global_clock_equivalent",False)
