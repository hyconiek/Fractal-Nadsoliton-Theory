#!/usr/bin/env python3
import math
import numpy as np
from scipy.optimize import brentq

N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
lam,Q=np.linalg.eigh(A)
src=0
idx=[j for j in range(N) if j!=src]
a=A[idx,src]
A2=A@A
A3=A2@A
b=A2[idx,src]
c=A3[idx,src]
V=float(np.sum(a*a))
T=float(np.sum(a*b))
r0=a*a/V

qW=-a*b/(6*V)+(a*a)*T/(6*V*V)
B=2*qW/math.sqrt(V)

u4=b*b/4-a*c/3
QU=float(np.sum(u4))
qU=u4/V-(a*a)*QU/(V*V)
KU=qU/V

def probs_unitary(t):
    U=(Q*np.exp(-1j*lam*t))@Q.T
    return np.abs(U[:,src])**2

def probs_wave(t):
    C=(Q*np.cos(np.sqrt(np.clip(lam,0,None))*t))@Q.T
    return np.abs(C[:,src])**2

def SD(p):
    off=p[idx]
    S=float(off.sum())
    r=off/S
    return S,float(np.linalg.norm(r-r0)),r

def grid(kind, ts):
    rows=[]
    for t in ts:
        rows.append(SD(probs_unitary(t) if kind=="U" else probs_wave(t))[:2])
    rows=np.asarray(rows)
    gamma=np.gradient(np.log(rows[:,1]),np.log(rows[:,0]))
    return rows,gamma

gU,gamU=grid("U",np.logspace(-3,0,1000))
gW,gamW=grid("W",np.logspace(-2,0,1000))

assert np.ptp(b/a) > 8.0
assert np.linalg.norm(B) > 0.08
assert np.linalg.norm(KU) > 0.10

print("PASS")
print("spectrum",np.linalg.eigvalsh(A))
print("V",repr(V))
print("B_l2",repr(float(np.linalg.norm(B))))
print("KU_l2",repr(float(np.linalg.norm(KU))))
print("ratio_spread",repr(float(np.ptp(b/a))))
print("heat_vs_squared_profile_l2",
      repr(float(np.linalg.norm(W[idx,src]/W[idx,src].sum()-r0))))
