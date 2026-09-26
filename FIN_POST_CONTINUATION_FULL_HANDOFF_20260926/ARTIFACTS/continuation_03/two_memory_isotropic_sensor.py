#!/usr/bin/env python3
import math, numpy as np
from scipy.optimize import brentq
Q=12
th=2*np.pi*np.arange(Q)/Q
# strict kernel
W=np.array([[0.0 if i==j else math.cos(.18575*min(abs(i-j),Q-abs(i-j))+.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)])
A=np.diag(W.sum(1))-W
lam=np.fft.fft(A[0]).real[:7]
amp={3:math.sqrt(lam[3]/6),4:math.sqrt(lam[4]/6),5:math.sqrt(lam[5]/6),6:math.sqrt(lam[6]/12)}
# x=b/a solves equality of hidden k1 and k2 response amplitudes for f=a c4+b c5.
poly=lambda x: 8*x**3-20*x**2+15*x-18
x=brentq(poly,2,3)
t=x*amp[4]/amp[5]
H=np.column_stack([np.cos(th),np.sin(th),np.cos(2*th),np.sin(2*th)])
R=np.zeros((Q,4))
for r in range(Q):
    d=2*np.pi*r/Q
    f=amp[4]*np.cos(4*(th+d))+t*amp[5]*np.cos(5*(th+d))
    for q in range(4):
        R[r,q]=(1/Q)*sum(H[i,q]*(f[j]-f[i])**4 for i in range(Q) for j in range(Q))
sv=np.linalg.svd(R,compute_uv=False)
a=amp[4];b=t*amp[5]
S1=3*a*b*(9*a*a+10*b*b)
S2=1.5*b*b*(15*a*a+8*b*b)
print('lambda4',repr(lam[4]))
print('lambda5',repr(lam[5]))
print('x=b_over_a',repr(x))
print('t_feature_ratio',repr(t))
print('S1',repr(S1),'S2',repr(S2),'difference',repr(S1-S2))
print('singular_values',*[repr(float(v)) for v in sv])
print('condition_number',repr(float(sv[0]/sv[-1])))
print('max_DFT_model_residual',repr(float(np.max(np.abs(R[:,0]-S1*np.cos(2*np.pi*np.arange(Q)/Q))))))
print('rank',np.linalg.matrix_rank(R,tol=1e-10))
assert abs(S1-S2)<1e-11
assert sv[-1]>0 and sv[0]/sv[-1]-1<1e-12
assert np.linalg.matrix_rank(R,tol=1e-10)==4
