import math, numpy as np
N=12
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/(1+min(abs(i-j),N-abs(i-j))**1.8) for j in range(N)] for i in range(N)])
A=np.diag(W.sum(1))-W
L=np.fft.fft(A[0]).real[:7]
l5=float(L[5]); g=3.7183448981203875
coef=l5*l5*(g*l5-6)/144
g0=6/l5
gamma=1-g*l5/12
print('lambda5',repr(l5))
print('g_eq',repr(g))
print('gamma5',repr(gamma))
print('coefficient',repr(coef))
print('sign_change_g',repr(g0))
assert abs(coef-0.093453912275861)<5e-10
