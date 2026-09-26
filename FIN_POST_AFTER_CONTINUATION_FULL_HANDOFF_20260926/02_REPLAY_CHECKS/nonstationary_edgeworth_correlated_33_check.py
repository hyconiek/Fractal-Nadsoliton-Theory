#!/usr/bin/env python3
import math
import numpy as np

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

# Retained normalized strict feature columns:
cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(L[k]/6)*np.cos(2*np.pi*k*j/N),
        np.sqrt(L[k]/6)*np.sin(2*np.pi*k*j/N)
    ]
cols += [np.sqrt(L[6]/12)*(-1.)**j]
X=np.column_stack(cols)

# Hidden Euclidean-orthonormal Fourier basis k=1,2.
hidden=[]
for k in (1,2):
    hidden += [
        np.sqrt(2/N)*np.cos(2*np.pi*k*j/N),
        np.sqrt(2/N)*np.sin(2*np.pi*k*j/N)
    ]
Y=np.column_stack(hidden)

T=[X.T@np.diag(Y[:,a])@X for a in range(4)]
G=np.array([[np.trace(T[a]@T[b]) for b in range(4)] for a in range(4)])
Geig=np.linalg.eigvalsh(G)

assert Geig.min()>2.0
assert np.linalg.cond(G)<1.21

# Exact isolating entries from Fourier multiplication identities.
lam3,lam4,lam6=L[3],L[4],L[6]
pred=[
    math.sqrt(lam3*lam4)*math.sqrt(6)/12,
    math.sqrt(lam3*lam4)*math.sqrt(6)/12,
    math.sqrt(lam4*lam6)*math.sqrt(3)/6,
   -math.sqrt(lam4*lam6)*math.sqrt(3)/6
]
obs=[
    T[0][0,2],
    T[1][0,3],
    T[2][2,6],
    T[3][3,6]
]
assert np.max(np.abs(np.array(pred)-np.array(obs)))<1e-12

# Dual detector matrices Q_b.
Ginv=np.linalg.inv(G)
Q=[sum(Ginv[b,c]*T[c] for c in range(4)) for b in range(4)]

dual=np.array([[np.trace(T[a]@Q[b]) for b in range(4)] for a in range(4)])
assert np.max(np.abs(dual-np.eye(4)))<1e-12

# Synthetic reconstruction of m and M.
rng=np.random.default_rng(20260926)
m=rng.normal(size=4)
M=rng.normal(size=(4,7))
delta=.07
preps=[np.zeros(7)]+[delta*np.eye(7)[k] for k in range(7)]

# By theorem L_hid f_b = m_b + M_b x.
R=np.array([[m[b]+M[b]@x for b in range(4)] for x in preps])

mhat=R[0].copy()
Mhat=np.column_stack([(R[k+1]-R[0])/delta for k in range(7)])
assert np.max(np.abs(mhat-m))<1e-13
assert np.max(np.abs(Mhat-M))<1e-13

# Direct operator check on f_b=x^T Q_b x:
# Hessian=2Q_b and 1/2 sum_a h_a T_a:Hess = h_b.
for _ in range(10):
    x=rng.normal(size=7)
    h=m+M@x
    direct=[]
    for b in range(4):
        val=.5*sum(h[a]*np.sum(T[a]*(2*Q[b])) for a in range(4))
        direct.append(val)
    assert np.max(np.abs(np.array(direct)-h))<1e-12

print("PASS")
print("Gram=")
print(G)
print("Gram_eigenvalues=",Geig.tolist())
print("Gram_condition_number=",float(np.linalg.cond(G)))
print("isolating_entries=",obs)
print("dual_error=",float(np.max(np.abs(dual-np.eye(4)))))
print("reconstruction_error_m=",float(np.max(np.abs(mhat-m))))
print("reconstruction_error_M=",float(np.max(np.abs(Mhat-M))))
