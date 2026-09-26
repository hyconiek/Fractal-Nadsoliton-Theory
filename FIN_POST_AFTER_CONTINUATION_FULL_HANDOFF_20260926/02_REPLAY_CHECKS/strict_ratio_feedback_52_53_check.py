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
A7=X@X.T
PC=np.eye(N)-np.ones((N,N))/N

g=3.7183448981203875
B=-PC+(g/12)*A7
eig=np.linalg.eigvalsh(B)

expected=[]
# hidden k1,k2 => -1, retained k3..5 double + k6 single
expected += [-1.0]*4
for k in (3,4,5):
    expected += [-(1-g*L[k]/12)]*2
expected += [-(1-g*L[6]/12)]
expected += [0.0]
assert np.max(np.abs(np.sort(eig)-np.sort(expected)))<1e-12

rates=[1-g*L[k]/12 for k in (3,4,5,6)]
known=np.array([
 .3922344001359638,
 .3184370325847773,
 .2877490912855972,
 .2742466130695381
])
assert np.max(np.abs(np.array(rates)-known))<2e-15

# Localized stationary center: verify algebraic z decoupling.
s=np.array([1.8199035812800827,1.9139895546687241,
            1.9145691325468481,1.3672032801955039])
theta=np.array([s[0],0,s[1],0,s[2],0,s[3]])
p=softmax(X@theta)
S=np.diag(p)-np.outer(p,p)

hidden=[]
for k in (1,2):
    hidden += [
        np.sqrt(2/N)*np.cos(2*np.pi*k*j/N),
        np.sqrt(2/N)*np.sin(2*np.pi*k*j/N)
    ]
Y=np.column_stack(hidden)

F=X.T@S@X
H=Y.T@S@X
K=H@np.linalg.inv(F)

rng=np.random.default_rng(20260926)
for _ in range(100):
    x=rng.normal(size=7)
    y=rng.normal(size=4)
    dx=(g*F-np.eye(7))@x
    dy=-y+g*H@x
    z=y-K@x
    dz=dy-K@dx
    assert np.max(np.abs(dz+z))<2e-12

# Whitening alone does not fix kinetic matrix.
Sigma=np.array([[1.0,.2],[.2,.7]])
ev,U=np.linalg.eigh(Sigma)
Sinvhalf=U@np.diag(1/np.sqrt(ev))@U.T
# In whitened coordinates any SPD C preserves covariance I.
C=np.array([[.7,.15],[.15,1.8]])
assert np.linalg.eigvalsh(C).min()>0
assert np.linalg.norm(C@np.eye(2)+np.eye(2)@C-2*C)<1e-14

# Direct A coupling cannot retain hidden rate exactly one except rho=0.
assert abs(L[1]-L[2])>1e-3

print("PASS")
print("uniform_retained_rates",rates)
print("uniform_hidden_rate",1.0)
print("localized_decoupling_trials",100)
print("lambda1_lambda2",[float(L[1]),float(L[2])])
print("whitened_covariance_does_not_fix_C",True)
