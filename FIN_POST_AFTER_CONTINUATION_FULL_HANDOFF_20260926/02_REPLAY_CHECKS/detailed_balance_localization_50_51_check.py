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
  1.8199035812800827,
  1.9139895546687241,
  1.9145691325468481,
  1.3672032801955039
])
theta=np.array([s[0],0,s[1],0,s[2],0,s[3]])
h=X@theta
p=softmax(h)

gamma=1.0
eta=.2
B=np.ones((N,N))*gamma/12 + eta*W
np.fill_diagonal(B,0.0)

def make_Q(kind):
    Q=np.zeros((N,N))
    for i in range(N):
        for k in range(N):
            if i==k: continue
            x=h[k]-h[i]
            if kind=="half":
                f=np.exp(x/2)
            elif kind=="barker":
                f=2/(1+np.exp(-x))
            elif kind=="metro":
                f=min(1.,np.exp(x))
            elif kind=="smooth":
                f=np.exp(x/2+.01*x*x)
            Q[i,k]=B[i,k]*f
        Q[i,i]=-sum(Q[i,k] for k in range(N) if k!=i)
    return Q

Qs={k:make_Q(k) for k in ("half","barker","metro","smooth")}
P=np.diag(p)

for name,Q in Qs.items():
    assert np.linalg.norm(P@Q-Q.T@P)<1e-13
    assert np.max(np.abs(Q@np.ones(N)))<1e-11
    off=Q.copy()
    np.fill_diagonal(off,np.nan)
    assert np.nanmin(off)>0

SqrtP=np.diag(np.sqrt(p))
Inv=np.diag(1/np.sqrt(p))
spectra={}
for name,Q in Qs.items():
    H=SqrtP@Q@Inv
    spectra[name]=np.linalg.eigvalsh(H)

assert np.max(np.abs(spectra["half"]-spectra["barker"]))>1e-2
assert np.max(np.abs(spectra["half"]-spectra["metro"]))>1e-2

PC=np.eye(N)-np.ones((N,N))/N
Q0=-gamma*PC-eta*A
Qu=B.copy()
for i in range(N):
    Qu[i,i]=-sum(Qu[i,k] for k in range(N) if k!=i)
assert np.max(np.abs(Qu-Q0))<1e-12

for x in (1e-4,2e-4,5e-4):
    fh=np.exp(x/2)
    fb=2/(1+np.exp(-x))
    fs=np.exp(x/2+.01*x*x)
    assert abs((fh-1)/x-.5)<2e-4
    assert abs((fb-1)/x-.5)<2e-4
    assert abs((fs-1)/x-.5)<2e-4

for x,y in [(.3,-.1),(.2,.4),(-.25,-.35)]:
    f=lambda z: np.exp(z/2)
    assert abs(f(x+y)-f(x)*f(y))<1e-14

def barker(z): return 2/(1+np.exp(-z))
def metro(z): return min(1.,np.exp(z))
assert abs(barker(.3+.2)-barker(.3)*barker(.2))>1e-3
assert abs(metro(.3-.2)-metro(.3)*metro(-.2))>1e-3

print("PASS")
print("localized_p_max",float(p.max()))
print("spectral_difference_half_vs_barker",
      float(np.max(np.abs(spectra["half"]-spectra["barker"]))))
print("spectral_difference_half_vs_metropolis",
      float(np.max(np.abs(spectra["half"]-spectra["metro"]))))
print("half_density_cocycle",True)
print("barker_cocycle",False)
print("metropolis_cocycle",False)
