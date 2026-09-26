#!/usr/bin/env python3
import numpy as np, math
from scipy.optimize import root
from scipy.special import softmax

N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
lam=np.fft.fft(A[0]).real[:7]
j=np.arange(N)

cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/N),
        np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/N)
    ]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols)
l3,l6=lam[3],lam[6]
a3=np.sqrt(l3/6)
a6=np.sqrt(l6/12)

# R7P-095 midpoint.
r=.5*(0.41421132290+0.41421132293)
K=a6*r
gp=12*K/(l6*np.tanh(K))

th=np.zeros(7); th[6]=r
p=softmax(X@th)
S=np.diag(p)-np.outer(p,p)
H=np.eye(7)/gp-X.T@S@X

assert abs(H[0,0])<5e-12
other=np.delete(np.linalg.eigvalsh(H),0)
assert np.linalg.eigvalsh(H)[1]>1e-3

# Cubic zero and reduced quartic.
Xc=X-np.ones((N,1))*(p@X)
u=Xc[:,0]
D3=-np.sum(p*u**3)
assert abs(D3)<1e-12

stable=list(range(1,7))
Hs=H[np.ix_(stable,stable)]
t=-np.sum(p[:,None]*(u[:,None]**2)*Xc,axis=0)[stable]
m2=np.sum(p*u*u)
D4=-(np.sum(p*u**4)-3*m2*m2)
D4red=D4-3*t@np.linalg.solve(Hs,t)
beta=D4red/24
assert beta<0

# Exact derivative formulas.
tt=np.tanh(K); sech2=1/np.cosh(K)**2
dgdK=12/l6*(tt-K*sech2)/(tt*tt)
dldK=(l6*(K*sech2-tt)/(K*K)-l3*sech2)/12
ell=dldK/dgdK
assert ell<0

A_theta=np.sqrt(abs(ell)/(-4*beta))
A_J=a3*A_theta

# Compare with MP7-039 point.
gD3=5.171841831942818
JD3=.038434753978898215
Jpred=A_J*np.sqrt(gp-gD3)
assert abs(Jpred/JD3-1)<0.02

# Two-harmonic stationary equations and radial fold.
def xy(J,K):
    d=np.cosh(J)+np.exp(-2*K)
    return np.sinh(J)/d,(np.cosh(J)-np.exp(-2*K))/d

def derivs(J,K):
    e=np.exp(-2*K); ch=np.cosh(J); sh=np.sinh(J); d=ch+e
    x=sh/d; y=(ch-e)/d
    xJ=(1+e*ch)/d**2
    xK=2*e*sh/d**2
    yJ=xK
    yK=4*e*ch/d**2
    return x,y,xJ,xK,yJ,yK

def fold_eq(v):
    J,K,g=v
    x,y,xJ,xK,yJ,yK=derivs(J,K)
    F1=J-g*l3*x/6
    F2=K-g*l6*y/12
    M=np.array([[1-g*l3*xJ/6,-g*l3*xK/6],
                [-g*l6*yJ/12,1-g*l6*yK/12]])
    return [F1,F2,np.linalg.det(M)]

fr=root(fold_eq,[.7185,.4909,4.6212],tol=1e-13)
assert np.linalg.norm(fold_eq(fr.x),np.inf)<1e-12

print("PASS")
print("pure_k6_uniform_bifurcation_gain",float(12/l6))
print("pure_k6_k3_pitchfork_gain",float(gp))
print("reduced_quartic_beta",float(beta))
print("soft_slope_dlambda_dg",float(ell))
print("subcritical_J_coefficient",float(A_J))
print("MP7_039_J_prediction",float(Jpred))
print("two_harmonic_fold",fr.x.tolist())
