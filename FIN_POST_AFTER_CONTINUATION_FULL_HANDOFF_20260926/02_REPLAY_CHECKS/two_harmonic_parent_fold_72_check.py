#!/usr/bin/env python3
import numpy as np
from scipy.optimize import root

l3=1.96140686197644
l6=2.3421820411463

def xy(J,K):
    q=np.exp(-2*K)
    C=np.cosh(J);S=np.sinh(J);D=C+q
    return S/D,(C-q)/D

def jac2(J,K,g):
    q=np.exp(-2*K)
    C=np.cosh(J);S=np.sinh(J);D=C+q
    xJ=(1+C*q)/(D*D)
    xK=2*q*S/(D*D)
    yJ=xK
    yK=4*q*C/(D*D)
    return np.array([
        [1-g*l3*xJ/6,-g*l3*xK/6],
        [-g*l6*yJ/12,1-g*l6*yK/12]
    ])

def F(z):
    J,K,g=z
    x,y=xy(J,K)
    return np.r_[
        J-g*l3*x/6,
        K-g*l6*y/12,
        np.linalg.det(jac2(J,K,g))
    ]

sol=root(F,[.72,.49,4.62],tol=1e-13)
assert np.linalg.norm(F(sol.x))<1e-12

target=np.array([0.718505246005630,
                 0.490908854614046,
                 4.621196599489125])
assert np.linalg.norm(sol.x-target)<2e-12

print("PASS")
print("two_harmonic_fold",sol.x.tolist())
