#!/usr/bin/env python3
import math
import numpy as np

# Numerical replay of the arbitrary-split identity.
rng=np.random.default_rng(20260926)

def effective_fraction(a,b,u0,u1,ba,bb):
    um=(b*u0+a*u1)/(a+b)
    phia=ba*(u0+um)/2
    phib=bb*(um+u1)/2
    return (a*phia+b*phib)/(a+b)

# Constant beta must pass arbitrary splits exactly.
for _ in range(1000):
    a,b=rng.uniform(.01,5,2)
    u0,u1=rng.normal(size=2)
    beta=rng.normal()
    eff=effective_fraction(a,b,u0,u1,beta,beta)
    target=beta*(u0+u1)/2
    assert abs(eff-target)<1e-12

# Nonconstant beta fails generically.
fails=0
for _ in range(100):
    a,b=rng.uniform(.1,3,2)
    u0,u1=rng.normal(size=2)
    ba=1+a
    bb=1+b
    bL=1+a+b
    eff=effective_fraction(a,b,u0,u1,ba,bb)
    target=bL*(u0+u1)/2
    if abs(eff-target)>1e-6:
        fails+=1
assert fails>90

# Strict-kernel nonnearest fraction from the repository model.
N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
row=W[0]
total=row.sum()
nearest=row[1]+row[-1]
nonnearest=(total-nearest)/total
assert abs(nonnearest-0.4338569990558586)<1e-12

print("PASS")
print("constant_beta_random_split_checks",1000)
print("nonconstant_beta_generic_failures",fails)
print("strict_non_nearest_exit_fraction",float(nonnearest))
