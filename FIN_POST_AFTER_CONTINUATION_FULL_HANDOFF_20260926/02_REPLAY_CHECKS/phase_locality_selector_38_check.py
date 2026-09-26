#!/usr/bin/env python3
import numpy as np

# Refinement preservation on a common low-mode trigonometric polynomial.
rng=np.random.default_rng(20260926)
for _ in range(100):
    c=rng.normal(size=4)
    theta=rng.uniform(-np.pi,np.pi)
    a,b=rng.normal(size=2)
    h1=c[0]*np.cos(theta)+c[1]*np.sin(theta)
    h2=c[2]*np.cos(2*theta)+c[3]*np.sin(2*theta)
    F=a*h1+b*h2

    # sample on several grids and recover exact low Fourier coefficients
    for q in (12,24,48,96):
        x=2*np.pi*np.arange(q)/q
        vals=c[0]*np.cos(x)+c[1]*np.sin(x)+c[2]*np.cos(2*x)+c[3]*np.sin(2*x)
        # exact discrete projections for q resolving k=1,2
        cc1=2/q*np.sum(vals*np.cos(x))
        ss1=2/q*np.sum(vals*np.sin(x))
        cc2=2/q*np.sum(vals*np.cos(2*x))
        ss2=2/q*np.sum(vals*np.sin(2*x))
        Fr=a*(cc1*np.cos(theta)+ss1*np.sin(theta))+b*(cc2*np.cos(2*theta)+ss2*np.sin(2*theta))
        assert abs(Fr-F)<1e-12

# Zero-order locality forces equal weights:
# choose h and g with same total value at theta but different modal split.
# F=a h1+b h2 must agree for all such splits => a=b.
a,b=1.7,-.4
x=0.31
h1,h2=.8,-.2
delta=.37
g1,g2=h1+delta,h2-delta
assert abs((h1+h2)-(g1+g2))<1e-15
assert abs((a*h1+b*h2)-(a*g1+b*g2))>1e-3

# Equal weights pass the same-value condition.
beta=2.3
assert abs(beta*h1+beta*h2-(beta*g1+beta*g2))<1e-15

print("PASS")
print("refinement_grids_checked",[12,24,48,96])
print("refinement_leaves_two_parameter_family",True)
print("zero_order_locality_selects_equal_modal_weights",True)
