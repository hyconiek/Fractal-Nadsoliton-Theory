#!/usr/bin/env python3
import numpy as np

# Moving-frame deterministic system for k=1,2 at constant omega.
omega=0.37
A=np.array([
    [-1, omega, 0, 0],
    [-omega,-1,0,0],
    [0,0,-1,2*omega],
    [0,0,-2*omega,-1]
],float)

eig=np.linalg.eigvals(A)
target=np.array([-1+1j*omega,-1-1j*omega,-1+2j*omega,-1-2j*omega])
# unordered spectral check
for z in target:
    assert np.min(np.abs(eig-z))<1e-12

# Closure obstruction: same u=u1+u2, different v1/v2 -> different du.
x1=np.array([.4,.0,.6,.0])    # u=1
x2=np.array([.4,1.0,.6,.0])   # u=1, different v1
c=np.array([1,0,1,0.])
assert abs(c@x1-c@x2)<1e-15
du1=c@(A@x1)
du2=c@(A@x2)
assert abs(du1-du2)>1e-3

# Same static stiffness, different transfer functions at s=0.
K=2.3
eta=.8
M=.7
b=np.array([.4,.2])
lam=np.array([1.1,3.2])

def chiR(s): return 1/(K+eta*s)
def chiI(s): return 1/(K+M*s*s)
def chiM(s): return 1/(K+eta*s+np.sum(b*s/(s+lam)))

assert abs(chiR(0)-1/K)<1e-15
assert abs(chiI(0)-1/K)<1e-15
assert abs(chiM(0)-1/K)<1e-15

# But dynamic responses differ.
s=1.7
vals=np.array([chiR(s),chiI(s),chiM(s)])
assert np.ptp(vals)>1e-2

# Adiabatic steady quadrature for each k: v=-k omega u, hence
# du=-(1+k^2 omega^2)u.
for k in (1,2):
    u=.7
    v=-k*omega*u
    du=-u+k*omega*v
    assert abs(du+(1+k*k*omega*omega)*u)<1e-12

print("PASS")
print("moving_frame_eigenvalues",eig.tolist())
print("same_u_derivative_split",float(du2-du1))
print("static_susceptibility",1/K)
print("dynamic_transfer_values_at_s=1.7",vals.tolist())
