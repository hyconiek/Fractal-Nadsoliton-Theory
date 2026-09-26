#!/usr/bin/env python3
import numpy as np

omega=.37

A=np.array([
    [-1, omega, 0, 0],
    [-omega,-1,0,0],
    [0,0,-1,2*omega],
    [0,0,-2*omega,-1]
],float)
c=np.array([1.,0.,1.,0.])

# Exact spectrum.
eig=np.linalg.eigvals(A)
target=np.array([
    -1+1j*omega,-1-1j*omega,
    -1+2j*omega,-1-2j*omega
])
for z in target:
    assert np.min(np.abs(eig-z))<1e-12

# Observability rank of scalar point-evaluation output.
O=np.vstack([c, c@A, c@A@A, c@A@A@A])
rank=np.linalg.matrix_rank(O,tol=1e-12)
assert rank==4

# Equal-input transfer numerator/denominator has no pole cancellation.
# H(s) = (s+1)/D1 + (s+1)/D2.
# Verify at generic complex points against state-space resolvent.
b=c.copy()
def H_state(s):
    return c @ np.linalg.solve(s*np.eye(4)-A,b)
def H_formula(s):
    D1=(s+1)**2+omega**2
    D2=(s+1)**2+4*omega**2
    return (s+1)/D1+(s+1)/D2

for s in (0.2+0.1j,1.3,2.1+0.7j):
    assert abs(H_state(s)-H_formula(s))<1e-12

# Exact variable-omega Volterra reconstruction against direct ODE integration
# for one k sector using a fine deterministic discretization.
k=2
T=1.2
n=40001
t=np.linspace(0,T,n)
dt=t[1]-t[0]
om=.22+.11*np.sin(.8*t)
u=np.zeros(n); v=np.zeros(n)
u[0]=.7; v[0]=-.31

# RK4 direct integration.
def f(tt,xx):
    oo=.22+.11*np.sin(.8*tt)
    uu,vv=xx
    return np.array([-uu+k*oo*vv,-vv-k*oo*uu])
for q in range(n-1):
    tt=t[q]; x=np.array([u[q],v[q]])
    h=dt
    k1=f(tt,x)
    k2=f(tt+h/2,x+h*k1/2)
    k3=f(tt+h/2,x+h*k2/2)
    k4=f(tt+h,x+h*k3)
    xn=x+h*(k1+2*k2+2*k3+k4)/6
    u[q+1],v[q+1]=xn

# Check exact eliminated relation at sampled points using trapezoidal convolution.
errs=[]
for q in np.linspace(1000,n-1,25,dtype=int):
    tq=t[q]
    s=t[:q+1]
    kernel=(k*k)*om[q]*om[:q+1]*np.exp(-(tq-s))
    integ=np.trapz(kernel*u[:q+1],s)
    rhs=-u[q]+k*om[q]*np.exp(-tq)*v[0]-integ
    direct=-u[q]+k*om[q]*v[q]
    errs.append(abs(rhs-direct))
assert max(errs)<2e-7

print("PASS")
print("constant_omega_eigenvalues",eig.tolist())
print("scalar_observability_rank",rank)
print("variable_omega_volterra_max_error",max(errs))
print("frame_kernel_strength_ratio_k2_to_k1",4.0)
