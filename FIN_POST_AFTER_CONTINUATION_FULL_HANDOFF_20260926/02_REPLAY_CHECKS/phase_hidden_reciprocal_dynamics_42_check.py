#!/usr/bin/env python3
import numpy as np
from scipy.integrate import solve_ivp

# One-cell invariant fixture.
I=1.3
m1=.8
m2=1.7
k1=.9
k2=1.2
g=.6

def energy(q):
    th,x1,y1,x2,y2=q
    u1=x1*np.cos(th)+y1*np.sin(th)
    u2=x2*np.cos(2*th)+y2*np.sin(2*th)
    u=u1+u2
    return .5*k1*(x1*x1+y1*y1)+.5*k2*(x2*x2+y2*y2)+.5*g*u*u

def gradE(q):
    th,x1,y1,x2,y2=q
    c1,s1=np.cos(th),np.sin(th)
    c2,s2=np.cos(2*th),np.sin(2*th)
    u1=x1*c1+y1*s1
    u2=x2*c2+y2*s2
    u=u1+u2
    v1=-x1*s1+y1*c1
    v2=-x2*s2+y2*c2
    return np.array([
        g*u*(v1+2*v2),
        k1*x1+g*u*c1,
        k1*y1+g*u*s1,
        k2*x2+g*u*c2,
        k2*y2+g*u*s2
    ])

# Check global phase invariance numerically.
rng=np.random.default_rng(20260926)
max_inv=0.
for _ in range(100):
    q=rng.normal(size=5)
    a=rng.normal()
    th,x1,y1,x2,y2=q
    z1=(x1+1j*y1)*np.exp(1j*a)
    z2=(x2+1j*y2)*np.exp(2j*a)
    qp=np.array([th+a,z1.real,z1.imag,z2.real,z2.imag])
    max_inv=max(max_inv,abs(energy(qp)-energy(q)))
assert max_inv<1e-12

# Inertial equations: state [q, qdot].
mass=np.array([I,m1,m1,m2,m2])
def rhs_inertial(t,Y):
    q=Y[:5]; v=Y[5:]
    return np.r_[v,-gradE(q)/mass]

Y0=np.array([.3,.4,-.2,.1,.5, .2,.1,.3,-.15,.05])
sol=solve_ivp(rhs_inertial,[0,15],Y0,rtol=2e-11,atol=1e-12,
              method="DOP853",dense_output=False,max_step=.02)

def total_energy(Y):
    q=Y[:5];v=Y[5:]
    return energy(q)+.5*np.sum(mass*v*v)

def noether(Y):
    th,x1,y1,x2,y2=Y[:5]
    thd,xd1,yd1,xd2,yd2=Y[5:]
    return I*thd + m1*(x1*yd1-y1*xd1) + 2*m2*(x2*yd2-y2*xd2)

Es=np.array([total_energy(sol.y[:,i]) for i in range(sol.y.shape[1])])
Qs=np.array([noether(sol.y[:,i]) for i in range(sol.y.shape[1])])
assert np.ptp(Es)<2e-9
assert np.ptp(Qs)<2e-9

# Gradient flow and exact energy monotonicity.
Gam=np.array([1.1,.7,.7,1.6,1.6])
def rhs_grad(t,q):
    return -gradE(q)/Gam

q0=Y0[:5]
sg=solve_ivp(rhs_grad,[0,12],q0,rtol=2e-11,atol=1e-12,
             method="DOP853",max_step=.02)
Eg=np.array([energy(sg.y[:,i]) for i in range(sg.y.shape[1])])
assert np.max(np.diff(Eg))<2e-11

# Moving-frame kinetic identity on random fixture.
th=.41; thd=-.27
w1=.3+.5j; wd1=-.2+.1j
w2=-.4+.2j; wd2=.13-.31j
z1=np.exp(1j*th)*w1
zd1=np.exp(1j*th)*(wd1+1j*thd*w1)
z2=np.exp(2j*th)*w2
zd2=np.exp(2j*th)*(wd2+2j*thd*w2)
lhs=.5*m1*abs(zd1)**2+.5*m2*abs(zd2)**2
rhs=.5*m1*abs(wd1+1j*thd*w1)**2+.5*m2*abs(wd2+2j*thd*w2)**2
assert abs(lhs-rhs)<1e-14

print("PASS")
print("phase_invariance_error",max_inv)
print("inertial_energy_drift",float(np.ptp(Es)))
print("noether_charge_drift",float(np.ptp(Qs)))
print("gradient_energy_drop",float(Eg[0]-Eg[-1]))
print("gradient_max_upstep",float(np.max(np.diff(Eg))))
print("independent_inertial_coefficients",3)
print("dimensionless_ratios_after_clock_quotient",2)
