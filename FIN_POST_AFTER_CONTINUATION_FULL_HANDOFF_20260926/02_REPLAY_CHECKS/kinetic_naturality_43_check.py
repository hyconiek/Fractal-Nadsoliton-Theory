#!/usr/bin/env python3
import numpy as np

# Arbitrary-split additive storage replay.
rng=np.random.default_rng(20260926)
mu=np.array([.8,1.3,2.1])

for _ in range(1000):
    a,b=rng.uniform(.01,5,2)
    lhs=mu*(a+b)
    rhs=mu*a+mu*b
    assert np.max(np.abs(lhs-rhs))<1e-14

# Independent ratios survive a common scale quotient.
mu_theta,mu1,mu2=mu
ratios=np.array([mu1/mu_theta,mu2/mu_theta])

scaled=7.3*mu
ratios_scaled=np.array([scaled[1]/scaled[0],scaled[2]/scaled[0]])
assert np.max(np.abs(ratios-ratios_scaled))<1e-14

# Same static curvatures but different relative masses -> different frequency ratios.
K=np.array([1.2,.9,1.7])
w=np.sqrt(K/mu)
mu_alt=np.array([mu_theta,2*mu1,.5*mu2])
w_alt=np.sqrt(K/mu_alt)
assert abs((w[1]/w[0])-(w_alt[1]/w_alt[0]))>1e-3
assert abs((w[2]/w[0])-(w_alt[2]/w_alt[0]))>1e-3

# D12 invariant block metric fixture:
# phase scalar + k1 plane + k2 plane
G=np.diag([mu_theta,mu1,mu1,mu2,mu2])

def R2(t):
    return np.array([[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]])

alpha=2*np.pi/12
Rep=np.block([
    [np.ones((1,1)),np.zeros((1,2)),np.zeros((1,2))],
    [np.zeros((2,1)),R2(alpha),np.zeros((2,2))],
    [np.zeros((2,1)),np.zeros((2,2)),R2(2*alpha)]
])
assert np.linalg.norm(Rep.T@G@Rep-G)<1e-12

print("PASS")
print("additive_split_checks",1000)
print("surviving_dimensionless_ratios",ratios.tolist())
print("frequency_ratios_original",(w/w[0]).tolist())
print("frequency_ratios_changed_metric",(w_alt/w_alt[0]).tolist())
print("invariant_metric_coefficients",3)
print("ratios_after_clock_quotient",2)
