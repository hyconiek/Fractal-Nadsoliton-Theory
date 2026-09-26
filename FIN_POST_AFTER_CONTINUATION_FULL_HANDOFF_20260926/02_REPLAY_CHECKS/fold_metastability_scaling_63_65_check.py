#!/usr/bin/env python3
import math, numpy as np

# Certified fold intervals.
glo=3.5156447068395917
ghi=3.5156447268395917
alo=-0.212571640084013874
ahi=-0.212571626020459008
blo=0.118982527450498775
bhi=0.118982667103223664

def C(a,b):
    return (2*b/3)*((-2*a/b)**1.5)

vals=[C(a,b) for a in (alo,ahi) for b in (blo,bhi)]
Clo,Chi=min(vals),max(vals)
assert 0.5357594<Clo<Chi<0.5357599

# Noise coefficient from the fold null equation Fv=v/g and ||v||=1.
siglo=math.sqrt(2*glo)
sighi=math.sqrt(2*ghi)
assert 2.6516578<siglo<sighi<2.6516580

# Scaling algebra: with delta=N^-2/3 Delta, y=N^-1/3 Y,
# t=N^1/3 tau, each transformed contribution is N-independent.
for N in (100,1000,10000):
    Delta=1.7
    Y=.8
    # original drift magnitude coefficient excluding dt
    delta=N**(-2/3)*Delta
    y=N**(-1/3)*Y
    a=.5*(alo+ahi); b=.5*(blo+bhi); g=.5*(glo+ghi)

    drift=-g*(a*delta+.5*b*y*y)
    # after multiplying dy equation by N^1/3 and dt=N^1/3 d tau
    scaled_drift=N**(2/3)*drift
    target=-g*(a*Delta+.5*b*Y*Y)
    assert abs(scaled_drift-target)<1e-12

    # noise: N^(1/3)*sqrt(2g/N)*sqrt(N^(1/3)) = sqrt(2g)
    scaled_noise=N**(1/3)*math.sqrt(2*g/N)*N**(1/6)
    assert abs(scaled_noise-math.sqrt(2*g))<1e-12

# Barrier variable in critical scaling is N-independent.
for N in (100,1000,10000):
    Delta=2.2
    delta=N**(-2/3)*Delta
    B=N*.5*(Clo+Chi)*delta**1.5
    target=.5*(Clo+Chi)*Delta**1.5
    assert abs(B-target)<1e-12

print("PASS")
print("barrier_coefficient_interval",[Clo,Chi])
print("soft_noise_amplitude_interval",[siglo,sighi])
print("critical_parameter_exponent",-2/3)
print("critical_amplitude_exponent",-1/3)
print("critical_time_exponent",1/3)
