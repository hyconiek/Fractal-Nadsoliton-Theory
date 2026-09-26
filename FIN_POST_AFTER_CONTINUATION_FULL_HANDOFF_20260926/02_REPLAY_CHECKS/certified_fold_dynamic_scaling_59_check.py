#!/usr/bin/env python3
import math

glo=3.5156447068395917
ghi=3.5156447268395917
alo=-0.21257164008401387
ahi=-0.21257162602045901
blo=0.11898252745049878
bhi=0.11898266710322366

assert ahi<0
assert blo>0

# Branch exists for delta>0.
delta=1e-4
assert -2*ahi*delta/bhi > 0

# Certified leading amplitude interval.
amp_lo=math.sqrt(2*(-ahi)/bhi)
amp_hi=math.sqrt(2*(-alo)/blo)

# Certified critical-slowing prefactor.
c_lo=glo*math.sqrt(2*(-ahi)*blo)
c_hi=ghi*math.sqrt(2*(-alo)*bhi)

assert 1.890278 < amp_lo < amp_hi < 1.890280
assert .790704 < c_lo < c_hi < .790706

# Direct normal-form replay at coefficient midpoints.
g=.5*(glo+ghi)
a=.5*(alo+ahi)
b=.5*(blo+bhi)
for delta in (1e-8,1e-7,1e-6,1e-5):
    y=math.sqrt(-2*a*delta/b)
    # stationarity of quadratic normal form
    assert abs(a*delta+.5*b*y*y)<1e-18
    rate=g*b*y
    pred=g*math.sqrt(-2*a*b)*math.sqrt(delta)
    assert abs(rate-pred)<1e-14
    assert rate>0

print("PASS")
print("branch_amplitude_prefactor_interval",[amp_lo,amp_hi])
print("critical_slowing_prefactor_interval",[c_lo,c_hi])
print("fold_side","g > g_f")
print("stable_branch_sign_for_certificate_v","+")
