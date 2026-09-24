#!/usr/bin/env python3
import math, json
rows=[]
gamma=2.0
I=1.7
for q in [12,18,24,37]:
    Delta=0.3
    def U(a): return Delta*(1-math.cos(q*a))
    a0=0.0; a1=2*math.pi/q
    dU=U(a1)-U(a0)
    c=dU/(gamma*I)
    assert abs(dU)<1e-12 and abs(c)<1e-12
    rows.append({'q':q,'deltaU_equal_wells':dU,'speed_from_balance':c})
# Held-out explicit bias: reversing it reverses drift sign.
q=18; F=0.01
period=2*math.pi/q
# U_bias = periodic - F*a => deltaU=-F*period
c_plus=(-F*period)/(gamma*I)
c_minus=(F*period)/(gamma*I)
assert c_plus*c_minus<0 and abs(c_plus+c_minus)<1e-15
# Show the PHA-002 all-orders upper bound decreases with q.
B=.5; rho=.2
M=math.exp(B*sum(math.cosh(k*rho) for k in range(1,6)))
def barrier(q): return 4*M*math.exp(-rho*q)/(1-math.exp(-rho*q))
bounds=[(q,barrier(q)) for q in [12,24,48,96]]
assert all(bounds[i+1][1]<bounds[i][1] for i in range(len(bounds)-1))
print(json.dumps({'equal_well_checks':rows,'bias_speed':{'F':F,'c':c_plus,'reversed_c':c_minus},'barrier_bounds':bounds},indent=2))
