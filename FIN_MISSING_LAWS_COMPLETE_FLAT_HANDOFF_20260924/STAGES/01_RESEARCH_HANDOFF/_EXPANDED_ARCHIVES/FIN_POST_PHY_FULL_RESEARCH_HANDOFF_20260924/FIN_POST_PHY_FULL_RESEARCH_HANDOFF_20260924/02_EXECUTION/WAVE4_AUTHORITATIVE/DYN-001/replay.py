#!/usr/bin/env python3
import math, json
mu=0.4; c=1.3; V=0.2; L=2*math.pi
# h(l)=(e^(mu l)-1)/l is increasing for l>0, so a degree-2 cycle gives this bound.
Gamma=2*c*(math.exp(mu*L)-1)/L
C=2*mu*V+Gamma
# One illustrative tail evaluation; not a physical calibration.
t=0.5; R=math.pi
p_bound=math.exp(-mu*R+C*t)
# Near-crossing regularization: raw intensity diverges, metric action stays finite.
v=0.7; T=0.1
rows=[]
for eps in [1e-2,1e-4,1e-6,1e-8]:
    # exact integral of c/sqrt((v t)^2+eps^2), symmetric about 0
    raw=(2*c/v)*math.asinh(v*T/eps)
    # numerical midpoint integral of lambda*(exp(mu ell)-1)
    n=200000
    dt=2*T/n
    s=0.0
    for k in range(n):
        x=-T+(k+0.5)*dt
        ell=math.sqrt((v*x)**2+eps**2)
        lam=c/ell
        s += lam*(math.exp(mu*ell)-1)*dt
    rows.append({'eps':eps,'raw_integrated_rate':raw,'metric_action_integral':s})
assert all(rows[i+1]['raw_integrated_rate']>rows[i]['raw_integrated_rate'] for i in range(len(rows)-1))
assert max(r['metric_action_integral'] for r in rows)<1.0
print(json.dumps({'Gamma_mu_bound':Gamma,'cone_exponent_rate':C,'illustrative_tail_bound':p_bound,'near_collision':rows},indent=2))
