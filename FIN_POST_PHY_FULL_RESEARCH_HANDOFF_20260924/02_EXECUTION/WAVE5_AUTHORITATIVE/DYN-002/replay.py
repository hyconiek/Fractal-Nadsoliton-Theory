#!/usr/bin/env python3
import math,json
beta=1.7; M0=.8; x=.31; k=1.2
# E=.5*k*x^2; f=x^3+0.2*x
def E(x): return .5*k*x*x
def f(x): return x**3+.2*x
def fp(x): return 3*x*x+.2
def fpp(x): return 6*x
def target(M): return -M*(k*x)*fp(x)+(M/beta)*fpp(x)
rows=[]
for eps in [0.1,0.05,0.02,0.01,0.005]:
 qp=M0/(beta*eps*eps)*math.exp(-beta*(E(x+eps)-E(x))/2)
 qm=M0/(beta*eps*eps)*math.exp(-beta*(E(x-eps)-E(x))/2)
 Lex=qp*(f(x+eps)-f(x))+qm*(f(x-eps)-f(x))
 rows.append({'eps':eps,'L_exact':Lex,'target':target(M0),'abs_error':abs(Lex-target(M0))})
assert rows[-1]['abs_error']<rows[0]['abs_error']/10
# Rate rescaling control.
assert abs(target(3*M0)-3*target(M0))<1e-14
gap=1.028201984333153e-5
tau=1/(M0*gap)
print(json.dumps({'generator_convergence':rows,'rate_rescaling_ratio':target(3*M0)/target(M0),'tau_perp_bound_for_M0':tau},indent=2))
