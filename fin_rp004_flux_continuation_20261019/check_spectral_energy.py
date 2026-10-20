#!/usr/bin/env python3
# Bounded independent diagnostic: continuum oscillator energy from BOTH q roots.
import math
from scipy.integrate import quad
O=math.sqrt(6.0)
q0=math.sqrt(2.0)
def M2(q):
    if abs(q)<1e-8: return 0.0
    F4=math.pi*q*(q*q+4)/(6*math.sinh(math.pi*q/2)) if q<30 else 0.0
    return 81*q*q*(q*q+4)*F4*F4/(1024*(1+q*q))
Gamma=M2(q0)/(2*q0)
def Econt(t):
    def integr(q):
        nu=math.sqrt(q*q+4)
        d=nu*nu-O*O
        if abs(d)<1e-8:
            b= -t*math.sin(O*t)/(2*O)
            bp= -(math.sin(O*t)+O*t*math.cos(O*t))/(2*O)
        else:
            b=-(math.cos(O*t)-math.cos(nu*t))/d
            bp=(O*math.sin(O*t)-nu*math.sin(nu*t))/d
        return M2(q)*(bp*bp+nu*nu*b*b)/(2*math.pi)
    parts=[(0,q0-0.1),(q0-0.1,q0),(q0,q0+0.1),(q0+0.1,20)]
    return sum(quad(integr,a,b,epsabs=1e-6,limit=200)[0] for a,b in parts)
print('Gamma',Gamma,'expected global E/t limit',O*Gamma/2)
for t in [40,80,160]:
    val=Econt(t)/t
    print('t',t,'Econt/t',val,'relative-to-limit',val/(O*Gamma/2))
