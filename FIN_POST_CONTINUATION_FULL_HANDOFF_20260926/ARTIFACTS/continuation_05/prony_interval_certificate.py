#!/usr/bin/env python3
import math

def cert(m,e):
    m0,m1,m2,m3=m; e0,e1,e2,e3=e
    D=m0*m2-m1*m1
    dD=abs(m2)*e0+abs(m0)*e2+2*abs(m1)*e1+e0*e2+e1*e1
    N1=m0*m3-m1*m2
    dN1=abs(m3)*e0+abs(m0)*e3+abs(m2)*e1+abs(m1)*e2+e0*e3+e1*e2
    N2=m1*m3-m2*m2
    dN2=abs(m3)*e1+abs(m1)*e3+2*abs(m2)*e2+e1*e3+e2*e2
    Dl,Du=D-dD,D+dD
    if Dl<=0 or N1-dN1<=0 or N2-dN2<=0:
        return D,Dl,None,None,None,False
    s1=((N1-dN1)/Du,(N1+dN1)/Dl)
    s2=((N2-dN2)/Du,(N2+dN2)/Dl)
    disc=(s1[0]**2-4*s2[1],s1[1]**2-4*s2[0])
    return D,Dl,s1,s2,disc,disc[0]>0

def moments(delta,a1=1.0,a2=1.0,lam=1.0):
    l1,l2=lam-delta,lam+delta
    return [a1*l1**n+a2*l2**n for n in range(4)]

rho=1e-6
print('relative_moment_box',rho)
print('delta\tD\tD_lower\tdisc_lower\tcertified')
for d in [.2,.1,.075,.06,.055,.05,.04,.02,.01,.005,.002,.001]:
    m=moments(d); e=[rho*max(1.0,abs(x)) for x in m]
    D,Dl,s1,s2,disc,ok=cert(m,e)
    print(f'{d:.6g}\t{D:.12g}\t{Dl:.12g}\t{float("nan") if disc is None else disc[0]:.12g}\t{ok}')

# Bisection threshold for this conservative interval certificate.
lo,hi=.001,.2
for _ in range(60):
    d=(lo+hi)/2;m=moments(d);e=[rho*max(1.0,abs(x)) for x in m]
    if cert(m,e)[-1]: hi=d
    else: lo=d
print('certification_threshold_delta',hi)
print('true_pole_separation_threshold',2*hi)
