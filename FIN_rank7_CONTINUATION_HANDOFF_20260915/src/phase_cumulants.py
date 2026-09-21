"""R7P-081--084: angular amplitude convention and exact cumulant/resonance formulas.

The imported handoff uses complex mode amplitudes z3,z4,z5 and real z6 with
||theta||^2 = 2 sum_{k=3}^5 |z_k|^2/lambda_k + z6^2/lambda_6.
With this convention
  theta_{k,c}+i*(-theta_{k,s}) = sqrt(2/lambda_k) z_k,
so h=X7 theta has harmonic amplitude |z_k|/sqrt(3) for k=3,4,5 and
z6/sqrt(12) for k=6.
"""
from __future__ import annotations
import cmath, itertools, math
import numpy as np
from .model import feature_spaces

N=12
ACTIVE=(3,4,5)

def theta_from_z(z3:complex,z4:complex,z5:complex,z6:float, L=None):
    if L is None: L=feature_spaces()[2]
    out=[]
    for k,z in zip(ACTIVE,(z3,z4,z5)):
        scale=math.sqrt(2.0/L[k])
        out.extend([scale*z.real, -scale*z.imag])
    out.append(float(z6)/math.sqrt(L[6]))
    return np.array(out,float)

def z_from_theta(theta,L=None):
    if L is None: L=feature_spaces()[2]
    theta=np.asarray(theta,float)
    zs=[]
    for idx,k in enumerate(ACTIVE):
        tc,ts=theta[2*idx:2*idx+2]
        zs.append(math.sqrt(L[k]/2.0)*(tc-1j*ts))
    z6=math.sqrt(L[6])*theta[6]
    return (*zs,z6)

def theta_radius_from_z(z3,z4,z5,z6,L=None):
    if L is None: L=feature_spaces()[2]
    return math.sqrt(sum(2*abs(z)**2/L[k] for k,z in zip(ACTIVE,(z3,z4,z5)))+float(z6)**2/L[6])

def field_from_z(z3,z4,z5,z6):
    j=np.arange(N)
    h=np.zeros(N,float)
    for k,z in zip(ACTIVE,(z3,z4,z5)):
        # Re[z exp(i omega j)]/sqrt(3)
        h += np.real(z*np.exp(2j*np.pi*k*j/N))/math.sqrt(3.0)
    h += float(z6)*((-1.0)**j)/math.sqrt(12.0)
    return h

def fourier_coefficients(z3,z4,z5,z6):
    """Complex Fourier coefficients c_k with h_j=sum c_k exp(2pi i k j/12)."""
    c={}
    for k,z in zip(ACTIVE,(z3,z4,z5)):
        c[k]=complex(z)/(2*math.sqrt(3.0))
        c[(-k)%N]=complex(z).conjugate()/(2*math.sqrt(3.0))
    c[6]=float(z6)/math.sqrt(12.0)
    return c

def resonance_moment(order,z3,z4,z5,z6):
    c=fourier_coefficients(z3,z4,z5,z6)
    ks=tuple(c)
    total=0j
    for tup in itertools.product(ks,repeat=order):
        if sum(tup)%N==0:
            prod=1+0j
            for k in tup: prod*=c[k]
            total+=prod
    return float(total.real)

def direct_moments(z3,z4,z5,z6):
    h=field_from_z(z3,z4,z5,z6)
    return {n:float(np.mean(h**n)) for n in (2,3,4)}

def cumulants(z3,z4,z5,z6):
    m=direct_moments(z3,z4,z5,z6)
    k2=m[2]; k3=m[3]; k4=m[4]-3*m[2]**2
    return {'kappa2':k2,'kappa3':k3,'kappa4':k4,
            'K4':k2/2+k3/6+k4/24, **{f'm{n}':v for n,v in m.items()}}

def analytic_moments(r3,r4,r5,z6,phi3,phi4,phi5):
    """Closed resonance formulas; r3,r4,r5>=0 and signed real z6."""
    a=r3/(2*math.sqrt(3.0)); b=r4/(2*math.sqrt(3.0)); c=r5/(2*math.sqrt(3.0)); d=z6/math.sqrt(12.0)
    m2=2*(a*a+b*b+c*c)+d*d
    m3=(6*a*a*d*math.cos(2*phi3)
        +12*a*b*c*math.cos(phi3+phi4+phi5)
        +2*b**3*math.cos(3*phi4))
    base=(6*(a**4+b**4+c**4)+d**4
          +24*(a*a*b*b+a*a*c*c+b*b*c*c)
          +12*d*d*(a*a+b*b+c*c))
    phase=(2*a**4*math.cos(4*phi3)
           +24*a*b*b*c*math.cos(phi3-2*phi4+phi5)
           +48*a*b*c*d*math.cos(-phi3+phi4+phi5)
           +8*a*c**3*math.cos(phi3-3*phi5)
           +24*b*c*c*d*math.cos(phi4-2*phi5))
    return {2:m2,3:m3,4:base+phase}

def cubic_locks(sign:int):
    if sign not in (-1,1): raise ValueError('sign must be +/-1')
    p3s=(0.0,math.pi) if sign==1 else (math.pi/2,3*math.pi/2)
    out=[]
    for p3 in p3s:
        for m in range(3):
            p4=2*math.pi*m/3
            p5=(-p3-p4)%(2*math.pi)
            out.append((p3%(2*math.pi),p4%(2*math.pi),p5))
    return out

def quartic_lock_cosines(sign:int, phases):
    p3,p4,p5=phases
    return np.array([
        math.cos(4*p3),
        math.cos(p3-2*p4+p5),
        sign*math.cos(-p3+p4+p5),
        math.cos(p3-3*p5),
        sign*math.cos(p4-2*p5),
    ])

def cubic_lock_cosines(sign:int, phases):
    p3,p4,p5=phases
    return np.array([sign*math.cos(2*p3),math.cos(p3+p4+p5),math.cos(3*p4)])

def _phase_term_data(r3,r4,r5,z6):
    a=r3/(2*math.sqrt(3.0)); b=r4/(2*math.sqrt(3.0)); c=r5/(2*math.sqrt(3.0)); d=z6/math.sqrt(12.0)
    cubic=[
        (6*a*a*d, np.array([2.,0.,0.])),
        (12*a*b*c, np.array([1.,1.,1.])),
        (2*b**3, np.array([0.,3.,0.])),
    ]
    quartic=[
        (2*a**4,np.array([4.,0.,0.])),
        (24*a*b*b*c,np.array([1.,-2.,1.])),
        (48*a*b*c*d,np.array([-1.,1.,1.])),
        (8*a*c**3,np.array([1.,0.,-3.])),
        (24*b*c*c*d,np.array([0.,1.,-2.])),
    ]
    m2=2*(a*a+b*b+c*c)+d*d
    base4=(6*(a**4+b**4+c**4)+d**4
          +24*(a*a*b*b+a*a*c*c+b*b*c*c)
          +12*d*d*(a*a+b*b+c*c))
    return m2,base4,cubic,quartic

def k4_phase_value_grad_hess(r3,r4,r5,z6,phases):
    ph=np.asarray(phases,float)
    m2,base4,cubic,quartic=_phase_term_data(r3,r4,r5,z6)
    m3=0.; m4=base4; grad=np.zeros(3); H=np.zeros((3,3))
    for coef,n in cubic:
        ang=float(n@ph); co=math.cos(ang); si=math.sin(ang)
        m3+=coef*co
        grad += (-coef/6.0)*n*si
        H += (-coef/6.0)*np.outer(n,n)*co
    for coef,n in quartic:
        ang=float(n@ph); co=math.cos(ang); si=math.sin(ang)
        m4+=coef*co
        grad += (-coef/24.0)*n*si
        H += (-coef/24.0)*np.outer(n,n)*co
    k4=m4-3*m2*m2
    K4=m2/2+m3/6+k4/24
    return float(K4),grad,H

def full_phase_value_grad_hess(r3,r4,r5,z6,phases):
    ph=np.asarray(phases,float); amps=np.array([r3,r4,r5],float); ks=np.array([3.,4.,5.])
    j=np.arange(N,dtype=float)
    angles=2*np.pi*j[:,None]*ks[None,:]/N + ph[None,:]
    h=(np.cos(angles)*amps[None,:]).sum(axis=1)/math.sqrt(3.0) + z6*((-1.0)**j)/math.sqrt(12.0)
    hp=-np.sin(angles)*amps[None,:]/math.sqrt(3.0)
    hpp=-np.cos(angles)*amps[None,:]/math.sqrt(3.0)
    m=float(np.max(h)); ew=np.exp(h-m); p=ew/ew.sum()
    K=m+math.log(float(ew.mean()))
    grad=p@hp
    centered=hp-grad
    H=centered.T@(p[:,None]*centered)
    H += np.diag(p@hpp)
    return float(K),grad,H
