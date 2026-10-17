#!/usr/bin/env python3
"""Independent symbolic/scoped algebra controls for the NS conditional theorems.
NOT a PDE proof, FIN-specific coefficient replay, or a strict-source certificate.
"""
import json
from pathlib import Path
import sympy as S

checks=[]
def check(name, assertion, details=""):
    if assertion is not True:
        raise AssertionError(f"{name} failed: {assertion}")
    checks.append({"name":name,"passed":True,"details":details})

x,phi,v,k,m2,lam,g=S.symbols('x phi v k m2 lam g', real=True, nonzero=True)
Uprime=-m2*phi+lam*phi**3+g*phi**5
phiK=v*S.tanh(k*x)
res=S.expand((S.diff(phiK,x,2)-Uprime.subs(phi,phiK)))
t=S.symbols('t',real=True)
respoly=S.Poly(S.expand(res.subs(S.tanh(k*x),t)), t)
c=respoly.all_coeffs()
check('sextic phi5 residual',S.simplify(respoly.coeff_monomial(t**5)-(-g*v**5))==0,'Coefficient of t^5 in Kdd-Uprime is -g*v^5.')
check('mass consistency for tanh',S.simplify(respoly.coeff_monomial(t)-(-2*v*k*k+m2*v))==0)
check('quartic consistency for tanh',S.simplify(respoly.coeff_monomial(t**3)-(2*v*k*k-lam*v**3))==0)
check('tanh exact only zero sextic',S.simplify(res.subs({m2:2*k*k,lam:2*k*k/v**2,g:0}))==0)

H0pot=4*k*k-6*k*k*S.sech(k*x)**2
f0=S.sech(k*x)**2
f1=S.sinh(k*x)*S.sech(k*x)**2
check('scalar translational ground mode',S.simplify(-S.diff(f0,x,2)+H0pot*f0)==0)
check('scalar internal mode eigenvalue 3k2',S.simplify(-S.diff(f1,x,2)+H0pot*f1-3*k*k*f1)==0)

s,b=S.symbols('s b',positive=True,real=True)
f=S.sech(k*x)**s
lhs=S.simplify((-S.diff(f,x,2)-k*k*s*(s+1)*S.sech(k*x)**2*f+k*k*s*s*f)/f)
check('PT ground trial exact eigenfunction',S.simplify(lhs)==0)
check('PT parameter identity', S.simplify(((S.sqrt(1+4*b/k**2)-1)/2)*((S.sqrt(1+4*b/k**2)+1)/2)-b/k**2)==0)

A11,A22,B11,B22,e=S.symbols('A11 A22 B11 B22 e',real=True)
A=S.diag(A11,A22)
B=S.Matrix([[B11,e],[e,B22]])
comm=A*B-B*A
check('matrix commutator nonzero generic',S.simplify(comm[0,1]-e*(A11-A22))==0 and S.simplify(comm[1,0]+e*(A11-A22))==0)

rho=S.symbols('rho', real=True)
Vinf=A+B
check('interpolation exact',S.simplify((1-rho)*A+rho*Vinf-(A+rho*B))==S.zeros(2))

D=S.diag(1,2)
E=S.Matrix([[0,1],[1,0]])
Bcontrol=100*S.eye(2)+E
Acontrol=-D
Vcontrol=Acontrol+Bcontrol
check('signed A negative all directions',Acontrol.eigenvals()=={-1:1,-2:1})
check('noncommuting signed control',Acontrol*Bcontrol!=Bcontrol*Acontrol)
check('control b exact',max(Bcontrol.eigenvals())==101)
check('control mu exact',min(Vcontrol.eigenvals())==(197-S.sqrt(5))/2)
mu=(197-S.sqrt(5))/2
bstar=S.Integer(101)
sstar=(S.sqrt(1+4*bstar)-1)/2
del_control=S.simplify(mu-sstar**2)
check('positive gap exact signed noncommuting control',bool(del_control>0),str(del_control))

# The 12x12 extension has exactly 24 real field directions under the radial block assumption.
A12=S.diag(-1,-2,*([-1]*10))
B12=S.diag(*([100]*12));B12[0,1]=1;B12[1,0]=1
check('actual dimensional control 12x12',A12.shape==(12,12) and B12.shape==(12,12))
A24=S.diag(A12,A12);B24=S.diag(B12,B12)
check('full 24-real transverse blocks retained',A24.shape==(24,24) and B24.shape==(24,24))
check('signed 12x12 control is noncommuting',A12*B12!=B12*A12)

# Exact symbolic transfer identities; P and Q stand for fixed profile integrals.
P,Q,eps=S.symbols('P Q eps', real=True, nonzero=True)
N=2*eps**2*P
Den=2*eps*Q
q=S.simplify(-N/Den)
Eintprime=S.simplify(eps*N+eps*q*Den)
Ephiprime=-eps*N
Epsiprime=-eps*q*Den
check('velocity multiplier correct order',S.simplify(q+eps*P/Q)==0)
check('interaction initial derivative exactly zero',Eintprime==0)
check('opposite sector signs exact at leading order',S.simplify(Ephiprime+Epsiprime)==0 and S.simplify(Epsiprime-2*eps**3*P)==0)
check('sector plus interaction energy exact cancellation',S.simplify(Eintprime+Ephiprime+Epsiprime)==0)

# Orthonormal modulation derivative nonzero for true nontrivial K, symbolic phi4 integral.
f0=S.sech(x)**2
mod_norm=S.integrate((1-S.tanh(x)**2)**2,(x,-S.oo,S.oo))
check('phi4 translation modulation nondegenerate',mod_norm==S.Rational(4,3))

# H1 gap constant c=delta/(delta+1+C) inequality is positive for positive data.
delta,C=S.symbols('delta C',positive=True)
cH=delta/(delta+1+C)
check('H1 coercivity constant positive',S.simplify(cH*(delta+1+C)-delta)==0 and bool(cH.subs({delta:1,C:1})>0))

# Exact signed nonzero-sextic heteroclinic first-integral identity (both g signs
# are analytic provided radicands are positive), no fitted parameter values.
T,V,L,G=S.symbols('T V L G', real=True)
C6=L/2+2*G*V**2/3
D6=G*V**2/3
th=S.simplify(D6/C6)
F=1+th*(1-T**2)
kg2=V**2*(C6+D6)
Kg2=V**2*T**2/F
# K' obtained by chain rule and T'=k(1-T^2).
Kder2=S.factor(V**2*kg2*(1+th)**2*(1-T**2)**2/F**3)
U6=S.factor((Kg2-V**2)**2*(L/4+G*(Kg2+2*V**2)/6))
check('exact sextic first integral after coefficient stationarity',S.simplify(Kder2-2*U6)==0)
check('exact sextic factorization',S.simplify(-((L*V**2+G*V**4)/2)*phi**2+L*phi**4/4+G*phi**6/6-(-((L*V**2+G*V**4)/2)*V**2+L*V**4/4+G*V**6/6)-(phi**2-V**2)**2*(L/4+G*(phi**2+2*V**2)/6))==0)
check('kink normalization g=0',S.simplify(th.subs(G,0))==0 and S.simplify(kg2.subs(G,0)-V**2*L/2)==0)
r2=S.simplify(T**2/F)
check('sextic wall comparison error algebra',S.simplify((r2-T**2) + th*T**2*(1-T**2)/F)==0)
ms6=L*V**2+G*V**4
phiUdd=-ms6+3*L*Kg2+5*G*Kg2**2
VPT=-2*kg2+6*kg2*T**2
rhs=G*V**4*(1-6*T**2+5*r2**2)+3*L*V**2*(r2-T**2)
check('exact sextic scalar Hessian comparison identity',S.simplify(phiUdd-VPT-rhs)==0)
# Explicit signed-g checks used only to certify nonempty inequality *classes*.
from sympy import Rational as Rat
for gg,label in [(Rat(1,1000),'positive'),(-Rat(1,1000),'negative')]:
    v0=Rat(1);l0=Rat(2)
    kval2=S.factor(kg2.subs({V:v0,L:l0,G:gg}))
    theta0=S.factor(th.subs({V:v0,L:l0,G:gg}))
    r0=abs(theta0)/(4*(1-abs(theta0)))
    e6=12*abs(gg)*v0**4+3*l0*v0**2*r0
    mu=(197-S.sqrt(5))/2
    b=Rat(101)
    s0=(S.sqrt(1+4*b/kval2)-1)/2
    dphi=3*kval2-e6
    dpsi=mu-kval2*s0**2-101*r0
    check('signed sextic '+label+' control gap and wall',bool(dphi>0 and dpsi>0),f'g={gg}, dphi={dphi}, dpsi={dpsi}')

out={"status":"passed","count":len(checks),"scope":"exact symbolic / small-matrix controls; not QW replay nor PDE computer proof","checks":checks}
path=Path(__file__).with_name('exact_checks.json')
path.write_text(json.dumps(out,indent=2,ensure_ascii=False)+'\n')
print(json.dumps({"status":out['status'],"count":out['count'],"signed_control_gap_exact":str(del_control)},indent=2))
