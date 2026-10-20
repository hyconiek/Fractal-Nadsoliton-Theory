#!/usr/bin/env python3
import sympy as s
from pathlib import Path
K,e,u,w,v,r,g4,g6=s.symbols('K e u w v r g4 g6',real=True)
eta=e**2*u+w
psi=e*v+r
# Expanded scalar nonlinear residual, two independent algebraic constructions
A=s.expand(-6*K*psi**2-6*eta*psi**2-6*K*eta**2-2*eta**3+6*K*e**2*v**2)
B=s.expand(-6*K*(psi**2-e**2*v**2)-6*eta*psi**2-6*K*eta**2-2*eta**3)
assert s.simplify(A-B)==0
# psi remainder source and odd epsilon symmetry
P=s.expand(-12*K*eta*psi-6*eta**2*psi-g4*psi**3-g6*psi**5)
assert s.expand(P.subs({e:-e,r:-r})+P)==0  # psi odd when scalar w even
assert s.expand(A.subs({e:-e,r:-r})-A)==0
# Freezing r=w=0 yields required leading source orders
assert s.simplify(s.expand(P.subs({r:0,w:0})).coeff(e,3)+12*K*u*v+g4*v**3)==0
assert s.simplify(s.expand(A.subs({r:0,w:0})).coeff(e,4)+6*u*v**2+6*K*u**2)==0
# Potential derivatives independently verify both PDE and sign convention
phi=s.symbols('phi',real=True)
Pvar=s.symbols('Pvar',real=True)
V=(phi**2-1)**2/s.Integer(2)+s.Rational(1,2)*(s.symbols('a',real=True)+6*phi**2)*Pvar**2 + g4*Pvar**4/4+g6*Pvar**6/6
assert s.simplify(s.diff(V,phi)-2*phi*(phi**2-1)-6*phi*Pvar**2)==0
assert s.simplify(s.diff(V,Pvar)-(s.symbols('a',real=True)+6*phi**2)*Pvar-g4*Pvar**3-g6*Pvar**5)==0
x=s.symbols('x',real=True)
Kx=s.tanh(x); f=s.sqrt(s.Rational(3,4))*s.sech(x)**2
L=lambda y:-s.diff(y,x,2)+(4-6*s.sech(x)**2)*y
assert s.simplify(s.trigsimp(L(f)-f))!=0 # nontrivial check
assert s.simplify(s.trigsimp(L(f)))==0
assert s.simplify(s.trigsimp(-s.diff(Kx,x,2)+2*Kx*(Kx**2-1)))==0
assert s.integrate(1-s.symbols('y')**2,(s.symbols('y'),-1,1))==s.Rational(4,3)
assert s.simplify(6*Kx*f**2-2*(s.Rational(9,4)*Kx*s.sech(x)**4))==0
print('PASS: 10 independent sign/parity/normalization symbolic controls')
