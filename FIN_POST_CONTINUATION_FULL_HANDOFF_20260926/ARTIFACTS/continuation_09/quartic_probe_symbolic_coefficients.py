#!/usr/bin/env python3
import sympy as sp
Q=12
c3,c4,c5,c6,g=sp.symbols('c3 c4 c5 c6 g', nonnegative=True, real=True)
l3,l4,l5,l6=sp.symbols('lambda3 lambda4 lambda5 lambda6', positive=True, real=True)
cs=[c3,c4,c5,c6]; ls=[l3,l4,l5,l6]; ks=[3,4,5,6]
vecs=[]
for k in ks:
    if k==6:
        v=sp.Matrix([(-1)**n/sp.sqrt(12) for n in range(Q)])
    else:
        v=sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.cos(2*sp.pi*k*n/Q) for n in range(Q)])
    vecs.append(v)
phi=sum((cs[i]*vecs[i] for i in range(4)), sp.zeros(Q,1))
Aphi=sum((ls[i]*cs[i]*vecs[i] for i in range(4)), sp.zeros(Q,1))
hidden=[]
for k in (1,2):
    hidden += [sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.cos(2*sp.pi*k*n/Q) for n in range(Q)]),
               sp.Matrix([sp.sqrt(sp.Rational(2,Q))*sp.sin(2*sp.pi*k*n/Q) for n in range(Q)])]
PH=sum((v*v.T for v in hidden), sp.zeros(Q,Q))
h=sp.simplify(PH*phi.multiply_elementwise(phi))
k=sp.simplify(PH*phi.multiply_elementwise(Aphi))
C=sp.Poly(sp.expand(-12*h.dot(h)/Q+2*g*h.dot(k)/Q),c3,c4,c5,c6)
for mon,coef in C.terms(): print(mon, sp.factor(coef))
