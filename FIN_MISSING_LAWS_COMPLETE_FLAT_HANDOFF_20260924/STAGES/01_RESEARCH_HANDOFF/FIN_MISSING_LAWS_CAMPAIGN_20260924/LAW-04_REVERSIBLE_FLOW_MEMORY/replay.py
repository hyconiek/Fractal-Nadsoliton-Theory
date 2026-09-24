#!/usr/bin/env python3
import sympy as sp, json
s,w,O,g=sp.symbols('s w O g', nonzero=True)
D=sp.Matrix([[s**2+w**2,g],[g,s**2+O**2]])
schur=sp.simplify(D[0,0]-D[0,1]*D[1,1]**-1*D[1,0])
x,y,px,py=sp.symbols('x y px py')
H=sp.Rational(1,2)*(px**2+py**2+w**2*x**2+O**2*y**2+2*g*x*y)
dH=sp.diff(H,x)*px+sp.diff(H,y)*py+sp.diff(H,px)*(-w**2*x-g*y)+sp.diff(H,py)*(-O**2*y-g*x)
out={'schur':str(schur),'expected':str(s**2+w**2-g**2/(s**2+O**2)),'schur_identity':bool(sp.simplify(schur-(s**2+w**2-g**2/(s**2+O**2)))==0),'Hamiltonian_derivative':str(sp.simplify(dH))}
print(json.dumps(out,indent=2,sort_keys=True))
