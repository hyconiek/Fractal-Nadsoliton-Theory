#!/usr/bin/env python3
"""Independent bounded algebra checks for the V-RP4-001 input scope. NOT proof of cone/bootstrap."""
import math,json,hashlib
from pathlib import Path
import sympy as sp
here=Path(__file__).resolve().parent
z=sp.symbols('z',real=True)
K=sp.tanh(z);sech2=1-sp.tanh(z)**2
f=sp.sqrt(sp.Rational(3,4))*sech2
S=sp.Rational(9,4)*K*sech2**2
V=(sp.Symbol('phi')**2-1)**2/sp.Integer(2)+(sp.Symbol('a')+6*sp.Symbol('phi')**2)*sp.Symbol('psi')**2/2+sp.Symbol('g4')*sp.Symbol('psi')**4/4+sp.Symbol('g6')*sp.Symbol('psi')**6/6
assert sp.expand(sp.diff(V,sp.Symbol('phi'))-(2*sp.Symbol('phi')*(sp.Symbol('phi')**2-1)+6*sp.Symbol('phi')*sp.Symbol('psi')**2))==0
assert sp.expand(sp.diff(V,sp.Symbol('psi'))-((sp.Symbol('a')+6*sp.Symbol('phi')**2)*sp.Symbol('psi')+sp.Symbol('g4')*sp.Symbol('psi')**3+sp.Symbol('g6')*sp.Symbol('psi')**5))==0
assert sp.simplify(-sp.diff(K,z,2)+2*K*(K**2-1))==0
Lphi=lambda e: -sp.diff(e,z,2)+(4-6*sech2)*e
assert sp.simplify(Lphi(f))==0
assert sp.simplify(3*K*f*f-S)==0
p,q=sp.symbols('p q',real=True,positive=True)
M=sp.Rational(9,4)*sp.I*q*(q*q+4)*p/(8*(1+sp.I*q)*(2+sp.I*q))
Gamma=sp.simplify(sp.expand_complex(M*sp.conjugate(M))/(2*q))
assert sp.simplify(Gamma-81*q*(q*q+4)*p**2/(2048*(1+q*q)))==0
# Verify exact Jost eigenfunction in coordinate y=tanh(x) and D_x=(1-y²)D_y.
y=sp.symbols('y',real=True);h=3*y*y-3*sp.I*q*y-(q*q+1)
D=lambda g: (1-y*y)*sp.diff(g,y)+sp.I*q*g
assert sp.simplify(-D(D(h))+(4-6*(1-y*y)-(q*q+4))*h)==0
vmin=1/math.sqrt(5); vmax=math.sqrt(3/7)
assert vmin>1/3
c_ref=625/(884736*math.sqrt(3));c_report=625/12386304
assert c_ref>c_report>0
for it in json.loads((here/'INPUTS.json').read_text())['actual_checked_input_files']:
 fp=Path(it['path'])
 assert fp.is_file() and hashlib.sha256(fp.read_bytes()).hexdigest()==it['sha256']
result={"bounded_checks_passed":True,"symbolic_potential_PDE":True,"source_S_harmonic":True,"Jost_eigen_equation":True,"Gamma_algebra":True,"velocity_range":[vmin,vmax],"fixed_source_coefficient_lower_conditional":c_ref,"reported_final_coefficient":c_report,"ratio":c_ref/c_report,"caveat":"NO retarded cone, uniform bootstrap or final flux theorem proof checked"}
(here/'VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
