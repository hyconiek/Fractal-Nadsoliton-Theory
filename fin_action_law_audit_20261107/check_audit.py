"""Exact projection countercontrol and analytic bookkeeping checks.
No nonlinear PDE simulation or independent theorem certification.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp
HERE=Path(__file__).resolve().parent
checks=[]
def check(name, condition):
    checks.append({'name':name,'pass':bool(condition)})
    if not condition:
        raise AssertionError(name)
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
    check('source integrity: '+row['path'],hashlib.sha256((HERE.parent/row['path']).read_bytes()).hexdigest()==row['sha256'])
P=sp.ones(2)/2
O=sp.diag(2,3)
I=sp.eye(2)
J=I.row_join(sp.I*O.inv()).col_join(I.row_join(-sp.I*O.inv()))
D=sp.diag(P,P)
C=sp.simplify(J*D*J.inv())
check('conjugated projection idempotent',C*C==C)
check('componentwise half-wave projection is different',C!=D)
s=sp.Matrix([1,1]);b=sp.Matrix([1,0])
correct=(s.T*P*O*b)[0]; wrong=(s.T*O*P*b)[0]
check('correct dual pairing retains operator order',correct==(s.T*O*b)[0]==2)
check('illegal commutation is detected',wrong==sp.Rational(5,2) and correct!=wrong)
r,t,nu=sp.symbols('r t nu',positive=True,real=True)
# Exact finite-primitive integration by parts for a nonconstant envelope.
F=(sp.exp(sp.I*nu*r)-1)/(sp.I*nu)
A=1+r
rhs=(A*F).subs(r,t)-sp.integrate(F,(r,0,t))
lhs=sp.integrate(A*sp.exp(sp.I*nu*r),(r,0,t))
check('finite primitive keeps endpoint',sp.simplify(lhs-rhs)==0)
check('dropping endpoint changes result',sp.simplify(lhs+sp.integrate(F,(r,0,t)))!=0)
x,y=sp.symbols('x y',nonzero=True)
Pz=(x-1/x)/(2*sp.I);Z=(x+1/x)/2
Rc=(y-1/y)/(2*sp.I);Cpos=(y+1/y)/2
check('self work has no DC',sp.expand(Pz*Z).coeff(x,0)==0 and sp.expand(Rc*Cpos).coeff(y,0)==0)
exponents=set()
for term in sp.Add.make_args(sp.expand(Pz*Cpos)):
    powers=term.as_powers_dict();exponents.add((int(powers.get(x,0)),int(powers.get(y,0))))
check('mixed work retains sum and difference carriers',exponents=={(1,1),(1,-1),(-1,1),(-1,-1)})
eps,S=sp.symbols('eps S',positive=True)
T=S/eps**2
check('source history exponent',sp.simplify(eps**3*T)==S*eps)
check('envelope total variation exponent',sp.simplify(sp.sqrt(T)*eps**3)==sp.sqrt(S)*eps**2)
check('Gronwall exponent stays finite at fixed S',sp.simplify(eps**2*T)==S)
result={'status':'PASS','scope':'Exact algebra and source integrity only; not independent proof review.',
'checks':checks,'check_count':len(checks),'projection_countercontrol':{'correct':str(correct),'illegal_reordering':str(wrong)},'sympy':sp.__version__}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':'PASS','check_count':len(checks),'countercontrol':result['projection_countercontrol']}))
