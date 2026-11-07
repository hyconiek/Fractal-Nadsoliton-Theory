"""Exact degeneracy witness for matching a nonlinear PDE from its Hessian."""
from pathlib import Path
import json
import sympy as s
phi,psi,a,g4,g6,d4,d6=s.symbols('phi psi a g4 g6 d4 d6')
V=(phi*phi-1)**2/2+(a+6*phi*phi)*psi*psi/2+g4*psi**4/4+g6*psi**6/6
change=d4*psi**4/4+d6*psi**6/6
checks=[]
def check(name,cond):
    checks.append({'name':name,'pass':bool(cond)})
    if not cond:raise AssertionError(name)
check('phi force',s.expand(s.diff(V,phi)-2*phi*(phi*phi-1)-6*phi*psi*psi)==0)
check('psi force',s.expand(s.diff(V,psi)-(a+6*phi*phi)*psi-g4*psi**3-g6*psi**5)==0)
check('same kink potential and first variation',change.subs(psi,0)==0 and all(s.diff(change,q).subs(psi,0)==0 for q in [phi,psi]))
check('same entire kink Hessian',s.hessian(change,[phi,psi]).subs(psi,0)==s.zeros(2))
check('different nonlinear equations',s.diff(change,psi)==d4*psi**3+d6*psi**5)
check('quartic matching target',s.diff(V,psi,4).subs(psi,0)==6*g4)
check('sixth-order matching target',s.diff(V,psi,6)==120*g6)
check('mixed matching target',s.diff(V,phi,2,psi,2)==12)
Path(__file__).with_name('bridge_verification.json').write_text(json.dumps({'status':'PASS','scope':'Spectral matching underdetermination only, not a universal FIN no-go result.','checks':checks},indent=2)+'\n')
print('8 exact bridge controls passed')
