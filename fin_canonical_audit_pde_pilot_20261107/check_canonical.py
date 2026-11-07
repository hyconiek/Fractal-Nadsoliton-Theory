"""Exact CM homological and preparation checks, independently reconstructed."""
from pathlib import Path
import json,hashlib
import sympy as s
import numpy as np
import canonical as C
H=Path(__file__).resolve().parent;checks=[]
def check(name,condition):
    checks.append({'name':name,'pass':bool(condition)})
    if not condition:raise AssertionError(name)
for row in json.loads((H/'INPUTS.json').read_text())['files']:
    check('input: '+row['path'],hashlib.sha256((H.parent/row['path']).read_bytes()).hexdigest()==row['sha256'])
z,p,c,r=C.coords;w,g4,d=C.w,C.g4,C.d
check('cubic homological sign',s.simplify(C.bracket(C.G3,C.H2)-C.H3)==0)
D,Pi,Q,Z,G4,blocks=C.construct()
check('19-dimensional space',len(C.basis)==19)
check('angle projector idempotent',(Pi*Pi-Pi).applyfunc(s.simplify)==s.zeros(19))
check('angle projector annihilates derivative',D*Pi==s.zeros(19) and Pi*D==s.zeros(19))
det=s.factor(s.prod(b.det() for b in blocks))
check('quartic uniform invertibility determinant',s.simplify(det-2**20*w**3*27*(w-3)**2)==0)
check('quartic homological identity',s.simplify(C.bracket(G4,C.H2)-(Q-Z))==0)
check('quartic zero average',s.simplify(C.average(G4))==0)
Jz=p*p+w*z*z;Jc=r*r+3*c*c
Az=27*g4/(560*w*w)+9*d*d*(8*w-9)/(12*w*w*(3-4*w))
Bz=9/(210*w)+8*d*d/(3-4*w)
Cz=s.Rational(27,5040)-15*d*d/27
check('all three coherent coefficients',s.simplify(Z-Az*Jz**2-Bz*Jz*Jc-Cz*Jc**2)==0)
for name,G in [('G3',C.G3),('G4',G4)]:
    check(name+' momentum reversal',s.simplify(G.subs({p:-p,r:-r},simultaneous=True)+G)==0)
    check(name+' psi sign symmetry',s.simplify(G.subs({z:-z,p:-p},simultaneous=True)-G)==0)
check('exceptional cubic fixed point',all(s.simplify(s.diff(C.G3,q).subs({p:0,c:0,r:0,w:s.Rational(3,2)}))==0 for q in C.coords))
check('quartic preserves zero shape',all(s.simplify(s.diff(G4,q).subs({c:0,r:0}))==0 for q in [c,r]))
check('shape reconstruction DC',s.simplify(C.A+4*w/(3-4*w)+2)==0)
check('shape reconstruction harmonic',s.simplify(C.A-4*w/(3-4*w)+6/(3-4*w))==0)
err=max(np.max(abs(C.pullback([.2,0,0,0],wv,1.,16)-C.pullback([.2,0,0,0],wv,1.,32))) for wv in [1.25,1.5,1.75])
check('auxiliary inverse integration refinement',err<1e-9)
out={'status':'PASS','checks':checks,'check_count':len(checks),'inverse_map_refinement':float(err),'scope':'Exact polynomial identities plus finite-dimensional flow refinement; not independent PDE proof review.'}
(H/'canonical_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'checks':len(checks),'refinement':err}))
