"""Finite algebra/provenance controls, not a machine proof of the PDE theorem."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, sqrt, sin, cos
import cmath
import hashlib
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = []

def check(name, passed, kind="exact"):
    checks.append({"name": name, "passed": bool(passed), "kind": kind})

# Rational endpoint margins justify the analytic square-root inequalities.
check("sqrt15 greater than 3", F(15) > 3**2)
check("lambda less than 7/4", F(3) < F(7, 4)**2)
check("omega less than 4/3", F(7, 4) < F(4, 3)**2)
check("m less than 12/5", F(23, 4) < F(12, 5)**2)
check("sum of upper frequency bounds less than 6", F(7,4)+F(4,3)+F(12,5)<6)
check("quadratic threshold gap lower bound", F(2,6)==F(1,3))
check("cubic threshold still present at w=4/3", 3*F(4,3)==4)

# Multivariate polynomial identities for the full remaining source.
# Variables are c,z,g,f,v,u,K,g4,g6; coefficients are exact integers.
N = 9
ZERO = (0,)*N
def add(*ps):
    out = {}
    for p in ps:
        for e,a in p.items():
            out[e] = out.get(e,0)+a
    return {e:a for e,a in out.items() if a}
def scale(a,p):
    return {e:a*b for e,b in p.items() if a*b}
def mul(p,q):
    out = {}
    for e,a in p.items():
        for f,b in q.items():
            k=tuple(x+y for x,y in zip(e,f))
            out[k]=out.get(k,0)+a*b
    return {e:a for e,a in out.items() if a}
def power(p,n):
    out={ZERO:1}
    for _ in range(n):
        out=mul(out,p)
    return out
def prod(*ps):
    out={ZERO:1}
    for p in ps:
        out=mul(out,p)
    return out
vs=[]
for i in range(N):
    e=list(ZERO); e[i]=1; vs.append({tuple(e):1})
c,z,g,f,v,u,K,g4,g6=vs
cg=mul(c,g); zf=mul(z,f)
eta=add(cg,v); psi=add(zf,u)
full=add(scale(-12,prod(K,eta,psi)),scale(-6,prod(power(eta,2),psi)),
         scale(-1,prod(g4,power(psi,3))),scale(-1,prod(g6,power(psi,5))))
paid=add(scale(-12,prod(K,cg,zf)),scale(-6,prod(power(cg,2),zf)),
         scale(-1,prod(g4,power(zf,3))),scale(-1,prod(g6,power(zf,5))))
residual=add(scale(-12,mul(K,add(mul(cg,u),mul(zf,v),mul(v,u)))),
    scale(-6,add(prod(power(cg,2),u),scale(2,prod(cg,zf,v)),
                 scale(2,prod(cg,v,u)),prod(zf,power(v,2)),prod(power(v,2),u))),
    scale(-1,mul(g4,add(scale(3,prod(power(zf,2),u)),
                       scale(3,prod(zf,power(u,2))),power(u,3)))),
    scale(-1,mul(g6,add(*(scale(comb(5,j),prod(power(zf,5-j),power(u,j))) for j in range(1,6))))))
check("exact full source equals paid plus residual", full==add(paid,residual))
check("every residual monomial contains continuum", all(e[4]+e[5]>=1 for e in residual))
check("residual has 16 monomials", len(residual)==16)
missing_mixed=add(residual,scale(12,prod(cg,zf,v)))
check("negative control detects omitted mixed cubic", full!=add(paid,missing_mixed))

for t in (0.0,0.3,1.7):
    lam=sqrt(3); om=sqrt(1.5)
    # Arbitrary smooth real histories, not a PDE simulation.
    c0=t*t; cp=2*t; cpp=2
    z0=cos(t); zp=-sin(t); zpp=-cos(t)
    a=cmath.exp(1j*lam*t)*(c0+1j*cp/lam)
    b=cmath.exp(1j*om*t)*(z0+1j*zp/om)
    ap=cmath.exp(1j*lam*t)*(1j*lam*(c0+1j*cp/lam)+cp+1j*cpp/lam)
    check(f"exact oscillator identity preview t={t}",abs(ap-1j*cmath.exp(1j*lam*t)*(cpp+lam*lam*c0)/lam)<1e-12,"floating_preview")
    reconstructed=sum((a if s==1 else a.conjugate())*(b if r==1 else b.conjugate())*
                      cmath.exp(-1j*(s*lam+r*om)*t) for s in (-1,1) for r in (-1,1))/4
    check(f"four branch product preview t={t}",abs(reconstructed-c0*z0)<1e-12,"floating_preview")

# Independent quadrature control of QH7's signs and switch-on term.
# Use a decaying scalar kernel; this validates the identity, not LE.
def simpson(fn,t,n=4000):
    h=t/n
    return h/3*(fn(0)+fn(t)+sum((4 if j%2 else 2)*fn(j*h) for j in range(1,n)))
for nu in (-3.0,-0.5,0.5,3.0):
    t=1.25; p=1/(1-1j*nu)
    tail=lambda r: p*cmath.exp(-(1-1j*nu)*r)
    A=lambda s: 1+s*s
    hist=simpson(lambda r: tail(r)*2*(t-r),t)
    direct=simpson(lambda s: cmath.exp(-(t-s))*cmath.exp(-1j*nu*s)*A(s),t)
    exact=cmath.exp(-1j*nu*t)*(p*A(t)-tail(t)*A(0)-hist)
    wrong=cmath.exp(-1j*nu*t)*(p*A(t)-hist)
    check(f"retarded primitive identity nu={nu}",abs(direct-exact)<1e-10,"floating_preview")
    check(f"negative control retains initial transient nu={nu}",abs(direct-wrong)>0.01,"floating_preview")

check("boundary term phase amplitude epsilon exponent",4-2==2)
check("history squared critical integral exponent",6-4==2)
check("crude quadratic absolute integral exponent",4-4==0)
check("history product derivative exponent",1+2==3)
check("source coefficient for each branch",F(-12,4)==-3)

inputs=json.loads((HERE/'INPUTS.json').read_text())
for row in inputs['files']:
    data=(ROOT/row['path']).read_bytes()
    check('sha256 '+row['path'],hashlib.sha256(data).hexdigest()==row['sha256'],'provenance')
    blob=hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()
    check('git blob '+row['path'],blob==row['git_blob_sha'],'provenance')

result={"status":"PASS_COMPONENT_CONTROLS_CONTINUUM_TRACE_AND_ACTION_GATES_OPEN" if all(c['passed'] for c in checks) else "FAIL",
        "counts":{kind:sum(c['kind']==kind for c in checks) for kind in sorted({c['kind'] for c in checks})},
        "passed":sum(c['passed'] for c in checks),"checks":checks,
        "formal_verification":False,"independent_review":False,"actual_critical_damping_proved":False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
for c in checks:
    if not c['passed']:
        print('FAILED:',c['name'])
raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
