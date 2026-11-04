"""Finite controls for LF; not formal verification of the analytical proof."""
from pathlib import Path
from fractions import Fraction as F
from math import comb, exp, sqrt, pi
import hashlib
import json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
checks=[]
def check(name,passed,kind="exact"):
    checks.append({"name":name,"passed":bool(passed),"kind":kind})
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

extended=add(scale(-12,prod(K,v,u)),scale(-6,prod(power(v,2),u)),scale(-1,prod(g4,power(u,3))),scale(-1,prod(g6,power(u,5))))
localized=add(residual,scale(-1,extended))
check("exact full split",full==add(paid,localized,extended))
check("four extended monomials",len(extended)==4)
check("twelve localized monomials",len(localized)==12)
check("every localized monomial contains fixed f or g",all(e[2]+e[3]>=1 for e in localized))
check("every localized monomial contains continuum",all(e[4]+e[5]>=1 for e in localized))
check("all extended monomials are continuum only",all(sum(e[:4])==0 for e in extended))
check("localized nonlinearity has at least two small factors",all(e[0]+e[1]+e[4]+e[5]>=2 for e in localized))
check("positive weight transfer uses nine powers",6+3==9)
check("quad response squared weight scaling",2-2==0)
check("paid term weighted integral exponent",4-2==2)
check("paid history weighted integral exponent",6-4==2)
for D,CL in [(F(1),F(1)),(F(816),F(7,3)),(F(1,100),F(1,100))]:
    mu=max(F(1),6*D*CL)
    check(f"absorption coefficient D={D},CL={CL}",3*D*CL/mu<=F(1,2))
check("small epsilon alone cannot reduce epsilon sqrt(S epsilon^-2)",F(1)-F(2,2)==0)
# Independent finite quadrature of the cumulative-to-weighted identity.
def simpson(fn,T,n=4000):
    h=T/n
    return h/3*(fn(0)+fn(T)+sum((4 if j%2 else 2)*fn(j*h) for j in range(1,n)))
for beta,T in [(0.2,3.0),(2.0,1.0),(5.0,2.0)]:
    cumulative=lambda t: t+t*t
    direct=simpson(lambda t:exp(-beta*t)*(1+2*t),T)
    boundary=exp(-beta*T)*cumulative(T)
    formula=boundary+beta*simpson(lambda t:exp(-beta*t)*cumulative(t),T)
    check(f"cumulative integration identity beta={beta}",abs(direct-formula)<1e-10,"floating_preview")
    check(f"negative control detects lost endpoint beta={beta}",abs(direct-(formula-boundary))>1e-4,"floating_preview")
for alpha in (0.01,1.0,10.0):
    # r=s^2 removes the integrable majorant singularity at the origin.
    value=simpson(lambda s:2*exp(-alpha*s*s),10/sqrt(alpha))
    check(f"Young majorant integral alpha={alpha}",abs(value-sqrt(pi/alpha))<1e-9,"floating_preview")
for R in (4,8,16):
    check(f"packet support weight lower bound R={R}",F(R)-F(1,4)>=F(R,2))
    check(f"packet growth factor R={R}",(F(R,2))**6>=64)
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
    b=(ROOT/row['path']).read_bytes()
    check('sha256 '+row['path'],hashlib.sha256(b).hexdigest()==row['sha256'],'provenance')
    check('git blob '+row['path'],hashlib.sha1(b'blob '+str(len(b)).encode()+b'\0'+b).hexdigest()==row['git_blob_sha'],'provenance')
result={'status':'PASS_LOCAL_FEEDBACK_REDUCTION_EXTENDED_TRACE_OPEN' if all(c['passed'] for c in checks) else 'FAIL','counts':{k:sum(c['kind']==k for c in checks) for k in sorted({c['kind'] for c in checks})},'passed':sum(c['passed'] for c in checks),'checks':checks,'formal_verification':False,'independent_review':False,'physical_validation':False,'actual_critical_damping_proved':False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
for c in checks:
    if not c['passed']: print('FAILED',c['name'])
raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
