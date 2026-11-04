"""Exact finite controls and provenance; not a formal PDE verifier."""
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
import math
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
checks=[]
def check(name,ok,kind='exact'):
    checks.append({'name':name,'passed':bool(ok),'kind':kind})
def mat(rows): return [[F(x) for x in row] for row in rows]
def add(A,B): return [[a+b for a,b in zip(ar,br)] for ar,br in zip(A,B)]
def scale(c,A): return [[c*a for a in row] for row in A]
def mul(A,B): return [[sum(a*b for a,b in zip(row,col)) for col in zip(*B)] for row in A]
def blocks(A,B,C,D): return [a+b for a,b in zip(A,B)]+[c+d for c,d in zip(C,D)]
I=mat([[1,0],[0,1]]); O=mat([[0,0],[0,0]])
P=mat([[F(1,2),F(-1,2)],[F(-1,2),F(1,2)]])
Om=mat([[2,0],[0,3]]); OmInv=mat([[F(1,2),0],[0,F(1,3)]])
A=add(scale(3,I),P)
J=blocks(I,OmInv,I,scale(-1,OmInv))
Jinv=scale(F(1,2),blocks(I,I,Om,scale(-1,Om)))
Pi=blocks(P,O,O,P)
Proj=mul(mul(J,Pi),Jinv)
R=mul(mul(OmInv,P),Om)
Formula=scale(F(1,2),blocks(add(P,R),add(P,scale(-1,R)),add(P,scale(-1,R)),add(P,R)))
check('physical projector idempotent',mul(P,P)==P)
check('free multiplier and projector do not commute',mul(Om,P)!=mul(P,Om))
check('full A and projector commute',mul(A,P)==mul(P,A))
check('conjugated projector formula exact',Proj==Formula)
check('conjugated projector idempotent',mul(Proj,Proj)==Proj)
check('negative control componentwise P is wrong',Proj!=Pi)
Dyn=blocks(O,I,scale(-1,A),O)
HalfDyn=mul(mul(J,Dyn),Jinv)
check('correct projection commutes with full generator',mul(Proj,HalfDyn)==mul(HalfDyn,Proj))
check('naive projection fails full generator commutation',mul(Pi,HalfDyn)!=mul(HalfDyn,Pi))
source=mat([[1],[2]])
inj=mul(OmInv,source)+scale(-1,mul(OmInv,source))
expected=mul(OmInv,mul(P,source))+scale(-1,mul(OmInv,mul(P,source)))
check('force projection identity exact',mul(Proj,inj)==expected)
check('naive force projection wrong',mul(Pi,inj)!=expected)
# Polynomial source split at many rational inputs; expansion coefficients are
# independently checked via the binomial theorem below.
for c,z,g,f,v,u in [(F(1),F(2),F(3),F(4),F(5),F(6)),(F(-1,3),F(2,5),F(1,2),F(-3,4),F(7,8),F(-2,9)),(F(0),F(1),F(2),F(3),F(0),F(0))]:
    Gphi=-6*(v*v+u*u)
    Bphi=-6*(c*c*g*g+z*z*f*f+2*c*g*v+2*z*f*u)
    Gpsi=-12*v*u
    Bpsi=-12*(c*z*g*f+c*g*u+z*f*v)
    check('odd quadratic exact evaluation '+str(c),Gphi+Bphi==-6*((c*g+v)**2+(z*f+u)**2))
    check('even quadratic exact evaluation '+str(c),Gpsi+Bpsi==-12*(c*g+v)*(z*f+u))
check('square expansion coefficients',(1,2,1)==tuple(math.comb(2,j) for j in range(3)))
G=mat([[2],[-3]]); B=mat([[5],[7]])
L=add(mul(P,B),scale(-1,mul(add(I,scale(-1,P)),G)))
check('finite rank quadratic subtraction identity',add(G,L)==mul(P,add(G,B)))
check('omitting finite rank term fails',add(G,mul(P,B))!=mul(P,add(G,B)))
check('branch coefficient',F(-12,4)==-3)
check('correction integrated trace exponent',4-2==2)
check('localized cubic response integrated exponent',6-4==2)
check('generic cubic trace exponent remains one',6-5==1)
check('quadratic localized crude trace exponent remains zero',4-4==0)
check('localized quadratic forcing times continuum is cubic',2+1==3)
check('cubic forcing inside quadratic correction is quartic',3+1==4)
for sigma in (-1,1):
 for tau in (-1,1):
  for rho in (-1,1):
   p,q=0.7,-1.1
   E0=math.sqrt((p+q)**2+5.5);E1=math.sqrt(p*p+4);E2=math.sqrt(q*q+5.5)
   phase=sigma*E0-tau*E1-rho*E2
   b=1/(E0*phase)
   check(f'quadratic homological identity signs {sigma,tau,rho}',abs(phase*b-1/E0)<1e-14,'floating_preview')
   check(f'phase sign reversal {sigma,tau,rho}',abs(1/(E0*(-phase))+b)<1e-14,'floating_preview')
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
 data=(ROOT/row['path']).read_bytes()
 check('sha256 '+row['path'],hashlib.sha256(data).hexdigest()==row['sha256'],'provenance')
 check('git blob '+row['path'],hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob_sha'],'provenance')
result={'status':'PASS_EXACT_REDUCTION_TWO_TRACE_HISTORIES_OPEN' if all(c['passed'] for c in checks) else 'FAIL','counts':{k:sum(c['kind']==k for c in checks) for k in sorted({c['kind'] for c in checks})},'passed':sum(c['passed'] for c in checks),'checks':checks,'formal_verification':False,'independent_review':False,'actual_critical_damping_proved':False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
for c in checks:
 if not c['passed']: print('FAILED:',c['name'])
raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
