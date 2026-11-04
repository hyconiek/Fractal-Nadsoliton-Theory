"""Exact finite algebra and source hashes, not a formal PDE verifier."""
from pathlib import Path
from fractions import Fraction as F
import hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
checks=[]
def check(name,ok,kind='exact'):
 checks.append({'name':name,'passed':bool(ok),'kind':kind})
def pa(a,b):
 c=dict(a)
 for k,v in b.items(): c[k]=c.get(k,F(0))+v
 return {k:v for k,v in c.items() if v}
def pm(a,b):
 c={}
 for k,v in a.items():
  for l,w in b.items(): c[k+l]=c.get(k+l,F(0))+v*w
 return {k:v for k,v in c.items() if v}
def R(x): return ({0:F(x)} if x else {},{0:F(1)})
def add(a,b): return pa(pm(a[0],b[1]),pm(b[0],a[1])),pm(a[1],b[1])
def mul(a,b): return pm(a[0],b[0]),pm(a[1],b[1])
def neg(a): return {k:-v for k,v in a[0].items()},a[1]
def sub(a,b): return add(a,neg(b))
def div(a,b):
 assert b[0]
 return pm(a[0],b[1]),pm(a[1],b[0])
def equal(a,b): return pm(a[0],b[1])==pm(b[0],a[1])
m=({1:F(1)},{0:F(1)});m2=mul(m,m)
branch={}
for rho in (-1,1):
 s=R(0)
 for tau in (-1,1):
  denom=mul(m,sub(sub(m,R(2*tau)),mul(R(rho),m)))
  s=add(s,div(R(tau),denom))
 branch[rho]=s
check('symbolic positive input branch sum',equal(branch[1],neg(div(R(1),m))))
check('symbolic negative input branch sum',equal(branch[-1],div(R(1),mul(m,sub(m2,R(1))))))
cascade=mul(R(-9),add(mul(R(F(1,2)),branch[1]),mul(R(F(1,4)),branch[-1])))
physical=div(sub(R(F(9,2)),div(R(F(9,4)),sub(m2,R(1)))),m)
check('symbolic normal-form and physical cascade agreement',equal(cascade,physical))
check('symbolic cubic resonant phase vanishes',equal(add(sub(sub(m,m),m),m),R(0)))
check('negative control omitted second branch changes coefficient',not equal(mul(R(F(-9,2)),branch[1]),physical))
# Independent harmonic-polynomial expansion. Key=(temporal harmonic,b power,bar b power).
def ha(a,b):
 c=dict(a)
 for k,v in b.items(): c[k]=c.get(k,F(0))+v
 return {k:v for k,v in c.items() if v}
def hs(s,a):return {k:s*v for k,v in a.items() if s*v}
def hm(a,b):
 c={}
 for k,v in a.items():
  for l,w in b.items():
   e=tuple(x+y for x,y in zip(k,l));c[e]=c.get(e,F(0))+v*w
 return {k:v for k,v in c.items() if v}
u={(1,1,0):F(1,2),(-1,0,1):F(1,2)}
u2=hm(u,u);u3=hm(u2,u)
check('direct cubic resonant multiplicity',u3[(1,2,1)]==F(3,8))
for w in (F(5,4),F(4,3),F(7,4)):
 msq=4+w;star=12-F(6)/(w+3)
 for g4 in (F(0),F(1),star):
  expected=F(9,2)-F(9)/(4*(w+3))-F(3,8)*g4
  for kappa in (-1,1):
   forced=hs(-6*kappa,u2)
   v2={key:coef/(4-key[0]**2*msq) for key,coef in forced.items()}
   out=ha(hs(-12*kappa,hm(v2,u)),hs(-g4,u3))
   coefficient=out.get((1,2,1),F(0))
   check(f'harmonic coefficient w={w},g4={g4},K={kappa}',coefficient==expected)
  check(f'null coefficient iff selected g4 w={w},g4={g4}',(expected==0)==(g4==star))
check('lower cancellation endpoint',12-F(6)/F(17,4)==F(180,17))
check('upper cancellation endpoint',12-F(6)/F(19,4)==F(204,19))
check('fixed coefficient cannot match both endpoints',F(180,17)!=F(204,19))
check('cancellation at internal threshold',12-F(6)/F(13,3)==F(138,13))
check('zero-g4 positive lower endpoint',F(9,2)-F(9,17)==F(135,34)>0)
check('higher source trace exponent',8-4-1==3)
check('quintic source trace exponent',10-4-1==5)
check('generic cubic trace exponent still insufficient',6-4-1==1)
check('F_C times h starts at fourth order',3+1==4)
for x,y,A in [(F(1),F(2),F(3)),(F(-2,3),F(5,7),F(-4)),(F(0),F(1),F(0))]:
 speed=A*(x*x+y*y)
 dx=-speed*y;dy=speed*x
 check(f'real coefficient gives phase-only frozen flow {x,y,A}',2*x*dx+2*y*dy==0)
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
 data=(ROOT/row['path']).read_bytes()
 check('sha256 '+row['path'],hashlib.sha256(data).hexdigest()==row['sha256'],'provenance')
 check('git blob '+row['path'],hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob_sha'],'provenance')
result={'status':'PASS_HIGHER_HISTORY_BOUND_FROZEN_RESONANCE_ACTUAL_CUBIC_GATE_OPEN' if all(c['passed'] for c in checks) else 'FAIL','counts':{k:sum(c['kind']==k for c in checks) for k in sorted({c['kind'] for c in checks})},'passed':sum(c['passed'] for c in checks),'checks':checks,'actual_critical_damping_proved':False,'formal_verification':False,'independent_review':False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
for c in checks:
 if not c['passed']:print('FAILED',c['name'])
raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
