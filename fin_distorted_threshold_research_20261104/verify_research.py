"""Exact operator algebra, finite quadrature and source hashes; not a PDE proof."""
from pathlib import Path
from fractions import Fraction as F
from math import comb,sqrt,tanh,cosh,sinh,pi
import cmath,hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
checks=[]
def check(name,ok,kind='exact'):
 checks.append({'name':name,'passed':bool(ok),'kind':kind})
def add(*ps):
 out={}
 for p in ps:
  for k,v in p.items():out[k]=out.get(k,F(0))+v
 return {k:v for k,v in out.items() if v}
def scale(a,p):return {k:a*v for k,v in p.items() if a*v}
def mul(p,q):
 out={}
 for (a,b),v in p.items():
  for (c,d),w in q.items():out[a+c,b+d]=out.get((a+c,b+d),F(0))+v*w
 return {k:v for k,v in out.items() if v}
def power(p,n):
 out={(0,0):F(1)}
 for _ in range(n):out=mul(out,p)
 return out
one={(0,0):F(1)};x={(1,0):F(1)};y={(0,1):F(1)}
s=add(one,scale(-1,power(y,2)))
def derivative(p):
 dx={(a-1,b):a*v for (a,b),v in p.items() if a}
 dy={(a,b-1):b*v for (a,b),v in p.items() if b}
 return add(dx,mul(s,dy))
def H(p):return add(scale(-1,derivative(derivative(p))),mul(add(scale(-2,one),scale(6,power(y,2))),p))
r0=scale(F(1,2),add(scale(3,power(y,2)),scale(-1,one)))
source=mul(y,power(r0,2))
U0=add(scale(F(5,8),y),scale(F(-3,8),power(y,3)),scale(F(-1,2),mul(x,s)))
check('exact threshold eigenfunction',H(r0)==scale(4,r0))
check('static odd response identity',H(U0)==source)
check('translation profile zero mode',H(s)=={})
check('x times translation identity',H(mul(x,s))==scale(4,mul(y,s)))
wrong=add(scale(F(5,8),y),scale(F(-3,8),power(y,3)))
check('negative control omission of x sech^2 fails',H(wrong)!=source)
# Integrate decaying even y-polynomials against dx=dy/(1-y^2).
def divide_by_s(p):
 assert all(a==0 for a,b in p)
 rem={b:v for (a,b),v in p.items()};out={}
 while rem and max(rem)>=2:
  n=max(rem);c=-rem[n];k=n-2
  out[k]=out.get(k,F(0))+c
  rem[k]=rem.get(k,F(0))-c;rem[n]=rem.get(n,F(0))+c
  rem={k:v for k,v in rem.items() if v}
 assert not rem,'nondecaying integrand'
 return out
def integral_dx(p):return sum(v*F(2,k+1) for k,v in divide_by_s(p).items() if k%2==0)
def integral_dy(p):
 assert all(a==0 for a,b in p)
 return sum(v*F(2,b+1) for (a,b),v in p.items() if b%2==0)
check('r0 orthogonal f',integral_dy(r0)==0)
check('r0^3 f overlap divided by sqrt3',F(1,2)*integral_dy(power(r0,3))==F(2,35))
check('fourth threshold defect',integral_dx(add(power(r0,4),scale(-1,one)))==F(-132,35))
check('second threshold defect',integral_dx(add(power(r0,2),scale(-1,one)))==-3)
check('projected direct cubic defect',integral_dx(add(power(r0,4),scale(-1,power(r0,2))))==F(-27,35))
# Arcsine moments give the normalized shape overlap exactly.
p=mul(y,source)
shape=sum(v*F(comb(b,b//2),2**b) for (a,b),v in p.items() if b%2==0)
check('shape overlap divided by pi sqrt(3/2)',shape==F(17,64))
for z in (F(21,4),F(16,3),F(23,4)):
 a=-F(1)/(4*(z-1))
 J=add(scale(-3-F(1)/(z-1),mul(y,s)),scale(F(9,4),mul(y,power(s,2))))
 lhs=add(source,scale(-1,add(H(scale(a,y)),scale(-4*z,scale(a,y)))))
 check(f'exact second-harmonic residual z={z}',lhs==J)
 q2=4*z-4
 factor=-3-F(1)/(z-1)+F(9,4)*(q2+4)/6
 expected=-3-F(1)/(z-1)+F(3,2)*z
 check(f'on-shell transform factor z={z}',factor==expected)
 check(f'on-shell factor lower bound z={z}',factor>=F(631,136))
 check(f'open odd channel z={z}',17<=q2<=19)
 check(f'shape denominator nonzero z={z}',3-4*z!=0)
 check(f'homogeneous static coefficient recovered z={z}',18*F(1,4)==F(9,2))
 check(f'homogeneous harmonic coefficient recovered z={z}',9*a==-F(9)/(4*(z-1)))
check('minimum transform factor exact',-3-F(4,17)+F(3,2)*F(21,4)==F(631,136))
# Independent quadrature of actual distorted Fourier integrals.
def simpson(fn,lo,hi,n=40000):
 h=(hi-lo)/n
 return h/3*(fn(lo)+fn(hi)+sum((4 if j%2 else 2)*fn(lo+j*h) for j in range(1,n)))
for z in (5.25,16/3,5.75):
 k=sqrt(4*z-4);den=(1+1j*k)*(2+1j*k)
 def integrand(t):
  K=tanh(t);ss=1/cosh(t)**2
  J=(-3-1/(z-1))*K*ss+2.25*K*ss*ss
  psi_bar=cmath.exp(-1j*k*t)*(3*K*K-1-k*k+3j*k*K)/den
  return psi_bar*J
 numeric=simpson(integrand,-20,20)
 base=1j*k*(k*k+4)*(pi*k/sinh(pi*k/2))/(8*den)
 exact=(-3-1/(z-1)+1.5*z)*base
 check(f'distorted transform quadrature z={z}',abs(numeric-exact)<1e-9,'floating_preview')
 check(f'nonzero spectral density preview z={z}',abs(exact)**2/(2*k)>0,'floating_preview')
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
 data=(ROOT/row['path']).read_bytes()
 check('sha256 '+row['path'],hashlib.sha256(data).hexdigest()==row['sha256'],'provenance')
 check('git blob '+row['path'],hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob_sha'],'provenance')
result={'status':'PASS_FULL_POTENTIAL_REFERENCE_ACTUAL_PACKET_HISTORY_OPEN' if all(c['passed'] for c in checks) else 'FAIL','counts':{k:sum(c['kind']==k for c in checks) for k in sorted({c['kind'] for c in checks})},'passed':sum(c['passed'] for c in checks),'checks':checks,'actual_critical_damping_proved':False,'formal_verification':False,'independent_review':False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
for c in checks:
 if not c['passed']:print('FAILED',c['name'])
raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
