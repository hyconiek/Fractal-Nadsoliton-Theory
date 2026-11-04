"""Finite exact and quadrature controls; not a nonlinear PDE simulation."""
from pathlib import Path
from fractions import Fraction as F
from math import sqrt,pi,exp,tanh,asinh
import cmath,csv,hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
checks=[];rows=[]
def check(name,ok,kind='exact'):
 checks.append({'name':name,'passed':bool(ok),'kind':kind})
def simpson(fn,a,b,n=6000):
 h=(b-a)/n
 return h/3*(fn(a)+fn(b)+sum((4 if j%2 else 2)*fn(a+j*h) for j in range(1,n)))
# Normalized exp(-q^2) moments follow M_2n=(2n-1)M_2n-2/2.
M={0:F(1)}
for n in (1,2,3):M[2*n]=F(2*n-1,2)*M[2*n-2]
check('Gaussian spectral second moment',M[2]==F(1,2))
check('Gaussian spectral fourth moment',M[4]==F(3,4))
check('Gaussian second derivative norm',M[4]-2*M[2]+M[0]==F(3,4))
check('local source error exponent',F(-1,2)+F(-5,2)==-3)
check('trace leading power in packet width',-1+2==1)
check('trace error power in packet width',F(-1,2)+F(-5,2)+2==-1)
check('small packet trace epsilon exponent',2-1==1)
check('small packet trace remainder epsilon exponent',2+1==3)
# Numerical Gaussian moments and exact-envelope Fourier integral.
for power,target in [(0,1.0),(2,0.5),(4,0.75)]:
 val=simpson(lambda q:q**power*exp(-q*q)/sqrt(pi),-10,10)
 check(f'Gaussian moment quadrature {power}',abs(val-target)<1e-11,'floating_preview')
for m in (sqrt(5.25),sqrt(16/3),sqrt(5.75)):
 for tau in (0.5,2.0):
  A=pi**(-0.25)*(1+1j*tau/m)**(-0.5)
  val=simpson(lambda q:pi**(-0.25)/sqrt(2*pi)*cmath.exp(-q*q/2-1j*tau*q*q/(2*m)),-10,10)
  check(f'Schrodinger envelope quadrature m={m},tau={tau}',abs(val-A)<1e-11,'floating_preview')
  limit=(1+1j*tau/(2*m))**(-0.5)
  check(f'frozen carrier limit differs m={m},tau={tau}',2-2*limit.real>0,'floating_preview')
  for R in (4,8,16):
   C=simpson(lambda q:exp(-q*q)/sqrt(pi)*cmath.exp(-1j*tau*q*q/(sqrt(m*m+q*q/(R*R))+m)),-10,10)
   bound=3*tau/(32*m**3*R**2)
   check(f'autocorrelation remainder m={m},tau={tau},R={R}',abs(C-limit)<=bound+1e-11,'floating_preview')
# Direct distorted-wave packet quadrature, with stable rationalized dispersion.
m=sqrt(16/3)
for tau in (0.0,0.5,2.0):
 for x in (0.0,0.75,2.0):
  K=tanh(x);r0=(3*K*K-1)/2
  A=pi**(-0.25)*(1+1j*tau/m)**(-0.5)
  for R in (8,16,32):
   def integrand(q):
    k=q/R
    symbol=(3*K*K-1-k*k-3j*k*K)/sqrt((1+k*k)*(4+k*k))
    phase=q*x/R-tau*q*q/(sqrt(m*m+k*k)+m)
    return pi**(-0.25)/sqrt(2*pi)*exp(-q*q/2)*cmath.exp(1j*phase)*symbol
   val=simpson(integrand,-10,10)
   err=abs(val-r0*A)
   rows.append({'mass':m,'tau':tau,'x':x,'R':R,'scaled_field_real':val.real,'scaled_field_imag':val.imag,'error_scaled_field':err,'R_squared_error':R*R*err})
   check(f'local packet preview tau={tau},x={x},R={R}',R*R*err<10,'floating_preview')
  a,b,c=rows[-3:]
  # A finite convergence preview; the analytical proof covers the uniform estimate.
  check(f'local error decreases tau={tau},x={x}',c['error_scaled_field']<b['error_scaled_field']<a['error_scaled_field'],'floating_preview')
cq=2*asinh(1)
for m in (sqrt(5.25),sqrt(16/3),sqrt(5.75)):
 S=2.0
 val=simpson(lambda t:abs(pi**(-0.25)*(1+1j*t/m)**(-0.5))**2,0,S)
 expected=m/sqrt(pi)*asinh(S/m)
 check(f'trace coefficient integral m={m}',abs(val-expected)<1e-11,'floating_preview')
 for R in (4,8,16):
  realval=simpson(lambda tau:(cmath.exp(-1j*m*R*R*tau)*pi**(-0.25)*(1+1j*tau/m)**(-0.5)).real**2,0,S,n=24000)
  check(f'real carrier averaging m={m},R={R}',abs(realval-expected/2)*R*R<2,'floating_preview')
with (HERE/'packet_preview.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
 data=(ROOT/row['path']).read_bytes()
 check('sha256 '+row['path'],hashlib.sha256(data).hexdigest()==row['sha256'],'provenance')
 check('git blob '+row['path'],hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob_sha'],'provenance')
result={'status':'PASS_LINEAR_PACKET_ASYMPTOTICS_ACTUAL_NONLINEAR_HISTORY_OPEN' if all(c['passed'] for c in checks) else 'FAIL','counts':{k:sum(c['kind']==k for c in checks) for k in sorted({c['kind'] for c in checks})},'passed':sum(c['passed'] for c in checks),'checks':checks,'actual_critical_damping_proved':False,'formal_verification':False,'independent_review':False,'nonlinear_pde_simulation':False}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
for c in checks:
 if not c['passed']:print('FAILED:',c['name'])
raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
