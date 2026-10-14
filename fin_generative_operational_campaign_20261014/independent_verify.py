#!/usr/bin/env python3
"""Independent small verifier: no imports from source_model.py or run_operational.py.
Rebuilds FIN via Laplacian eigenspectrum -> X X^T, N2 Gibbs law, N2 dense generators,
and two-time joint laws by scipy.linalg.expm rather than producer expm_multiply.
Does not certify floating roundoff or actual measurement capability.
"""
import json,math
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.special import softmax
ROOT=Path(__file__).resolve().parent
q=12
x=np.arange(q);omega=.18575;phi=.1625;eta=1.8
W=np.array([[0 if a==b else math.cos(omega*min(abs(a-b),q-abs(a-b))+phi)/(1+min(abs(a-b),q-abs(a-b))**eta) for b in x] for a in x])
lap=np.diag(W.sum(axis=1))-W
lam=np.fft.fft(lap[0]).real
v=[]
for k in (3,4,5):
  v.append(np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*x/q))
  v.append(np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*x/q))
v.append(np.sqrt(lam[6]/12)*(-1.)**x)
X=np.array(v).T;FIN=X@X.T;tr=np.trace(FIN)
Potts=(np.eye(q)-np.ones((q,q))/q)*tr/11
flat=np.zeros((q,q))
for k in (3,4,5):
  vec_cos=np.cos(2*np.pi*k*x/q);vec_sin=np.sin(2*np.pi*k*x/q)
  flat+=(tr/7)/6*(np.outer(vec_cos,vec_cos)+np.outer(vec_sin,vec_sin))
flat+=(tr/7)/12*np.outer((-1.)**x,(-1.)**x)
models={'FIN_A7':FIN,'FULL_POTTS_TRACE':Potts,'FLAT_P7_TRACE':flat}
source=json.loads((ROOT/'RG-E02/SHAPE_PREDICTIONS.json').read_text())
kin=json.loads((ROOT/'RG-E04/MULTITIME_PREDICTION_DESIGN.json').read_text())['N2_two_time']['joint_13x13_by_model_kinetic']
cert=json.loads((ROOT/'RG-E03/STATIC_INTERVAL_CERTIFICATE.json').read_text())
checks=[]
def check(name,ok,details=None):
  checks.append({'name':name,'pass':bool(ok),'detail':details})

def n2_L(A,kind,g=3):
 L=np.zeros((q*q,q*q))
 for ix in range(q*q):
  a,b=divmod(ix,q)
  for site in (0,1):
   old,other=(a,b) if site==0 else (b,a)
   field=g/2*A[:,other]
   hb=softmax(field)
   for new in range(q):
    if new==old:continue
    jx=q*new+b if site==0 else q*a+new
    if kind=='heat_bath':rate=hb[new]
    else:
      delta=field[new]-field[old]
      if kind=='metropolis':rate=(1./q)*min(1.,math.exp(delta))
      else:rate=(2./q)/(1+math.exp(-delta))
    L[ix,jx]+=rate;L[ix,ix]-=rate
 return L

def pairp(A,g=3):
 p=np.exp(g/2*A[0]-max(g/2*A[0]));return p/p.sum()

for key,A in models.items():
 check(key+': symmetric',np.max(abs(A-A.T))<1e-13)
 check(key+': equal_trace',abs(np.trace(A)-tr)<1e-13)
 reference=np.array(source['observable_laws'][key]['N2_diff_g3'])
 check(key+': N2 reproduction',np.max(abs(pairp(A)-reference))<1e-13)
 if key=='FIN_A7':
  check('FIN: source Fourier spectrum',np.max(abs(np.fft.fft(A[0]).real[[3,4,5,6]]-np.array([1.9614068619764458,2.199568849333211,2.298606272079097,2.3421820411462977])))<1e-12)
 d=np.fromfunction(lambda i,j:(j-i)%q,(q,q),dtype=int).ravel()
 p0=np.ones(q*q)/q**2
 g0hb=n2_L(A,'heat_bath',g=0)
 for method in ['heat_bath','metropolis','barker']:
  L=n2_L(A,method,g=3)
  pstat=np.exp(1.5*A.ravel()-np.max(1.5*A));pstat/=pstat.sum()
  check(key+'/'+method+': generator normalization',np.max(abs(L.sum(axis=1)))<1e-13)
  check(key+'/'+method+': detailed balance',np.max(abs(pstat[:,None]*L-(pstat[:,None]*L).T))<1e-13)
  check(key+'/'+method+': g0 calibration',np.max(abs(n2_L(A,method,0)-g0hb))<1e-13)
  p1=expm(.4*L.T)@p0;T=expm(.6*L.T)
  actual=np.zeros((q+1,q+1))
  for first in range(q):
    pfirst=T@(p1*(d==first))
    for second in range(q):actual[first,second]=pfirst[d==second].sum()
  recorded=next(row['joint_13x13'] for row in kin[key] if row['kinetic']==method)
  diff=np.max(np.abs(actual-np.asarray(recorded)))
  check(key+'/'+method+': independent dense two-time law',diff<2e-12,{'max_abs':float(diff)})
  check(key+'/'+method+': nonnegative joint',min(actual.ravel())>=-1e-14 and abs(actual.sum()-1)<1e-13)

for i,x in enumerate(models):
 for y in list(models)[i+1:]:
  observed=float(np.abs(pairp(models[x])-pairp(models[y])).sum()/2)
  lower,upper=[float(v) for v in cert['pairwise_tv_public_outward'][f'{x}_vs_{y}']]
  check('interval_static_'+x+'_'+y,lower<=observed<=upper,{'double_tv':observed,'iv_public_interval':[lower,upper]})
positive=sum(c['pass'] for c in checks)
result={'verifier':'independent dense N2 generator + X7 X7^T FIN + source interval brackets',
        'checks_total':len(checks),'checks_passed':positive,'all_passed':positive==len(checks),
        'limits':['SciPy matrix exponentials double precision, not fully outward certified',
                  'interval assumptions inherited from mpmath.iv, separately authored static certificate',
                  'no physical records or time-custody proof'], 'checks':checks}
(ROOT/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'passed':positive,'total':len(checks),'all_passed':result['all_passed']}))
if not result['all_passed']:raise SystemExit(1)
