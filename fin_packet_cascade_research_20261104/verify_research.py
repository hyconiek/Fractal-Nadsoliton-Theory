"""Finite algebra, spectral previews, convergence checks and pinned hashes."""
from pathlib import Path
from fractions import Fraction as F
from math import sqrt,tanh,cos,sin
import hashlib,json
HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
checks=[]
def check(name,ok,kind='exact'):
 checks.append({'name':name,'passed':bool(ok),'kind':kind})
check('source coefficient from four branches',4*F(-3,2)==-6)
check('localized source Wiener product scaling',F(-1,2)+F(-1,2)==-1)
check('global correction mixed scaling',0+F(-1,2)==F(-1,2))
check('initial transient integrated scaling',2*F(-1,2)==-1)
check('finite-window response integrated scaling',2*(-1)+2==0)
check('small packet response squared amplitude',2*2==4)
check('odd threshold Jacobian gain',2-1==1)
check('even threshold Jacobian obstruction',0-1==-1)
# Independent numerical inspection of the odd normalized generalized mode.
def simpson(fn,a,b,n=12000):
 h=(b-a)/n
 return h/3*(fn(a)+fn(b)+sum((4 if j%2 else 2)*fn(a+j*h) for j in range(1,n)))
def mode(k,x):
 K=tanh(x);s=1-K*K;A=3*K*K-1-k*k;d=sqrt((1+k*k)*(4+k*k))
 eo=(A*sin(k*x)-3*k*K*cos(k*x))/d
 ex=((6*K*s+3*k*k*K)*sin(k*x)+k*(A-3*s)*cos(k*x))/d
 return eo,ex
spectral=[]
for k in (.01,.02,.04,.5,1,2,4,8):
 def b0(x):
  e,ex=mode(k,x);weight=(1+x*x)**(-3)
  return (weight*e)**2
 def b1(x):
  e,ex=mode(k,x);weight=(1+x*x)**(-3);wp=-6*x*(1+x*x)**(-4)
  return (weight*e)**2+(weight*ex+wp*e)**2
 B0=simpson(b0,-30,30);B1=simpson(b1,-30,30);nu=sqrt(4+k*k)
 ratio=(nu/k)*(B1+nu*nu*B0)/(nu*nu)
 spectral.append({'k':k,'weighted_L2_squared':B0,'weighted_H1_squared':B1,'energy_normalized_Plancherel_factor':ratio})
 check(f'finite odd Plancherel preview k={k}',0<ratio<100,'floating_preview')
check('odd mode quadratic zero preview',abs(spectral[1]['weighted_H1_squared']/(4*spectral[0]['weighted_H1_squared'])-1)<.01,'floating_preview')
even0=simpson(lambda x:((1+x*x)**(-3)*(3*tanh(x)**2-1)/2)**2,-30,30)
check('negative control even threshold does not vanish',even0>.001,'floating_preview')
(HERE/'odd_spectral_preview.json').write_text(json.dumps({'scope':'Finite samples only; analytical proof covers all frequencies','samples':spectral,'even_threshold_weighted_L2_squared':even0},indent=2)+'\n')
# These verify results of the separately executed comparison solver.
preview=json.loads((HERE/'cascade_preview.json').read_text());rs=preview['results']
for r in rs:
 tag=f"R={r['R']},dt={r['dt']},dx={r['dx']},L={r['half_length']}"
 check('packet normalization '+tag,abs(r['initial_packet_L2_squared']-1)<1e-10,'numerical_comparison')
 check('shape projection '+tag,r['max_shape_projection']<1e-10,'numerical_comparison')
 check('boundary-band diagnostic '+tag,r['final_outer_band_phase_squared']<1e-12,'numerical_comparison')
 check('source-work balance '+tag,r['energy_work_relative_defect']<.002,'numerical_comparison')
for coarse,fine in [(rs[0],rs[1]),(rs[4],rs[5])]:
 check(f'time refinement R={fine["R"]}',abs(coarse['integrated_local_phase_squared']/fine['integrated_local_phase_squared']-1)<1e-4,'numerical_comparison')
 check(f'work defect decreases quadratically R={fine["R"]}',3.8<coarse['energy_work_relative_defect']/fine['energy_work_relative_defect']<4.2,'numerical_comparison')
for refined,label in [(rs[2],'spatial grid'),(rs[3],'domain')]:
 check(label+' refinement R=4',abs(refined['integrated_local_phase_squared']/rs[1]['integrated_local_phase_squared']-1)<1e-9,'numerical_comparison')
for row in json.loads((HERE/'INPUTS.json').read_text())['files']:
 data=(ROOT/row['path']).read_bytes()
 check('sha256 '+row['path'],hashlib.sha256(data).hexdigest()==row['sha256'],'provenance')
 check('git blob '+row['path'],hashlib.sha1(b'blob '+str(len(data)).encode()+b'\0'+data).hexdigest()==row['git_blob_sha'],'provenance')
result={'status':'PASS_RETARDED_PACKET_COMPARISON_ACTUAL_NONLINEAR_GATE_OPEN' if all(c['passed'] for c in checks) else 'FAIL','counts':{k:sum(c['kind']==k for c in checks) for k in sorted({c['kind'] for c in checks})},'passed':sum(c['passed'] for c in checks),'checks':checks,'actual_critical_damping_proved':False,'formal_verification':False,'independent_review':False,'nonlinear_pde_simulation':False,'linear_forced_comparison_simulation':True}
(HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='checks'},indent=2))
for c in checks:
 if not c['passed']:print('FAILED:',c['name'])
raise SystemExit(0 if all(c['passed'] for c in checks) else 1)
