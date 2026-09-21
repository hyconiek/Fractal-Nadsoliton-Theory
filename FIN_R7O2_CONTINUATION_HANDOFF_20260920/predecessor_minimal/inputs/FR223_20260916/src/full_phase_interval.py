"""R7P-092 direct validated isolation of full log-mgf phase roots."""
from __future__ import annotations
import json,math
from pathlib import Path
import numpy as np, mpmath as mp
from scipy.optimize import root
from .phase_cumulants import full_phase_value_grad_hess
FIXTURE=(0.1131879146,0.1698528641,0.2269339093,-0.3380663037)
N=12

def wrap(x): return np.mod(np.asarray(x,float),2*math.pi)
def tdist(a,b):
 d=np.abs(np.asarray(a)-np.asarray(b)); d=np.minimum(d,2*math.pi-d); return float(np.linalg.norm(d))

def solve_from_quartic(cat):
 out=[]
 for qr in cat['roots']:
  x0=np.asarray(qr['phase'],float)
  sol=root(lambda x:full_phase_value_grad_hess(*FIXTURE,x)[1],x0,jac=lambda x:full_phase_value_grad_hess(*FIXTURE,x)[2],tol=1e-12)
  x=wrap(sol.x); v,g,H=full_phase_value_grad_hess(*FIXTURE,x); eig=np.linalg.eigvalsh(H)
  out.append({'quartic_id':qr['id'],'phase':x.tolist(),'residual_norm':float(np.linalg.norm(g)),'hessian_eigenvalues':eig.tolist(),
              'negative_index':int(np.sum(eig<0)),'displacement_from_quartic':tdist(x,x0)})
 return out

def _bounds(x): return float(x.a),float(x.b)
def _FH(vars):
 mp.iv.dps=50
 r3,r4,r5,z6=FIXTURE; amps=[r3,r4,r5]; ks=[3,4,5]
 h=[]; hp=[]; hpp=[]
 for jj in range(N):
  hv=mp.iv.mpf(repr(z6*((-1.0)**jj)/math.sqrt(12.0))); row=[]; row2=[]
  for a,k,ph in zip(amps,ks,vars):
   ang=mp.iv.mpf(repr(2*math.pi*k*jj/N))+ph; aa=mp.iv.mpf(repr(a/math.sqrt(3.0)))
   hv += aa*mp.iv.cos(ang); row.append(-aa*mp.iv.sin(ang)); row2.append(-aa*mp.iv.cos(ang))
  h.append(hv); hp.append(row); hpp.append(row2)
 ew=[mp.iv.exp(x) for x in h]; Z=sum(ew,mp.iv.mpf('0')); p=[x/Z for x in ew]
 g=[sum((p[j]*hp[j][i] for j in range(N)),mp.iv.mpf('0')) for i in range(3)]
 H=[[mp.iv.mpf('0') for _ in range(3)] for __ in range(3)]
 for i in range(3):
  for k in range(3):
   eprod=sum((p[j]*hp[j][i]*hp[j][k] for j in range(N)),mp.iv.mpf('0'))
   H[i][k]=eprod-g[i]*g[k]
   if i==k: H[i][k]+=sum((p[j]*hpp[j][i] for j in range(N)),mp.iv.mpf('0'))
 return g,H

def certify(rec,radius=1e-5):
 c=np.asarray(rec['phase'],float); C=[mp.iv.mpf(repr(float(v))) for v in c]; X=[mp.iv.mpf([repr(float(v-radius)),repr(float(v+radius))]) for v in c]
 Fc,Jc=_FH(C); _,JX=_FH(X); H0=full_phase_value_grad_hess(*FIXTURE,c)[2]; A=np.linalg.inv(H0)
 zero=mp.iv.mpf('0'); one=mp.iv.mpf('1'); D=mp.iv.mpf([repr(-radius),repr(radius)])
 K=[]
 for i in range(3):
  b=zero
  for j in range(3): b += mp.iv.mpf(repr(float(A[i,j])))*Fc[j]
  ss=-b
  for j in range(3):
   Bij=(one if i==j else zero)
   for k in range(3): Bij -= mp.iv.mpf(repr(float(A[i,k])))*JX[k][j]
   ss += Bij*D
  K.append(_bounds(ss))
 maxoff=max(max(abs(a),abs(b)) for a,b in K)
 E=np.zeros((3,3))
 for i in range(3):
  for j in range(3):
   lo,hi=_bounds(JX[i][j]); E[i,j]=max(abs(lo-H0[i,j]),abs(hi-H0[i,j]))
 err=float(np.linalg.norm(E,'fro')); eig=np.linalg.eigvalsh(H0); margin=float(np.min(np.abs(eig))-err)
 return {'quartic_id':rec['quartic_id'],'radius':radius,'inclusion':maxoff<radius,'max_abs_krawczyk_offset':maxoff,
         'hessian_error_bound':err,'inertia_margin':margin,'negative_index':int(np.sum(eig<0)) if margin>0 else None}

def run(qcat_path,out_path,radius=1e-5):
 qcat=json.load(open(qcat_path)); roots=solve_from_quartic(qcat); cert=[certify(r,radius) for r in roots]
 out={'fixture':dict(zip(['r3','r4','r5','z6'],FIXTURE)),'count':len(roots),'roots':roots,'certificates':cert,
      'all_inclusions':all(c['inclusion'] for c in cert),'all_inertia_certified':all(c['inertia_margin']>0 for c in cert),
      'max_displacement':max(r['displacement_from_quartic'] for r in roots),'median_displacement':float(np.median([r['displacement_from_quartic'] for r in roots])),
      'index_changes':sum(r['negative_index']!=q['negative_index'] for r,q in zip(roots,qcat['roots'])),
      'min_abs_full_hessian_eigenvalue':min(min(abs(x) for x in r['hessian_eigenvalues']) for r in roots),
      'max_krawczyk_offset':max(c['max_abs_krawczyk_offset'] for c in cert),'min_inertia_margin':min(c['inertia_margin'] for c in cert)}
 Path(out_path).write_text(json.dumps(out,indent=2)+'\n'); return out
if __name__=='__main__':
 import sys
 o=run(sys.argv[1],sys.argv[2]); print(json.dumps({k:v for k,v in o.items() if k not in ('roots','certificates')},indent=2))
