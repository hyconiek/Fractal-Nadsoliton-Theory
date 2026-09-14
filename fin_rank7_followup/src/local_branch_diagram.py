"""R7P-032 numerical continuation between certified local landmarks."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
from scipy.optimize import root
from .model import feature_spaces
from .derivatives import dual_all

def solve(theta0,g,C):
 sol=root(lambda x:dual_all(x,g,C)[1],theta0,jac=lambda x:dual_all(x,g,C)[2],tol=1e-11)
 v,gr,H,*_=dual_all(sol.x,g,C)
 return sol.x,float(v),float(np.linalg.norm(gr)),np.linalg.eigvalsh(H),bool(sol.success)

def run(certdir,out_path):
 certdir=Path(certdir)
 eq=json.load(open(certdir/'R7P-026_equal_energy_event.json'))
 sad=json.load(open(certdir/'R7P-029_barrier_saddle.json'))
 fold=json.load(open(certdir/'R7P-031_simple_fold.json'))
 _,_,L,_,C,_=feature_spaces()
 loc0=np.array([float(x) for x in eq['root_center'][:4]])
 sad0=np.array([(float(a)+float(b))/2 for a,b in sad['saddle_root_box']])
 gc=float(eq['root_center'][4]); gf=sum(map(float,fold['fold_gain_interval']))/2
 # Continue downward from crossing; stop before certified fold singularity.
 gs=np.linspace(gc, gf+2e-4, 81)
 branches={}
 for name,th0 in [('localized',loc0),('saddle',sad0)]:
  th=th0.copy(); rows=[]
  for g in gs:
   th,val,res,eig,ok=solve(th,float(g),C)
   rows.append({'g':float(g),'theta':th.tolist(),'Phi':val,'residual':res,'H4_eigenvalues':eig.tolist(),'solver_success':ok})
  branches[name]=rows
 out={'certified_nodes':{'fold_gain_interval':fold['fold_gain_interval'],'equal_energy_root_box':eq['root_box'],
                          'saddle_at_rational_gain_box':sad['saddle_root_box'],
                          'uniform_spinodal':12/float(L[6])},
      'numerical_continuation':branches,
      'continuation_gain_range':[float(gs[-1]),float(gs[0])],
      'joins_certified':False,
      'missing_join':'No overlapping interval continuation tubes were built between the isolated fold and equal-energy boxes; numerical branches are retained only as guidance.',
      'scope':'local reflection-even branch diagram; not exhaustive hysteresis or physical time'}
 Path(out_path).write_text(json.dumps(out,indent=2)+'\n');return out
if __name__=='__main__':
 import sys; o=run(sys.argv[1],sys.argv[2]); print(json.dumps({k:v for k,v in o.items() if k!='numerical_continuation'},indent=2))
