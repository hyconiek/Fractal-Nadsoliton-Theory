"""R7P-089 numerical quartic phase-root discovery and catalog."""
from __future__ import annotations
import json,math
from pathlib import Path
import numpy as np
from scipy.optimize import root
from .phase_cumulants import k4_phase_value_grad_hess,theta_from_z,z_from_theta
from .model import feature_spaces,d12_actions
FIXTURE=(0.1131879146,0.1698528641,0.2269339093,-0.3380663037)

def wrap(x): return np.mod(np.asarray(x,float),2*math.pi)
def torus_delta(a,b):
    d=np.abs(np.asarray(a)-np.asarray(b)); return np.minimum(d,2*math.pi-d)
def torus_dist(a,b): return float(np.linalg.norm(torus_delta(a,b)))

def discover(nstarts=4000,seed=89089,tol=1e-7):
    rng=np.random.default_rng(seed); roots=[]; failures=[]
    args=FIXTURE
    for n in range(nstarts):
        x0=rng.uniform(0,2*math.pi,3)
        sol=root(lambda x:k4_phase_value_grad_hess(*args,x)[1],x0,
                 jac=lambda x:k4_phase_value_grad_hess(*args,x)[2],method='hybr',tol=1e-11)
        x=wrap(sol.x); val,g,H=k4_phase_value_grad_hess(*args,x); res=float(np.linalg.norm(g))
        if res>1e-8:
            failures.append({'start':x0.tolist(),'residual':res,'solver_success':bool(sol.success)})
            continue
        if not any(torus_dist(x,y)<tol for y in roots): roots.append(x)
    roots=sorted(roots,key=lambda x:tuple(np.round(x,12)))
    records=[]
    for i,x in enumerate(roots):
        val,g,H=k4_phase_value_grad_hess(*args,x); eig=np.linalg.eigvalsh(H)
        records.append({'id':i,'phase':x.tolist(),'K4':val,'residual_norm':float(np.linalg.norm(g)),
                        'hessian_eigenvalues':eig.tolist(),'negative_index':int(np.sum(eig<0))})
    # Symmetry graph within the fixed negative-z6 sign. D12 actions with transformed z6<0.
    _,_,L,X,_,_=feature_spaces(); acts=d12_actions(X); r3,r4,r5,z6=args
    for rec in records:
        ph=np.array(rec['phase']); th=theta_from_z(r3*np.exp(1j*ph[0]),r4*np.exp(1j*ph[1]),r5*np.exp(1j*ph[2]),z6,L=L)
        targets=set(); stab=0
        for key,(P,T) in acts.items():
            zs=z_from_theta(T@th,L=L)
            if zs[3]>=0: continue
            ph2=wrap([np.angle(zs[0]),np.angle(zs[1]),np.angle(zs[2])])
            ds=[torus_dist(ph2,r['phase']) for r in records]; j=int(np.argmin(ds))
            if ds[j]<5e-7: targets.add(j)
            if torus_dist(ph2,ph)<5e-7: stab+=1
        rec['fixed_sign_symmetry_orbit_ids']=sorted(targets); rec['fixed_sign_stabilizer_size']=stab
    return {'fixture':{'r3':args[0],'r4':args[1],'r5':args[2],'z6':args[3]},'seed':seed,'starts':nstarts,
            'dedup_tolerance':tol,'root_count':len(records),'failure_count':len(failures),'roots':records,'failures':failures}

def main(path):
    out=discover(); Path(path).write_text(json.dumps(out,indent=2)+'\n')
    from collections import Counter
    print(json.dumps({'root_count':out['root_count'],'failure_count':out['failure_count'],
      'index_counts':dict(Counter(str(r['negative_index']) for r in out['roots'])),
      'min_abs_hessian_eigenvalue':min(min(abs(x) for x in r['hessian_eigenvalues']) for r in out['roots'])},indent=2))
if __name__=='__main__':
    import sys; main(sys.argv[1] if len(sys.argv)>1 else 'R7P-089_quartic_roots.json')
