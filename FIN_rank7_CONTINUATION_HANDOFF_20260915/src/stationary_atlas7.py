from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.stats import qmc
from model import feature_spaces, dual7, d12_actions
ROOT=Path(__file__).resolve().parents[1]

W,A,L,X7,C4,A7=feature_spaces(); ACTIONS=d12_actions(X7)
ROWMAX=float(np.max(np.linalg.norm(X7,axis=1)))

def solve(theta0,g):
    theta0=np.asarray(theta0,float)
    fun=lambda t: dual7(t,g,X7)[1]
    jac=lambda t: dual7(t,g,X7)[2]
    sol=root(fun,theta0,jac=jac,method='hybr',options={'xtol':1e-11,'maxfev':1000})
    th=np.asarray(sol.x,float); val,gr,H,p=dual7(th,g,X7)
    return {'success':bool(sol.success and np.linalg.norm(gr,np.inf)<1e-8),
            'theta':th,'phi':float(val),'residual_inf':float(np.linalg.norm(gr,np.inf)),
            'H_eigs':np.linalg.eigvalsh(H),'index':int(np.sum(np.linalg.eigvalsh(H)<-1e-7)),
            'stabilizer':None,'message':str(sol.message)}

def orbit_distance(a,b):
    best=1e99
    for P,T in ACTIONS.values(): best=min(best,float(np.linalg.norm(T@a-b)))
    return best

def stabilizer(theta,tol=2e-6):
    return sum(np.linalg.norm(T@theta-theta)<tol for P,T in ACTIONS.values())

def dedup(records,tol=2e-5):
    reps=[]
    for r in sorted((x for x in records if x['success']), key=lambda x:(x['phi'],np.linalg.norm(x['theta']))):
        if not any(orbit_distance(r['theta'],q['theta'])<tol for q in reps):
            r=dict(r); r['stabilizer']=stabilizer(r['theta']); reps.append(r)
    return reps

def seeds(g,batch):
    R=g*ROWMAX
    S=[np.zeros(7)]
    # coordinate axes at three scales, both signs
    for k in range(7):
        for frac in (.2,.55,.9):
            for sg in (-1,1):
                v=np.zeros(7); v[k]=sg*frac*R; S.append(v)
    # symmetric low-mode combinations
    for k in range(0,6,2):
        v=np.zeros(7); v[k]=.55*R; v[6]=.35*R; S += [v,-v]
    n=128 if batch==1 else 256
    sob=qmc.Sobol(d=7,scramble=True,seed=9137+batch).random_base2(int(math.log2(n)))
    # uniform in ball direction with deterministic radius profile
    Z=2*sob-1; norms=np.linalg.norm(Z,axis=1); norms[norms==0]=1
    rad=((np.arange(len(Z))+.5)/len(Z))**(1/7)*R
    S += list(Z/norms[:,None]*rad[:,None])
    rng=np.random.default_rng(15000+batch)
    for _ in range(96 if batch==1 else 160):
        z=rng.normal(size=7); z/=np.linalg.norm(z); z*=R*rng.random()**(1/7); S.append(z)
    return S

def run_gain(g):
    batches=[]; cumulative=[]; old=0
    for b in (1,2):
        rec=[solve(s,g) for s in seeds(g,b)]
        cumulative.extend(rec); reps=dedup(cumulative)
        batches.append({'batch':b,'seed_count':len(rec),'solve_successes':sum(x['success'] for x in rec),
                        'failed_solves':sum(not x['success'] for x in rec),'orbit_count':len(reps),
                        'new_orbits_vs_previous':len(reps)-old})
        old=len(reps)
    reps=dedup(cumulative)
    out=[]
    for r in reps:
        out.append({'theta':[float(x) for x in r['theta']], 'phi':r['phi'], 'residual_inf':r['residual_inf'],
                    'H_eigenvalues':[float(x) for x in r['H_eigs']], 'H7_index':r['index'],
                    'stabilizer_size':int(r['stabilizer']),'orbit_size':int(24//r['stabilizer'])})
    return {'g':g,'stationary_ball_radius':g*ROWMAX,'batches':batches,'orbits':out,
            'claim_scope':'numerical discovery saturation only; complement not covered'}

def main():
    res={'id':'R7P-037','status':'NUMERICAL_SATURATION_NOT_EXHAUSTION','gains':[run_gain(g) for g in (3.7,3.7183449,4.0,5.0)],
         'dedup_tolerance':2e-5,'group':'D12, 24 actions','nonconclusion':'No stationary exhaustion theorem.'}
    p=ROOT/'results/R7P-037_full7_stationary_atlas.json'; p.write_text(json.dumps(res,indent=2)+'\n')
    print(json.dumps({'counts':[(x['g'],len(x['orbits']),x['batches'][-1]['new_orbits_vs_previous']) for x in res['gains']]},indent=2))
if __name__=='__main__': main()
