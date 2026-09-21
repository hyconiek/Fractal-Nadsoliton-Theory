from pathlib import Path
import sys,json,math
import numpy as np
from scipy.stats import qmc
from scipy.optimize import differential_evolution
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916';sys.path.insert(0,str(H/'src'))
import off_face
TAU=67/250
bounds=np.array([[1/900,1],[1/128,1],[1/9,1],[1e-6,1]],float)

def evalpt(x):
    M,p=off_face.covariance_direct_compact(x); eig=np.linalg.eigvalsh(M)
    return float(eig[-2]),eig,p

def map_u(u):
    # log-uniform navigation across many decades, exact hull endpoints retained separately
    lo=np.log(bounds[:,0]); hi=np.log(bounds[:,1]); return np.exp(lo+u*(hi-lo))

sam=qmc.Sobol(d=4,scramble=False).random_base2(14)
pts=map_u(sam)
# Add coordinate faces/endpoints using 2048 shared 3d Sobol points each.
faces=[]; base=qmc.Sobol(d=3,scramble=False).random_base2(11)
for ax in range(4):
  others=[i for i in range(4) if i!=ax]
  for side in [0,1]:
    arr=np.empty((len(base),4)); arr[:,ax]=bounds[ax,side]
    for j,i in enumerate(others):
      arr[:,i]=np.exp(np.log(bounds[i,0])+base[:,j]*(np.log(bounds[i,1])-np.log(bounds[i,0])))
    faces.append(arr)
pts=np.vstack([pts]+faces)
best=(-1,None,None)
for x in pts:
    l,e,p=evalpt(x)
    if l>best[0]:best=(l,x.copy(),e.copy())

def obj(x):return -evalpt(x)[0]
de=[]
for seed in [18001,18002]:
    res=differential_evolution(obj,[tuple(z) for z in bounds],seed=seed,maxiter=100,popsize=12,tol=1e-9,polish=True,workers=1,updating='immediate')
    l,e,p=evalpt(res.x);de.append({'seed':seed,'lambda2':l,'tau_gap':l-TAU,'x':res.x.tolist(),'eigenvalues':e.tolist(),'nfev':int(res.nfev),'success':bool(res.success)})
    if l>best[0]:best=(l,res.x.copy(),e.copy())
out={'task':'R7N-018','scientific_status':'NUMERICAL_NEW','threshold_tau0':TAU,'sample_count':int(len(pts)),
     'sample_design':'16384 deterministic Sobol log-hull + 8 coordinate faces x 2048 + 2 differential-evolution runs',
     'best_lambda2':best[0],'best_tau_gap':best[0]-TAU,'best_point_rsty':best[1].tolist(),'best_eigenvalues':best[2].tolist(),
     'DE_runs':de,'interpretation':'Navigation only. A negative sampled/optimized tau gap is not a global proof; a positive stable candidate would require certification.'}
json.dump(out,open(ROOT/'results/R7N-018_target_p_navigation.json','w'),indent=2);print(json.dumps(out,indent=2))
