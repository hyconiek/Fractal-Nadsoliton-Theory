import sys, math, json, numpy as np
from scipy.optimize import differential_evolution
try:
    from . import off_face, boundary_ising as bi
    from .intervals import QI
except ImportError:
    import off_face, boundary_ising as bi
    from intervals import QI

L,sigma,_,_,C4=off_face.constants()
sg=bi.interval_signs(); T=QI(sg['t'][0],sg['t'][1]); R=(QI(1)-T)/(QI(1)+T); rstar=(float(R.lo)+float(R.hi))/2
rmin=1/900; smin=1/128; tmin=1/9; ymin=1/5_000_000
boxes=[
('FR8',1/6500,1/8192,1/4600,1/100000),('FR9',1/6500,1/2432,1/8192,1/100000),('FR15',1/6500,1/2432,1/8192,1/14000),('FR16',1/5000,1/8192,1/5600,1/100000),
('FR10',1/8192,1/4800,1/4800,1/3072),('FR11',1/7600,1/4800,1/4800,1/81920),
('FR12',1/4096,1/8192,1/8192,1/1024),('FR13',1/6200,1/6800,1/6800,1/400),
('FR14',1/6400,1/8192,1/4600,1/100000)]

def rec(x):
    r,s,t,logy=x; y=10**logy
    p=off_face.p_from_aligned_compact(r,s,t,y)
    mu=p@C4; Y=C4-mu; M=Y.T@(p[:,None]*Y)
    eig=np.linalg.eigvalsh(M); gap=float(eig[-2]-sigma)
    e=float(p[1::2].sum()); u=1-s; v=1-t
    reasons=[]
    for name,rx,ru,rv,re in boxes:
        if abs(r-rstar)<=rx+1e-12 and u<=ru+1e-12 and v<=rv+1e-12 and e<=re+1e-12: reasons.append(name)
    return gap,e,p,eig,reasons,y

def obj(x):
    gap,e,p,eig,reasons,y=rec(x)
    if reasons: return 100 + max(0,gap) + len(reasons)*.01
    return -gap

bounds=[(rmin,1),(smin,1),(tmin,1),(math.log10(ymin),0)]
rows=[]
for seed in (150927,150928,150929,150930):
    res=differential_evolution(obj,bounds,seed=seed,maxiter=350,popsize=18,tol=1e-10,polish=True,workers=1,updating='immediate')
    gap,e,p,eig,reasons,y=rec(res.x)
    r,s,t,logy=res.x
    rows.append(dict(seed=seed,success=bool(res.success),message=str(res.message),fun=float(res.fun),r=float(r),s=float(s),t=float(t),y=float(y),
                     J3=float(-.5*math.log(r)),J4=float(-(2/3)*math.log(s)),J5=float(-2*math.log(t)),J6=float(-.5*math.log(y)),
                     u=float(1-s),v=float(1-t),e=float(e),gap=gap,lambda2=float(eig[-2]),eigenvalues=eig.tolist(),p=p.tolist(),excluded_by=reasons,nfev=int(res.nfev)))
    print(seed,gap,'r,u,v,e,y',r,1-s,1-t,e,y,'reasons',reasons)
print('best',max(rows,key=lambda z:z['gap']))
open(str(__import__('pathlib').Path(__file__).resolve().parents[1]/'results/FR17_residual_adversarial_search.json'),'w').write(json.dumps({'id':'FR17-compact-residual-search','status':'NUMERICAL_SEARCH_ONLY','sigma':sigma,'domain':{'r':[rmin,1],'s':[smin,1],'t':[tmin,1],'y':[ymin,1]},'local_boxes':[dict(name=n,rx=rx,ru=ru,rv=rv,re=re) for n,rx,ru,rv,re in boxes],'rstar':rstar,'interpretation':'All saved maximizers have negative physical gap; the best lies numerically on an artificial certified-box wall and is navigation evidence only, not a theorem.', 'solver':{'method':'scipy differential_evolution','maxiter':350,'popsize':18,'tol':1e-10,'polish':True},'runs':rows,'best':max(rows,key=lambda z:z['gap'])},indent=2)+'\n')
