#!/usr/bin/env python3
import numpy as np, math
from scipy.special import softmax
from scipy.optimize import root

N=12
lam={3:1.96140686197644,4:2.19956884933321,
     5:2.2986062720790903,6:2.3421820411463}
j=np.arange(N)
cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/N),
        np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/N)
    ]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols)
Xp=np.linalg.pinv(X)

def grad(th,g):
    p=softmax(X@th)
    return th/g-X.T@p

def hess(th,g):
    p=softmax(X@th);mu=p@X
    return np.eye(7)/g-(X.T@(p[:,None]*X)-np.outer(mu,mu))

def phi(th,g):
    h=X@th;m=h.max()
    return .5*th@th/g-(m+np.log(np.exp(h-m).mean()))

C=[0,2,4,6]
def fold_eq(z):
    th=np.zeros(7);th[C]=z[:4];g=z[4]
    H=hess(th,g)[np.ix_(C,C)]
    return np.r_[grad(th,g)[C],np.linalg.det(H)]

# Upper fold.
zu=np.array([0.064551183235,0.005557251344,
             0.013192678416,0.418079565490,
             5.172231474684])
ru=root(fold_eq,zu,tol=1e-12)
assert np.linalg.norm(fold_eq(ru.x))<1e-10
assert abs(ru.x[4]-5.172231474684087)<2e-10

# Second signed-C4 fold.
zs=np.array([-1.330174026970,-0.651860029116,
              0.720774033303,1.069731466438,
              4.395526393543])
rs=root(fold_eq,zs,tol=1e-12)
assert np.linalg.norm(fold_eq(rs.x))<1e-10
assert abs(rs.x[4]-4.395526393543094)<2e-10

# g=5 branch roots from continuation.
roots={
"main_min":np.array([2.806464835949445,0,2.966968600794535,0,
                     3.018621651906571,0,2.154961049804607]),
"main_saddle":np.array([0.145841033987508,0,0.211819914533706,0,
                        0.245786571160258,0,0.207404551957227]),
"second_idx1":np.array([2.432480800965357,0,0.612210735571551,
                        1.060380098949041,0.644762114474649,
                       -1.116760741065632,1.871849732748578]),
"second_idx2":np.array([0.271673875058584,0,0.083546704080915,
                        0.144707136273071,0.121216308472001,
                       -0.209952804979437,0.396436330674618])
}
for th in roots.values():
    assert np.linalg.norm(grad(th,5.0))<5e-12

indices={k:int(np.sum(np.linalg.eigvalsh(hess(v,5.0))<-1e-8))
         for k,v in roots.items()}
assert indices=={"main_min":0,"main_saddle":1,
                 "second_idx1":1,"second_idx2":2}

atlas={
"main_min":np.array([0,2.806464835949382,2.9669686007946012,0,
                     0,-3.018621651906552,-2.154961049804644]),
"main_saddle":np.array([0,0.145841033987496,-0.105909957266844,
                        0.183441427013655,-0.212857414533923,
                        0.122893285579990,-0.207404551957113]),
"second_idx1":np.array([0,-2.432480800965327,-1.224421471143071,0,
                        0,-1.289524228949267,-1.871849732748546]),
"second_idx2":np.array([-0.271673875058560,0,-0.167093408161801,0,
                        0.242432616944010,0,0.396436330674561])
}

def transform(th,a,eps):
    h=X@th
    hp=np.array([h[(eps*q+a)%N] for q in range(N)])
    return Xp@hp

dist={}
for name,th in roots.items():
    best=1e9
    for a in range(N):
        for eps in (1,-1):
            best=min(best,np.linalg.norm(transform(th,a,eps)-atlas[name]))
    dist[name]=best
    assert best<3e-12

energies={k:phi(v,5.0) for k,v in roots.items()}

print("PASS")
print("upper_fold",ru.x.tolist())
print("second_fold",rs.x.tolist())
print("g5_indices",indices)
print("g5_energies",energies)
print("D12_match_distances",dist)
