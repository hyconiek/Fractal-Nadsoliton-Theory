#!/usr/bin/env python3
import numpy as np, math
from scipy.special import softmax
from scipy.optimize import root, brentq
from scipy.linalg import expm

# --- strict rank-7 data
Nlab=12
lam={3:1.96140686197644,4:2.19956884933321,
     5:2.2986062720790903,6:2.3421820411463}
j=np.arange(Nlab)
cols=[]
for k in (3,4,5):
    cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Nlab),
             np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Nlab)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols)
A=X@X.T

def grad(th,g):
    p=softmax(X@th)
    return th/g-X.T@p

def hess(th,g):
    p=softmax(X@th);mu=p@X
    return np.eye(7)/g-(X.T@(p[:,None]*X)-np.outer(mu,mu))

def phi(th,g):
    h=X@th;m=h.max()
    return .5*th@th/g-(m+np.log(np.exp(h-m).mean()))

# --- finite-N drift expansion check
rng=np.random.default_rng(42)
p=rng.dirichlet(np.ones(12)*2)
g=5.0
def loo_drift(p,g,N):
    f=g*A@p
    z=np.zeros(12)
    for i in range(12):
        z+=p[i]*softmax(f-(g/N)*A[:,i])
    return z-p
def approx(p,g,N):
    q=softmax(g*A@p)
    S=np.diag(q)-np.outer(q,q)
    return q-p-(g/N)*S@(A@p)
scaled=[]
for N in [20,40,80,160]:
    e=np.linalg.norm(loo_drift(p,g,N)-approx(p,g,N))
    scaled.append(e*N*N)
assert max(scaled)-min(scaled)<5e-4

# --- localized branch odd block
C=[0,2,4,6]; O=[1,3,5]
th5=np.array([2.806464835949445,0,2.966968600794535,0,
              3.018621651906571,0,2.154961049804607])
def solveC(g,z):
    def F(x):
        t=np.zeros(7);t[C]=x
        return grad(t,g)[C]
    s=root(F,z,tol=1e-12)
    t=np.zeros(7);t[C]=s.x
    return t
z=th5[C].copy()
min_scaled=1e9
for gg in np.linspace(5,3.51565,400):
    t=solveC(float(gg),z);z=t[C]
    eo=np.linalg.eigvalsh(hess(t,gg)[np.ix_(O,O)])
    min_scaled=min(min_scaled,float(gg*eo[0]))
assert min_scaled>0.63

# high-g replay
z=th5[C].copy()
for gg in [10,20,50,100,1000]:
    t=solveC(float(gg),z);z=t[C]
    eo=np.linalg.eigvalsh(hess(t,gg)[np.ix_(O,O)])
    assert eo[0]>0

# --- main saddle zero-energy crossing and k6 balance
ths5=np.array([0.145841033987508,0,0.211819914533706,0,
               0.245786571160258,0,0.207404551957227])
zs=ths5[C].copy()
def saddle(gg):
    return solveC(gg,zs)
l6=lam[6]; A6=np.sqrt(l6/12); g6=12/l6
def k6J(gg):
    c=gg*l6/12
    return brentq(lambda J:J-c*np.tanh(J),1e-12,5)
def k6th(gg):
    t=np.zeros(7);t[6]=k6J(gg)/A6;return t

gx=brentq(lambda gg:phi(saddle(gg),gg),5.15,5.151,xtol=1e-13)
assert abs(gx-5.150374750802444)<2e-10

gbal=brentq(lambda gg:-phi(k6th(gg),gg)-phi(saddle(gg),gg),
            5.13,5.15,xtol=1e-13)
assert abs(gbal-5.145228719489144)<2e-10

# --- weak projection counterexample
Jm=np.array([[1.,0.],[1.,0.],[0.,1.]])
Pm=np.array([[.5,.5,0.],[0.,0.,1.]])
Lf=np.array([[1.,0.,-1.],[0.,3.,-3.],[-2.,-2.,4.]])
Lc=Pm@Lf@Jm
assert np.linalg.norm(Pm@Jm-np.eye(2))<1e-14
assert np.linalg.norm(Pm@Lf@Jm-Lc)<1e-14
second=Pm@Lf@Lf@Jm-Lc@Lc
assert np.linalg.norm(second)>0.5
semidiff=np.linalg.norm(Pm@expm(-0.5*Lf)@Jm-expm(-0.5*Lc),2)
assert semidiff>0.05

# --- exact mod-3 lumpability for arbitrary positive rates
rates={1:.07,2:.11,3:2.3,4:.4,5:.13,6:.02}
Q=np.zeros((12,12))
for n in range(12):
    for d in range(1,6):
        k=rates[d]
        for m in {(n+d)%12,(n-d)%12}:
            Q[n,m]+=k
            Q[n,n]-=k
    m=(n+6)%12
    Q[n,m]+=rates[6];Q[n,n]-=rates[6]
# embedding block-constant observables
J3=np.zeros((12,3))
for n in range(12): J3[n,n%3]=1
kap=sum(rates[d] for d in [1,2,4,5])
Q3=kap*np.array([[-2,1,1],[1,-2,1],[1,1,-2]],float)
assert np.linalg.norm(Q@J3-J3@Q3)<1e-13

print("PASS")
print("drift_N2_scaled",scaled)
print("localized_min_scaled_odd_margin",min_scaled)
print("k6_escape_crossover_g",gx)
print("k6_balance_g",gbal)
print("projection_second_order_defect",second.tolist())
print("projection_semigroup_diff_t0p5",semidiff)
print("mod3_kappa",kap)
