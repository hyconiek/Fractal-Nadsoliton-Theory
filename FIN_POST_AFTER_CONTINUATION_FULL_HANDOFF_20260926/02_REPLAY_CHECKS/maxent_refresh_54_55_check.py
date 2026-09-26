#!/usr/bin/env python3
import numpy as np
from scipy.special import softmax

rng=np.random.default_rng(20260926)

# ---------- Report 54 ----------
Q=12
q=rng.uniform(.1,2,size=Q); q/=q.sum()
Kreset=np.ones((Q,1))*q[None,:]

def entropy(x):
    x=np.asarray(x)
    return -float(np.sum(x[x>0]*np.log(x[x>0])))

Hq=entropy(q)
Hcond=float(sum(q[i]*entropy(Kreset[i]) for i in range(Q)))
assert abs(Hcond-Hq)<1e-14
assert np.max(np.abs(q@Kreset-q))<1e-14

# Construct a different q-reversible kernel by lazy reset.
a=.37
Klazy=a*np.eye(Q)+(1-a)*Kreset
assert np.max(np.abs(q@Klazy-q))<1e-14
Hlazy=float(sum(q[i]*entropy(Klazy[i]) for i in range(Q)))
assert Hlazy < Hq-1e-5

# ---------- Report 55 ----------
# random zero-mean feature map
d=4
X=rng.normal(size=(Q,d))
X-=X.mean(axis=0)
g=.8
p=rng.uniform(.2,2,size=Q); p/=p.sum()
mu=X.T@p
h=g*X@mu
qt=softmax(h)

# Variational check against many random distributions.
def target_obj(r):
    return entropy(r)+h@r
best=target_obj(qt)
for _ in range(10000):
    r=rng.dirichlet(np.ones(Q))
    assert target_obj(r) <= best+1e-12

# Vg and directional derivative identity.
u=np.ones(Q)/Q
def V(r):
    return float(np.sum(r*np.log(r/u))-.5*g*np.linalg.norm(X.T@r)**2)

def KL(a,b):
    return float(np.sum(a*np.log(a/b)))

rhs=qt-p
grad=np.log(p/u)+1-g*X@mu
dV=grad@rhs
J=KL(p,qt)+KL(qt,p)
assert abs(dV+J)<1e-12
assert J>0

# finite difference replay
eps=1e-7
fd=(V(p+eps*rhs)-V(p-eps*rhs))/(2*eps)
assert abs(fd-dV)<2e-8

print("PASS")
print("target_entropy",Hq)
print("lazy_kernel_conditional_entropy",Hlazy)
print("maxent_target_object",best)
print("Jeffreys_dissipation",J)
print("analytic_dVdt",dV)
print("finite_difference_dVdt",fd)
