#!/usr/bin/env python3
import numpy as np
from scipy.special import softmax

rng=np.random.default_rng(20260926)
Q=8
d=3
X=rng.normal(size=(Q,d))
X-=X.mean(axis=0)
g=.2

# Find a weakly coupled interior fixed point (uniform is exact because X columns are mean zero).
p=np.ones(Q)/Q
q=softmax(g*X@(X.T@p))
assert np.max(np.abs(q-p))<1e-14

def logmean(a,b):
    if abs(a-b)<1e-13:
        return .5*(a+b)
    return (a-b)/(np.log(a)-np.log(b))

# Non-equilibrium random point for nonlinear Onsager identity.
pn=rng.uniform(.2,2,size=Q); pn/=pn.sum()
qn=softmax(g*X@(X.T@pn))
r=pn/qn
phi=np.log(r)
W=np.zeros((Q,Q))
for i in range(Q):
    for j in range(i+1,Q):
        w=qn[i]*qn[j]*logmean(r[i],r[j])
        W[i,j]=W[j,i]=w
L=np.diag(W.sum(axis=1))-W
assert np.max(np.abs(-L@phi-(qn-pn)))<1e-12

def KL(a,b): return np.sum(a*np.log(a/b))
J=KL(pn,qn)+KL(qn,pn)
assert abs(phi@L@phi-J)<1e-12

# Equilibrium mobility is S.
S=np.diag(p)-np.outer(p,p)
Weq=np.outer(p,p); np.fill_diagonal(Weq,0)
Leq=np.diag(Weq.sum(axis=1))-Weq
assert np.linalg.norm(Leq-S)<1e-14

# Label-space linearization / FDT.
A7=X@X.T
H=np.diag(1/p)-g*A7
B=-S@H
D=2*S

# Work in tangent orthonormal basis.
U=np.linalg.qr(np.column_stack([np.ones(Q),rng.normal(size=(Q,Q-1))]))[0][:,1:]
# safer construct tangent basis via QR on projector
evals,evecs=np.linalg.eigh(S)
U=evecs[:,evals>1e-12]
Ht=U.T@H@U
St=U.T@S@U
Bt=U.T@B@U
assert np.linalg.eigvalsh(Ht).min()>0

C=np.linalg.inv(Ht)
res=Bt@C+C@Bt.T+2*St
assert np.linalg.norm(res)<1e-11

# Similarity to symmetric operator => real spectrum.
sqrtS=np.linalg.cholesky(St)
# use symmetric S^1/2 from eig
se,sv=np.linalg.eigh(St)
Shalf=sv@np.diag(np.sqrt(se))@sv.T
sym=-Shalf@Ht@Shalf
be=np.linalg.eigvals(Bt)
assert np.max(np.abs(be.imag))<1e-12
assert np.max(np.linalg.eigvalsh(sym))<0

# Artificially soften one Hessian eigenvalue and verify dynamic rate tends to zero.
he,hv=np.linalg.eigh(Ht)
mins=[]
for eps in (1e-1,1e-2,1e-3,1e-4):
    Hsoft=hv@np.diag(np.r_[eps,he[1:]])@hv.T
    Bsoft=-St@Hsoft
    rates=-np.linalg.eigvals(Bsoft).real
    mins.append(rates.min())
assert mins[-1] < mins[0]/500

print("PASS")
print("nonlinear_onsager_error",float(np.max(np.abs(-L@phi-(qn-pn)))))
print("Jeffreys_identity_error",float(abs(phi@L@phi-J)))
print("equilibrium_mobility_error",float(np.linalg.norm(Leq-S)))
print("FDT_Lyapunov_error",float(np.linalg.norm(res)))
print("soft_mode_rates",mins)
