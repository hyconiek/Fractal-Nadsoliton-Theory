#!/usr/bin/env python3
"""Independent numerical diagnostics, NOT a proof or certified interval algebra.
Builds A and the fold without importing previous producer's theory_core.
"""
import json,math
import numpy as np
from scipy.optimize import root
from scipy.linalg import null_space,solve_continuous_lyapunov
from pathlib import Path

Q=12
j=np.arange(Q)
W=np.zeros((Q,Q))
for i in range(Q):
    for k in range(Q):
        if i!=k:
            d=min(abs(i-k),Q-abs(i-k))
            W[i,k]=math.cos(.18575*d+.1625)/(1+d**1.8)
L=np.diag(W.sum(axis=1))-W
ll=np.fft.fft(L[0]).real
X=np.column_stack([u for k in (3,4,5) for u in (np.sqrt(ll[k]/6)*np.cos(2*math.pi*k*j/Q),np.sqrt(ll[k]/6)*np.sin(2*math.pi*k*j/Q))]+[np.sqrt(ll[6]/12)*(-1.)**j])
A=X@X.T
Xe=X[:,[0,2,4,6]]
def sm(t):
    x=t-np.max(t); y=np.exp(x);return y/y.sum()
def fun(h):
    theta=h[:4];g=h[4]; p=sm(Xe@theta);C=np.diag(p)-np.outer(p,p); M=Xe.T@C@Xe
    return np.r_[theta/g-Xe.T@p,np.linalg.det(np.eye(4)/g-M)]
a0=np.array([1.3643114281782411,1.4331027014285429,1.4080457012550378,1.0056820773201627,3.5156447168395917]);res=root(fun,a0,tol=1e-12)
assert np.linalg.norm(fun(res.x))<1e-10
th=np.zeros(7);th[[0,2,4,6]]=res.x[:4];g=float(res.x[4]);p=sm(X@th)
Sigma=np.diag(p)-np.outer(p,p);M=X.T@Sigma@X
H=np.eye(7)/g-M
evalh,vs=np.linalg.eigh(H);v=vs[:,0]
if v[0]<0:v=-v
rf=Sigma@X@v;ell=g*v@X.T
J=g*Sigma@A-np.eye(Q)
T=null_space(np.ones((1,Q)))
E=null_space(np.vstack([np.ones(Q),ell]))
B=E.T@J@E
P=solve_continuous_lyapunov(B.T,-np.eye(10))
QL=np.eye(Q)-np.outer(rf,ell)
# r_f has Euclidean norm not necessarily one
scores={
 "model":"q12_loo_exact_decimal_formula_binary64_evaluation", "fold_g_nominal":g,
 "fold_equation_residual":float(np.linalg.norm(fun(res.x))),
 "fold_theta7":th.tolist(),"fold_v7":v.tolist(),"pf":p.tolist(),
 "full_H7_eigenvalues":evalh.tolist(),
 "J_fast_eigenvalues_real":sorted([float(np.real(t)) for t in np.linalg.eigvals(B)]),
 "B_fast_max_real_eigenvalue":float(max(np.linalg.eigvals(B).real)),
 "P_eigen_minmax":[float(np.linalg.eigvalsh(P)[0]),float(np.linalg.eigvalsh(P)[-1])],
 "norm_J_r":float(np.linalg.norm(J@rf)),
 "norm_ell_J_on_tangent":float(np.linalg.norm(ell@J@T)),
 "ell_r_minus_1":float(ell@rf-1),
 "Q_squared_minus_Q":float(np.linalg.norm(QL@QL-QL)),
 "Q_r":float(np.linalg.norm(QL@rf)),
 "Lyapunov_residual":float(np.linalg.norm(B.T@P+P@B+np.eye(10))),
 "A_symmetric_max":float(np.max(np.abs(A-A.T))),
 "a_nominal":float(-v@th/g**2),
 "b_nominal":float(-np.sum(p*((X@v)-(p@(X@v)))**3)),
 "soft_noise_variance_nominal":float(2*g),
 "precision_scope":"binary64 numerical consistency only; local source interval fold taken as separate premise"
}
# full occupation covariance and soft channel comparison at fold
D=2*Sigma
scores['soft_noise_direct']=float(g*g*v@X.T@D@X@v)
scores['soft_noise_direct_minus_2g']=scores['soft_noise_direct']-2*g
folder=Path(__file__).resolve().parent
(folder/'SF-001'/'NUMERICAL_CHART.json').write_text(json.dumps(scores,indent=2)+'\n')
print(json.dumps({k:scores[k] for k in ('fold_equation_residual','B_fast_max_real_eigenvalue','P_eigen_minmax','norm_J_r','norm_ell_J_on_tangent','Lyapunov_residual','soft_noise_direct_minus_2g')},indent=2))
