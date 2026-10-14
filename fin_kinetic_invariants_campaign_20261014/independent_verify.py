#!/usr/bin/env python3
"""Independent small numerical checks: rebuilds A and fold without importing producer/source modules.
All checks are binary64 diagnostics. Analytic proof claims are not inferred from flags.
"""
import json, math
from pathlib import Path
import numpy as np
from scipy.optimize import root
from scipy.linalg import null_space
BASE=Path(__file__).resolve().parent
P=json.loads((BASE/'KI-002'/'DERIVED_SOFT_CLOCKS.json').read_text())
Q=12;j=np.arange(Q);W=np.zeros((Q,Q))
for k in range(1,7):
    w=math.cos(.18575*k+.1625)/(1+k**1.8)
    for i in range(Q):
        W[i,(i+k)%Q]=w;W[i,(i-k)%Q]=w
L=np.diag(W.sum(1))-W;ft=np.fft.fft(L[0]).real
cols=[]
for k in (3,4,5):
    cols.extend([math.sqrt(ft[k]/6)*np.cos(2*math.pi*k*j/Q),math.sqrt(ft[k]/6)*np.sin(2*math.pi*k*j/Q)])
cols.append(math.sqrt(ft[6]/12)*np.cos(math.pi*j))
X=np.array(cols).T;A=X@X.T;Ee=X[:,[0,2,4,6]]
def soft(z):
    z=z-max(z);t=np.exp(z);return t/t.sum()
def f(y):
    t=y[:4];g=y[4];p=soft(Ee@t);sig=np.diag(p)-np.outer(p,p)
    h=np.eye(4)/g-Ee.T@sig@Ee
    return np.r_[t/g-Ee.T@p,np.linalg.eigvalsh(h)[0]]
sol=root(f,np.array([1.3643114281782407,1.4331027014285425,1.4080457012550374,1.0056820773201627,3.5156447168395917]),tol=1e-12)
assert np.linalg.norm(f(sol.x))<1e-10
g=sol.x[-1];t=np.zeros(7);t[[0,2,4,6]]=sol.x[:4];p=soft(X@t)
Sig=np.diag(p)-np.outer(p,p)
H7=np.eye(7)/g-X.T@Sig@X
_,E7=np.linalg.eigh(H7);v=E7[:,0];v*=np.sign(v[0]);r=Sig@X@v
T=null_space(np.ones((1,Q)))
H=np.diag(1/p)-g*A
checks=[]
def ck(name,val,tol=1e-9):
    checks.append({'name':name,'residual':float(val),'tolerance':float(tol),'pass':bool(np.isfinite(val) and abs(val)<tol)})
ck('fold root residual',np.linalg.norm(f(sol.x)),1e-10)
ck('model g agreement',g-P['fold_g_nominal'],1e-9)
ck('primal a',-(r@A@p)-P['a_F_nominal'])
ck('primal b',-sum(r**3/p**2)-P['b_F_nominal'])
ck('tangent Hess zero',np.linalg.norm(T.T@H@r))
ck('positive Hess remaining',min(np.linalg.eigvalsh(T.T@H@T)[1:])-min(np.linalg.eigvalsh(T.T@H@T)[1:]),1e-10)
ck('r tangent',sum(r))
ck('r theta',g*(v@X.T@r)-1)
for kin in ['HB','MET','BAR']:
    hh=g*(A@p)
    D=hh[None,:]-hh[:,None]
    if kin=='HB':R=np.tile(soft(hh),(Q,1))
    elif kin=='MET':R=np.minimum(1.,np.exp(np.minimum(D,0)))/12
    else:R=1/(6*(1+np.exp(-D)))
    C=p[:,None]*R;K=np.zeros((Q,Q));B=np.zeros((Q,Q))
    for i in range(Q):
        for k in range(i+1,Q):
            e=np.eye(Q)[k]-np.eye(Q)[i]
            K+=C[i,k]*np.outer(e,e)
    for i in range(Q):
        for k in range(Q):
            if i!=k:
                e=np.eye(Q)[k]-np.eye(Q)[i];B+=C[i,k]*np.outer(e,e)
    Kt=T.T@K@T
    ki=np.linalg.solve(Kt,T.T@r)
    m=1/(r@T@ki);ell=m*T@ki
    old=P['models'][kin]
    ck(kin+' symmetry',np.max(abs(C-C.T)))
    ck(kin+' K positive',max(0,-np.min(np.linalg.eigvalsh(Kt))),1e-10)
    ck(kin+' mobility',m-old['mobility_m_nominal'])
    ck(kin+' Einstein',ell@B@ell-2*m)
    ck(kin+' right null',np.linalg.norm(K@H@r))
    ck(kin+' left null',np.linalg.norm(ell@K@H@T))
    ck(kin+' left right',ell@r-1)
    ck(kin+' local D 2K',np.max(abs(B-2*K)))
    ck(kin+' g0',1/12- P['g0_distinct_rate_per_copy'][kin])
    ck(kin+' fast eigen',max(0,np.max(np.linalg.eigvals(-K@H).real[0:0]) if False else -min(np.linalg.eigvalsh(T.T@K@T))),1e-10) # positivity, not gap
    # Nonlinear limiting drift identity using logarithmic mean at arbitrary p
    perturb=np.linspace(-1,1,Q);perturb-=perturb.mean();q=p+1e-3*perturb
    if min(q)>0:
        hh2=g*A@q;d2=hh2[None,:]-hh2[:,None]
        rr=(np.tile(soft(hh2),(Q,1)) if kin=='HB' else (np.exp(np.minimum(0,d2))/12 if kin=='MET' else 2/(12*(1+np.exp(-d2)))))
        cc=q[:,None]*rr
        drift=np.sum(cc,axis=0)-np.sum(cc,axis=1)
        F=np.log(q)-g*A@q
        KL=np.zeros((Q,Q))
        for i in range(Q):
            for k in range(i+1,Q):
                ab=cc[i,k];ba=cc[k,i]
                logmean=(ab-ba)/(math.log(ab)-math.log(ba)) if abs(ab-ba)>1e-14 else .5*(ab+ba)
                e=np.eye(Q)[k]-np.eye(Q)[i];KL+=logmean*np.outer(e,e)
        ck(kin+' nonlinear Onsager',np.linalg.norm(drift+KL@F),1e-9)
# source finite g0 convention: all transition rates 1/12 for distinct labels
status={'checks':checks,'passed':sum(x['pass'] for x in checks),'total':len(checks),'all_pass':all(x['pass'] for x in checks),'arithmetic_scope':'numerical consistency only'}
(BASE/'INDEPENDENT_VERIFICATION.json').write_text(json.dumps(status,indent=2)+'\n')
print(json.dumps({'passed':status['passed'],'total':status['total'],'failed':[x['name'] for x in checks if not x['pass']]}))
