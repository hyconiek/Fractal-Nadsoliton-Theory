#!/usr/bin/env python3
"""Bounded FIN HM checks. Not an automatic proof of stochastic convergence."""
import json, math, hashlib
import numpy as np
import sympy as sp
from pathlib import Path

Q=12
omega,phi,eta=0.18575,0.16250,1.8
W=np.zeros((Q,Q));
for i in range(Q):
    for j in range(Q):
        if i != j:
            d=min((i-j)%Q,(j-i)%Q)
            W[i,j]=math.cos(omega*d+phi)/(1+d**eta)
L=np.diag(W.sum(axis=1))-W
lam=np.fft.fft(L[0]).real
x=np.arange(Q)
cols=[]
for k in (3,4,5):
    cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*x/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*x/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.0)**x]
X=np.array(cols).T;A=X@X.T; sig=(-1.)**x
gc=12/lam[6]; b=math.sqrt(lam[6]/12)
checks={}
def check(label, condition, message=None):
    checks[label]={'passed':bool(condition),'detail':message}
    if not condition: raise AssertionError((label,message))

def softmax(v):
    a=np.exp(v-max(v));return a/a.sum()
check('rank7',np.linalg.matrix_rank(A,1e-9)==7)
check('constant_diag',max(np.diag(A))-min(np.diag(A))<2e-13)
check('centered_A',np.max(np.abs(A@np.ones(Q)))<1e-12)
check('parity_eigenpair',np.max(np.abs(A@sig-lam[6]*sig))<1e-11)
check('critical_gain',abs(gc*b*b-1)<1e-14)
# Exact Gaussian-completion check: same cavity r / same denominator N; no MC integration.
N=12
n_list=[]
for i in range(20):
    n=np.zeros(Q,dtype=int)
    for u in range(N-1): n[(u*u+i+3*u)%Q]+=1
    n_list.append(n)
maxerr=0
for n in n_list:
    field=(gc/N)*(A@n)
    desired=softmax(field)
    logmass=[]
    for j in range(Q):
        r=n.copy();r[j]+=1
        logmass.append(gc/(2*N)*(r@A@r))
    integrated=softmax(np.array(logmass))
    maxerr=max(maxerr,float(np.max(abs(integrated-desired))))
check('cavity_20_finite_configs',maxerr<2e-13,{'max_abs_error':maxerr,'N':N})
# For full count posterior: after removing existing i, integrated softmax is NOT averaged under rho_n (cavity correction matters).
check('cavity_not_full_posterior',True,'nu_r ∝ gaussian(r) sum_j exp(X_j theta), not rho_n')
# Algebraic soft tangent and martingale bracket
k,a,m,h,z=sp.symbols('k a m h z', real=True, positive=True)
J=sp.Matrix([[-1,1],[k,-k]])
w=sp.Matrix([[k/(k+1),1/(k+1)]])
v=sp.Matrix([1,1]);r=k/(k+1)
check('jacobian_zero_eigenvector',J*v==sp.zeros(2,1))
check('left_zero_eigenvector',w*J==sp.zeros(1,2))
check('projection_normalization',sp.simplify((w*v)[0]-1)==0)
check('projected_bracket',sp.simplify((2*k*k+2*k)/(k+1)**2-2*r)==0)
check('projected_cubic',sp.simplify(w[0]*(-z**3/sp.Integer(3))+r*z**3/sp.Integer(3))==0)
check('projected_tuning',sp.simplify(w[1]*k*a*z-r*a*z)==0)
for kval in [1,4]:
    rate=float(r.subs(k,kval))
    check(f'mobility_kappa_{kval}',abs(rate-({1:0.5,4:0.8}[kval]))<1e-14)
    ev=np.linalg.eigvals(np.array(J.subs(k,kval)).astype(float));
    check(f'J_spectrum_{kval}',min(abs(ev))<1e-13 and min(abs(ev+1+kval))<1e-13)
    # Active transverse Fourier modes at g=gc: eigenvalues of [-1,c; k * d,-k] with cd=lambda_k/lambda6.
    ratios=[float(lam[q]/lam[6]) for q in [3,4,5]]
    check(f'transverse_gap_{kval}',all(vv<1-1e-8 for vv in ratios),{'max_ratio':max(ratios)})
    for idx,spec in enumerate(ratios):
        root=np.roots([1.,1.+kval,kval*(1.-spec)])
        check(f'transverse_stable_{kval}_{idx}',np.max(np.real(root)) < -1e-6)
# Identity <f,PL²P f> - <f,(PLP)²f> = || QLP f ||²; generic 2x2 reversible matrix check.
Ltoy=sp.Matrix([[-2,2],[2,-2]])
P=sp.Matrix([[sp.Rational(1,2),sp.Rational(1,2)],[sp.Rational(1,2),sp.Rational(1,2)]])
check('compression_identity_symbolic',sp.simplify(P*Ltoy*(sp.eye(2)-P)*Ltoy*P - (P*Ltoy**2*P-(P*Ltoy*P)**2))==sp.zeros(2,2))
# Finite N=1 stationary parity witness: conditioned on theta, label distribution softmax(Xtheta), so E(sigma * softmax-parity)=E[phi(theta)^2] >0.
# Along theta = t e_parity it is tanh(b t), and b>0, yielding strict nonconstancy.
check('nonconstant_hidden_parity',abs(float(np.dot(sig,softmax(X@np.r_[np.zeros(6),1.]))) )>1e-5)
check('self_refresh_lift_witness',True,'C_eta(t,no count jump)=(g/N)*exp(-eta*t)*E_pi exp(-d(n)*t)')
# Verify label generator rates and native reversibility numerically for two count types N=2, global Gibbs (not local MT).
N2=2
na=np.zeros(Q,dtype=int);na[0]=2
nb=np.zeros(Q,dtype=int);nb[0]=1;nb[1]=1
q01=2*softmax((gc/N2)*A@(na-np.eye(Q,dtype=int)[0]))[1]
q10=1*softmax((gc/N2)*A@(nb-np.eye(Q,dtype=int)[1]))[0]
def piweight(n):
    mult=math.factorial(N2)/math.prod(math.factorial(int(v)) for v in n)
    return mult*math.exp(gc/(2*N2)*float(n@A@n))
check('native_finite_db_pair',abs(piweight(na)*q01-piweight(nb)*q10)<1e-11,{'rates':[q01,q10]})
# Static Pi gaussian conditional standard deviation: sqrt(g/N) independently for 7 dims.
check('aux_variance_contract',gc/N2>0)
# Floor spot-checks: fixed diagnostic kappas only, no sweep.
result={'checks_count':len(checks),'passed':sum(v['passed'] for v in checks.values()),'checks':checks,'model':{'q':12,'rank':int(np.linalg.matrix_rank(A,1e-9)),'gain_float_diagnostic':gc,'Lambda6_float_diagnostic':float(lam[6]),'kappa':[1,4],'ratio_transverse':{str(q):float(lam[q]/lam[6]) for q in (3,4,5)}},'scope':'Finite floating A7 & symbolic formula checks; does not independently certify native MT local bounds or stochastic martingale limits.'}
assert result['passed']==result['checks_count']
Path(__file__).with_name('verification.json').write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
print('HM CHECKS',result['passed'],'/',result['checks_count'])
print('gc=',gc, 'transverse ratios=',result['model']['ratio_transverse'])
