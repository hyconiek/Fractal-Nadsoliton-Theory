#!/usr/bin/env python3
import math
import numpy as np
from scipy.special import softmax

Q=12
g=3.7183448981203875
j=np.arange(Q)
W=np.array([[0.0 if i==k else math.cos(0.18575*min(abs(i-k),Q-abs(i-k))+0.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)],float)
Aker=np.diag(W.sum(1))-W
lam=np.fft.fft(Aker[0]).real[:7]
cols=[]
for k in (3,4,5):
    cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),
             np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.0)**j]
X=np.column_stack(cols)
Lambda=X.T@X
F=Lambda/Q
R=X@np.linalg.inv(Lambda)
u=np.ones(Q)/Q


def q_exact(eps,x):
    return softmax(eps*g*(X@x))

def exact_L(eps,x,c,deg):
    # Formal analytic continuation in eps; p remains positive for tested eps.
    r=R@x
    p=u+eps*r
    q=q_exact(eps,x)
    s=float(c@x)
    y=X@c
    out=0.0
    for i in range(Q):
        for jj in range(Q):
            ds=eps*(y[jj]-y[i])
            df=(s+ds)**deg-s**deg
            out += (p[i]*q[jj]/eps**2)*df
    return out

def predicted(x,c,deg):
    r=R@x
    a=X@x
    m=(g/Q)*(Lambda@x)
    q1=(g/Q)*a
    a2=float(np.mean(a*a))
    q2=(g*g/(2*Q))*(a*a-a2)
    a3=float(np.mean(a*a*a))
    q3=(g**3/(6*Q))*(a*a*a-3*a2*a-a3)
    b0=m-x
    b1=X.T@q2
    b2=X.T@q3
    D0=2*F
    D1=X.T@((r+q1)[:,None]*X)
    D2=X.T@(q2[:,None]*X)-np.outer(x,m)-np.outer(m,x)
    y=X@c
    cF=float(c@F@c)
    C31=float(np.dot(q1-r,y**3)+3*(c@(m-x))*cF)
    C40=float(2*np.mean(y**4)+6*cF*cF)
    s=float(c@x)
    d1=deg*s**(deg-1) if deg>=1 else 0.0
    d2=deg*(deg-1)*s**(deg-2) if deg>=2 else 0.0
    d3=deg*(deg-1)*(deg-2)*s**(deg-3) if deg>=3 else 0.0
    d4=deg*(deg-1)*(deg-2)*(deg-3)*s**(deg-4) if deg>=4 else 0.0
    L0=d1*(c@b0)+0.5*d2*(c@D0@c)
    L1=d1*(c@b1)+0.5*d2*(c@D1@c)
    L2=d1*(c@b2)+0.5*d2*(c@D2@c)+(d3/6)*C31+(d4/24)*C40
    return L0,L1,L2

rng=np.random.default_rng(20260925)
max_l1=max_l2=0.0
for case in range(8):
    x=rng.normal(size=7); x*=0.35/np.linalg.norm(x)
    c=rng.normal(size=7); c/=np.linalg.norm(c)
    deg=1+(case%4)
    L0,L1,L2=predicted(x,c,deg)
    rows=[]
    for e in (0.04,0.02,0.01,0.005):
        lp=exact_L(e,x,c,deg); lm=exact_L(-e,x,c,deg)
        l1_est=(lp-lm)/(2*e)
        l2_est=((lp+lm)/2-L0)/(e*e)
        rows.append((e,l1_est,l2_est))
    # Richardson cancellation of leading e^2 error using final two scales.
    _,l1a,l2a=rows[-2]; _,l1b,l2b=rows[-1]
    l1R=(4*l1b-l1a)/3
    l2R=(4*l2b-l2a)/3
    e1=abs(l1R-L1); e2=abs(l2R-L2)
    max_l1=max(max_l1,e1); max_l2=max(max_l2,e2)
    print(f'case={case} deg={deg} L0={L0:.16e} L1={L1:.16e} L2={L2:.16e}')
    for rr in rows: print('  eps=%g L1est=%.16e L2est=%.16e'%rr)
    print(f'  richardson L1={l1R:.16e} err={e1:.3e} L2={l2R:.16e} err={e2:.3e}')
print('max_richardson_error_L1',repr(max_l1))
print('max_richardson_error_L2',repr(max_l2))
assert max_l1 < 2e-8
assert max_l2 < 2e-7
