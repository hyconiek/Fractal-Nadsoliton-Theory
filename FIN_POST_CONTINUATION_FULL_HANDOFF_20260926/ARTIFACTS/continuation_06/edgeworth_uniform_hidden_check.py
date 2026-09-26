#!/usr/bin/env python3
import math, numpy as np
from scipy.special import softmax
Q=12; g=3.7183448981203875; jj=np.arange(Q)
W=np.array([[0.0 if i==k else math.cos(0.18575*min(abs(i-k),Q-abs(i-k))+0.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)],float)
Aker=np.diag(W.sum(1))-W; lam=np.fft.fft(Aker[0]).real[:7]
cols=[]
for k in (3,4,5): cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
X=np.column_stack(cols); Lambda=X.T@X; F=Lambda/Q; R=X@np.linalg.inv(Lambda); u=np.ones(Q)/Q
Y=[]
for k in (1,2): Y += [np.sqrt(2/Q)*np.cos(2*np.pi*k*jj/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*jj/Q)]
Y=np.column_stack(Y)
Ts=np.array([X.T@(Y[:,a,None]*X) for a in range(4)])

def exact_L(eps,x,z,c,deg):
    d=R@x+Y@z; p=u+eps*d; q=softmax(eps*g*(X@x)); s=float(c@x); y=X@c; out=0.
    for i in range(Q):
      for j in range(Q):
        ds=eps*(y[j]-y[i]); out += p[i]*q[j]/eps**2*((s+ds)**deg-s**deg)
    return out

def pred(x,z,c,deg):
    r=R@x; yz=Y@z; a=X@x; m=g/Q*(Lambda@x); q1=g/Q*a
    a2=np.mean(a*a); q2=g*g/(2*Q)*(a*a-a2)
    a3=np.mean(a*a*a); q3=g**3/(6*Q)*(a**3-3*a2*a-a3)
    b0=m-x; b1=X.T@q2; b2=X.T@q3
    D0=2*F
    D1=X.T@((r+q1)[:,None]*X)+np.einsum('a,aij->ij',z,Ts)
    D2=X.T@(q2[:,None]*X)-np.outer(x,m)-np.outer(m,x)
    y=X@c; cF=c@F@c
    C31=np.dot(q1-r-yz,y**3)+3*(c@(m-x))*cF
    C40=2*np.mean(y**4)+6*cF*cF
    s=c@x
    d1=deg*s**(deg-1); d2=deg*(deg-1)*s**(deg-2) if deg>=2 else 0
    d3=deg*(deg-1)*(deg-2)*s**(deg-3) if deg>=3 else 0
    d4=deg*(deg-1)*(deg-2)*(deg-3)*s**(deg-4) if deg>=4 else 0
    return (d1*(c@b0)+.5*d2*(c@D0@c),
            d1*(c@b1)+.5*d2*(c@D1@c),
            d1*(c@b2)+.5*d2*(c@D2@c)+d3*C31/6+d4*C40/24)

rng=np.random.default_rng(91577); mx1=mx2=0
for case in range(8):
    x=rng.normal(size=7); x*=.25/np.linalg.norm(x)
    z=rng.normal(size=4); z*=.3/np.linalg.norm(z)
    c=rng.normal(size=7); c/=np.linalg.norm(c); deg=1+case%4
    L0,L1,L2=pred(x,z,c,deg)
    vals=[]
    for e in (.04,.02,.01,.005):
       lp=exact_L(e,x,z,c,deg); lm=exact_L(-e,x,z,c,deg)
       vals.append((e,(lp-lm)/(2*e),((lp+lm)/2-L0)/e**2))
    _,a1,a2=vals[-2]; _,b1,b2=vals[-1]
    r1=(4*b1-a1)/3; r2=(4*b2-a2)/3
    er1=abs(r1-L1);er2=abs(r2-L2);mx1=max(mx1,er1);mx2=max(mx2,er2)
    print(case,deg,'L1',L1,'L1R',r1,'err',er1,'L2',L2,'L2R',r2,'err',er2)
print('maxerr',mx1,mx2)
assert mx1<5e-8 and mx2<5e-7
