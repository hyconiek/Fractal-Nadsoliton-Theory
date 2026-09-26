#!/usr/bin/env python3
import math,numpy as np
from scipy.special import softmax
from scipy.optimize import root
Q=12;g=3.7183448981203875;jj=np.arange(Q)
W=np.array([[0.0 if i==j else math.cos(.18575*min(abs(i-j),Q-abs(i-j))+.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)])
L=np.diag(W.sum(1))-W;lam=np.fft.fft(L[0]).real[:7]
cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
X=np.column_stack(cols);Gamma=X.T@X;F=Gamma/Q;R=X@np.linalg.inv(Gamma);u=np.ones(Q)/Q
P0=np.ones((Q,Q))/Q;PV=X@np.linalg.inv(Gamma)@X.T;PH=np.eye(Q)-P0-PV

def pme(mu):
    th0=Q*np.linalg.solve(Gamma,mu)
    def fun(th): return X.T@softmax(X@th)-mu
    sol=root(fun,th0,tol=1e-12)
    if not sol.success and np.max(abs(fun(sol.x)))>1e-11: raise RuntimeError(sol.message)
    return softmax(X@sol.x)

def exact_L(eps,x,c,deg):
    mu=eps*x;p=pme(mu);q=softmax(eps*g*(X@x));s=c@x;y=X@c;out=0.
    for i in range(Q):
      for j in range(Q):
        ds=eps*(y[j]-y[i]);out+=p[i]*q[j]/eps**2*((s+ds)**deg-s**deg)
    return out

def pred(x,c,deg):
    r=R@x;a=X@x;m=g/Q*(Gamma@x);q1=g/Q*a
    a2=np.mean(a*a);q2=g*g/(2*Q)*(a*a-a2)
    a3=np.mean(a**3);q3=g**3/(6*Q)*(a**3-3*a2*a-a3)
    w=(Q/2)*PH@(r*r)
    b0=m-x;b1=X.T@q2;b2=X.T@q3
    D0=2*F;D1=X.T@((r+q1)[:,None]*X)
    D2=X.T@((q2+w)[:,None]*X)-np.outer(x,m)-np.outer(m,x)
    y=X@c;cF=c@F@c
    C31=np.dot(q1-r,y**3)+3*(c@(m-x))*cF
    C40=2*np.mean(y**4)+6*cF*cF
    s=c@x;d1=deg*s**(deg-1);d2=deg*(deg-1)*s**(deg-2) if deg>=2 else 0;d3=deg*(deg-1)*(deg-2)*s**(deg-3) if deg>=3 else 0;d4=deg*(deg-1)*(deg-2)*(deg-3)*s**(deg-4) if deg>=4 else 0
    return d1*c@b0+.5*d2*c@D0@c, d1*c@b1+.5*d2*c@D1@c, d1*c@b2+.5*d2*c@D2@c+d3*C31/6+d4*C40/24
rng=np.random.default_rng(1509);mx1=mx2=0
for case in range(8):
 x=rng.normal(size=7);x*=.22/np.linalg.norm(x);c=rng.normal(size=7);c/=np.linalg.norm(c);deg=1+case%4
 L0,L1,L2=pred(x,c,deg);vals=[]
 for e in (.03,.015,.0075,.00375):
  lp=exact_L(e,x,c,deg);lm=exact_L(-e,x,c,deg);vals.append((e,(lp-lm)/(2*e),((lp+lm)/2-L0)/e**2))
 _,a1,a2=vals[-2];_,b1,b2=vals[-1];r1=(4*b1-a1)/3;r2=(4*b2-a2)/3
 er1=abs(r1-L1);er2=abs(r2-L2);mx1=max(mx1,er1);mx2=max(mx2,er2)
 print(case,deg,'L1err',er1,'L2err',er2,'L2',L2,'L2R',r2)
print('maxerr',mx1,mx2)
assert mx1<5e-8 and mx2<5e-7
