#!/usr/bin/env python3
import math, numpy as np
from scipy.special import softmax
Q=12; g=3.7183448981203875; jj=np.arange(Q)
W=np.array([[0.0 if i==k else math.cos(.18575*min(abs(i-k),Q-abs(i-k))+.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)],float)
Aker=np.diag(W.sum(1))-W; lam=np.fft.fft(Aker[0]).real[:7]
cols=[]
for k in (3,4,5): cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
X=np.column_stack(cols); Gamma=X.T@X; R=X@np.linalg.inv(Gamma)
Y=[]
for k in (1,2): Y += [np.sqrt(2/Q)*np.cos(2*np.pi*k*jj/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*jj/Q)]
Y=np.column_stack(Y)

def ps(vals):
    th=np.zeros(7);th[[0,2,4,6]]=vals;return softmax(X@th)
states={
'uniform':np.ones(Q)/Q,
'saddle':ps([0.9409570673214491,1.0014394023775437,0.9621088536282347,0.6864149839950513]),
'localized':ps([1.8199035812800828,1.913989554668724,1.914569132546848,1.367203280195504])}

def setup(pstar):
    mu=X.T@pstar
    qstar=softmax(g*X@mu)
    S=np.diag(pstar)-np.outer(pstar,pstar)
    F=X.T@S@X; H=Y.T@S@X; G=Y.T@S@Y
    K=H@np.linalg.inv(F); Sig=G-H@np.linalg.solve(F,H.T)
    Xi=X-mu[None,:]
    T=np.array([Xi.T@(Y[:,a,None]*Xi) for a in range(4)])
    return mu,qstar,F,K,Sig,Xi,T

def exact_L(eps,pstar,mu,x,z,c,deg,K):
    d=(R+Y@K)@x+Y@z
    p=pstar+eps*d
    q=softmax(g*X@(mu+eps*x))
    s=float(c@x); y=X@c; out=0.
    for i in range(Q):
      for j in range(Q):
        ds=eps*(y[j]-y[i]);out += p[i]*q[j]/eps**2*((s+ds)**deg-s**deg)
    return out

def pred(pstar,mu,F,K,Xi,T,x,z,c,deg):
    dx=(R+Y@K)@x; yz=Y@z
    aa=g*(X@x); mean=np.dot(pstar,aa); at=aa-mean
    kap2=np.dot(pstar,at*at); kap3=np.dot(pstar,at**3)
    q1=pstar*at
    q2=.5*pstar*(at*at-kap2)
    q3=(pstar/6)*(at**3-3*kap2*at-kap3)
    m=X.T@q1
    b0=m-x; b1=X.T@q2; b2=X.T@q3
    D0=2*F
    D1=Xi.T@((dx+q1+yz)[:,None]*Xi)
    D2=Xi.T@(q2[:,None]*Xi)-np.outer(x,m)-np.outer(m,x)
    yy=Xi@c; cF=c@F@c
    C31=np.dot(q1-dx-yz,yy**3)+3*(c@(m-x))*cF
    C40=2*np.dot(pstar,yy**4)+6*cF*cF
    s=c@x; d1=deg*s**(deg-1); d2=deg*(deg-1)*s**(deg-2) if deg>=2 else 0
    d3=deg*(deg-1)*(deg-2)*s**(deg-3) if deg>=3 else 0
    d4=deg*(deg-1)*(deg-2)*(deg-3)*s**(deg-4) if deg>=4 else 0
    return (d1*c@b0+.5*d2*c@D0@c,
            d1*c@b1+.5*d2*c@D1@c,
            d1*c@b2+.5*d2*c@D2@c+d3*C31/6+d4*C40/24)

rng=np.random.default_rng(70251)
for name,pstar in states.items():
    mu,qstar,F,K,Sig,Xi,T=setup(pstar)
    print('STATE',name,'stationarity_residual',np.max(np.abs(qstar-pstar)),'Sig_eigs',np.linalg.eigvalsh(Sig))
    assert np.max(np.abs(qstar-pstar))<3e-10
    mx1=mx2=0.
    for case in range(6):
      x=rng.normal(size=7);x*=.18/np.linalg.norm(x)
      z=rng.normal(size=4);z*=.22/np.linalg.norm(z)
      c=rng.normal(size=7);c/=np.linalg.norm(c);deg=1+case%4
      L0,L1,L2=pred(pstar,mu,F,K,Xi,T,x,z,c,deg)
      vals=[]
      for e in (.03,.015,.0075,.00375):
        lp=exact_L(e,pstar,mu,x,z,c,deg,K);lm=exact_L(-e,pstar,mu,x,z,c,deg,K)
        vals.append((e,(lp-lm)/(2*e),((lp+lm)/2-L0)/e**2))
      _,a1,a2=vals[-2];_,b1,b2=vals[-1]
      r1=(4*b1-a1)/3;r2=(4*b2-a2)/3
      mx1=max(mx1,abs(r1-L1));mx2=max(mx2,abs(r2-L2))
    print(' max_errors',mx1,mx2)
    assert mx1<2e-7 and mx2<3e-6
