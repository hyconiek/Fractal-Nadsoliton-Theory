#!/usr/bin/env python3
import math,numpy as np
from scipy.special import softmax
Q=12;g=3.7183448981203875;jj=np.arange(Q)
W=np.array([[0.0 if i==j else math.cos(.18575*min(abs(i-j),Q-abs(i-j))+.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)])
L=np.diag(W.sum(1))-W;lam=np.fft.fft(L[0]).real[:7]
cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
X=np.column_stack(cols);Gamma=X.T@X;R=X@np.linalg.inv(Gamma)
Y=[]
for k in (1,2):Y += [np.sqrt(2/Q)*np.cos(2*np.pi*k*jj/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*jj/Q)]
Y=np.column_stack(Y)
def ps(vals):
 th=np.zeros(7);th[[0,2,4,6]]=vals;return softmax(X@th)
states={'uniform':np.ones(Q)/Q,'saddle':ps([0.9409570673214491,1.0014394023775437,0.9621088536282347,0.6864149839950513]),'localized':ps([1.8199035812800828,1.913989554668724,1.914569132546848,1.367203280195504])}

def setup(p):
 S=np.diag(p)-np.outer(p,p);F=X.T@S@X;H=Y.T@S@X;G=Y.T@S@Y;K=H@np.linalg.inv(F);Sig=G-H@np.linalg.solve(F,H.T)
 return S,F,K,Sig

def zeta_stationary(p,x):
 S,F,K,Sig=setup(p);d0=(R+Y@K)@x
 out=np.zeros(4)
 for i in range(Q):
  yi=Y[i]
  v=Sig@yi
  svar=yi@Sig@yi
  out += .5*v*((d0[i]*d0[i]+svar)/(p[i]*p[i])-1/p[i])
 return out

def zeta_me(p,x):
 S,F,K,Sig=setup(p);eta1=np.linalg.solve(F,x);aa=X@eta1;at=aa-p@aa;k2=p@(at*at);q2=.5*p*(at*at-k2)
 eta2=-np.linalg.solve(F,X.T@q2);d2=q2+S@X@eta2
 assert np.linalg.norm(X.T@d2)<1e-10 and abs(d2.sum())<1e-10
 return Y.T@d2
rng=np.random.default_rng(6122)
for name,p in states.items():
 print('STATE',name)
 for tag,x in [('zero',np.zeros(7)),('x1',rng.normal(size=7)),('x2',rng.normal(size=7))]:
  if np.linalg.norm(x):x=x*.2/np.linalg.norm(x)
  zs=zeta_stationary(p,x);zm=zeta_me(p,x)
  print(tag,'zeta_stat',zs.tolist(),'zeta_ME',zm.tolist(),'diffnorm',np.linalg.norm(zs-zm))
