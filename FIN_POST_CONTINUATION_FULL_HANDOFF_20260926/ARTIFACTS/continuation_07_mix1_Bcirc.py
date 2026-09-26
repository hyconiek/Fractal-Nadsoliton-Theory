#!/usr/bin/env python3
from functools import lru_cache
import math,numpy as np
Q=12;jj=np.arange(Q)
W=np.array([[0.0 if i==j else math.cos(.18575*min(abs(i-j),Q-abs(i-j))+.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)],float)
L=np.diag(W.sum(1))-W;lam=np.fft.fft(L[0]).real[:7]
cols=[]
for k in (3,4,5): cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*jj/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*jj/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
X=np.column_stack(cols);A=X@X.T
# deterministic zero-mean circulant correction B
b=np.array([math.sin(.53*d)+.4*math.cos(1.17*d)+.2*math.cos(2.31*d) for d in range(Q)])
b=b-b.mean(); B=np.array([[b[(l-i)%Q] for i in range(Q)] for l in range(Q)])
PV=X@np.linalg.inv(X.T@X)@X.T;P0=np.ones((Q,Q))/Q;PH=np.eye(Q)-P0-PV;u=1/Q

def pmul(a,b):
 out=np.zeros(3)
 for r in range(3):
  for s in range(3-r):out[r+s]+=a[r]*b[s]
 return out
@lru_cache(None)
def state(d):
 d=np.asarray(d,float);v=PV@d;w=6*PH@(v*v);Ad=A@d
 return v,w,Ad
@lru_cache(None)
def rates(d,closed,loo):
 d=np.asarray(d,float);v,w,Ad=state(tuple(d));p1=v if closed else d;p2=w if closed else np.zeros(Q)
 R=[np.zeros((Q,Q,3)) for _ in range(3)]
 for i in range(Q):
  if loo:
   t=Ad-B[:,i]
  else:
   t=Ad
  t2=t*t;mt2=np.mean(t2)
  for j in range(Q):
   R[0][i,j,0]=u*u
   R[1][i,j,0]=p1[i]*u
   R[1][i,j,1]=u*t[j]/Q
   R[2][i,j,0]=p2[i]*u
   R[2][i,j,1]=p1[i]*t[j]/Q
   R[2][i,j,2]=u*(t2[j]-mt2)/(2*Q)
 return R

def step(d,i,j):
 x=list(d);x[i]-=1;x[j]+=1;return tuple(x)
def coefficient(phi,loo):
 @lru_cache(None)
 def f0(d):return float(phi@np.asarray(d,float))**4
 @lru_cache(None)
 def F1(d,closed):
  R=rates(d,closed,loo);base=f0(d);out=[np.zeros(3) for _ in range(3)]
  for i in range(Q):
   for j in range(Q):
    if i==j:continue
    df=f0(step(d,i,j))-base
    for r in range(3):out[r]+=R[r][i,j]*df
  return tuple(out)
 @lru_cache(None)
 def F2(d,closed):
  R=rates(d,closed,loo);base=F1(d,closed);out=[np.zeros(3) for _ in range(3)]
  for i in range(Q):
   for j in range(Q):
    if i==j:continue
    ch=F1(step(d,i,j),closed);D=[ch[r]-base[r] for r in range(3)]
    out[0]+=pmul(R[0][i,j],D[0]);out[1]+=pmul(R[0][i,j],D[1])+pmul(R[1][i,j],D[0]);out[2]+=pmul(R[0][i,j],D[2])+pmul(R[1][i,j],D[1])+pmul(R[2][i,j],D[0])
  return tuple(out)
 def C1(closed):
  d=(0,)*Q;R=rates(d,closed,loo);base=F2(d,closed);out=np.zeros(3)
  for i in range(Q):
   for j in range(Q):
    if i==j:continue
    ch=F2(step(d,i,j),closed);D=[ch[r]-base[r] for r in range(3)]
    out += pmul(R[0][i,j],D[2])+pmul(R[1][i,j],D[1])+pmul(R[2][i,j],D[0])
  return out
 return C1(False)-C1(True)

cvec=np.array([.2,-.3,.4,.1,-.5,.6,.2]);cvec=cvec/np.linalg.norm(cvec);phi=X@cvec
cc=coefficient(phi,True);print("mix1_Bcirc",cc.tolist(),flush=True)
