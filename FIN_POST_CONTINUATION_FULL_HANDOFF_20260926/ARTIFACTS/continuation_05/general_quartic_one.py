from functools import lru_cache
import math, numpy as np
Q=12;j=np.arange(Q)
W=np.array([[0.0 if a==b else math.cos(0.18575*min(abs(a-b),Q-abs(a-b))+0.1625)/(1+min(abs(a-b),Q-abs(a-b))**1.8) for b in range(Q)] for a in range(Q)],float)
A=np.diag(W.sum(1))-W;lam=np.fft.fft(A[0]).real[:7]
cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols);A7=X@X.T;PV=X@np.linalg.inv(X.T@X)@X.T;P0=np.ones((Q,Q))/Q;PH=np.eye(Q)-P0-PV;u=1/Q

def pmul(a,b):
 out=np.zeros(3)
 for r in range(3):
  for s in range(3-r):out[r+s]+=a[r]*b[s]
 return out

def coefficient(phi):
 @lru_cache(None)
 def coeffs(d):
  d=np.asarray(d,float);v=PV@d;w=6*PH@(v*v);aa=A7@d;return v,w,aa,float(np.mean(aa*aa))
 @lru_cache(None)
 def rates(d,closed):
  v,w,aa,maa=coeffs(d);p1=v if closed else np.asarray(d,float);p2=w if closed else np.zeros(Q)
  R=[np.zeros((Q,Q,3)) for _ in range(3)];R[0][:,:,0]=u*u
  R[1][:,:,0]=p1[:,None]*u;R[1][:,:,1]=u*(aa[None,:]/12)
  R[2][:,:,0]=p2[:,None]*u;R[2][:,:,1]=p1[:,None]*(aa[None,:]/12);R[2][:,:,2]=u*((aa*aa-maa)[None,:]/24)
  return R
 def step(d,a,b):
  x=list(d);x[a]-=1;x[b]+=1;return tuple(x)
 @lru_cache(None)
 def f0(d):return float(phi@np.asarray(d,float))**4
 @lru_cache(None)
 def F1(d,closed):
  R=rates(d,closed);base=f0(d);out=[np.zeros(3) for _ in range(3)]
  for a in range(Q):
   for b in range(Q):
    if a==b:continue
    df=f0(step(d,a,b))-base
    for r in range(3):out[r]+=R[r][a,b]*df
  return tuple(out)
 @lru_cache(None)
 def F2(d,closed):
  R=rates(d,closed);base=F1(d,closed);out=[np.zeros(3) for _ in range(3)]
  for a in range(Q):
   for b in range(Q):
    if a==b:continue
    ch=F1(step(d,a,b),closed);D=[ch[r]-base[r] for r in range(3)]
    out[0]+=pmul(R[0][a,b],D[0]);out[1]+=pmul(R[0][a,b],D[1])+pmul(R[1][a,b],D[0]);out[2]+=pmul(R[0][a,b],D[2])+pmul(R[1][a,b],D[1])+pmul(R[2][a,b],D[0])
  return tuple(out)
 def C1(closed):
  d=(0,)*Q;R=rates(d,closed);base=F2(d,closed);out=np.zeros(3)
  for a in range(Q):
   for b in range(Q):
    if a==b:continue
    ch=F2(step(d,a,b),closed);D=[ch[r]-base[r] for r in range(3)]
    out+=pmul(R[0][a,b],D[2])+pmul(R[1][a,b],D[1])+pmul(R[2][a,b],D[0])
  return out
 return C1(False)-C1(True)

def pred(phi):
 h=PH@(phi*phi); k=PH@(phi*(A7@phi))
 return np.array([-12*np.mean(h*h),2*np.mean(h*k),0.])

tests={'mixed_all':X@np.array([.3,-.2,.4,.1,-.5,.25,.6])}
for name,f in tests.items():
 got=coefficient(f);pr=pred(f);res=got-pr
 print(name,'got',got.tolist(),'pred',pr.tolist(),'res',res.tolist(),'max',float(np.max(np.abs(res))))
