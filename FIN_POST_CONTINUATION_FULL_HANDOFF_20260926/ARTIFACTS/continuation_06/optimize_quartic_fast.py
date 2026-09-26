import numpy as np, math
from scipy.optimize import minimize
Q=12;j=np.arange(Q);g=3.7183448981203875
W=np.array([[0.0 if a==b else math.cos(.18575*min(abs(a-b),Q-abs(a-b))+.1625)/(1+min(abs(a-b),Q-abs(a-b))**1.8) for b in range(Q)] for a in range(Q)])
L=np.diag(W.sum(1))-W;lam=np.fft.fft(L[0]).real[:7]
E=[];sec=[]
for k in (3,4,5):E += [np.sqrt(2/Q)*np.cos(2*np.pi*k*j/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*j/Q)];sec += [k,k]
E += [(-1.)**j/np.sqrt(Q)];sec += [6]
E=np.column_stack(E);Ad=np.array([lam[k] for k in sec]);PH=np.eye(Q)-np.ones((Q,Q))/Q-E@E.T

def val(z):
 z=np.asarray(z,float);z=z/np.linalg.norm(z);phi=E@z;psi=E@(Ad*z);h=PH@(phi*phi);k=PH@(phi*psi);return -12*np.mean(h*h)+2*g*np.mean(h*k)
def batch(Z):
 # Z rows normalized
 Phi=Z@E.T; Psi=(Z*Ad)@E.T; H=(Phi*Phi)@PH.T;K=(Phi*Psi)@PH.T
 return -12*np.mean(H*H,axis=1)+2*g*np.mean(H*K,axis=1)
rng=np.random.default_rng(777)
keepmax=[];keepmin=[]
for b in range(40):
 Z=rng.normal(size=(50000,7));Z/=np.linalg.norm(Z,axis=1)[:,None];v=batch(Z)
 im=np.argpartition(v,-20)[-20:]; ii=np.argpartition(v,20)[:20]
 keepmax.extend([(float(v[i]),Z[i].copy()) for i in im]);keepmin.extend([(float(v[i]),Z[i].copy()) for i in ii])
keepmax=sorted(keepmax,reverse=True,key=lambda x:x[0])[:100];keepmin=sorted(keepmin,key=lambda x:x[0])[:100]
print('randommax',keepmax[0][0],'randommin',keepmin[0][0])
# unconstrained normalized objective; numerical gradient by scipy BFGS is okay for 200 starts
for name,seeds,sign in [('max',keepmax,-1),('min',keepmin,1)]:
 sols=[]
 for _,x in seeds:
  r=minimize(lambda y:sign*val(y),x,method='BFGS',options={'gtol':1e-10,'maxiter':500})
  z=r.x/np.linalg.norm(r.x);sols.append((val(z),z,r.success))
 sols.sort(key=lambda t:t[0],reverse=(name=='max'))
 v,z,s=sols[0]
 print(name,'best',repr(float(v)),'success',s,'z',z.tolist(),'sector_norms',[float(np.hypot(z[0],z[1])),float(np.hypot(z[2],z[3])),float(np.hypot(z[4],z[5])),abs(float(z[6]))])
 # show distinct values rounded
 print(name,'distinct',sorted(set(round(float(x[0]),12) for x in sols),reverse=(name=='max'))[:10])
print('pure5',val(np.array([0,0,0,0,1,0,0.])))
