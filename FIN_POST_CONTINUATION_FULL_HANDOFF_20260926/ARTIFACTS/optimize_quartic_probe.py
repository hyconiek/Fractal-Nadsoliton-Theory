import numpy as np, math
from scipy.optimize import minimize
Q=12; j=np.arange(Q)
# strict lambda
W=np.array([[0.0 if a==b else math.cos(.18575*min(abs(a-b),Q-abs(a-b))+.1625)/(1+min(abs(a-b),Q-abs(a-b))**1.8) for b in range(Q)] for a in range(Q)])
L=np.diag(W.sum(1))-W; lam=np.fft.fft(L[0]).real[:7]
E=[]; sectors=[]
for k in (3,4,5):
 E += [np.sqrt(2/Q)*np.cos(2*np.pi*k*j/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*j/Q)];sectors += [k,k]
E += [(-1.)**j/np.sqrt(Q)]; sectors += [6]
E=np.column_stack(E); Ad=np.array([lam[k] for k in sectors]);
PV=E@E.T;PH=np.eye(Q)-np.ones((Q,Q))/Q-PV
g=3.7183448981203875

def parts(z):
 z=np.asarray(z,float); z=z/np.linalg.norm(z); phi=E@z
 h=PH@(phi*phi);k=PH@(phi*(E@(Ad*z)))
 c0=-12*np.mean(h*h);c1=2*g*np.mean(h*k)
 # Fourier hidden energy split using DFT
 hh=np.fft.fft(h)/Q
 n1=2*abs(hh[1])**2; n2=2*abs(hh[2])**2
 return c0+c1,c0,c1,n1,n2,z

def obj(z,mode):
 c,*_=parts(z)
 return {'max':-c,'min':c,'abs':-abs(c)}[mode]
rng=np.random.default_rng(20260925)
for mode in ['max','min','abs']:
 best=None
 for _ in range(80):
  x=rng.normal(size=7); x/=np.linalg.norm(x)
  r=minimize(lambda y:obj(y,mode),x,method='BFGS',options={'gtol':1e-11,'maxiter':1000})
  tup=parts(r.x)
  score={'max':tup[0],'min':-tup[0],'abs':abs(tup[0])}[mode]
  if best is None or score>best[0]: best=(score,tup,r.fun,r.success)
 score,tup,_,success=best;c,c0,c1,n1,n2,z=tup
 print(mode,'C',repr(c),'c0',repr(c0),'c1',repr(c1),'hidden1',repr(n1),'hidden2',repr(n2),'z',z.tolist(),'success',success)
# pure strict real-sector unit E basis directions
for idx,label in enumerate(['3c','3s','4c','4s','5c','5s','6']):
 z=np.zeros(7);z[idx]=1;print('pure',label,parts(z)[:5])
