import numpy as np, math, json, time
from scipy.special import gammaln, logsumexp, softmax
from scipy.optimize import minimize
G=5.145228719489142;Q=12;j=np.arange(Q)
W=np.array([[0.0 if i==k else math.cos(0.18575*min(abs(i-k),Q-abs(i-k))+0.1625)/(1+min(abs(i-k),Q-abs(i-k))**1.8) for k in range(Q)] for i in range(Q)])
Lap=np.diag(W.sum(1))-W; lam=np.fft.fft(Lap[0]).real[:7];cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/Q),np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/Q)]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols);A=X@X.T;THR=2.3812285731142

def phigr(th):
 h=X@th;p=softmax(h);return float(th@th/(2*G)-logsumexp(h)+math.log(12)), th/G-X.T@p
mins=[]
for a in range(12):
 r=minimize(lambda z:phigr(z)[0],G*X[a],jac=lambda z:phigr(z)[1],method='BFGS',options={'gtol':1e-12,'maxiter':2000});mins.append(r.x)
mins=np.array(mins)
def comps(n,k,p=()):
 if k==1: yield p+(n,)
 else:
  for i in range(n+1): yield from comps(n-i,k-1,p+(i,))
def lowstates(N,K=6):return np.array([(N-D,)+m for D in range(K+1) for m in comps(D,11)],dtype=np.int16)
def logw(s,N):return gammaln(N+1)-gammaln(s+1).sum(1)+(G/(2*N))*np.einsum('bi,ij,bj->b',s,A,s)
def corelab(s,N):
 th=(G/N)*(s@X);d=((th[:,None,:]-mins[None,:,:])**2).sum(2);oo=np.argsort(d,axis=1);gap=d[np.arange(len(s)),oo[:,1]]-d[np.arange(len(s)),oo[:,0]];lab=oo[:,0].astype(np.int8);lab[gap<=THR]=-1;return lab,gap
def descend_label(s,N):
 th0=(G/N)*(s@X);r=minimize(lambda z:phigr(z)[0],th0,jac=lambda z:phigr(z)[1],method='BFGS',options={'gtol':1e-12,'maxiter':4000});p=softmax(X@r.x);mx=p.max();return int(p.argmax()) if mx>.9 else -1, float(mx), float(np.linalg.norm(r.jac,np.inf) if hasattr(r,'jac') else np.nan)
out={}
for N,unmass in [(11,0.0006154901411593956),(12,0.00034544952969208946)]:
 allst=np.array(list(comps(N,12)),dtype=np.int16); logZ=logsumexp(logw(allst,N)); ls=lowstates(N); lw=logw(ls,N); lab,gap=corelab(ls,N)
 need=np.where(lab!=0)[0]; changes=[]
 for ii in need:
  dl,mx,gr=descend_label(ls[ii],N);changes.append((int(ii),int(lab[ii]),dl,mx,gr));lab[ii]=dl
 abs_low=float(np.exp(logsumexp(lw[lab==0])-logZ)); den_hi=1/12; den_lo=(1-unmass)/12
 out[str(N)]={'need_descent':len(need),'after_j0':int((lab==0).sum()),'after_unlocalized':int((lab<0).sum()),'after_other':int(((lab>=0)&(lab!=0)).sum()),'abs_low_mass_j0':abs_low,'tail_upper':1-abs_low/den_hi,'tail_lower':max(0,1-abs_low/den_lo),'changes_sample':changes[:50]}
 print(N,out[str(N)])
open('/mnt/data/defect_research_inputs/DEFECT_N11_N12_EXACT_LOW.json','w').write(json.dumps(out,indent=2))
