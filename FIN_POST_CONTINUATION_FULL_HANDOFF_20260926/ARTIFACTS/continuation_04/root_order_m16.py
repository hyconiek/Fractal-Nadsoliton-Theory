import math,functools,numpy as np
from scipy.special import softmax,logsumexp
N=12;alpha=0.7341335684322641
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/(1+min(abs(i-j),N-abs(i-j))**1.8) for j in range(N)] for i in range(N)])
A=np.diag(W.sum(1))-W;L=np.fft.fft(A[0]).real[:7];jj=np.arange(N);cols=[]
for k in (3,4,5):cols += [np.sqrt(L[k]/6)*np.cos(2*np.pi*k*jj/N),np.sqrt(L[k]/6)*np.sin(2*np.pi*k*jj/N)]
cols += [np.sqrt(L[6]/12)*(-1.)**jj];X=np.column_stack(cols);th=np.zeros(7);th[[0,2,4,6]]=[1.8199035812800828,1.913989554668724,1.914569132546848,1.367203280195504]
p0=softmax(X@th);mus=np.array([X.T@np.roll(p0,a) for a in range(12)]);V=mus[[1,2,4,5,7,8,10,11]];s=2.0
@functools.lru_cache(None)
def splits(c):
 h=sum(c)//2;out=[]
 def rec(i,rem,a):
  if i==8:
   if rem==0:
    aa=tuple(a);bb=tuple(x-y for x,y in zip(c,aa))
    if aa<=bb:out.append((aa,bb))
   return
  tail=sum(c[i+1:])
  for x in range(max(0,rem-tail),min(c[i],rem)+1):rec(i+1,rem-x,a+[x])
 rec(0,h,[]);return tuple(out)
def delta(a,b):
 na=sum(a);nb=sum(b);ma=np.array(a)@V/na;mb=np.array(b)@V/nb;return na*nb/(na+nb)*float((ma-mb)@(ma-mb))
@functools.lru_cache(None)
def logZ(c,beta):
 if sum(c)==1:return 0.
 xs=[]
 for a,b in splits(c):
  d=delta(a,b);la=logZ(a,beta*s)
  if a!=b:xs.append(-beta*d+la+logZ(b,beta*s))
  else:xs.append(-beta*d+np.logaddexp(2*la,logZ(a,2*beta*s))-math.log(2))
 return float(logsumexp(xs))
root=(2,)*8; beta=alpha/2
# ordered root split weights use child logZ at 2 beta=alpha
rows=[]
def gen(i,rem,a):
 if i==8:
  if rem==0:
   c=tuple(a);b=tuple(2-x for x in c);rows.append((c,b))
  return
 tail=(7-i)*2
 for x in range(max(0,rem-tail),min(2,rem)+1):gen(i+1,rem-x,a+[x])
gen(0,8,[])
lw=[];q=[]
for a,b in rows:
 lw.append(logZ(a,alpha)+logZ(b,alpha));q.append(sum((x-1)**2 for x in a))
lw=np.array(lw);p=np.exp(lw-logsumexp(lw));q=np.array(q);h=np.bincount(q,weights=p)
print('states',len(rows),'hist',[(i,float(v)) for i,v in enumerate(h) if v>1e-10])
peaks=[]
for i in range(1,len(h)-1):
 if h[i]>h[i-1] and h[i]>=h[i+1] and h[i]>1e-4:peaks.append((i,float(h[i])))
print('peaks',sorted(peaks,key=lambda x:x[1],reverse=True))
