import math, functools, numpy as np
from scipy.special import softmax
N=12
W=np.array([[0.0 if i==j else math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/(1+min(abs(i-j),N-abs(i-j))**1.8) for j in range(N)] for i in range(N)])
A=np.diag(W.sum(1))-W; L=np.fft.fft(A[0]).real[:7]; jj=np.arange(N)
cols=[]
for k in (3,4,5): cols += [np.sqrt(L[k]/6)*np.cos(2*np.pi*k*jj/N),np.sqrt(L[k]/6)*np.sin(2*np.pi*k*jj/N)]
cols += [np.sqrt(L[6]/12)*(-1.)**jj]
X=np.column_stack(cols); th=np.zeros(7); th[[0,2,4,6]]=[1.8199035812800828,1.913989554668724,1.914569132546848,1.367203280195504]
p0=softmax(X@th)
mus=np.array([X.T@np.roll(p0,a) for a in range(12)])
V=mus[[0,3,6,9]]; s=2.0
@functools.lru_cache(None)
def splits(c):
    h=sum(c)//2; out=[]
    def rec(i,rem,a):
        if i==len(c):
            if rem==0:
                aa=tuple(a);bb=tuple(x-y for x,y in zip(c,aa))
                if aa<=bb:out.append((aa,bb))
            return
        tail=sum(c[i+1:]); lo=max(0,rem-tail);hi=min(c[i],rem)
        for x in range(lo,hi+1):rec(i+1,rem-x,a+[x])
    rec(0,h,[]);return tuple(out)
def delta(a,b):
    na=sum(a);nb=sum(b);ma=np.array(a)@V/na;mb=np.array(b)@V/nb
    return na*nb/(na+nb)*float((ma-mb)@(ma-mb))
def lse(xs):
    m=max(xs);return m+math.log(sum(math.exp(x-m) for x in xs))
@functools.lru_cache(None)
def minE(c):
    if sum(c)==1:return 0.,1
    best=1e300;deg=0
    for a,b in splits(c):
        ea,ga=minE(a);eb,gb=minE(b);E=delta(a,b)+s*(ea+eb);g=ga*gb if a!=b else ga*(ga+1)//2
        if E<best-1e-11:best,deg=E,g
        elif abs(E-best)<1e-11:deg+=g
    return best,deg
@functools.lru_cache(None)
def logZ(c,beta):
    if sum(c)==1:return 0.
    xs=[]
    for a,b in splits(c):
        d=delta(a,b);la=logZ(a,beta*s)
        if a!=b: xs.append(-beta*d+la+logZ(b,beta*s))
        else: xs.append(-beta*d+lse([2*la,logZ(a,2*beta*s)])-math.log(2))
    return lse(xs)
def pmin(c,beta):
    e,g=minE(c);return g*math.exp(-beta*e-logZ(c,float(beta)))
for n in (1,2,4,8):
    c=(n,n,n,n);lo,hi=0.,8.
    while pmin(c,hi)<.95:hi*=2
    for _ in range(48):
        mid=(lo+hi)/2
        if pmin(c,mid)<.95:lo=mid
        else:hi=mid
    beta=(lo+hi)/2
    print('n',n,'beta95',beta,'n_beta95',n*beta,'Pmin_alpha5',pmin(c,5/n))
