import sys,math,functools,numpy as np
from scipy.special import softmax
N=int(sys.argv[1]) if len(sys.argv)>1 else 120
mode=int(sys.argv[2]) if len(sys.argv)>2 else 2
eps=float(sys.argv[3]) if len(sys.argv)>3 else 1e-5
Q=12;jj=np.arange(Q);th=2*np.pi*jj/Q
W=np.array([[0.0 if i==j else math.cos(.18575*min(abs(i-j),Q-abs(i-j))+.1625)/(1+min(abs(i-j),Q-abs(i-j))**1.8) for j in range(Q)] for i in range(Q)])
A=np.diag(W.sum(1))-W; lam=np.fft.fft(A[0]).real[:7]
cols=[]
for k in (3,4,5):cols += [np.sqrt(lam[k]/6)*np.cos(k*th),np.sqrt(lam[k]/6)*np.sin(k*th)]
cols += [np.sqrt(lam[6]/12)*(-1.)**jj]
X=np.column_stack(cols);A7=X@X.T;u=np.ones(Q)/Q;g=3.7183448981203875
# translated k5 probes, exact Z12 translations
phis=np.array([np.sqrt(lam[5]/6)*np.cos(5*(th+2*np.pi*r/Q)) for r in range(Q)])
hidden=np.cos(mode*th) # unnormalized declared perturbation

def run(sign):
    base=u+sign*eps*hidden
    @functools.lru_cache(None)
    def rec(depth,off):
        o=np.asarray(off,dtype=float); p=base+o/N; mu=X.T@p
        if depth==0:
            x=np.sqrt(N)*(phis@p)
            return x**4
        q=softmax(g*(A7@p)); cur=rec(depth-1,off); out=np.zeros(Q)
        for i in range(Q):
            if p[i]==0:continue
            for k in range(Q):
                if i==k:continue
                rate=N*p[i]*q[k]
                oo=list(off);oo[i]-=1;oo[k]+=1
                out += rate*(rec(depth-1,tuple(oo))-cur)
        return out
    z=(0,)*Q
    ans=rec(3,z)
    print('cache',sign,rec.cache_info(),file=sys.stderr)
    return ans
plus=run(1);minus=run(-1)
der=N*(plus-minus)/(2*eps)
print('N',N,'mode',mode,'eps',eps)
for r,x in enumerate(der):print(r,repr(float(x)))
# DFT real coefficients in translation index
for m in range(7):
 c=2/Q*np.sum(der*np.cos(2*np.pi*m*np.arange(Q)/Q)) if m else np.mean(der)
 s=2/Q*np.sum(der*np.sin(2*np.pi*m*np.arange(Q)/Q)) if m else 0.0
 print('harm',m,'cos',repr(float(c)),'sin',repr(float(s)))
