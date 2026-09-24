import numpy as np

def op(x,f):
    h=(np.roll(x,-1)-x)%(2*np.pi); hm=np.roll(h,1); v=(h+hm)/2
    L=((np.roll(f,-1)-f)/h + (np.roll(f,1)-f)/hm)/v
    return L,h,v

def f(x): return np.sin(3*x)+.3*np.cos(2*x)+.1*np.sin(5*x)
def f2(x): return -9*np.sin(3*x)-1.2*np.cos(2*x)-2.5*np.sin(5*x)
for eps in [0,.35,.7]:
    errs=[]
    for n in [64,128,256,512,1024,2048]:
        s=2*np.pi*np.arange(n)/n; x=s+eps*np.sin(s)
        L,h,v=op(x,f(x)); errs.append(np.max(np.abs(L-f2(x))))
        assert abs(v@L)<1e-10
    assert errs[-1] < errs[0]/100
rng=np.random.default_rng(20260924)
errs=[]
for n in [128,256,512,1024,2048,4096]:
    vals=[]
    for _ in range(12):
        x=np.sort(rng.uniform(0,2*np.pi,n)); L,h,v=op(x,f(x)); vals.append(np.max(np.abs(L-f2(x))))
        assert abs(v@L)<1e-10
    errs.append(np.median(vals))
assert errs[-1] < errs[0]/10
print('NL-02 PASS replay',errs)
