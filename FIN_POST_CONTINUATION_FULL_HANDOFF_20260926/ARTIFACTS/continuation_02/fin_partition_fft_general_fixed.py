#!/usr/bin/env python3
import numpy as np, math, ast, sys, time, gc
from scipy import fft
V=np.array(ast.literal_eval(open('/mnt/data/v12_repr.txt').read()),float)
class SA:
    __slots__=('a','s')
    def __init__(self,a,s):self.a=a;self.s=float(s)

def geom(root,axes,dep,m):
    shape=tuple(root[i]+1 for i in axes)
    sm=np.zeros(shape,dtype=np.int16)
    for pos,i in enumerate(axes):
        sh=[1]*len(axes);sh[pos]=root[i]+1
        sm += np.arange(root[i]+1,dtype=np.int16).reshape(sh)
    cd=m-sm;valid=(cd>=0)&(cd<=root[dep])
    norm=np.zeros(shape,float)
    for d in range(V.shape[1]):
        t=(m*V[dep,d])*np.ones(shape,float)
        for pos,i in enumerate(axes):
            sh=[1]*len(axes);sh[pos]=root[i]+1
            t += np.arange(root[i]+1,dtype=float).reshape(sh)*(V[i,d]-V[dep,d])
        norm+=t*t;del t
    return valid,norm,cd

def base(root,axes,dep,beta0,M,L):
    shape=tuple(root[i]+1 for i in axes);nv=np.sum(V*V,axis=1);out={}
    for q in [1<<k for k in range(L+1)]:
        beta=beta0*M*q;logs=-beta*nv/2;s=float(np.max(logs));a=np.zeros(shape,float)
        a[(0,)*len(axes)]=math.exp(logs[dep]-s)
        for pos,i in enumerate(axes):
            idx=[0]*len(axes);idx[pos]=1;a[tuple(idx)]=math.exp(logs[i]-s)
        out[q]=SA(a,s)
    return out

def solve(root,beta0,workers=-1,verbose=False):
    root=tuple(map(int,root));M=sum(root);L=int(round(math.log2(M)));assert 2**L==M
    dep=min(range(len(root)),key=lambda i:root[i]);axes=[i for i in range(len(root)) if i!=dep]
    if not all(root[i]+1>root[dep] for i in axes):raise ValueError('alias-free modulus condition failed')
    shape=tuple(root[i]+1 for i in axes)
    prev=base(root,axes,dep,beta0,M,L)
    if verbose:print('M',M,'dep',dep,'shape',shape,'points',np.prod(shape),flush=True)
    for lev in range(1,L+1):
        m=1<<lev;valid,norm,cd=geom(root,axes,dep,m);new={};maxq=1<<(L-lev)
        for q in [1<<k for k in range(int(math.log2(maxq))+1)]:
            t0=time.time();c=prev[q];d=prev[2*q]
            F=fft.fftn(c.a,workers=workers);conv=fft.ifftn(F*F,workers=workers,overwrite_x=True).real;del F
            mn=float(conv.min());
            if mn<-1e-7:raise RuntimeError(('negative FFT convolution',mn))
            conv[conv<0]=0
            diag=np.zeros(shape,float)
            src=[];tgt=[]
            for i in axes:
                src.append(slice(0,root[i]//2+1));tgt.append(slice(0,root[i]+1,2))
            # Equal-child correction exists only when *all* parent counts are even,
            # including the dependent coordinate cd=m-sum(explicit axes).
            block=d.a[tuple(src)]
            dep_even=(cd[tuple(tgt)]%2)==0
            view=diag[tuple(tgt)]
            view[dep_even]=block[dep_even]
            diag[tuple(tgt)]=view
            S=max(2*c.s,d.s);raw=math.exp(2*c.s-S)*conv+math.exp(d.s-S)*diag;del conv,diag
            beta=beta0*(M/m)*q;fac=beta*norm/(2*m)-math.log(2)
            pos=valid&(raw>0);adj=np.full(shape,-np.inf,float);adj[pos]=np.log(raw[pos])+fac[pos]
            ma=float(np.max(adj[pos]));arr=np.zeros(shape,float);arr[pos]=np.exp(adj[pos]-ma)
            new[q]=SA(arr,S+ma);del raw,adj
            if verbose:print(' lev',lev,'m',m,'q',q,'sec',time.time()-t0,flush=True)
        prev=new;del valid,norm,cd;gc.collect()
    ridx=tuple(root[i] for i in axes);r=prev[1];logY=r.s+math.log(r.a[ridx])
    T=np.array(root)@V;logZ=logY+beta0*float(T@T)/(2*M)
    return logZ,dep,shape
if __name__=='__main__':
    root=tuple(int(x) for x in sys.argv[1].split(','));beta=float(sys.argv[2]);workers=int(sys.argv[3]) if len(sys.argv)>3 else -1
    t=time.time();z,dep,shape=solve(root,beta,workers,True);print('RESULT root',root,'beta',beta,'logZ',repr(z),'dep',dep,'shape',shape,'sec',time.time()-t)
