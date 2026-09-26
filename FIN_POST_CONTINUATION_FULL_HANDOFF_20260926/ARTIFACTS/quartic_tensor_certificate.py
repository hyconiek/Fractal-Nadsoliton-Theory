from functools import lru_cache
from itertools import combinations_with_replacement
import numpy as np, math
Q=12; D=7; u=1/Q; j=np.arange(Q)
# strict operator
W=np.array([[0.0 if a==b else math.cos(0.18575*min(abs(a-b),Q-abs(a-b))+0.1625)/(1+min(abs(a-b),Q-abs(a-b))**1.8) for b in range(Q)] for a in range(Q)],float)
L=np.diag(W.sum(1))-W; lam=np.fft.fft(L[0]).real[:7]
E=[]; sectors=[]
for k in (3,4,5):
    E += [np.sqrt(2/Q)*np.cos(2*np.pi*k*j/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*j/Q)]; sectors += [k,k]
E += [((-1.)**j)/np.sqrt(Q)]; sectors += [6]
E=np.column_stack(E)
A=E@np.diag([lam[k] for k in sectors])@E.T
PV=E@E.T; P0=np.ones((Q,Q))/Q; PH=np.eye(Q)-P0-PV

# degree-4 monomials alpha tuple counts length D
alphas=[]
def rec(rem,pos,cur):
    if pos==D-1:
        alphas.append(tuple(cur+[rem])); return
    for v in range(rem+1): rec(rem-v,pos+1,cur+[v])
rec(4,0,[])
M=len(alphas); idx={a:i for i,a in enumerate(alphas)}
from math import factorial
mult=np.array([factorial(4)/np.prod([factorial(x) for x in a]) for a in alphas],float)

def poly4_linear(x):
    # (z.x)^4
    out=np.empty(M)
    for n,a in enumerate(alphas):
        v=mult[n]
        for i,p in enumerate(a):
            if p: v*=x[i]**p
        out[n]=v
    return out

def pmul(a,b):
    out=np.zeros((3,M))
    for r in range(3):
        for s in range(3-r): out[r+s]+=a[r]*b[s]
    return out

def step(d,a,b):
    x=list(d); x[a]-=1; x[b]+=1; return tuple(x)

@lru_cache(None)
def coeffs(d):
    dv=np.asarray(d,float); v=PV@dv; w=6*PH@(v*v); aa=A@dv
    return v,w,aa,float(np.mean(aa*aa))
@lru_cache(None)
def rates(d,closed):
    v,w,aa,maa=coeffs(d); p1=v if closed else np.asarray(d,float); p2=w if closed else np.zeros(Q)
    R=[np.zeros((Q,Q,3)) for _ in range(3)]
    R[0][:,:,0]=u*u
    R[1][:,:,0]=p1[:,None]*u; R[1][:,:,1]=u*(aa[None,:]/Q)
    R[2][:,:,0]=p2[:,None]*u; R[2][:,:,1]=p1[:,None]*(aa[None,:]/Q); R[2][:,:,2]=u*((aa*aa-maa)[None,:]/(2*Q))
    return R
@lru_cache(None)
def f0(d):
    x=E.T@np.asarray(d,float); return poly4_linear(x)
@lru_cache(None)
def F1(d,closed):
    R=rates(d,closed); base=f0(d); out=np.zeros((3,M))
    for a in range(Q):
        for b in range(Q):
            if a==b: continue
            df=f0(step(d,a,b))-base
            for r in range(3): out[r]+=R[r][a,b,:,None]*df
    return out
@lru_cache(None)
def F2(d,closed):
    R=rates(d,closed); base=F1(d,closed); out=np.zeros((3,M))
    for a in range(Q):
        for b in range(Q):
            if a==b: continue
            Dp=F1(step(d,a,b),closed)-base
            out += pmul(R[a,b] if False else np.zeros((3,M)), np.zeros((3,M))) if False else 0
            for r in range(3):
                for s in range(3-r): out[r+s]+=R[r][a,b,s]*Dp[s] if False else 0
            # explicit convolution scalar g-polys with polynomial vectors
            for rr in range(3):
                for ss in range(3-rr): out[rr+ss]+=R[rr][a,b,ss]*Dp[ss]
    return out

def C1(closed):
    d=(0,)*Q; R=rates(d,closed); base=F2(d,closed); out=np.zeros((3,M))
    for a in range(Q):
        for b in range(Q):
            if a==b: continue
            Dp=F2(step(d,a,b),closed)-base
            for rr in range(3):
                for ss in range(3-rr): out[rr+ss]+=R[rr][a,b,ss]*Dp[ss]
    return out

# WAIT: rate structure R[r] is coefficient in h^r and inside vector is g^0,g^1,g^2.
# F recurrences above incorrectly conflated h and g. Reimplement from scalar reference:
@lru_cache(None)
def F1b(d,closed):
    R=rates(d,closed); base=f0(d); out=[np.zeros((3,M)) for _ in range(3)]
    for a in range(Q):
      for b in range(Q):
       if a==b: continue
       df=f0(step(d,a,b))-base
       for r in range(3): out[r]+=R[r][a,b,:,None]*df
    return tuple(out)
@lru_cache(None)
def F2b(d,closed):
    R=rates(d,closed); base=F1b(d,closed); out=[np.zeros((3,M)) for _ in range(3)]
    for a in range(Q):
      for b in range(Q):
       if a==b: continue
       ch=F1b(step(d,a,b),closed); Dp=[ch[r]-base[r] for r in range(3)]
       # polynomial convolution in g, truncated degree2
       for hr, pairs in enumerate([[(0,0)],[(0,1),(1,0)],[(0,2),(1,1),(2,0)]]):
         for ri,di in pairs:
           for gr in range(3):
             for gd in range(3-gr): out[hr][gr+gd]+=R[ri][a,b,gr]*Dp[di][gd]
    return tuple(out)
def C1b(closed):
    d=(0,)*Q; R=rates(d,closed); base=F2b(d,closed); out=np.zeros((3,M))
    for a in range(Q):
      for b in range(Q):
       if a==b: continue
       ch=F2b(step(d,a,b),closed); Dp=[ch[r]-base[r] for r in range(3)]
       # h^2 term only: (r0,D2)+(r1,D1)+(r2,D0)
       for ri,di in [(0,2),(1,1),(2,0)]:
         for gr in range(3):
           for gd in range(3-gr): out[gr+gd]+=R[ri][a,b,gr]*Dp[di][gd]
    return out
actual=C1b(False)-C1b(True)

# polynomial helpers degree2->degree4
qalphas=[]
def rec2(rem,pos,cur):
    if pos==D-1: qalphas.append(tuple(cur+[rem])); return
    for v in range(rem+1): rec2(rem-v,pos+1,cur+[v])
rec2(2,0,[]); qidx={a:i for i,a in enumerate(qalphas)}; QM=len(qalphas)
mult2=np.array([factorial(2)/np.prod([factorial(x) for x in a]) for a in qalphas],float)
def quad_from_rows(B):
    # vector of nodewise quadratics: B[:,r,s] z_r z_s represented conventional polynomial coeffs
    out=np.zeros((Q,QM))
    for node in range(Q):
      for n,a in enumerate(qalphas):
       inds=[]
       for ii,p in enumerate(a): inds += [ii]*p
       if len(inds)==2:
         r,s=inds; out[node,n]=B[node,r,s]*(2 if r!=s else 1)
    return out

def mul_quad(q1,q2):
    out=np.zeros(M)
    for i,a in enumerate(qalphas):
      if q1[i]==0: continue
      for jj,b in enumerate(qalphas):
       if q2[jj]==0: continue
       c=tuple(a[k]+b[k] for k in range(D)); out[idx[c]]+=q1[i]*q2[jj]
    return out
# phi_node = E z. phi^2 tensor node r,s = E_nr E_ns
Bphi=np.einsum('nr,ns->nrs',E,E)
# phi*(A phi): A on E coord diagonal lamsector
Ad=np.array([lam[k] for k in sectors]); BA=np.einsum('nr,ns,s->nrs',E,E,Ad)
# symmetrize because z_r z_s polynomial
BA=(BA+np.swapaxes(BA,1,2))/2
Hnode=np.einsum('nm,mrs->nrs',PH,Bphi)
Knode=np.einsum('nm,mrs->nrs',PH,BA)
Hq=quad_from_rows(Hnode); Kq=quad_from_rows(Knode)
pred=np.zeros((3,M))
for n in range(Q):
    pred[0] += (-12/Q)*mul_quad(Hq[n],Hq[n])
    pred[1] += (2/Q)*mul_quad(Hq[n],Kq[n])
res=actual-pred
print('M',M,'states_F1',F1b.cache_info(),'states_F2',F2b.cache_info())
for gdeg in range(3):
    imax=int(np.argmax(np.abs(res[gdeg]))); print('gdeg',gdeg,'max_abs',repr(float(np.max(np.abs(res[gdeg])))),'monomial',alphas[imax],'actual',repr(float(actual[gdeg,imax])),'pred',repr(float(pred[gdeg,imax])))
print('g2_actual_max',repr(float(np.max(np.abs(actual[2])))))
# random eval cross-check
rng=np.random.default_rng(123)
for t in range(3):
    z=rng.normal(size=D)
    powers=np.array([np.prod([z[i]**a[i] for i in range(D)]) for a in alphas])
    print('eval',t,(actual@powers).tolist(),(pred@powers).tolist(),((actual-pred)@powers).tolist())
