from functools import lru_cache
import numpy as np, math, random
Q=12; j=np.arange(Q); u=1/Q
# orthonormal retained real Fourier basis (Euclidean norm 1)
E=[]; sectors=[]
for k in (3,4,5):
    E += [np.sqrt(2/Q)*np.cos(2*np.pi*k*j/Q), np.sqrt(2/Q)*np.sin(2*np.pi*k*j/Q)]
    sectors += [k,k]
E += [((-1.)**j)/np.sqrt(Q)]; sectors += [6]
E=np.column_stack(E)
PV=E@E.T; P0=np.ones((Q,Q))/Q; PH=np.eye(Q)-P0-PV

def pmul(a,b):
    out=np.zeros(3)
    for r in range(3):
        for s in range(3-r):out[r+s]+=a[r]*b[s]
    return out

def coefficient(phi,A):
    @lru_cache(None)
    def coeffs(d):
        d=np.asarray(d,float); v=PV@d; w=6*PH@(v*v); aa=A@d
        return v,w,aa,float(np.mean(aa*aa))
    @lru_cache(None)
    def rates(d,closed):
        v,w,aa,maa=coeffs(d); p1=v if closed else np.asarray(d,float); p2=w if closed else np.zeros(Q)
        R=[np.zeros((Q,Q,3)) for _ in range(3)]
        R[0][:,:,0]=u*u
        R[1][:,:,0]=p1[:,None]*u; R[1][:,:,1]=u*(aa[None,:]/Q)
        R[2][:,:,0]=p2[:,None]*u; R[2][:,:,1]=p1[:,None]*(aa[None,:]/Q); R[2][:,:,2]=u*((aa*aa-maa)[None,:]/(2*Q))
        # Since q=softmax(h g A d): q1=(A d)/Q, q2=((Ad)^2-mean)/2Q.
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

def pred(phi,A):
    h=PH@(phi*phi); k=PH@(phi*(A@phi))
    return np.array([-12*np.mean(h*h),2*np.mean(h*k),0.])

rng=np.random.default_rng(20260925)
worst=0
for case in range(12):
    z=rng.normal(size=7)
    phi=E@z
    # arbitrary diagonal self-adjoint retained operator, repeated cos/sin per sector
    vals={3:rng.uniform(.2,3),4:rng.uniform(.2,3),5:rng.uniform(.2,3),6:rng.uniform(.2,3)}
    diag=np.array([vals[k] for k in sectors])
    A=E@np.diag(diag)@E.T
    got=coefficient(phi,A); pr=pred(phi,A); err=np.max(np.abs(got-pr));worst=max(worst,err)
    print(case, 'err',repr(float(err)),'g2',repr(float(got[2])), 'got', got.tolist(), 'pred',pr.tolist())
print('worst',repr(float(worst)))
