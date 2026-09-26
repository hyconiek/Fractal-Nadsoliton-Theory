import numpy as np, math
from itertools import product
from collections import defaultdict
from math import factorial
Q=12; D=7; u=1/Q; j=np.arange(Q)
W=np.array([[0.0 if a==b else math.cos(0.18575*min(abs(a-b),Q-abs(a-b))+0.1625)/(1+min(abs(a-b),Q-abs(a-b))**1.8) for b in range(Q)] for a in range(Q)],float)
L=np.diag(W.sum(1))-W; lam=np.fft.fft(L[0]).real[:7]
E=[]; sectors=[]
for k in (3,4,5):
    E += [np.sqrt(2/Q)*np.cos(2*np.pi*k*j/Q),np.sqrt(2/Q)*np.sin(2*np.pi*k*j/Q)]; sectors += [k,k]
E += [((-1.)**j)/np.sqrt(Q)]; sectors += [6]
E=np.column_stack(E); Adiag=np.array([lam[k] for k in sectors]); A=E@np.diag(Adiag)@E.T
PV=E@E.T; P0=np.ones((Q,Q))/Q; PH=np.eye(Q)-P0-PV
zero=(0,)*Q

def step(d,a,b):
    x=list(d); x[a]-=1; x[b]+=1; return tuple(x)

def rate_poly(d,closed,hr,a,b):
    # returns coeffs in g degree 0..2 for h^hr rate p_i q_j
    dv=np.asarray(d,float); v=PV@dv; w=6*PH@(v*v); aa=A@dv; maa=float(np.mean(aa*aa))
    p1=v if closed else dv; p2=w if closed else np.zeros(Q)
    z=np.zeros(3)
    if hr==0: z[0]=u*u
    elif hr==1:
        z[0]=p1[a]*u; z[1]=u*aa[b]/Q
    elif hr==2:
        z[0]=p2[a]*u; z[1]=p1[a]*aa[b]/Q; z[2]=u*(aa[b]*aa[b]-maa)/(2*Q)
    return z

def apply_adjoint(measure,closed,hr):
    # measure: dict state -> g-poly len3; returns after left operator G_hr
    out=defaultdict(lambda:np.zeros(3))
    for d,wp in measure.items():
        for a in range(Q):
            for b in range(Q):
                if a==b: continue
                rp=rate_poly(d,closed,hr,a,b)
                # convolution wp * rp truncated g2
                cp=np.zeros(3)
                for x in range(3):
                    for y in range(3-x): cp[x+y]+=wp[x]*rp[y]
                nd=step(d,a,b)
                out[nd]+=cp; out[d]-=cp
    return dict(out)

def final_measure(closed):
    total=defaultdict(lambda:np.zeros(3))
    tuples=[t for t in product(range(3),repeat=3) if sum(t)==2]
    for t in tuples:
        meas={zero:np.array([1.,0.,0.])}
        # adjoint applies leftmost first for (G_a G_b G_c P)(0)
        for hr in t: meas=apply_adjoint(meas,closed,hr)
        for d,w in meas.items(): total[d]+=w
    return dict(total)
full=final_measure(False); closed=final_measure(True)
weights=defaultdict(lambda:np.zeros(3))
for d,w in full.items(): weights[d]+=w
for d,w in closed.items(): weights[d]-=w
print('state_counts',len(full),len(closed),len(weights))
# drop numerical tiny weights for speed? retain all.
# Build homogeneous degree4 monomial coefficients in retained coordinates.
alphas=[]
def rec(rem,pos,cur):
    if pos==D-1: alphas.append(tuple(cur+[rem])); return
    for v in range(rem+1): rec(rem-v,pos+1,cur+[v])
rec(4,0,[]); M=len(alphas)
mult=np.array([factorial(4)/np.prod([factorial(x) for x in a]) for a in alphas],float)
actual=np.zeros((3,M))
for d,w in weights.items():
    if np.max(np.abs(w))<1e-18: continue
    x=E.T@np.asarray(d,float)
    vals=np.empty(M)
    for n,a in enumerate(alphas):
        v=mult[n]
        for ii,p in enumerate(a):
            if p:v*=x[ii]**p
        vals[n]=v
    actual += w[:,None]*vals[None,:]
# candidate polynomial by sampling monomial coefficients with quadratic dict algebra.
qalphas=[]
def rec2(rem,pos,cur):
    if pos==D-1: qalphas.append(tuple(cur+[rem])); return
    for v in range(rem+1): rec2(rem-v,pos+1,cur+[v])
rec2(2,0,[]); idx4={a:i for i,a in enumerate(alphas)}
def quad_node(mat,node):
    # z^T mat[node] z, mat node,D,D symmetric
    q={}
    for r in range(D):
        for s in range(r,D):
            c=mat[node,r,s]*(2 if r!=s else 1)
            if c:
                a=[0]*D;a[r]+=1;a[s]+=1;q[tuple(a)]=c
    return q
def qmul(q1,q2):
    out=np.zeros(M)
    for a,ca in q1.items():
        for b,cb in q2.items():
            c=tuple(a[i]+b[i] for i in range(D));out[idx4[c]]+=ca*cb
    return out
Bphi=np.einsum('nr,ns->nrs',E,E)
BA=np.einsum('nr,ns,s->nrs',E,E,Adiag);BA=(BA+np.swapaxes(BA,1,2))/2
H=np.einsum('nm,mrs->nrs',PH,Bphi);K=np.einsum('nm,mrs->nrs',PH,BA)
pred=np.zeros((3,M))
for n in range(Q):
    hq=quad_node(H,n); kq=quad_node(K,n)
    pred[0]+=(-12/Q)*qmul(hq,hq); pred[1]+=(2/Q)*qmul(hq,kq)
res=actual-pred
for gd in range(3):
    imax=np.argmax(np.abs(res[gd]));print('gdeg',gd,'max_abs',format(np.max(np.abs(res[gd])),'.17g'),'alpha',alphas[int(imax)],'actual',format(actual[gd,imax],'.17g'),'pred',format(pred[gd,imax],'.17g'))
print('g2_actual_max',format(np.max(np.abs(actual[2])),'.17g'))
# evaluate previously used probes
for name,z in [('mix',np.array([.3,-.2,.4,.1,-.5,.25,.6])),('34',np.array([1,0,.37,0,0,0,0]))]:
    powers=np.array([np.prod([z[i]**a[i] for i in range(D)]) for a in alphas])
    print(name,'actual',(actual@powers).tolist(),'pred',(pred@powers).tolist(),'res',(res@powers).tolist())
np.savez('/mnt/data/FIN_NEXT_RESEARCH_20260925/quartic_tensor_certificate.npz',actual=actual,pred=pred,res=res,alphas=np.array(alphas),lam=lam)
