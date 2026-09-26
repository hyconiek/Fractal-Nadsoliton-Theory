#!/usr/bin/env python3
import numpy as np, math
from itertools import permutations
from scipy.special import softmax
from mpmath import iv, mp

mp.dps=70
iv.dps=55
N=12

L={
3:iv.mpf(["1.96140686197643","1.96140686197645"]),
4:iv.mpf(["2.1995688493332","2.19956884933322"]),
5:iv.mpf(["2.29860627207908","2.2986062720791"]),
6:iv.mpf(["2.34218204114629","2.34218204114631"])
}
LM={3:1.96140686197644,4:2.19956884933321,
    5:2.29860627207909,6:2.34218204114630}

x0=np.array([0.064551183234886,
             0.005557251343608,
             0.013192678416272,
             0.418079565489821,
             5.172231474684087])
rad=2e-8
DET_SCALE=1e6

def det_iv(M):
    n=len(M); total=iv.mpf(0)
    for perm in permutations(range(n)):
        inv=sum(1 for i in range(n) for j in range(i+1,n)
                if perm[i]>perm[j])
        t=iv.mpf(-1 if inv%2 else 1)
        for i in range(n):
            t*=M[i][perm[i]]
        total+=t
    return total

def minor(M,ri,cj):
    return [[M[i][j] for j in range(len(M)) if j!=cj]
            for i in range(len(M)) if i!=ri]

def adj_iv(M):
    n=len(M); A=[[None]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            c=det_iv(minor(M,j,i))
            if (i+j)%2: c=-c
            A[i][j]=c
    return A

def feats_iv():
    rows=[]
    for j in range(N):
        row=[]
        for k in (3,4,5):
            sc=iv.sqrt(L[k]/6)
            c=mp.cos(2*mp.pi*k*j/12)
            row.append(sc*iv.mpf([str(c),str(c)]))
        row.append(iv.sqrt(L[6]/12)*
                   iv.mpf([str((-1)**j),str((-1)**j)]))
        rows.append(row)
    return rows

def eval_GJ(box):
    s=box[:4];g=box[4]
    ft=feats_iv()
    h=[sum(s[a]*ft[j][a] for a in range(4)) for j in range(N)]
    w=[iv.exp(z) for z in h]; Z=sum(w); p=[z/Z for z in w]
    m=[sum(p[j]*ft[j][a] for j in range(N)) for a in range(4)]
    c=[[ft[j][a]-m[a] for a in range(4)] for j in range(N)]
    C=[[sum(p[j]*c[j][a]*c[j][b] for j in range(N))
        for b in range(4)] for a in range(4)]
    H=[[(iv.mpf(1)/g if a==b else iv.mpf(0))-C[a][b]
        for b in range(4)] for a in range(4)]
    F=[s[a]/g-m[a] for a in range(4)]
    dh=det_iv(H)

    # third derivative of Phi = - third centered moment
    T=[[[-sum(p[j]*c[j][a]*c[j][b]*c[j][k] for j in range(N))
          for k in range(4)] for b in range(4)] for a in range(4)]

    Adj=adj_iv(H)
    J=[[None]*5 for _ in range(5)]
    for a in range(4):
        for b in range(4): J[a][b]=H[a][b]
        J[a][4]=-s[a]/(g*g)
    for k in range(4):
        z=iv.mpf(0)
        for a in range(4):
            for b in range(4):
                z+=Adj[b][a]*T[a][b][k]
        J[4][k]=DET_SCALE*z
    J[4][4]=DET_SCALE*(-sum(Adj[a][a] for a in range(4))/(g*g))
    G=F+[DET_SCALE*dh]
    return G,J,H,p,ft

def Gnum(x):
    s=x[:4];g=x[4]
    j=np.arange(N)
    ft=np.column_stack([
        np.sqrt(LM[3]/6)*np.cos(2*np.pi*3*j/N),
        np.sqrt(LM[4]/6)*np.cos(2*np.pi*4*j/N),
        np.sqrt(LM[5]/6)*np.cos(2*np.pi*5*j/N),
        np.sqrt(LM[6]/12)*(-1.)**j])
    p=softmax(ft@s);m=p@ft
    H=np.eye(4)/g-(ft.T@(p[:,None]*ft)-np.outer(m,m))
    return np.r_[s/g-m,DET_SCALE*np.linalg.det(H)]

def fdjac(f,x):
    y=f(x);J=np.zeros((len(y),len(x)))
    for k in range(len(x)):
        h=1e-7*max(1,abs(x[k]))
        e=np.zeros_like(x);e[k]=h
        J[:,k]=(f(x+e)-f(x-e))/(2*h)
    return J

R=np.linalg.inv(fdjac(Gnum,x0))
XI=[iv.mpf([str(v-rad),str(v+rad)]) for v in x0]
X0=[iv.mpf(str(v)) for v in x0]

G0,_,_,_,_=eval_GJ(X0)
_,JX,HX,pX,ftX=eval_GJ(XI)

Kr=[]
for i in range(5):
    base=X0[i]-sum(iv.mpf(str(R[i,j]))*G0[j] for j in range(5))
    corr=iv.mpf(0)
    for j in range(5):
        coeff=(iv.mpf(1) if i==j else iv.mpf(0))-sum(
            iv.mpf(str(R[i,k]))*JX[k][j] for k in range(5))
        corr+=coeff*(XI[j]-X0[j])
    Kr.append(base+corr)

for i in range(5):
    assert float(Kr[i].a)>float(XI[i].a)
    assert float(Kr[i].b)<float(XI[i].b)

detJG=det_iv(JX)
assert float(detJG.b)<0

# C4 complement inertia via midpoint eigenbasis + interval Gershgorin.
j=np.arange(N)
fc=np.column_stack([
    np.sqrt(LM[3]/6)*np.cos(2*np.pi*3*j/N),
    np.sqrt(LM[4]/6)*np.cos(2*np.pi*4*j/N),
    np.sqrt(LM[5]/6)*np.cos(2*np.pi*5*j/N),
    np.sqrt(LM[6]/12)*(-1.)**j])
p=softmax(fc@x0[:4]);m=p@fc
Hm=np.eye(4)/x0[4]-(fc.T@(p[:,None]*fc)-np.outer(m,m))
ev,U=np.linalg.eigh(Hm)
zero=int(np.argmin(np.abs(ev)))
keep=[i for i in range(4) if i!=zero]
Q=U[:,keep]

M=[[iv.mpf(0) for _ in keep] for __ in keep]
for a in range(3):
    for b in range(3):
        z=iv.mpf(0)
        for i in range(4):
            for k in range(4):
                z+=iv.mpf(str(Q[i,a]))*HX[i][k]*iv.mpf(str(Q[k,b]))
        M[a][b]=z

disc=[]
for a in range(3):
    rr=sum(max(abs(float(M[a][b].a)),abs(float(M[a][b].b)))
           for b in range(3) if b!=a)
    disc.append([float(M[a][a].a)-rr,float(M[a][a].b)+rr])

assert disc[0][1]<0
assert disc[1][0]>0 and disc[2][0]>0

# Odd block and its interval Gershgorin spectrum.
odd=[]
for jj in range(N):
    row=[]
    for k in (3,4,5):
        sc=iv.sqrt(L[k]/6)
        s=mp.sin(2*mp.pi*k*jj/12)
        row.append(sc*iv.mpf([str(s),str(s)]))
    odd.append(row)

om=[sum(pX[jj]*odd[jj][a] for jj in range(N)) for a in range(3)]
Ho=[[None]*3 for _ in range(3)]
for a in range(3):
    for b in range(3):
        cv=sum(pX[jj]*(odd[jj][a]-om[a])*(odd[jj][b]-om[b])
               for jj in range(N))
        Ho[a][b]=(iv.mpf(1)/XI[4] if a==b else iv.mpf(0))-cv

fo=np.column_stack([
    np.sqrt(LM[3]/6)*np.sin(2*np.pi*3*j/N),
    np.sqrt(LM[4]/6)*np.sin(2*np.pi*4*j/N),
    np.sqrt(LM[5]/6)*np.sin(2*np.pi*5*j/N)])
po=softmax(fc@x0[:4]);mo=po@fo
Hom=np.eye(3)/x0[4]-(fo.T@(po[:,None]*fo)-np.outer(mo,mo))
eo,Uo=np.linalg.eigh(Hom)

Mo=[[iv.mpf(0) for _ in range(3)] for __ in range(3)]
for a in range(3):
    for b in range(3):
        z=iv.mpf(0)
        for i in range(3):
            for k in range(3):
                z+=iv.mpf(str(Uo[i,a]))*Ho[i][k]*iv.mpf(str(Uo[k,b]))
        Mo[a][b]=z

odisc=[]
for a in range(3):
    rr=sum(max(abs(float(Mo[a][b].a)),abs(float(Mo[a][b].b)))
           for b in range(3) if b!=a)
    odisc.append([float(Mo[a][a].a)-rr,float(Mo[a][a].b)+rr])
assert all(d[0]>0 for d in odisc)

print("PASS")
print("krawczyk_image",[[float(z.a),float(z.b)] for z in Kr])
print("augmented_jacobian_det_interval",
      [float(detJG.a),float(detJG.b)])
print("C4_nonzero_eigenvalue_gershgorin",disc)
print("odd_eigenvalue_gershgorin",odisc)
