#!/usr/bin/env python3
import math
import numpy as np
from scipy.integrate import solve_ivp

N=12

def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],dtype=float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
j=np.arange(N)

hidden=[]
for k in (1,2):
    hidden += [
        np.sqrt(2/N)*np.cos(2*np.pi*k*j/N),
        np.sqrt(2/N)*np.sin(2*np.pi*k*j/N)
    ]
Y=np.column_stack(hidden)

def bridge(u):
    D=np.diag(u)
    return .5*(D@A+A@D)-.5*np.diag(A@u)

Js=[bridge(Y[:,a]) for a in range(4)]

# structural checks
for J in Js:
    assert np.linalg.norm(J-J.T)<1e-12
    assert np.linalg.norm(J@np.ones(N))<1e-12

# D12 equivariance
def Pmat(a,eps):
    P=np.zeros((N,N))
    for x in range(N):
        P[x,(eps*x+a)%N]=1
    return P

err=0.0
for a in range(N):
    for eps in (1,-1):
        P=Pmat(a,eps)
        for k in range(4):
            err=max(err,np.linalg.norm(bridge(P@Y[:,k])-P@Js[k]@P.T))
assert err<1e-12

JG=np.array([[np.trace(Js[a].T@Js[b]) for b in range(4)] for a in range(4)])
assert np.linalg.eigvalsh(JG).min()>1.0

# Profile Jacobians for source 0.
src=0
idx=[x for x in range(N) if x!=src]
w=W[src,idx]
pH=w/w.sum()
p2=w*w/(w@w)

JH=np.zeros((11,4))
JUW=np.zeros((11,4))
for a in range(4):
    g=.5*(Y[src,a]+Y[idx,a])
    JH[:,a]=pH*(g-pH@g)
    JUW[:,a]=2*p2*(g-p2@g)

assert np.linalg.matrix_rank(JH)==4
assert np.linalg.matrix_rank(JUW)==4

# Dynamic hidden fixture.
h=np.array([.7,-.4,.5,.2])
u=Y@h
Jh=bridge(u)
eps=.05
Ah=A+eps*Jh

# Heat profile
deltaW=np.zeros_like(W)
for i in range(N):
    for k in range(N):
        if i!=k:
            deltaW[i,k]=.5*W[i,k]*(u[i]+u[k])
Wh=W+eps*deltaW
r0H=Wh[idx,src]/Wh[idx,src].sum()

def heat_SD(t):
    e=np.zeros(N);e[src]=1.
    def rhs(s,p):
        return -(A+eps*np.exp(-s)*Jh)@p
    sol=solve_ivp(rhs,[0,t],e,method="DOP853",rtol=1e-12,atol=1e-14,t_eval=[t])
    off=sol.y[idx,-1]
    S=off.sum()
    r=off/S
    return S,np.linalg.norm(r-r0H)

# Unitary profile
a0=Ah[idx,src]
r0UW=a0*a0/np.sum(a0*a0)

def unitary_SD(t):
    e=np.zeros(N,dtype=complex);e[src]=1.
    def rhs(s,q):
        return -1j*(A+eps*np.exp(-s)*Jh)@q
    sol=solve_ivp(rhs,[0,t],e,method="DOP853",rtol=1e-12,atol=1e-14,t_eval=[t])
    q=sol.y[:,-1]
    off=np.abs(q[idx])**2
    S=off.sum()
    r=off/S
    return S,np.linalg.norm(r-r0UW)

def wave_SD(t):
    e=np.zeros(N);e[src]=1.
    y0=np.r_[e,np.zeros(N)]
    def rhs(s,y):
        q=y[:N];v=y[N:]
        return np.r_[v,-(A+eps*np.exp(-s)*Jh)@q]
    sol=solve_ivp(rhs,[0,t],y0,method="DOP853",rtol=1e-12,atol=1e-14,t_eval=[t])
    q=sol.y[:N,-1]
    off=q[idx]**2
    S=off.sum()
    r=off/S
    return S,np.linalg.norm(r-r0UW)

def slope(fun,ts):
    a=np.array([fun(t) for t in ts])
    return float(np.polyfit(np.log(a[:,0]),np.log(a[:,1]),1)[0])

gH=slope(heat_SD,np.logspace(-5,-4,8))
gU=slope(unitary_SD,np.logspace(-4,-3.3,8))
gW=slope(wave_SD,np.logspace(-3,-2.3,8))

assert abs(gH-1)<5e-3
assert abs(gU-.5)<3e-2
assert abs(gW-.25)<4e-2

print("PASS")
print("D12_covariance_error",err)
print("bridge_Gram_eigenvalues",np.linalg.eigvalsh(JG).tolist())
print("heat_profile_singular_values",np.linalg.svd(JH,compute_uv=False).tolist())
print("unitary_wave_profile_singular_values",np.linalg.svd(JUW,compute_uv=False).tolist())
print("dynamic_hidden_gamma_heat",gH)
print("dynamic_hidden_gamma_unitary",gU)
print("dynamic_hidden_gamma_wave",gW)
