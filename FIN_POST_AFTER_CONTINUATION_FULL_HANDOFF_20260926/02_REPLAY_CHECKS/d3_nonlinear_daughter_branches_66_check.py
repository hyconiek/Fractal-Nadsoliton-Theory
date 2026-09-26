#!/usr/bin/env python3
import numpy as np, math
from scipy.optimize import root
from scipy.special import softmax
from scipy.linalg import null_space

N=12
def strict_kernel():
    return np.array([[0.0 if i==j else
        math.cos(0.18575*min(abs(i-j),N-abs(i-j))+0.1625)/
        (1+min(abs(i-j),N-abs(i-j))**1.8)
        for j in range(N)] for i in range(N)],float)

W=strict_kernel()
A=np.diag(W.sum(axis=1))-W
lam=np.fft.fft(A[0]).real[:7]
j=np.arange(N)

cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(lam[k]/6)*np.cos(2*np.pi*k*j/N),
        np.sqrt(lam[k]/6)*np.sin(2*np.pi*k*j/N)
    ]
cols += [np.sqrt(lam[6]/12)*(-1.)**j]
X=np.column_stack(cols)

l3,l4,l5,l6=lam[3],lam[4],lam[5],lam[6]

def crossing_eq(v):
    J,K,g=v
    den=np.cosh(J)+np.exp(-2*K)
    x=np.sinh(J)/den
    y=(np.cosh(J)-np.exp(-2*K))/den
    D=(1/g-l4/12)*(1/g-l5/12)-l4*l5*x*x/144
    return np.array([J-g*l3*x/6,K-g*l6*y/12,D])

rr=root(crossing_eq,[.038434754,.1852720614,5.17184183194],tol=1e-13)
assert np.linalg.norm(crossing_eq(rr.x),np.inf)<1e-13
J,K,g=rr.x

theta=np.zeros(7)
theta[0]=J/np.sqrt(l3/6)
theta[6]=K/np.sqrt(l6/12)
p=softmax(X@theta)
S=np.diag(p)-np.outer(p,p)
H=np.eye(7)/g-X.T@S@X

# Explicit critical basis.
Bc=H[np.ix_([2,4],[2,4])]
Bs=H[np.ix_([3,5],[3,5])]
ec=np.linalg.eigh(Bc)[1][:,0]
es=np.linalg.eigh(Bs)[1][:,0]

# Orient to the report convention.
if ec[0]<0: ec=-ec
if es[0]<0: es=-es

ex=np.zeros(7); ex[[2,4]]=ec
ey=np.zeros(7); ey[[3,5]]=es

assert np.linalg.norm(H@ex)<1e-12
assert np.linalg.norm(H@ey)<1e-12
assert abs(ex@ey)<1e-14

Xc=X-np.ones((N,1))*(p@X)
ux=Xc@ex
uy=Xc@ey

D3xxx=-np.sum(p*ux**3)
D3xyy=-np.sum(p*ux*uy**2)
D3xxy=-np.sum(p*ux**2*uy)
D3yyy=-np.sum(p*uy**3)

assert abs(D3xxx+D3xyy)<1e-12
assert abs(D3xxy)<1e-12
assert abs(D3yyy)<1e-12
assert abs(D3xxx)>1e-3

# Quartic Lyapunov-Schmidt coefficient.
C=np.column_stack([ex,ey])
Sb=null_space(C.T)
Hs=Sb.T@H@Sb

def tpair(u,v):
    return -np.sum(p[:,None]*(u*v)[:,None]*Xc,axis=0)

txx=tpair(ux,ux)
tyy=tpair(uy,uy)
txy=tpair(ux,uy)
invHs=np.linalg.inv(Hs)

def sinv(t1,t2):
    a=Sb.T@t1; b=Sb.T@t2
    return a@invHs@b

Ex2=np.sum(p*ux**2)
Ey2=np.sum(p*uy**2)
Exy=np.sum(p*ux*uy)
D4xxxx=-(np.sum(p*ux**4)-3*Ex2**2)
D4xxyy=-(np.sum(p*ux**2*uy**2)-Ex2*Ey2-2*Exy**2)

Rxxxx=D4xxxx-3*sinv(txx,txx)
Rxxyy=D4xxyy-(sinv(txx,tyy)+2*sinv(txy,txy))

assert Rxxxx>1.0
assert abs(Rxxyy-Rxxxx/3)<1e-10

c=D3xxx/6
beta=Rxxxx/24

# Eigenvalue parameter derivative from MP7-039 determinant derivative.
delta_prime=0.001695886924714291
other=np.trace(Bc)
ell=delta_prime/other
assert ell>0

# Direct branch replay for +/-1e-4.
def F(th,gg):
    pp=softmax(X@th)
    return th/gg-X.T@pp

def HH(th,gg):
    pp=softmax(X@th)
    SS=np.diag(pp)-np.outer(pp,pp)
    return np.eye(7)/gg-X.T@SS@X

def base(gg):
    z0=theta[[0,6]]
    def f(z):
        th=np.zeros(7); th[[0,6]]=z
        return F(th,gg)[[0,6]]
    z=root(f,z0,tol=1e-13)
    th=np.zeros(7); th[[0,6]]=z.x
    return th

Alead=abs(ell/(3*c))
rows=[]
for mupar in (1e-4,-1e-4):
    gg=g+mupar
    thb=base(gg)
    phis=[0,2*np.pi/3,4*np.pi/3] if mupar>0 else [np.pi/3,np.pi,5*np.pi/3]
    for ph in phis:
        guess=thb+Alead*abs(mupar)*(np.cos(ph)*ex+np.sin(ph)*ey)
        z=root(lambda th:F(th,gg),guess,tol=1e-12)
        assert np.linalg.norm(F(z.x,gg),np.inf)<1e-10
        ev=np.linalg.eigvalsh(HH(z.x,gg))
        assert np.sum(ev<-1e-8)==2
        d=z.x-thb
        rows.append((mupar,float(d@ex),float(d@ey),float(np.linalg.norm(d)),2))

print("PASS")
print("crossing",rr.x.tolist())
print("critical_cubic_c",float(c))
print("reduced_quartic_beta",float(beta))
print("soft_eigenvalue_slope",float(ell))
print("daughter_amplitude_coefficient",float(Alead))
print("daughter_energy_mu3_coefficient",float(ell**3/(54*c**2)))
print("daughter_replay",rows)
