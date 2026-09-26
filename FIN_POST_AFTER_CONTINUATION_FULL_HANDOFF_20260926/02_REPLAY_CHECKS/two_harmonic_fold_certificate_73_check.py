#!/usr/bin/env python3
import numpy as np, sympy as sp
from mpmath import iv, mp

mp.dps=70
iv.dps=55

# Accepted spectral intervals.
L3=iv.mpf(["1.96140686197643","1.96140686197645"])
L6=iv.mpf(["2.34218204114629","2.34218204114631"])
L4=iv.mpf(["2.1995688493332","2.19956884933322"])
L5=iv.mpf(["2.29860627207908","2.2986062720791"])

x0=np.array([0.7185052460056293,
             0.4909088546140455,
             4.621196599489125])
rad=np.array([2e-7,2e-7,2e-7])

J,K,g,l3,l6=sp.symbols("J K g l3 l6")
eJ=sp.exp(J); emJ=sp.exp(-J); q=sp.exp(-2*K)
C=(eJ+emJ)/2; S=(eJ-emJ)/2; D=C+q
xx=S/D
yy=(C-q)/D
F1=J-g*l3*xx/6
F2=K-g*l6*yy/12
A11=sp.diff(F1,J);A12=sp.diff(F1,K)
A21=sp.diff(F2,J);A22=sp.diff(F2,K)
detA=sp.simplify(A11*A22-A12*A21)
G=sp.Matrix([F1,F2,detA])
JG=G.jacobian([J,K,g])

mods=[{"exp":iv.exp}]
Gf=[sp.lambdify((J,K,g,l3,l6),z,modules=mods) for z in G]
Jf=[[sp.lambdify((J,K,g,l3,l6),JG[i,j],modules=mods)
     for j in range(3)] for i in range(3)]

# Numeric preconditioner.
GJn=sp.lambdify((J,K,g,l3,l6),JG,"numpy")
l3m=1.96140686197644
l6m=2.3421820411463
R=np.linalg.inv(np.array(GJn(*x0,l3m,l6m),dtype=float))

XI=[iv.mpf([str(x0[i]-rad[i]),str(x0[i]+rad[i])])
    for i in range(3)]
X0=[iv.mpf(str(v)) for v in x0]

F0=[f(*X0,L3,L6) for f in Gf]
JX=[[Jf[i][j](*XI,L3,L6) for j in range(3)] for i in range(3)]

Kr=[]
for i in range(3):
    base=X0[i]-sum(iv.mpf(str(R[i,j]))*F0[j] for j in range(3))
    corr=iv.mpf(0)
    for j in range(3):
        coeff=(iv.mpf(1) if i==j else iv.mpf(0))-sum(
            iv.mpf(str(R[i,k]))*JX[k][j] for k in range(3))
        corr += coeff*(XI[j]-X0[j])
    Kr.append(base+corr)

for i in range(3):
    assert float(Kr[i].a)>float(XI[i].a)
    assert float(Kr[i].b)<float(XI[i].b)

# Fold transversality using unnormalized null vectors.
HF1=sp.hessian(F1,(J,K))
HF2=sp.hessian(F2,(J,K))
HF1f=[[sp.lambdify((J,K,g,l3,l6),HF1[i,j],modules=mods)
       for j in range(2)] for i in range(2)]
HF2f=[[sp.lambdify((J,K,g,l3,l6),HF2[i,j],modules=mods)
       for j in range(2)] for i in range(2)]
Fgf=[sp.lambdify((J,K,g,l3,l6),sp.diff(z,g),modules=mods)
     for z in (F1,F2)]

a,b=JX[0][0],JX[0][1]
c,d=JX[1][0],JX[1][1]
vr=[-b,a]
wl=[-c,a]
FG=[f(*XI,L3,L6) for f in Fgf]

af=wl[0]*FG[0]+wl[1]*FG[1]

def quad(M,v):
    return sum(v[i]*M[i][j]*v[j] for i in range(2) for j in range(2))

H1=[[HF1f[i][j](*XI,L3,L6) for j in range(2)] for i in range(2)]
H2=[[HF2f[i][j](*XI,L3,L6) for j in range(2)] for i in range(2)]
bv=[quad(H1,vr),quad(H2,vr)]
bf=wl[0]*bv[0]+wl[1]*bv[1]

assert float(af.b)<0
assert float(bf.a)>0

# Full-X7 transverse signs by direct interval categorical covariance.
Jv,Kv,gv=XI
cos3=[1,0,-1,0,1,0,-1,0,1,0,-1,0]
alt=[1,-1,1,-1,1,-1,1,-1,1,-1,1,-1]
hh=[Jv*cos3[j]+Kv*alt[j] for j in range(12)]
ww=[iv.exp(z) for z in hh]
Z=sum(ww)
p=[z/Z for z in ww]

def trig(k,j,kind):
    val=(mp.cos(2*mp.pi*k*j/12) if kind=="c"
         else mp.sin(2*mp.pi*k*j/12))
    return iv.mpf([str(val),str(val)])

def mean(v):
    return sum(p[j]*v[j] for j in range(12))

def cov(v,w):
    mv=mean(v);mw=mean(w)
    return mean([(v[j]-mv)*(w[j]-mw) for j in range(12)])

s4=iv.sqrt(L4/6);s5=iv.sqrt(L5/6);s3=iv.sqrt(L3/6)
f4=[s4*trig(4,j,"c") for j in range(12)]
f5=[s5*trig(5,j,"c") for j in range(12)]
f3s=[s3*trig(3,j,"s") for j in range(12)]

H44=1/gv-cov(f4,f4)
H55=1/gv-cov(f5,f5)
H45=-cov(f4,f5)
tr45=H44+H55
det45=H44*H55-H45*H45
H3s=1/gv-cov(f3s,f3s)

assert float(tr45.a)>0
assert float(det45.b)<0
assert float(H3s.a)>0

print("PASS")
print("krawczyk_image",
      [[float(z.a),float(z.b)] for z in Kr])
print("fold_a_interval",[float(af.a),float(af.b)])
print("fold_b_interval",[float(bf.a),float(bf.b)])
print("H45_trace_interval",[float(tr45.a),float(tr45.b)])
print("H45_det_interval",[float(det45.a),float(det45.b)])
print("H3sin_interval",[float(H3s.a),float(H3s.b)])
