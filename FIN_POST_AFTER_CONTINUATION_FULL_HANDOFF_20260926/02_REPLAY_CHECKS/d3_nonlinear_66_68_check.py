#!/usr/bin/env python3
import math, numpy as np
from scipy.special import softmax
from scipy.optimize import root
from mpmath import iv, mp

mp.dps=60
iv.dps=50
N=12

LINT={
3:(1.96140686197643,1.96140686197645),
4:(2.1995688493332,2.19956884933322),
5:(2.29860627207908,2.2986062720791),
6:(2.34218204114629,2.34218204114631),
}
lm={k:.5*(a+b) for k,(a,b) in LINT.items()}

Jbox=(0.038434553978898216,0.03843495397889821)
Kbox=(0.18527186140188867,0.18527226140188868)
gbox=(5.171841631942819,5.171842031942818)

J=iv.mpf([str(Jbox[0]),str(Jbox[1])])
K=iv.mpf([str(Kbox[0]),str(Kbox[1])])
gg=iv.mpf([str(gbox[0]),str(gbox[1])])
l4=iv.mpf([str(LINT[4][0]),str(LINT[4][1])])
l5=iv.mpf([str(LINT[5][0]),str(LINT[5][1])])

cos3=[1,0,-1,0,1,0,-1,0,1,0,-1,0]
alt=[1,-1,1,-1,1,-1,1,-1,1,-1,1,-1]
hh=[J*cos3[j]+K*alt[j] for j in range(N)]
ww=[iv.exp(x) for x in hh]
ZZ=sum(ww)
pp=[w/ZZ for w in ww]

def trig(k,j,kind):
    val=(mp.cos(2*mp.pi*k*j/12) if kind=="c"
         else mp.sin(2*mp.pi*k*j/12))
    return iv.mpf([str(val),str(val)])

s4=iv.sqrt(l4/6)
s5=iv.sqrt(l5/6)
f4=[s4*trig(4,j,"c") for j in range(N)]
f5=[s5*trig(5,j,"c") for j in range(N)]

def mean(vals):
    return sum(pp[j]*vals[j] for j in range(N))

def cov(f,h):
    mf=mean(f);mh=mean(h)
    return mean([(f[j]-mf)*(h[j]-mh) for j in range(N)])

A=1/gg-cov(f4,f4)
D=1/gg-cov(f5,f5)
C=-cov(f4,f5)

u=-C/A
nn=iv.sqrt(u*u+1)
aa=u/nn
bb=1/nn

Y=[aa*f4[j]+bb*f5[j] for j in range(N)]
my=mean(Y)
t=-mean([(Y[j]-my)**3 for j in range(N)])

tlo=float(t.a);thi=float(t.b)
assert thi<0
assert tlo < -0.01544 < thi

ddlo=0.0016957351021672582
ddhi=0.0016960387698923191
tr=A+D
trlo=float(tr.a);trhi=float(tr.b)
alpha_lo=ddlo/trhi
alpha_hi=ddhi/trlo
assert alpha_lo>0

l3=iv.mpf([str(LINT[3][0]),str(LINT[3][1])])
l6=iv.mpf([str(LINT[6][0]),str(LINT[6][1])])
s3=iv.sqrt(l3/6);s6=iv.sqrt(l6/12)
f3c=[s3*trig(3,j,"c") for j in range(N)]
f3s=[s3*trig(3,j,"s") for j in range(N)]
f6=[s6*alt[j] for j in range(N)]

H3s=1/gg-cov(f3s,f3s)
H33=1/gg-cov(f3c,f3c)
H66=1/gg-cov(f6,f6)
H36=-cov(f3c,f6)
det36=H33*H66-H36*H36

assert float(H3s.a)>0
assert float(det36.b)<0
assert float(tr.a)>0

tabs_min=abs(thi)
tabs_max=abs(tlo)
Cr_lo=2*alpha_lo/tabs_max
Cr_hi=2*alpha_hi/tabs_min
assert 18.4<Cr_lo<Cr_hi<18.7

j=np.arange(N)
cols=[]
for k in (3,4,5):
    cols += [
        np.sqrt(lm[k]/6)*np.cos(2*np.pi*k*j/N),
        np.sqrt(lm[k]/6)*np.sin(2*np.pi*k*j/N)
    ]
cols += [np.sqrt(lm[6]/12)*(-1.)**j]
X=np.column_stack(cols)

def xy(J,K):
    q=np.exp(-2*K);C=np.cosh(J);S=np.sinh(J);D=C+q
    return S/D,(C-q)/D

def crossfun(z):
    J,K,g=z;x,y=xy(J,K)
    det=(1/g-lm[4]/12)*(1/g-lm[5]/12)-lm[4]*lm[5]*x*x/144
    return np.array([J-g*lm[3]*x/6,K-g*lm[6]*y/12,det])

J0,K0,g0=root(crossfun,[.03843475,.18527206,5.17184183],tol=1e-13).x

def base(g):
    return root(lambda z:[
        z[0]-g*lm[3]*xy(*z)[0]/6,
        z[1]-g*lm[6]*xy(*z)[1]/12], [J0,K0],tol=1e-13).x

def grad(th,g):
    p=softmax(X@th)
    return th/g-X.T@p

th0=np.zeros(7)
th0[0]=J0/np.sqrt(lm[3]/6)
th0[6]=K0/np.sqrt(lm[6]/12)
p0=softmax(X@th0)
mu=p0@X
H=np.eye(7)/g0-(X.T@(p0[:,None]*X)-np.outer(mu,mu))
evc,Uc=np.linalg.eigh(H[np.ix_([2,4],[2,4])])
vec=Uc[:,0]
if vec[0]<0: vec=-vec

e1=np.zeros(7);e1[2]=vec[0];e1[4]=vec[1]
e2=np.zeros(7);e2[3]=vec[0];e2[5]=-vec[1]

amp_ratios=[]
indices=[]
for delta in (-1e-5,1e-5):
    g1=g0+delta
    J1,K1=base(g1)
    tb=np.zeros(7)
    tb[0]=J1/np.sqrt(lm[3]/6)
    tb[6]=K1/np.sqrt(lm[6]/12)
    pb=softmax(X@tb)
    mub=pb@X
    Hb=np.eye(7)/g1-(X.T@(pb[:,None]*X)-np.outer(mub,mub))
    lamsoft=np.linalg.eigvalsh(Hb[np.ix_([2,4],[2,4])])[0]
    angle=(0 if delta>0 else np.pi/3)
    edir=np.cos(angle)*e1+np.sin(angle)*e2
    tmid=.5*(tlo+thi)
    c3=(1 if delta>0 else -1)
    r0=-2*lamsoft/(tmid*c3)
    zz=root(lambda th:grad(th,g1),tb+r0*edir,tol=1e-12)
    assert np.linalg.norm(grad(zz.x,g1))<1e-10
    d=zz.x-tb
    rr=np.sqrt((d@e1)**2+(d@e2)**2)
    amp_ratios.append(rr/abs(delta))
    pz=softmax(X@zz.x);muz=pz@X
    Hz=np.eye(7)/g1-(X.T@(pz[:,None]*X)-np.outer(muz,muz))
    indices.append(int(np.sum(np.linalg.eigvalsh(Hz)<-1e-8)))

assert all(17.5<x<19.5 for x in amp_ratios)
assert indices==[2,2]

print("PASS")
print("cubic_t_interval",[tlo,thi])
print("soft_eigenvalue_slope_interval",[alpha_lo,alpha_hi])
print("amplitude_coefficient_interval",[Cr_lo,Cr_hi])
print("H3sin_lower",float(H3s.a))
print("H36_det_interval",[float(det36.a),float(det36.b)])
print("daughter_amp_ratio_replay",amp_ratios)
print("daughter_full_indices",indices)
print("parent_orbit_size",4)
print("daughter_orbit_size",12)
