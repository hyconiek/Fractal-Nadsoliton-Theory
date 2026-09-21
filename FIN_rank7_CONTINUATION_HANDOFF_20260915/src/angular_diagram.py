"""R7P-096 locked angular branch landmarks and local coexistence certificate."""
from __future__ import annotations
import json,math
from pathlib import Path
import numpy as np, mpmath as mp, sympy as sp
from scipy.optimize import root
from scipy.special import softmax,logsumexp
from .model import feature_spaces

PH=np.array([-math.pi/2,2*math.pi/3,-math.pi/6]); SIGN=-1

def basis_and_Y():
 _,_,L,X,_,_=feature_spaces();B=np.zeros((7,4))
 for q,p in enumerate(PH):B[2*q,q]=math.cos(p);B[2*q+1,q]=-math.sin(p)
 B[6,3]=SIGN;return L,B,X@B

def full_stats(a,Y):
 h=Y@a;p=softmax(h);mu=p@Y;D=Y-mu;H=D.T@(p[:,None]*D);K=logsumexp(h)-math.log(12);return K,mu,H

def full_F5(x,r,Y):
 a=x[:4];lam=x[4];K,g,H=full_stats(a,Y);return np.r_[g-lam*a,a@a-r*r]
def full_J5(x,r,Y):
 a=x[:4];lam=x[4];K,g,H=full_stats(a,Y);J=np.zeros((5,5));J[:4,:4]=H-lam*np.eye(4);J[:4,4]=-a;J[4,:4]=2*a;return J

def solve_full_coex(Y,L):
 z=np.array([.1131879146,.1698528641,.2269339093]);a=np.array([math.sqrt(2/L[k])*zz for k,zz in zip([3,4,5],z)]+[.3380663037/math.sqrt(L[6])])
 def F(v):
  aa=v[:4];lam=v[4];r=v[5];K,g,H=full_stats(aa,Y);Kp=logsumexp(Y[:,3]*r)-math.log(12)
  return np.r_[g-lam*aa,aa@aa-r*r,K-Kp]
 sol=root(F,np.r_[a,.1975,np.linalg.norm(a)],tol=1e-12);return sol.x,float(np.linalg.norm(F(sol.x)))

def full_fold(Y,coex):
 x=coex[:5].copy()
 for r in np.linspace(coex[5],.3465,80):
  x=root(lambda z:full_F5(z,r,Y),x,jac=lambda z:full_J5(z,r,Y),tol=1e-12).x
 _,_,vt=np.linalg.svd(full_J5(x,r,Y));v=vt[-1]
 def A(w):
  xx=w[:5];rr=w[5];vv=w[6:];return np.r_[full_F5(xx,rr,Y),full_J5(xx,rr,Y)@vv,vv@vv-1]
 sol=root(A,np.r_[x,r,v],tol=1e-11);return sol.x,float(np.linalg.norm(A(sol.x)))

def full_saddle_at_coex(Y,r):
 rng=np.random.default_rng(96096);sols=[]
 def F(x):return full_F5(x,r,Y)
 def J(x):return full_J5(x,r,Y)
 seeds=[np.r_[[0,0,0,r],.19]]
 for _ in range(300):
  a=np.abs(rng.normal(size=4));a=a/np.linalg.norm(a)*r;seeds.append(np.r_[a,.195])
 for x0 in seeds:
  s=root(F,x0,jac=J,tol=1e-11);x=s.x
  if np.linalg.norm(F(x))>1e-8 or np.min(x[:4])<-1e-7:continue
  if not any(np.linalg.norm(x[:4]-q[:4])<1e-7 for q in sols):sols.append(x)
 vals=[]
 for x in sols:
  K=full_stats(x[:4],Y)[0];vals.append((K,x))
 vals.sort(key=lambda q:q[0],reverse=True)
 # Two equal maxima are first; choose next lower nontrivial state as barrier saddle.
 kmax=vals[0][0]; candidates=[q for q in vals if kmax-q[0]>1e-9]
 return kmax,candidates[0][0],candidates[0][1],len(sols)

def quartic_functions(L):
 a=sp.symbols('a0:4',real=True); rmods=[sp.Float(math.sqrt(L[k]/2),30)*a[i] for i,k in enumerate([3,4,5])];z6=-sp.Float(math.sqrt(L[6]),30)*a[3]
 hs=[]
 for jj in range(12):
  v=0
  for kk,(mode,p) in enumerate(zip([3,4,5],PH)):v+=rmods[kk]/sp.sqrt(3)*sp.cos(2*sp.pi*mode*jj/12+sp.Float(float(p),30))
  v+=z6/sp.sqrt(12)*((-1)**jj);hs.append(v)
 m2=sp.expand(sum(x*x for x in hs)/12);m3=sp.expand(sum(x**3 for x in hs)/12);m4=sp.expand(sum(x**4 for x in hs)/12);K=sp.expand(m2/2+m3/6+(m4-3*m2*m2)/24)
 g=[sp.diff(K,x) for x in a];H=sp.hessian(K,a)
 return sp.lambdify(a,K,'numpy'),sp.lambdify(a,g,'numpy'),sp.lambdify(a,H,'numpy')

def quartic_landmarks(L):
 K,g,H=quartic_functions(L)
 def stats(a):return float(K(*a)),np.array(g(*a),float),np.array(H(*a),float)
 def coF(v):
  a=v[:4];lam=v[4];r=v[5];kv,gg,hh=stats(a);kp=K(0,0,0,r);return np.r_[gg-lam*a,a@a-r*r,kv-kp]
 seed=np.array([.1143,.1614,.2110,.2215,.1975,.36428]);co=root(coF,seed,tol=1e-12).x
 def F5(x,r):
  a=x[:4];lam=x[4];kv,gg,hh=stats(a);return np.r_[gg-lam*a,a@a-r*r]
 def J5(x,r):
  a=x[:4];lam=x[4];kv,gg,hh=stats(a);J=np.zeros((5,5));J[:4,:4]=hh-lam*np.eye(4);J[:4,4]=-a;J[4,:4]=2*a;return J
 x=co[:5].copy()
 for rr in np.linspace(co[5],.3468,80):x=root(lambda z:F5(z,rr),x,jac=lambda z:J5(z,rr),tol=1e-12).x
 _,_,vt=np.linalg.svd(J5(x,rr));v=vt[-1]
 def aug(w):return np.r_[F5(w[:5],w[5]),J5(w[:5],w[5])@w[6:],w[6:]@w[6:]-1]
 fold=root(aug,np.r_[x,rr,v],tol=1e-11).x
 return {'coexistence':co.tolist(),'fold':fold[:6].tolist()}

def certify_full_coex(Y,center,radius=1e-8):
 mp.iv.dps=50; n=6; C=[mp.iv.mpf(repr(float(x))) for x in center];Xv=[mp.iv.mpf([repr(float(x-radius)),repr(float(x+radius))]) for x in center]
 Ym=[[mp.iv.mpf(repr(float(Y[j,i]))) for i in range(4)] for j in range(12)]
 def FJ(v):
  a=v[:4];lam=v[4];rr=v[5];ew=[]
  for j in range(12): ew.append(mp.iv.exp(sum(Ym[j][i]*a[i] for i in range(4))))
  Z=sum(ew,mp.iv.mpf('0'));N=[sum(ew[j]*Ym[j][i] for j in range(12)) for i in range(4)];S=[[sum(ew[k]*Ym[k][i]*Ym[k][j] for k in range(12)) for j in range(4)] for i in range(4)]
  ep=[mp.iv.exp(Ym[j][3]*rr) for j in range(12)];Zp=sum(ep,mp.iv.mpf('0'));Np=sum(ep[j]*Ym[j][3] for j in range(12))
  F=[N[i]-lam*a[i]*Z for i in range(4)]+[sum(x*x for x in a)-rr*rr,Z-Zp]
  J=[[mp.iv.mpf('0') for _ in range(6)] for __ in range(6)]
  for i in range(4):
   for j in range(4):J[i][j]=S[i][j]-lam*((Z if i==j else 0)+a[i]*N[j])
   J[i][4]=-a[i]*Z
  for j in range(4):J[4][j]=2*a[j];J[4][5]=-2*rr;J[5][j]=N[j]
  J[5][5]=-Np
  return F,J
 Fc,Jc=FJ(C);_,JX=FJ(Xv)
 # midpoint Jacobian from tiny intervals
 J0=np.array([[(float(x.a)+float(x.b))/2 for x in row] for row in Jc]);A=np.linalg.inv(J0);D=mp.iv.mpf([repr(-radius),repr(radius)]);one=mp.iv.mpf('1');zero=mp.iv.mpf('0');mx=0.;boxes=[]
 for i in range(6):
  ss=zero
  for j in range(6):ss-=mp.iv.mpf(repr(float(A[i,j])))*Fc[j]
  for j in range(6):
   Bij=one if i==j else zero
   for k in range(6):Bij-=mp.iv.mpf(repr(float(A[i,k])))*JX[k][j]
   ss+=Bij*D
  lo,hi=float(ss.a),float(ss.b);mx=max(mx,abs(lo),abs(hi));boxes.append([lo,hi])
 return {'radius':radius,'strict_inclusion':mx<radius,'max_abs_krawczyk_offset':mx,'offset_intervals':boxes,'scope':'fixed nominal strict spectral tuple'}

def run(out):
 L,B,Y=basis_and_Y();co,res=solve_full_coex(Y,L);fold,fres=full_fold(Y,co);kmax,ksad,sad,nstat=full_saddle_at_coex(Y,co[5]);q=quartic_landmarks(L);cert=certify_full_coex(Y,co)
 data={'full':{'coexistence':co.tolist(),'coexistence_residual':res,'fold':fold[:6].tolist(),'fold_augmented_residual':fres,'coexistence_stationary_count_in_positive_locked_chart':nstat,'angular_saddle':sad.tolist(),'logmgf_barrier':kmax-ksad},'quartic':q,'full_coexistence_certificate':cert,'scope':'locked reflection-fixed angular branches on fixed theta-radius sphere; not radial g transition or global rank7 theorem'}
 Path(out).write_text(json.dumps(data,indent=2)+'\n');return data
if __name__=='__main__':
 import sys;o=run(sys.argv[1]);print(json.dumps(o,indent=2))
