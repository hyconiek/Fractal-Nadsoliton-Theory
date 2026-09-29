#!/usr/bin/env python3
import math, sys, numpy as np
Q=12; g=3.0
L=[0,0,0,1.9614068619764458,2.199568849333211,2.298606272079097,2.3421820411462977]
def first_row():
    out=np.zeros(Q)
    for k in range(1,7):
        out += (1 if k==6 else 2)*L[k]*np.cos(2*np.pi*k*np.arange(Q)/Q)/Q
    return out
def sm(x):
    x=x-x.max(); z=np.exp(x); return z/z.sum()
target=sm((g/2)*first_row())
rows=np.loadtxt(sys.argv[1],delimiter=",")
def near(curve,t):
    ii=np.where(curve>0)[0]
    z=ii[np.argmin(np.abs(np.log(curve[ii])-math.log(t)))]
    return curve[z],z,abs(math.log(curve[z])-math.log(t))
def assign(cost):
    dp={0:(0.0,())}
    for j in range(Q):
        nd={}
        for mask,(v,path) in dp.items():
            for i in range(cost.shape[1]):
                if mask>>i&1: continue
                nm=mask|(1<<i); nv=v+cost[j,i]
                if nm not in nd or nv<nd[nm][0]: nd[nm]=(nv,path+(i,))
        dp=nd
    return min(dp.values(),key=lambda z:z[0])[1]
def ev(c):
    nr=np.zeros((Q,len(rows))); code=np.zeros((Q,len(rows)),int); cost=np.zeros_like(nr)
    for j in range(Q):
        for i in range(len(rows)):
            r,k,d=near(rows[i],c*target[j]); nr[j,i]=r; code[j,i]=k; cost[j,i]=d*d
    ch=assign(cost); rates=np.array([nr[j,ch[j]] for j in range(Q)])
    q=rates/rates.sum(); tv=.5*np.abs(q-target).sum()
    return tv,q,ch,[int(code[j,ch[j]]) for j in range(Q)]
best=None
for c in np.geomspace(.03,50,80):
    x=ev(c)
    if best is None or x[0]<best[0]: best=(x[0],c,*x[1:])
for c in best[1]*np.exp(np.linspace(-.25,.25,101)):
    x=ev(c)
    if x[0]<best[0]: best=(x[0],c,*x[1:])
print("TV",best[0]); print("scale",best[1]); print("programmed",best[2].tolist())
print("channels",list(best[3])); print("codes",best[4])
