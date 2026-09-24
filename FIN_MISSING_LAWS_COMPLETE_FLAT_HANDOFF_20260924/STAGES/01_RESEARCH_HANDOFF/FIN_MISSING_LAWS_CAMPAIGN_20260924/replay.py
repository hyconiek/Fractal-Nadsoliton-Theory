#!/usr/bin/env python3
from pathlib import Path
import json, numpy as np, math
root=Path(__file__).resolve().parent
# Verify the numerical headline fields without rewriting artifacts.
w=np.array([0.4699856726450201,0.1920435516901028,0.09142861427792495,0.0470291687456504,0.02413122336363006,0.011070817321442113])
def A(q):
 x=np.zeros((q,q))
 for i in range(q):
  for j in range(q):
   if i!=j:
    d=min((j-i)%q,(i-j)%q)
    if 1<=d<=6:x[i,j]=-w[d-1]
  x[i,i]=-x[i].sum()
 return x
def C(q):
 t=2*np.pi*np.arange(q)/q
 return np.vstack([np.cos(t),np.sin(t),np.cos(2*t),np.sin(2*t)])
def E(a):
 c=C(len(a));return a-a@c.T@np.linalg.inv(c@a@c.T)@c@a
assert sum(np.linalg.eigvalsh(E(A(12)))>1e-10)==7
assert sum(np.linalg.eigvalsh(E(A(18)))>1e-10)==13
M=np.array([[-1,1,0],[1,-2,1],[4,-3,0]])
assert round(np.linalg.det(M))==1
assert np.all(M@np.array([3,4,5])==np.array([1,0,0]))
print(json.dumps({'status':'PASS','rank12':7,'rank18':13,'phase_basis_det':1}))
