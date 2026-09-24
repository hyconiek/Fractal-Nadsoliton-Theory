#!/usr/bin/env python3
import numpy as np, json
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
def cond(a):
 c=C(len(a));return a-a@c.T@np.linalg.inv(c@a@c.T)@c@a
x=A(12); e=cond(x); vals,U=np.linalg.eigh(x); ix=np.argsort(vals)[-7:]; top=(U[:,ix]*vals[ix])@U[:,ix].T
xp=x.copy(); d=.2*w[0]; xp[0,0]+=d;xp[1,1]+=d;xp[0,1]-=d;xp[1,0]-=d
ep=cond(xp); vp,Up=np.linalg.eigh(xp); ix=np.argsort(vp)[-7:]; tp=(Up[:,ix]*vp[ix])@Up[:,ix].T
out={'rank12':int(sum(np.linalg.eigvalsh(e)>1e-10)),'Aeff_top7_fro':float(np.linalg.norm(e-top)),'q_ranks':{str(q):int(sum(np.linalg.eigvalsh(cond(A(q)))>1e-10)) for q in [12,18,24,30]},'perturbed_relative_difference':float(np.linalg.norm(ep-tp)/np.linalg.norm(ep)),'perturbed_conditioned_constraint':float(np.linalg.norm(C(12)@ep)),'perturbed_top7_constraint':float(np.linalg.norm(C(12)@tp))}
print(json.dumps(out,indent=2,sort_keys=True))
