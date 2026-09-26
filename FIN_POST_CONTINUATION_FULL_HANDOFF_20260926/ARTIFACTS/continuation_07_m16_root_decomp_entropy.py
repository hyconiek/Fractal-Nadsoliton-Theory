import math,functools,numpy as np
from scipy.special import logsumexp
# strict 8-phase vectors from validated replay
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.44591493030279611 if False else -0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]],float)
s=2.0
@functools.lru_cache(None)
def splits(c):
 h=sum(c)//2;out=[]
 def rec(i,rem,a):
  if i==len(c):
   if rem==0:
    aa=tuple(a);bb=tuple(x-y for x,y in zip(c,aa))
    if aa<=bb:out.append((aa,bb))
   return
  tail=sum(c[i+1:])
  for x in range(max(0,rem-tail),min(c[i],rem)+1):rec(i+1,rem-x,a+[x])
 rec(0,h,[]);return tuple(out)
def delta(a,b):
 na=sum(a);nb=sum(b); ma=np.array(a)@V/na;mb=np.array(b)@V/nb
 return na*nb/(na+nb)*float((ma-mb)@(ma-mb))
def lse(xs): return float(logsumexp(np.array(xs,float)))
def make_logZ():
 @functools.lru_cache(None)
 def logZ(c,beta):
  if sum(c)==1:return 0.
  xs=[]
  for a,b in splits(c):
   d=delta(a,b);la=logZ(a,beta*s)
   if a!=b: xs.append(-beta*d+la+logZ(b,beta*s))
   else: xs.append(-beta*d+np.logaddexp(2*la,logZ(a,2*beta*s))-math.log(2))
  return lse(xs)
 return logZ
def root_branches(alpha):
 logZ=make_logZ(); root=(2,)*8; beta=alpha/2
 out=[]; ids=[]
 for a,b in splits(root):
  d=delta(a,b); la=logZ(a,2*beta)
  if a!=b:
   out.append(-beta*d+la+logZ(b,2*beta));ids.append(('ab',a,b))
  else:
   out.append(-beta*d+2*la-math.log(2));ids.append(('eqprod',a))
   out.append(-beta*d+logZ(a,4*beta)-math.log(2));ids.append(('eqdiag',a))
 return ids,np.array(out),logZ(root,beta)
alphas=[.724,.729,.734,.739,.744]; L=[]; ids0=None; Z=[]
for a in alphas:
 ids,l,z=root_branches(a);print('point',a,len(l),z);L.append(l);Z.append(z)
 if ids0 is None:ids0=ids
 else:assert ids==ids0
L=np.stack(L);Z=np.array(Z);h=.005;a0=.734
lp=(L[0]-8*L[1]+8*L[3]-L[4])/(12*h)
lpp=(-L[4]+16*L[3]-30*L[2]+16*L[1]-L[0])/(12*h*h)
ln=logsumexp(L[2]);p=np.exp(L[2]-ln);ent=float(-p@np.log(p+1e-300)); print('root_entropy',ent,'root_Neff',math.exp(ent),'pmax',float(p.max()))
mean=float(p@lp);between=float(p@((lp-mean)**2));within=float(p@lpp);total=between+within
zpp=(-Z[4]+16*Z[3]-30*Z[2]+16*Z[1]-Z[0])/(12*h*h)
print('CH_between',a0*a0*between,'CH_within',a0*a0*within,'CH_total',a0*a0*total,'frac_between',between/total)
print('CH_direct5',a0*a0*zpp,'resid',a0*a0*(zpp-total))
