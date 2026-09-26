import math
import numpy as np
from scipy.special import logsumexp
alpha=.7459; Q=8
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]])
G=V@V.T
def dec(h):
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def norm(c):
 a=np.array(c,float);return float(a@G@a)
D={}
with open('/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02/level16_07459.txt') as f:
 for line in f:
  h,z1,z2,z4=line.split();c=dec(h)
  if max(c)<=4:
   D[c]=float(z1)-(alpha/2)*norm(c)/32
root=(4,)*8
cs=[];lw=[];q=[];delta=[]
for c,y in D.items():
 if sum(c)!=16: continue
 b=tuple(4-x for x in c)
 if b not in D:continue
 cs.append(c);lw.append(y+D[b]);q.append(sum((x-2)**2 for x in c));delta.append(norm(c)/8)
lw=np.array(lw);p=np.exp(lw-logsumexp(lw));q=np.array(q);delta=np.array(delta)
H=-float(np.sum(p*np.log(p)))
print('states',len(cs),'q_mean',float(p@q),'q_sd',float(np.sqrt(p@((q-p@q)**2))))
print('entropy',H,'Neff',math.exp(H),'maxp',p.max(),'delta_mean',float(p@delta))
h=np.bincount(q.astype(int),weights=p,minlength=int(max(q))+1)
print('hist_nonzero',[(i,float(v)) for i,v in enumerate(h) if v>1e-8])
mx=max(q); print('qmax',mx,'tail_ge_3/4max',float(p[q>=.75*mx].sum()),'low_le_1/4max',float(p[q<=.25*mx].sum()))
ix=np.argsort(p)[-10:][::-1]
for i in ix:print(cs[i],p[i],q[i],delta[i])
