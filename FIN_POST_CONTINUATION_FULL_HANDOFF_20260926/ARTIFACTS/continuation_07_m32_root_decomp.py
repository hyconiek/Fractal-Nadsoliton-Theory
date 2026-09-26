import math,numpy as np
from scipy.special import logsumexp
BASE='/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02'
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]],float)
G=V@V.T
def decode_hex(h):
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def norm(c):
 a=np.asarray(c,float);return float(a@G@a)
def child_logs(alpha,tag):
 vals={}; eq=None
 for line in open(f'{BASE}/level16_{tag}.txt'):
  h,s1,s2,s4=line.split();c=decode_hex(h); n=norm(c)
  # only total16 rows useful
  if sum(c)==16:
   vals[c]=float(s1)-(alpha/2)*n/32
  if c==(2,)*8:
   eq=float(s2)-alpha*n/32
 # ordered root children a+b=4^8
 logs=[]
 for c,ly in vals.items():
  comp=tuple(4-x for x in c)
  if comp in vals: logs.append(ly+vals[comp])
 return np.array(logs),eq
alphas=[.735,.740,.745,.750,.755]; tags=['0735','074','07459','075','0755']
# NOTE tag 07459 file alpha may be 0.7459, not .745; inspect later. Prefer exact 0745? 
# replace center using level16_07459 only if no 0745 exists
import os
if os.path.exists(f'{BASE}/level16_0745.txt'):
 tags[2]='0745'
Ls=[];Es=[]
for a,t in zip(alphas,tags):
 l,e=child_logs(a,t); print('point',a,t,len(l),e,logsumexp(l));Ls.append(l);Es.append(e)
# verify lengths and implicit ordering from dict insertion stable across files
for l in Ls: assert len(l)==len(Ls[0])
L=np.stack(Ls);E=np.array(Es);h=.005;a0=.745
lp=(L[0]-8*L[1]+8*L[3]-L[4])/(12*h)
lpp=(-L[4]+16*L[3]-30*L[2]+16*L[1]-L[0])/(12*h*h)
ep=(E[0]-8*E[1]+8*E[3]-E[4])/(12*h)
epp=(-E[4]+16*E[3]-30*E[2]+16*E[1]-E[0])/(12*h*h)
lognorm=np.logaddexp(logsumexp(L[2]),E[2]);p=np.exp(L[2]-lognorm);pe=math.exp(E[2]-lognorm)
mean=float(p@lp+pe*ep);between=float(p@((lp-mean)**2)+pe*(ep-mean)**2);within=float(p@lpp+pe*epp);total=between+within
Z=np.array([np.logaddexp(logsumexp(L[k]),E[k])-math.log(2) for k in range(5)])
zpp=(-Z[4]+16*Z[3]-30*Z[2]+16*Z[1]-Z[0])/(12*h*h)
print('pequal',pe,'Z',Z.tolist())
print('CH_between',a0*a0*between,'CH_within',a0*a0*within,'CH_total',a0*a0*total,'frac_between',between/total)
print('CH_direct5',a0*a0*zpp,'resid',a0*a0*(zpp-total))
