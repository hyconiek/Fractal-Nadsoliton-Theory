# Dump M64 root ordered branch log weights at one alpha.
import math,sys,numpy as np
from scipy import fft
from scipy.special import logsumexp
alpha=float(sys.argv[1]); tag=sys.argv[2]; out=sys.argv[3]
BASE='/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02'
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]],float)
G=V@V.T;shape=(9,)*8;NGRID=9**8
def decode_hex(h):
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def flat5(c):
 x=0
 for v in c:x=x*5+v
 return x
def norm_one(c):
 a=np.asarray(c,float);return float(a@G@a)
rows=[];z2=np.full(5**8,np.nan);z4half=None; beta16=alpha/2
with open(f'{BASE}/level16_{tag}.txt') as f:
 for line in f:
  h,s1,s2,s4=line.split();c=decode_hex(h);n=norm_one(c);rows.append((c,float(s1)-beta16*n/32))
  if s2!='nan':z2[flat5(c)]=float(s2)-alpha*n/32
  if c==(2,)*8 and s4!='nan':z4half=float(s4)-2*alpha*n/32
M=max(v for _,v in rows);A=np.zeros(shape)
for c,v in rows:A[c]=math.exp(v-M)
F=fft.rfftn(A,workers=5);F*=F;C=fft.irfftn(F,s=shape,workers=5);del A,F
flat=np.arange(NGRID,dtype=np.int64);tmp=flat.copy();digits=[];total=np.zeros(NGRID,dtype=np.uint8)
for _ in range(8):
 d=(tmp%9).astype(np.uint8);digits.append(d);total+=d;tmp//=9
digits=digits[::-1];mask=(total==32);idx=flat[mask];coords=[d[idx] for d in digits];del flat,tmp,total,mask
# norm
normv=np.zeros(idx.size,float)
for i in range(8):
 ci=coords[i].astype(float);normv+=G[i,i]*ci*ci
 for j in range(i):normv+=2*G[i,j]*ci*coords[j]
cv=C.ravel()[idx];del C;lord=np.log(cv)+2*M
all_even=np.ones(idx.size,dtype=bool)
for ci in coords:all_even&=((ci&1)==0)
hidx=np.zeros(all_even.sum(),dtype=np.int64)
for ci in coords:hidx=hidx*5+(ci[all_even]//2)
lord[all_even]=np.logaddexp(lord[all_even],z2[hidx])
beta32=alpha/4;logY32=beta32*normv/64-math.log(2)+lord
comp9=np.zeros(idx.size,dtype=np.int64)
for ci in coords:comp9=comp9*9+(8-ci)
pos=np.searchsorted(idx,comp9);assert np.all(idx[pos]==comp9)
logw=logY32+logY32[pos]
# equal correction
valid=np.flatnonzero(np.isfinite(z2));vv=valid.copy();digs5=[]
for _ in range(8):digs5.append(vv%5);vv//=5
digs5=digs5[::-1];comp=np.zeros(valid.size,dtype=np.int64)
for d in digs5:comp=comp*5+(4-d)
pos2=np.searchsorted(valid,comp);lord2=float(logsumexp(z2[valid]+z2[pos2]));n4=norm_one((4,)*8);beta32q2=alpha/2
leq=beta32q2*n4/64-math.log(2)+float(np.logaddexp(lord2,z4half))
np.save(out,logw)
open(out+'.eq','w').write(repr(leq)+'\n')
print(alpha,logw.shape,leq,float(logsumexp(logw)))
