import math,sys,numpy as np
from scipy import fft
from scipy.special import logsumexp
BASE='/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02'
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],
[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],
[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],
[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],
[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],
[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],
[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],
[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]],float)
G=V@V.T
shape=(9,)*8; NGRID=9**8

def decode_hex(h):
 x=int(h,16); return tuple((x>>(4*i))&15 for i in range(8))
def flat5(c):
 x=0
 for v in c:x=x*5+v
 return x
def norm_one(c):
 a=np.asarray(c,float);return float(a@G@a)
def norm_arr(coords):
 out=np.zeros(len(coords[0]))
 for i in range(8):
  ci=coords[i].astype(float);out+=G[i,i]*ci*ci
  for j in range(i):out+=2*G[i,j]*ci*coords[j]
 return out
# static indices total32 once
flat=np.arange(NGRID,dtype=np.int64); tmp=flat.copy();digits=[];total=np.zeros(NGRID,dtype=np.uint8)
for _ in range(8):
 d=(tmp%9).astype(np.uint8);digits.append(d);total+=d;tmp//=9
digits=digits[::-1];mask=(total==32);idx=flat[mask];coords=[d[idx] for d in digits]
comp9=np.zeros(idx.size,dtype=np.int64)
for ci in coords: comp9=comp9*9+(8-ci)
posc=np.searchsorted(idx,comp9);assert np.all(idx[posc]==comp9)
normv=norm_arr(coords)
all_even=np.ones(idx.size,dtype=bool)
for ci in coords:all_even &= ((ci&1)==0)
hidx=np.zeros(all_even.sum(),dtype=np.int64)
for ci in coords:hidx=hidx*5+(ci[all_even]//2)
# base5 valid/complements static

def branch_logs(alpha,tag):
 path=f'{BASE}/level16_{tag}.txt'
 rows=[];z2=np.full(5**8,np.nan);z4half=None; beta16=alpha/2
 with open(path) as f:
  for line in f:
   h,s1,s2,s4=line.split();c=decode_hex(h);n=norm_one(c);rows.append((c,float(s1)-beta16*n/32))
   if s2!='nan':z2[flat5(c)]=float(s2)-alpha*n/32
   if c==(2,)*8 and s4!='nan':z4half=float(s4)-2*alpha*n/32
 M=max(v for _,v in rows);A=np.zeros(shape)
 for c,v in rows:A[c]=math.exp(v-M)
 F=fft.rfftn(A,workers=5);F*=F;C=fft.irfftn(F,s=shape,workers=5)
 cv=C.ravel()[idx]; lord=np.log(cv)+2*M
 lord[all_even]=np.logaddexp(lord[all_even],z2[hidx])
 beta32=alpha/4; logY32=beta32*normv/64-math.log(2)+lord
 ordered=logY32+logY32[posc] # common root -log2 omitted
 # equal correction branch at root
 valid=np.flatnonzero(np.isfinite(z2));vv=valid.copy();digs5=[]
 for _ in range(8):digs5.append(vv%5);vv//=5
 digs5=digs5[::-1];comp=np.zeros(valid.size,dtype=np.int64)
 for d in digs5:comp=comp*5+(4-d)
 pos2=np.searchsorted(valid,comp);assert np.all(valid[pos2]==comp)
 lord2=float(logsumexp(z2[valid]+z2[comp]));n4=norm_one((4,)*8); beta32q2=alpha/2
 ly32q2=beta32q2*n4/64-math.log(2)+float(np.logaddexp(lord2,z4half))
 return ordered,ly32q2
alphas=[0.8235,0.8240,0.8245,0.8250,0.8255]
tags=['08235','08240','08245','0825','08255']
Ls=[]; Es=[]
for a,t in zip(alphas,tags):
 print('build',a,flush=True);l,e=branch_logs(a,t);Ls.append(l);Es.append(e)
L=np.stack(Ls,axis=0);E=np.array(Es)
h=0.0005;a0=0.8245
# 5point derivatives for log branch
lp=(L[0]-8*L[1]+8*L[3]-L[4])/(12*h)
lpp=(-L[4]+16*L[3]-30*L[2]+16*L[1]-L[0])/(12*h*h)
ep=(E[0]-8*E[1]+8*E[3]-E[4])/(12*h)
epp=(-E[4]+16*E[3]-30*E[2]+16*E[1]-E[0])/(12*h*h)
# root branch probs, including equal correction; common 1/2 cancels
lognorm=np.logaddexp(logsumexp(L[2]),E[2]);p=np.exp(L[2]-lognorm);pe=math.exp(E[2]-lognorm)
meanp=float(p@lp + pe*ep)
between=float(p@((lp-meanp)**2)+pe*(ep-meanp)**2)
within=float(p@lpp+pe*epp)
total=between+within
print('pequal',pe)
print('d2_between',between,'d2_within',within,'d2_total',total)
print('CH_between',a0*a0*between,'CH_within',a0*a0*within,'CH_total',a0*a0*total)
print('fractions',between/total,within/total)
# independent total finite diff logZ from branch normalization values
Z=np.array([np.logaddexp(logsumexp(L[k]),E[k])-math.log(2) for k in range(5)])
zpp=(-Z[4]+16*Z[3]-30*Z[2]+16*Z[1]-Z[0])/(12*h*h)
print('CH_direct5',a0*a0*zpp,'resid',a0*a0*(zpp-total))
