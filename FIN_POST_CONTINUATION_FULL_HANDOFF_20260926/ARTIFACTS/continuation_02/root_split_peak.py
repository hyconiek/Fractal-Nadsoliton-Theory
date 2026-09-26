import math, numpy as np
from scipy import fft
from scipy.special import logsumexp
alpha=0.8245
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
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def flat5(c):
 x=0
 for v in c:x=x*5+v
 return x
def norm_one(c):
 a=np.asarray(c,float);return float(a@G@a)
rows=[];z2=np.full(5**8,np.nan);z4half=None; beta16=alpha/2
with open('/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02/level16_08245.txt') as f:
 for line in f:
  h,s1,s2,s4=line.split();c=decode_hex(h);n=norm_one(c); ly1=float(s1)-beta16*n/32;rows.append((c,ly1))
  if s2!='nan':z2[flat5(c)]=float(s2)-alpha*n/32
  if c==(2,)*8 and s4!='nan':z4half=float(s4)-(2*alpha)*n/32
M=max(x[1] for x in rows);A=np.zeros(shape,float)
for c,ly in rows:A[c]=math.exp(ly-M)
F=fft.rfftn(A,workers=5);F*=F;C=fft.irfftn(F,s=shape,workers=5);del F,A
flat=np.arange(NGRID,dtype=np.int64);tmp=flat.copy();digits=[];total=np.zeros(NGRID,dtype=np.uint8)
for _ in range(8):
 d=(tmp%9).astype(np.uint8);digits.append(d);total+=d;tmp//=9
digits=digits[::-1];mask=(total==32);idx=flat[mask];coords=[d[idx] for d in digits];del flat,tmp,total,mask
normv=np.zeros(idx.size,float)
for i in range(8):
 ci=coords[i].astype(float);normv+=G[i,i]*ci*ci
 for j in range(i):normv+=2*G[i,j]*ci*coords[j]
cv=C.ravel()[idx];del C
lord=np.log(np.maximum(cv,1e-300))+2*M
all_even=np.ones(idx.size,dtype=bool)
for ci in coords:all_even &= ((ci&1)==0)
if all_even.any():
 hidx=np.zeros(all_even.sum(),dtype=np.int64)
 for ci in coords:hidx=hidx*5+(ci[all_even]//2)
 ly2=z2[hidx];lord[all_even]=np.logaddexp(lord[all_even],ly2)
beta32=alpha/4;logY32=beta32*normv/64-math.log(2)+lord
# complement position
comp9=np.zeros(idx.size,dtype=np.int64)
for ci in coords:comp9=comp9*9+(8-ci)
posc=np.searchsorted(idx,comp9)
logw=logY32+logY32[posc]
lognorm=logsumexp(logw);p=np.exp(logw-lognorm)
# root Ward split cost; symmetric root T=0, m=64: delta=(norm(a)+norm(b))/32=norm(a)/16
root_delta=normv/16.0
# composition imbalance relative to 4 each
imb=np.zeros(idx.size,float)
for ci in coords:imb+=(ci.astype(float)-4.0)**2
# effective number / entropy of ordered root splits
H=-float(np.sum(p*np.log(np.maximum(p,1e-300))));Neff=math.exp(H);ipr=float(np.sum(p*p))
print('ordered_states',idx.size,'entropy',H,'Neff',Neff,'IPR',ipr,'invIPR',1/ipr)
for x,name in [(root_delta,'root_delta'),(imb,'imbalance')]:
 mean=float(np.sum(p*x));var=float(np.sum(p*(x-mean)**2));
 qs=[];order=np.argsort(x);cs=np.cumsum(p[order])
 for q in [.05,.25,.5,.75,.95]:qs.append(float(x[order[np.searchsorted(cs,q)]]))
 print(name,'mean',mean,'sd',math.sqrt(var),'q',qs)
# coarse histogram root_delta with probability bins
lo,hi=np.quantile(root_delta,[0,1]); bins=np.linspace(root_delta.min(),root_delta.max(),81); hist=np.zeros(len(bins)-1)
b=np.minimum(np.searchsorted(bins,root_delta,side='right')-1,len(hist)-1); ok=(b>=0)&(b<len(hist)); np.add.at(hist,b[ok],p[ok])
peaks=[]
for i in range(1,len(hist)-1):
 if hist[i]>=hist[i-1] and hist[i]>=hist[i+1]:peaks.append((hist[i],(bins[i]+bins[i+1])/2))
print('hist_peaks_top',sorted(peaks,reverse=True)[:8])
# top ordered root split weights; divide intuition only, complements appear twice
for ii in np.argsort(p)[-12:][::-1]:
 c=tuple(int(x[ii]) for x in coords)
 print('top',p[ii],c,'delta',root_delta[ii],'imb',imb[ii])
