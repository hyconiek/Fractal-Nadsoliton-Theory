# Reconstruct M32 homogeneous array exactly as fft_m64_vec_arg.py, then analyse M64 root split weights.
import math,sys,time,numpy as np
from scipy import fft
from scipy.special import logsumexp
alpha=float(sys.argv[1]) if len(sys.argv)>1 else 0.8245
lev=sys.argv[2] if len(sys.argv)>2 else f'/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_02/level16_{str(alpha).replace("0.","0").replace(".","")}.txt'
# hard use provided existing file path from argv recommended
shape=(9,)*8; NGRID=9**8
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

def decode_hex(h):
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def flat5(c):
 x=0
 for v in c:x=x*5+v
 return x
def norm_arr(coords):
 out=np.zeros(len(coords[0]))
 for i in range(8):
  ci=coords[i].astype(float);out+=G[i,i]*ci*ci
  for j in range(i):out+=2*G[i,j]*ci*coords[j]
 return out
def norm_one(c):
 a=np.asarray(c,float);return float(a@G@a)
rows=[];z2=np.full(5**8,np.nan);z4half=None; beta16=alpha/2
with open(lev) as f:
 for line in f:
  h,s1,s2,s4=line.split();c=decode_hex(h);n=norm_one(c);rows.append((c,float(s1)-beta16*n/32))
  if s2!='nan':z2[flat5(c)]=float(s2)-alpha*n/32
  if c==(2,)*8 and s4!='nan':z4half=float(s4)-2*alpha*n/32
M=max(v for _,v in rows);A=np.zeros(shape)
for c,v in rows:A[c]=math.exp(v-M)
F=fft.rfftn(A,workers=5);F*=F;C=fft.irfftn(F,s=shape,workers=5); del A,F
flat=np.arange(NGRID,dtype=np.int64);tmp=flat.copy();digits=[];total=np.zeros(NGRID,dtype=np.uint8)
for _ in range(8):
 d=(tmp%9).astype(np.uint8);digits.append(d);total+=d;tmp//=9
digits=digits[::-1];mask=(total==32);idx=flat[mask];coords=[d[idx] for d in digits];del flat,tmp,total,mask
normv=norm_arr(coords);cv=C.ravel()[idx];lord=np.log(cv)+2*M
all_even=np.ones(idx.size,dtype=bool)
for ci in coords:all_even&=((ci&1)==0)
if all_even.any():
 hidx=np.zeros(all_even.sum(),dtype=np.int64)
 for ci in coords:hidx=hidx*5+(ci[all_even]//2)
 lord[all_even]=np.logaddexp(lord[all_even],z2[hidx])
beta32=alpha/4;logY32=beta32*normv/64-math.log(2)+lord
# complement under root 8^8
comp9=np.zeros(idx.size,dtype=np.int64)
for ci in coords:comp9=comp9*9+(8-ci)
pos=np.searchsorted(idx,comp9);assert np.all(idx[pos]==comp9)
logw=logY32+logY32[pos]
lognorm=logsumexp(logw); p=np.exp(logw-lognorm)
# root Ward energy Delta=normv/16 since total vector is zero and children are 32/32
Delta=normv/16
mean=float(p@Delta);var=float(p@((Delta-mean)**2));sk=float(p@((Delta-mean)**3))/var**1.5;kurt=float(p@((Delta-mean)**4))/var**2
ent=float(-p@np.log(p+1e-300));neff=math.exp(ent)
print('alpha',alpha,'states',len(p),'rootDelta_mean',mean,'sd',math.sqrt(var),'skew',sk,'kurt',kurt,'entropy',ent,'Neff',neff)
# quantiles
ordr=np.argsort(Delta);cdf=np.cumsum(p[ordr])
for q in [.01,.05,.1,.25,.5,.75,.9,.95,.99]:
 k=np.searchsorted(cdf,q);print('q',q,'Delta',Delta[ordr[k]])
# weighted histogram 80 bins; find local maxima robust after simple smoothing
lo,hi=Delta.min(),Delta.max();hist,edges=np.histogram(Delta,bins=80,range=(lo,hi),weights=p);cent=(edges[:-1]+edges[1:])/2
sm=np.convolve(hist,[.25,.5,.25],mode='same');peaks=[]
for i in range(1,len(sm)-1):
 if sm[i]>sm[i-1] and sm[i]>=sm[i+1]:peaks.append((sm[i],cent[i],hist[i]))
print('hist_peaks',sorted(peaks,reverse=True)[:10])
# top individual splits, collapse complements by only idx<=comp index
keep=np.arange(len(idx))<=pos
top=np.argsort(p*keep)[-20:][::-1]
for rank,k in enumerate(top[:12],1):
 c=tuple(int(ci[k]) for ci in coords)
 print('top',rank,'p_ordered',p[k],'pair_p',p[k]*(1 if idx[k]==comp9[k] else 2),'Delta',Delta[k],'c',c)
