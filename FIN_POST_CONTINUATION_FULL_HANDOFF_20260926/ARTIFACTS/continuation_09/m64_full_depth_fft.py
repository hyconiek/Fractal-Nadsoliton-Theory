import math, numpy as np
from scipy import fft
from scipy.special import logsumexp
alpha=0.8245
path='/mnt/data/FIN_NEXT_RESEARCH_20260925/continuation_09/level16_depth_08245.txt'
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]],float)
G=V@V.T; shape=(9,)*8; NGRID=9**8

def decode_hex(h):
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def flat5(c):
 x=0
 for v in c:x=x*5+v
 return x
def normc(c):
 a=np.asarray(c,float);return float(a@G@a)
# parse raw16 records
rows=[]; q2={}; q4={}
with open(path) as f:
 for line in f:
  p=line.split(); c=decode_hex(p[0]); vals=list(map(lambda x:float(x) if x!='nan' else np.nan,p[1:]))
  s1=np.array(vals[0:7]); rows.append((c,s1))
  if np.isfinite(vals[7]): q2[c]=np.array(vals[7:14])
  if np.isfinite(vals[14]): q4[c]=np.array(vals[14:21])
print('rows',len(rows),'q2',len(q2),'q4',len(q4),flush=True)
# valid total32 index geometry
flat=np.arange(NGRID,dtype=np.int64); tmp=flat.copy(); digits=[]; total=np.zeros(NGRID,dtype=np.uint8)
for _ in range(8):
 d=(tmp%9).astype(np.uint8);digits.append(d);total+=d;tmp//=9
digits=digits[::-1];mask=(total==32);idx=flat[mask];coords=[d[idx] for d in digits]
comp9=np.zeros(idx.size,dtype=np.int64)
for ci in coords:comp9=comp9*9+(8-ci)
posc=np.searchsorted(idx,comp9);assert np.all(idx[posc]==comp9)
normv=np.zeros(idx.size,float)
for i in range(8):
 ci=coords[i].astype(float);normv+=G[i,i]*ci*ci
 for j in range(i):normv+=2*G[i,j]*ci*coords[j]
all_even=np.ones(idx.size,dtype=bool)
for ci in coords:all_even&=((ci&1)==0)
# build q1 Y16 arrays
logys=[]; ds=[]; levs=[]
for c,s in rows:
 n=normc(c);logys.append(s[0]-alpha*n/64);ds.append(s[1]-n/64);levs.append(s[2:6])
logys=np.array(logys);ds=np.array(ds);levs=np.array(levs);M=float(logys.max())
A=np.zeros(shape,float);Ad=np.zeros(shape,float);Ad2=np.zeros(shape,float);AL=[np.zeros(shape,float) for _ in range(4)]
for (c,_),ly,d,lv in zip(rows,logys,ds,levs):
 z=math.exp(ly-M);A[c]=z;Ad[c]=z*d;Ad2[c]=z*d*d
 for k in range(4):AL[k][c]=z*lv[k]
F0=fft.rfftn(A,workers=5);del A

def conv_extract(B):
 FB=fft.rfftn(B,workers=5); del B
 C=fft.irfftn(FB*F0,s=shape,workers=5); del FB
 out=C.ravel()[idx].copy(); del C
 return out
S0=conv_extract(np.exp(0)*np.zeros(shape)) if False else None
# S0 from F0^2
C=fft.irfftn(F0*F0,s=shape,workers=5);S0=C.ravel()[idx].copy();del C
S1=2*conv_extract(Ad)
Ssq=2*conv_extract(Ad2)
# +2 conv(Ad,Ad)
# rebuild Ad cheaply
Ad=np.zeros(shape,float)
for (c,_),ly,d in zip(rows,logys,ds):Ad[c]=math.exp(ly-M)*d
Fd=fft.rfftn(Ad,workers=5);del Ad
C=fft.irfftn(Fd*Fd,s=shape,workers=5);Ssq+=2*C.ravel()[idx];del C,Fd
Lnum=[]
for k in range(4):Lnum.append(2*conv_extract(AL[k]))
del F0
# q2 equal branch per even state
# make dict transformed y2 keyed half composition c16; arrays aligned only even idx
Eq=np.zeros(idx.size);Dq=np.zeros(idx.size);Lq=np.zeros((4,idx.size))
if all_even.any():
 inds=np.flatnonzero(all_even)
 for t,ii in enumerate(inds):
  c=tuple(int(coords[d][ii]//2) for d in range(8)); s=q2[c]; n=normc(c)
  ly=s[0]-alpha*n/32; Eq[ii]=math.exp(ly-2*M); Dq[ii]=s[1]-n/32; Lq[:,ii]=s[2:6]
Den=S0+Eq
mean=(S1+Eq*Dq)/Den
root0=(Ssq+Eq*Dq*Dq)/Den-mean*mean
lev32=np.empty((4,idx.size));lev32[0]=root0
for k in range(3):lev32[k+1]=(Lnum[k]+Eq*Lq[k])/Den
# total check all levels and d2 using q2 lev3 omitted contributes to depth4? for M32 max depth root+M16 4 entries; M16 lev3 zero in outputs generally.
# logY32 + d1 transform
logY32=np.log(Den)+2*M - math.log(2)+alpha*normv/256
dY32=mean+normv/256
# q2 M32 homogeneous state (4,...,4): ordered q2 children + q4 correction
valid=list(q2.keys()); logs=[]; slopes=[]; childlev=[]
for c in valid:
 b=tuple(4-x for x in c)
 if b not in q2: continue
 s=q2[c];t=q2[b];nc=normc(c);nb=normc(b)
 ly1=s[0]-alpha*nc/32;ly2=t[0]-alpha*nb/32
 logs.append(ly1+ly2);slopes.append((s[1]-nc/32)+(t[1]-nb/32));childlev.append(s[2:6]+t[2:6])
# q4 diagonal correction c=(2,..,2)
c4=(2,)*8;s4=q4[c4];n4h=normc(c4);logs.append(s4[0]-alpha*n4h/16);slopes.append(s4[1]-n4h/16);childlev.append(s4[2:6])
logs=np.array(logs);slopes=np.array(slopes);childlev=np.array(childlev);lm=logsumexp(logs);pp=np.exp(logs-lm);md=float(pp@slopes);r0q=float(pp@((slopes-md)**2));levq=np.zeros(4);levq[0]=r0q
for k in range(3):levq[k+1]=float(pp@childlev[:,k])
# common transform for M32 q2
nroot4=normc((4,)*8); logY32q=lm-math.log(2)+alpha*nroot4/128; dY32q=md+nroot4/128
# M64 root ordered branches all idx plus equal q2 branch
lw=logY32+logY32[posc]; slope=dY32+dY32[posc]; childL=lev32+lev32[:,posc]
leq=logY32q; lognorm=np.logaddexp(logsumexp(lw),leq); p=np.exp(lw-lognorm);pe=math.exp(leq-lognorm)
mu=float(p@slope+pe*dY32q); root=float(p@((slope-mu)**2)+pe*(dY32q-mu)**2)
lev64=np.zeros(5);lev64[0]=root
for k in range(4):lev64[k+1]=float(p@childL[k]+pe*levq[k])
CH=alpha*alpha*lev64
print('M32 validation homogeneous index',flush=True)
# homogeneous c=(4,..4) locate index
flat_h=0
for x in (4,)*8:flat_h=flat_h*9+x
pos=np.searchsorted(idx,flat_h);print('logY32_hom',logY32[pos],'dY',dY32[pos],'lev',lev32[:,pos])
print('M32q2 logY',logY32q,'dY',dY32q,'lev',levq)
print('M64 peq',pe,'d1',mu)
for k,x in enumerate(CH):print('depth',k,'CH',repr(float(x)),'fraction',repr(float(x/CH.sum())))
print('sumCH',repr(float(CH.sum())))
