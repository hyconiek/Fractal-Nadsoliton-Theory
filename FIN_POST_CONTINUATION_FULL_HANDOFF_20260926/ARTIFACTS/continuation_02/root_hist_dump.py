import math,sys,time
import numpy as np
from scipy import fft
from scipy.special import logsumexp
alpha=float(sys.argv[1]); infile=sys.argv[2]
shape=(9,)*8;NGRID=9**8
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

def dec(h):
 x=int(h,16);return tuple((x>>(4*i))&15 for i in range(8))
def norm(c):
 a=np.asarray(c,float);return float(a@G@a)
def flat5(c):
 x=0
 for v in c:x=x*5+v
 return x
rows=[];z2=np.full(5**8,np.nan);z4half=None;beta16=alpha/2
with open(infile) as f:
 for line in f:
  h,s1,s2,s4=line.split();c=dec(h);n=norm(c);rows.append((c,float(s1)-beta16*n/32))
  if s2!='nan':z2[flat5(c)]=float(s2)-alpha*n/32
  if c==(2,)*8 and s4!='nan':z4half=float(s4)-(2*alpha)*n/32
M=max(v for _,v in rows);A=np.zeros(shape)
for c,v in rows:A[c]=math.exp(v-M)
F=fft.rfftn(A,workers=5);F*=F;C=fft.irfftn(F,s=shape,workers=5);del A,F
flat=np.arange(NGRID,dtype=np.int64);tmp=flat.copy();digits=[];tot=np.zeros(NGRID,dtype=np.uint8)
for _ in range(8):
 d=(tmp%9).astype(np.uint8);digits.append(d);tot+=d;tmp//=9
digits=digits[::-1];mask=(tot==32);idx=flat[mask];coords=[d[idx] for d in digits];del flat,tmp,tot,mask
normv=np.zeros(idx.size)
for i in range(8):
 ci=coords[i].astype(float);normv+=G[i,i]*ci*ci
 for j in range(i):normv+=2*G[i,j]*ci*coords[j]
cv=C.ravel()[idx];del C
lord=np.log(cv)+2*M
all_even=np.ones(idx.size,dtype=bool)
for ci in coords:all_even&=((ci&1)==0)
if all_even.any():
 hidx=np.zeros(all_even.sum(),dtype=np.int64)
 for ci in coords:hidx=hidx*5+(ci[all_even]//2)
 lord[all_even]=np.logaddexp(lord[all_even],z2[hidx])
beta32=alpha/4
logY32=beta32*normv/64-math.log(2)+lord
comp9=np.zeros(idx.size,dtype=np.int64)
for ci in coords:comp9=comp9*9+(8-ci)
pos=np.searchsorted(idx,comp9);assert np.all(idx[pos]==comp9)
logw=logY32+logY32[pos];lognorm=float(logsumexp(logw));p=np.exp(logw-lognorm)
# ordered root distribution stats
H=-float(np.sum(p*np.log(p)))
Neff=math.exp(H)
q2=np.zeros(idx.size)
for ci in coords:q2+=(ci.astype(float)-4.)**2
root_delta=normv/16.0 # root total feature sum is zero
for name,x in [('q2',q2),('root_delta',root_delta)]:
 mu=float(np.sum(p*x));var=float(np.sum(p*(x-mu)**2));print(name,'mean',mu,'sd',math.sqrt(var))
print('ordered_entropy',H,'effective_ordered_splits',Neff,'maxprob',float(p.max()))
# exact histogram of root compositional segregation q2
q2i=q2.astype(int)
h=np.bincount(q2i,weights=p,minlength=129)
loc=[]
for q in range(1,128):
    if h[q]>=h[q-1] and h[q]>=h[q+1] and h[q]>1e-6: loc.append((q,float(h[q])))
print('q2_hist_modes',sorted(loc,key=lambda x:x[1],reverse=True)[:12])
print('q2_tail_ge96',float(h[96:].sum()),'q2_low_le32',float(h[:33].sum()))
# equal composition probability in unordered root recursion
# exact diagonal weight uses Y32(q=2), but use value from existing derivation
# recompute q2 root child quickly using z2
valid=np.flatnonzero(np.isfinite(z2));vv=valid.copy();d5=[]
for _ in range(8):d5.append(vv%5);vv//=5
d5=d5[::-1];comp=np.zeros(valid.size,dtype=np.int64)
for d in d5:comp=comp*5+(4-d)
pos2=np.searchsorted(valid,comp); lord2=float(logsumexp(z2[valid]+z2[comp]));n4=norm((4,)*8)
ly32q2=(alpha/2)*n4/64-math.log(2)+float(np.logaddexp(lord2,z4half))
peq=math.exp(ly32q2-float(np.logaddexp(lognorm,ly32q2)))
print('equal_child_probability',peq)
# top ordered compositions
sel=np.argpartition(p,-12)[-12:];sel=sel[np.argsort(p[sel])[::-1]]
for k in sel:
 c=tuple(int(ci[k]) for ci in coords)
 print('top',c,'p_ordered',float(p[k]),'q2',float(q2[k]),'delta',float(root_delta[k]))
print('FULL_HIST_BEGIN')
for q,val in enumerate(h):
    if val>1e-8: print(q,float(val))
print('FULL_HIST_END')
