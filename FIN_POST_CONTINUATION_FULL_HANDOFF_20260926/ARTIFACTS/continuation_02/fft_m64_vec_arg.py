import math,time,sys
import numpy as np
from scipy import fft
from scipy.special import logsumexp
alpha=float(sys.argv[1]) if len(sys.argv)>1 else 0.745
shape=(9,)*8
NGRID=9**8
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
    x=int(h,16); return tuple((x>>(4*i))&15 for i in range(8))
def flat9(c):
    x=0
    for v in c:x=x*9+v
    return x
def flat5(c):
    x=0
    for v in c:x=x*5+v
    return x

def norm_one(c):
    a=np.asarray(c,float);return float(a@G@a)

# load level16 logs and fill dense homogeneous input A
rows=[]; logy1=[]
z2=np.full(5**8,np.nan); z4_all2=None
beta16=alpha/2
with open(sys.argv[2] if len(sys.argv)>2 else 'level16_a0745.txt') as f:
    for line in f:
        h,s1,s2,s4=line.split(); c=decode_hex(h); zz1=float(s1); n=norm_one(c); ly1=zz1-beta16*n/32
        rows.append((c,ly1));
        if s2!='nan':
            ly2=float(s2)-alpha*n/32; z2[flat5(c)]=ly2
        if c==(2,)*8 and s4!='nan': z4_all2=float(s4)-(2*alpha)*n/32
print('loaded',len(rows),'z2',np.isfinite(z2).sum(),'z4half',z4_all2,flush=True)
M=max(x[1] for x in rows)
A=np.zeros(shape,float)
for c,ly in rows:A[c]=math.exp(ly-M)
print('A ready',A.nbytes/1e6,'M',M,flush=True)
t=time.time();F=fft.rfftn(A,workers=5);print('rfft',time.time()-t,flush=True);del A
F*=F
t=time.time();C=fft.irfftn(F,s=shape,workers=5);print('irfft',time.time()-t,'min',C.min(),flush=True);del F
# Build total-degree-32 flat indices without Python composition loops.
t=time.time(); flat=np.arange(NGRID,dtype=np.int64); tmp=flat.copy(); digits=[]; total=np.zeros(NGRID,dtype=np.uint8)
for _ in range(8):
    d=(tmp%9).astype(np.uint8);digits.append(d);total+=d;tmp//=9
digits=digits[::-1] # c0..c7
mask=(total==32); idx=flat[mask]; del flat,tmp,total,mask
coords=[d[idx] for d in digits]; del digits
print('idx32',idx.size,'build sec',time.time()-t,flush=True)
# norm on selected coordinates using Gram matrix
normv=np.zeros(idx.size,float)
for i in range(8):
    ci=coords[i].astype(float); normv += G[i,i]*ci*ci
    for j in range(i): normv += 2*G[i,j]*ci*coords[j]
cv=C.ravel()[idx]; del C
if np.any(cv<=0): print('WARNING nonpositive desired',np.count_nonzero(cv<=0),cv.min(),flush=True)
lord=np.log(cv)+2*M
# equal-child correction via q2 logY table
all_even=np.ones(idx.size,dtype=bool)
for ci in coords: all_even &= ((ci&1)==0)
if all_even.any():
    hidx=np.zeros(all_even.sum(),dtype=np.int64)
    for ci in coords: hidx=hidx*5+(ci[all_even]//2)
    ly2=z2[hidx]
    assert np.isfinite(ly2).all()
    lord[all_even]=np.logaddexp(lord[all_even],ly2)
beta32=alpha/4
logY32=beta32*normv/64-math.log(2)+lord
# M32 regression root c=4^8
c4idx=flat9((4,)*8); pos=np.searchsorted(idx,c4idx);assert idx[pos]==c4idx
z32=logY32[pos]+beta32*normv[pos]/64
print('M32 logZ',z32,'resid',z32+6.9887714979191902,flush=True)
# q2 value for M32 equal root child c=(4,...,4), using z2 homogeneous array.
valid=np.flatnonzero(np.isfinite(z2)); # base5 digits total16 only
# complement in base5: digits -> 4-digit; compute via digit decode
vv=valid.copy(); digs5=[]
for _ in range(8):digs5.append(vv%5);vv//=5
digs5=digs5[::-1]
comp=np.zeros(valid.size,dtype=np.int64)
for d in digs5:comp=comp*5+(4-d)
pos2=np.searchsorted(valid,comp);assert np.all(valid[pos2]==comp)
lord2=float(logsumexp(z2[valid]+z2[comp]))
# z2 already logY16 at child q2; equal correction is logY16 q4(all2)
lsum2=float(np.logaddexp(lord2,z4_all2))
n4=norm_one((4,)*8); beta32q2=alpha/2
ly32q2=beta32q2*n4/64-math.log(2)+lsum2
print('ly32q2',ly32q2,'lord2',lord2,'diag',z4_all2,flush=True)
# root ordered complement using vectorized base9 complement positions
comp9=np.zeros(idx.size,dtype=np.int64)
for ci in coords:comp9=comp9*9+(8-ci)
posc=np.searchsorted(idx,comp9);assert np.all(idx[posc]==comp9)
logordroot=float(logsumexp(logY32+logY32[posc]))
root=(8,)*8; nr=norm_one(root); beta64=alpha/8
lsr=float(np.logaddexp(logordroot,ly32q2))
ly64=beta64*nr/128-math.log(2)+lsr
z64=ly64+beta64*nr/128
print('ROOT ordered',logordroot,'equal',ly32q2,flush=True)
print('M64 alpha',alpha,'logZ',z64,flush=True)
open((sys.argv[3] if len(sys.argv)>3 else 'fft_m64_a0745_OUTPUT.txt'),'w').write(f'M32_validation_logZ {z32:.17g}\nM32_validation_residual {z32+6.9887714979191902:.17g}\nM64_alpha {alpha:.17g}\nM64_logZ {z64:.17g}\nroot_ordered_log {logordroot:.17g}\nroot_equal_logY32_q2 {ly32q2:.17g}\n')
