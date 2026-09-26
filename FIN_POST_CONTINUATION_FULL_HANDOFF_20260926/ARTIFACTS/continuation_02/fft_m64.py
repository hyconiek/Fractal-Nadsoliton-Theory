import math, time, os
import numpy as np
from scipy import fft
from scipy.special import logsumexp

alpha=0.745
Q=8
shape=(9,)*8
V=np.array([
[4.3368086899420177e-17,0.48943915401716565,-0.2573711701187591,0.44577994304914381,-0.44591493030279611,0.25744910504599255,-0.36769135668039343],
[-0.48943915401716553,7.9797279894933126e-17,-0.25737117011875937,-0.44577994304914353,0.25744910504599255,-0.44591493030279611,0.36769135668039338],
[0.48943915401716559,-9.384495385472472e-17,-0.25737117011875882,0.44577994304914398,-0.25744910504599244,-0.44591493030279605,0.36769135668039338],
[5.6725457664441592e-16,0.48943915401716559,-0.25737117011875904,-0.44577994304914381,0.44591493030279589,0.25744910504599283,-0.36769135668039343],
[-1.6653345369377348e-16,-0.48943915401716553,-0.25737117011875871,0.44577994304914398,0.44591493030279605,-0.25744910504599261,-0.36769135668039343],
[0.48943915401716559,-2.0643209364124004e-16,-0.25737117011875987,-0.4457799430491432,-0.25744910504599267,0.44591493030279611,0.36769135668039338],
[-0.48943915401716559,1.1307217996749377e-15,-0.25737117011875937,0.4457799430491437,0.257449105045992,0.4459149303027965,0.36769135668039338],
[5.4296844798074062e-16,-0.48943915401716559,-0.25737117011876004,-0.4457799430491432,-0.44591493030279561,-0.25744910504599328,-0.36769135668039338]
],float)

def dec(code): return tuple((code>>(4*i))&15 for i in range(8))
def norm(c):
    t=np.asarray(c,float)@V
    return float(t@t)

rows=[]
with open('level16_a0745.txt') as f:
    for line in f:
        h,z1,z2,z4=line.split()
        code=int(h,16); c=dec(code)
        rows.append((c,float(z1),float(z2) if z2!='nan' else None,float(z4) if z4!='nan' else None))
print('rows',len(rows),flush=True)
# logY q1 at m16, beta=alpha/2
beta16=alpha/2
logy=[]
for c,z1,_,_ in rows: logy.append(z1-beta16*norm(c)/(2*16))
M=max(logy)
print('maxlogY16',M,'min',min(logy),flush=True)
A=np.zeros(shape,dtype=np.float64)
for (c,_,_,_),ly in zip(rows,logy): A[c]=math.exp(ly-M)
print('array MB',A.nbytes/1e6,'sum',A.sum(),flush=True)
t=time.time(); F=fft.rfftn(A,workers=5); print('rfft sec',time.time()-t,'MB',F.nbytes/1e6,flush=True)
del A
F*=F
t=time.time(); C=fft.irfftn(F,s=shape,workers=5); print('irfft sec',time.time()-t,'MB',C.nbytes/1e6,flush=True)
del F
# validate alias/noise and build logY32 q1 on desired raw compositions total32 cap8
neg=float(C.min()); print('conv min',neg,'max',float(C.max()),flush=True)
# q2/q4 lookup for level16 raw
z2map={c:z2 for c,_,z2,_ in rows if z2 is not None}
z4map={c:z4 for c,_,_,z4 in rows if z4 is not None}

def gen(total,cap=8):
    a=[0]*8
    def rec(i,rem):
        if i==8:
            if rem==0: yield tuple(a)
            return
        tail=(7-i)*cap
        for x in range(max(0,rem-tail),min(cap,rem)+1):
            a[i]=x; yield from rec(i+1,rem-x)
    yield from rec(0,total)

beta32=alpha/4
logY32={}
small_neg=0
for idx,c in enumerate(gen(32,8)):
    cv=float(C[c])
    if cv<=0:
        small_neg+=1
        cv=max(cv,1e-300)
    lord=math.log(cv)+2*M
    if all((x%2)==0 for x in c):
        h=tuple(x//2 for x in c)
        z2=z2map[h]
        ly2=z2-alpha*norm(h)/(2*16) # beta m16 q2 = alpha
        lsum=np.logaddexp(lord,ly2)
    else: lsum=lord
    ly=beta32*norm(c)/(2*32)-math.log(2)+lsum
    logY32[c]=float(ly)
print('logY32 states',len(logY32),'nonposconv',small_neg,flush=True)
# Validate n=4 M32 root c=4^8 against direct -6.9887714979191902
c4=(4,)*8; ly4=logY32[c4]; z4root=ly4+beta32*norm(c4)/(2*32)
print('M32 validation logZ',z4root,'resid',z4root-(-6.9887714979191902),flush=True)
# compute Y32 q2 for c=4^8 directly from level16 q2 values
# beta parent m32 q2=alpha/2; child m16 q2=alpha; equal child correction q4 m16 beta=2alpha
terms=[]
for a in gen(16,4):
    b=tuple(4-x for x in a)
    za=z2map[a]; zb=z2map[b]
    lya=za-alpha*norm(a)/(32)
    lyb=zb-alpha*norm(b)/(32)
    terms.append(lya+lyb)
lord2=float(logsumexp(terms))
half=(2,)*8
z4h=z4map[half]
ly4h=z4h-(2*alpha)*norm(half)/(32)
lsum2=float(np.logaddexp(lord2,ly4h))
beta32q2=alpha/2
ly32q2=beta32q2*norm(c4)/(64)-math.log(2)+lsum2
print('Y32q2 log',ly32q2,'ordered',lord2,'diag',ly4h,flush=True)
# root M64 ordered convolution coefficient by direct complement sum across 2.306m states
terms_root=[]
# avoid storing list ~18MB okay, but stream logsumexp stable incrementally
lr=-math.inf
cnt=0
for c,ly in logY32.items():
    b=tuple(8-x for x in c)
    x=ly+logY32[b]
    lr=float(np.logaddexp(lr,x)); cnt+=1
print('root ordered log',lr,'count',cnt,flush=True)
root=(8,)*8; beta64=alpha/8
lsumroot=float(np.logaddexp(lr,ly32q2))
ly64=beta64*norm(root)/(128)-math.log(2)+lsumroot
z64=ly64+beta64*norm(root)/(128)
print('M64 alpha',alpha,'logZ',z64,'logY',ly64,'diag_logY32q2',ly32q2,flush=True)
# save compact outputs needed for derivative campaign maybe not full dict
with open('fft_m64_a0745_OUTPUT.txt','w') as f:
    f.write(f'M32_validation_logZ {z4root:.17g}\n')
    f.write(f'M32_validation_residual {z4root-(-6.9887714979191902):.17g}\n')
    f.write(f'M64_alpha {alpha:.17g}\nM64_logZ {z64:.17g}\n')
    f.write(f'root_ordered_log {lr:.17g}\nroot_equal_logY32_q2 {ly32q2:.17g}\n')
