import math, numpy as np
from scipy.special import logsumexp
base='/mnt/data/FIN_NEXT_RESEARCH_20260925'
paths=[f'{base}/b8235.npy',f'{base}/b8240.npy',f'{base}/b8245.npy',f'{base}/b8250.npy',f'{base}/b8255.npy']
L=[np.load(p,mmap_mode='r') for p in paths]
E=np.array([-26.625618501391713,-26.644855707483927,-26.664092118536949,-26.683327736256764,-26.702562562342486])
h=0.0005;a0=.8245
l0=np.asarray(L[2])
lognorm=np.logaddexp(logsumexp(l0),E[2])
p=np.exp(l0-lognorm);pe=math.exp(E[2]-lognorm)
# vectorized 5pt derivatives; create one array at a time
lp=(np.asarray(L[0])-8*np.asarray(L[1])+8*np.asarray(L[3])-np.asarray(L[4]))/(12*h)
lpp=(-np.asarray(L[4])+16*np.asarray(L[3])-30*l0+16*np.asarray(L[1])-np.asarray(L[0]))/(12*h*h)
ep=(E[0]-8*E[1]+8*E[3]-E[4])/(12*h)
epp=(-E[4]+16*E[3]-30*E[2]+16*E[1]-E[0])/(12*h*h)
mean=float(p@lp+pe*ep)
between=float(p@((lp-mean)**2)+pe*(ep-mean)**2)
within=float(p@lpp+pe*epp)
total=between+within
Z=np.array([np.logaddexp(logsumexp(np.asarray(L[k])),E[k])-math.log(2) for k in range(5)])
zpp=(-Z[4]+16*Z[3]-30*Z[2]+16*Z[1]-Z[0])/(12*h*h)
print('pequal',repr(pe))
print('alpha',a0)
print('dlog_branch_mean',repr(mean))
print('d2_between',repr(between))
print('d2_within',repr(within))
print('d2_total',repr(total))
print('CH_between',repr(a0*a0*between))
print('CH_within',repr(a0*a0*within))
print('CH_total_decomp',repr(a0*a0*total))
print('fraction_between',repr(between/total))
print('fraction_within',repr(within/total))
print('CH_direct5',repr(a0*a0*zpp))
print('residual',repr(a0*a0*(zpp-total)))
print('Z5',Z.tolist())
