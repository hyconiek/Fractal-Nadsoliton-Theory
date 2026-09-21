from pathlib import Path
from fractions import Fraction as F
import sys,itertools,numpy as np
ROOT=Path('/mnt/data/fin_rank7_next_campaign');H=ROOT/'inputs/FR223_20260916';sys.path[:0]=[str(ROOT/'src'),str(H/'src')]
from intervals import QI
import target_p_trace_tight_rounded as tr
import compression_interval_probe as old
import off_face
TAU=F(67,250)

def center_normal_and_c(cell,den=10**6):
    x=np.array([float((a+b)/2) for a,b in cell]); p=off_face.p_from_aligned_compact(*x)
    _,_,_,_,C4=off_face.constants(); mu=p@C4; Y=C4-mu; M=Y.T@(p[:,None]*Y); w,V=np.linalg.eigh(M); n0=V[:,-1]
    # orientation irrelevant; rational integer normal
    n=np.rint(n0*den).astype(np.int64)
    if not np.any(n): raise ValueError('zero normal')
    # exact-rational center near actual mean
    c=[F(format(float(v),'.16g')) for v in mu]
    return n,w,c

def h_upper_for_state(i,n,c):
    n=[F(int(x)) for x in n]; n2=sum(x*x for x in n)
    d=[old.OBS[i][r]-QI(c[r]) for r in range(4)]
    norm2=QI(0); nd=QI(0)
    for r in range(4): norm2 += d[r]*d[r]; nd += n[r]*d[r]
    h=norm2-(nd*nd)/n2
    # P is PSD analytically, so clamp only tiny negative lower enclosure; upper is rigorous.
    return max(F(0),h.hi)

def max_weighted_average(wb,h):
    # For a linear-fractional weighted average, an optimum over the box has
    # upper weights on a prefix of h sorted descending and lower weights on the rest.
    order=sorted(range(len(h)),key=lambda i:h[i],reverse=True)
    best=F(-1)
    for k in range(len(h)+1):
      hi=set(order[:k]); ws=[wb[i][1] if i in hi else wb[i][0] for i in range(len(wb))]
      D=sum(ws); val=sum(ws[i]*h[i] for i in range(len(wb)))/D
      if val>best:best=val
    return best

def certify(cell,den=10**6):
    n,eigs,c=center_normal_and_c(cell,den); wb=tr.weights(cell)
    h=[h_upper_for_state(i,n,c) for i in range(7)]
    u=max_weighted_average(wb,h)
    return {'ok':u<=TAU,'reason':'PROJECTED_TRACE' if u<=TAU else 'FAILED','upper':str(u),'normal_num':[int(x) for x in n],
            'normal_scale':den,'center_c':[str(x) for x in c],'center_eigs':eigs.tolist(),'h_upper':[str(x) for x in h]}
