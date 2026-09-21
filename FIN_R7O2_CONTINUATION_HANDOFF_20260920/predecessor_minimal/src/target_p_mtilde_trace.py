from fractions import Fraction as F
from pathlib import Path
import sys,itertools
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); sys.path[:0]=[str(ROOT/'src'),str(ROOT/'inputs/FR223_20260916/src')]
import target_p_trace_tight_rounded as tr
import target_p_trace_cover as b
from intervals import QI
from compression_interval_probe import OBS
TAU=F(67,250)
Q=10**60

def floorq(x): return F((x*Q).__floor__(),Q)
def ceilq(x): return F((x*Q).__ceil__(),Q)
# first-three-coordinate squared distances
D3=[[F(0)]*7 for _ in range(7)]
for i in range(7):
  for j in range(i+1,7):
    s=QI(0)
    for k in range(3): s += (OBS[i][k]-OBS[j][k])**2
    D3[i][j]=D3[j][i]=ceilq(s.hi)

def pair_up_subset(wb,inds,i,j):
    li,ui=wb[i];lj,uj=wb[j]
    R=sum((wb[k][0] for k in inds if k not in (i,j)),F(0))
    def fup(x,y):
      den=floorq((x+y+R)**2); num=ceilq(x*y)
      return ceilq(num/den)
    vals=[fup(li,lj),fup(li,uj),fup(ui,lj),fup(ui,uj)]
    for x in (li,ui):
      y=min(max(x+R,lj),uj); vals.append(fup(x,y))
    for y in (lj,uj):
      x=min(max(y+R,li),ui); vals.append(fup(x,y))
    return max(vals)

def cond_trace_up(wb,inds):
    s=F(0)
    for i,j in itertools.combinations(inds,2): s=ceilq(s+ceilq(pair_up_subset(wb,inds,i,j)*D3[i][j]))
    return s

def cond_mean_intervals(wb,inds):
    dlo=sum(wb[i][0] for i in inds); dhi=sum(wb[i][1] for i in inds)
    D=QI(dlo,dhi); out=[]
    for k in range(3):
      n=QI(0)
      for i in inds: n += QI(wb[i][0],wb[i][1])*OBS[i][k]
      out.append(n/D)
    return out

def q_interval(wb):
    A0=sum(wb[i][0] for i in range(4)); A1=sum(wb[i][1] for i in range(4))
    O0=sum(wb[i][0] for i in range(4,7)); O1=sum(wb[i][1] for i in range(4,7))
    qlo=floorq(A0/(A0+O1)); qhi=ceilq(A1/(A1+O0))
    if qlo<F(1,2): qlo=F(1,2) # accepted parity theorem q>=1/2 on physical domain
    return qlo,qhi

def upper(cell):
    wb=tr.weights(cell); E=(0,1,2,3); O=(4,5,6)
    te=cond_trace_up(wb,E); to=cond_trace_up(wb,O)
    qlo,qhi=q_interval(wb)
    # q*te + (1-q)*to, maximize affine bound over physical q interval.
    mix=max(qlo*te+(1-qlo)*to, qhi*te+(1-qhi)*to)
    me=cond_mean_intervals(wb,E); mo=cond_mean_intervals(wb,O)
    d2=F(0)
    for a,c in zip(me,mo):
      z=a-c; m=max(abs(z.lo),abs(z.hi)); d2=ceilq(d2+ceilq(m*m))
    # q>=1/2 => q(1-q) decreases in q; maximum at qlo.
    qe=ceilq(qlo*(1-qlo))
    l6=b.L[6]
    c=ceilq(l6.hi/(3*TAU)); eta_lo=floorq(1-c*qe)
    if eta_lo<=0: raise ValueError('nonpositive eta bound')
    f=ceilq(qe/eta_lo)
    return ceilq(mix+ceilq(f*d2)), {'te':str(te),'to':str(to),'q':[str(qlo),str(qhi)],'d2':str(d2),'qe':str(qe),'eta_lo':str(eta_lo),'mix':str(mix),'f':str(f)}
