from fractions import Fraction as F
from pathlib import Path
import sys,itertools,math
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import target_p_trace_cover as b
import target_p_trace_tight_rounded as tr
from intervals import QI
Q=10**60
def floorq(x):return F((x*Q).__floor__(),Q)
def ceilq(x):return F((x*Q).__ceil__(),Q)
# squared parallelogram area for each triple of 7 representative feature vectors
AREA={}
for ia,ib,ic in itertools.combinations(range(7),3):
    labels=[b.REPS[ia],b.REPS[ib],b.REPS[ic]];i,j,k=labels
    def dotdiff(a,c,d):
        out=QI(0)
        for arr,kk,den in zip(b.CV,b.KS,b.DEN):out += b.L[kk]/den*(arr[c]-arr[a])*(arr[d]-arr[a])
        return out
    g11=dotdiff(i,j,j);g22=dotdiff(i,k,k);g12=dotdiff(i,j,k);det=g11*g22-g12*g12
    AREA[(ia,ib,ic)]=ceilq(max(F(0),det.hi))

def f2max(yb,zb,C):
    ly,uy=yb;lz,uz=zb
    def f(y,z):return y*z/(C+y+z)**3
    vals=[f(ly,lz),f(ly,uz),f(uy,lz),f(uy,uz)]
    # interior y=z=C
    if ly<=C<=uy and lz<=C<=uz:vals.append(f(C,C))
    for y in (ly,uy):
        z=min(max((C+y)/2,lz),uz);vals.append(f(y,z))
    for z in (lz,uz):
        y=min(max((C+z)/2,ly),uy);vals.append(f(y,z))
    return max(vals)

def triple_prob_upper(wb,inds):
    inds=tuple(inds);R=sum((wb[k][0] for k in range(7) if k not in inds),F(0));best=F(0)
    for pos,k in enumerate(inds):
        others=[q for q in inds if q!=k]
        for x in wb[k]:
            v=x*f2max(wb[others[0]],wb[others[1]],R+x)
            if v>best:best=v
    return ceilq(best)

def e2_upper(cell):
    wb=tr.weights(cell);s=F(0)
    for inds,a in AREA.items():s += ceilq(triple_prob_upper(wb,inds)*a)
    return ceilq(s)
