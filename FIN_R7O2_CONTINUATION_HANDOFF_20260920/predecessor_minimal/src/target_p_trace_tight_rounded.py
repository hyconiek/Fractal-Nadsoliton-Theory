from fractions import Fraction as F
from pathlib import Path
import sys,math
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import target_p_trace_cover as b
Q=10**60

def floorq(x): return F((x*Q).__floor__(),Q)
def ceilq(x): return F((x*Q).__ceil__(),Q)
D2=[[ceilq(x) for x in row] for row in b.D2HI]

def weights(cell):
 out=[]
 for lo,hi in b.weight_bounds(cell):out.append((floorq(lo),ceilq(hi)))
 return out

def div_up(n,d): return ceilq(n/d)
def pair_up(wb,i,j):
 li,ui=wb[i];lj,uj=wb[j];R=sum((wb[k][0] for k in range(7) if k not in (i,j)),F(0))
 def fup(x,y):return div_up(ceilq(x*y), floorq((x+y+R)**2))
 vals=[fup(li,lj),fup(li,uj),fup(ui,lj),fup(ui,uj)]
 for x in (li,ui):
  y=min(max(x+R,lj),uj);vals.append(fup(x,y))
 for y in (lj,uj):
  x=min(max(y+R,li),ui);vals.append(fup(x,y))
 return max(vals)

def trace_e(cell):
 wb=weights(cell);tr=F(0)
 for i in range(7):
  for j in range(i+1,7): tr += ceilq(pair_up(wb,i,j)*D2[i][j])
 tr=ceilq(tr)
 A_lo=sum(wb[i][0] for i in range(4));A_hi=sum(wb[i][1] for i in range(4));O_lo=sum(wb[i][0] for i in range(4,7));O_hi=sum(wb[i][1] for i in range(4,7))
 elo=floorq(O_lo/(A_hi+O_lo)) if O_lo else F(0);ehi=ceilq(O_hi/(A_lo+O_hi))
 return tr,(elo,ehi)
