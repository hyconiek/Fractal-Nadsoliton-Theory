from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,math,time,heapq
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); SRC=ROOT/'inputs/FR223_20260916'; sys.path.insert(0,str(SRC/'src')); sys.path.insert(0,str(ROOT/'src'))
from generic_threshold_shifted import pow_frac_point,_constants,AL,AH
from intervals import QI,sqrt_interval
import boundary_ising as bi
TAU=F(67,250)
HULL=((F(1,900),F(1)),(F(1,128),F(1)),(F(1,9),F(1)),(F(1,1000000),F(1)))
# Canonical sigma-safe boxes (99 direct + 18 repair leaves)
reg=json.load(open(ROOT/'results/safe_union_v2_audited.json'))
MASKS=[]
for r in reg['direct_domains']: MASKS.append((r['name'],tuple(tuple(map(F,p)) for p in r['exact_decimal_box'])))
for r in reg['repaired_domains']:
    for leaf in r['leaves']: MASKS.append((r['name']+':'+leaf['path'],tuple(tuple(map(F,p)) for p in leaf['box'])))
# rstar strict interval
RSTAR=_constants(TAU)[2]
# 7 representative labels corresponding to aggregate weights:
# even 0,4,6,2; odd 3,5,1
REPS=[0,4,6,2,3,5,1]
L=bi.strict_intervals(); rt3=sqrt_interval(QI(3),40)
c3=[QI(x) for x in ([1,0,-1,0]*3)]
c4=[QI(1),QI(F(-1,2)),QI(F(-1,2))]*4
c5=[QI(1),-rt3/2,QI(F(1,2)),QI(0),QI(F(-1,2)),rt3/2,QI(-1),rt3/2,QI(F(-1,2)),QI(0),QI(F(1,2)),-rt3/2]
c6=[QI((-1)**j) for j in range(12)]
CV=[c3,c4,c5,c6]; DEN=[6,6,6,12]; KS=[3,4,5,6]
D2HI=[[F(0)]*7 for _ in range(7)]
for a,i in enumerate(REPS):
  for b,j in enumerate(REPS):
    if a>=b: continue
    d=QI(0)
    for arr,k,den in zip(CV,KS,DEN): d += L[k]/den*(arr[i]-arr[j])**2
    D2HI[a][b]=D2HI[b][a]=d.hi

def pwr_bounds(lo:F,hi:F,p:int,q:int):
    # positive exponent, 0<lo<=hi<=1
    return pow_frac_point(lo,p,q).lo,pow_frac_point(hi,p,q).hi

def weight_bounds(cell):
    (r0,r1),(s0,s1),(t0,t1),(y0,y1)=cell
    sr0,sr1=pwr_bounds(r0,r1,1,2)
    plus_lo,plus_hi=pwr_bounds(t0,t1,2*AH[1]+AH[0],AH[1]) # 2+AH lower at t0; overwritten upper below
    # upper for 2+sqrt3 uses smaller exponent AL
    plus_hi=pow_frac_point(t1,2*AL[1]+AL[0],AL[1]).hi
    minus_lo=pow_frac_point(t0,3,11).lo  # exponent 2-AL = 3/11 is larger -> lower
    minus_hi=pow_frac_point(t1,4,15).hi  # exponent 2-AH = 4/15 is smaller -> upper
    wb=[
      (F(1),F(1)),
      (2*s0*t0**3,2*s1*t1**3),
      (r0*t0**4,r1*t1**4),
      (2*r0*s0*t0,2*r1*s1*t1),
      (2*sr0*t0**2*y0,2*sr1*t1**2*y1),
      (2*sr0*s0*minus_lo*y0,2*sr1*s1*minus_hi*y1),
      (2*sr0*s0*plus_lo*y0,2*sr1*s1*plus_hi*y1),
    ]
    return wb

def trace_upper_and_e(cell):
    wb=weight_bounds(cell); dlo=sum(a for a,b in wb)
    num=F(0)
    for i in range(7):
      for j in range(i+1,7): num += wb[i][1]*wb[j][1]*D2HI[i][j]
    tr=num/(dlo*dlo)
    A_lo=sum(wb[i][0] for i in range(4));A_hi=sum(wb[i][1] for i in range(4))
    O_lo=sum(wb[i][0] for i in range(4,7));O_hi=sum(wb[i][1] for i in range(4,7))
    e_lo=O_lo/(A_hi+O_lo) if O_lo else F(0)
    e_hi=O_hi/(A_lo+O_hi)
    return tr,(e_lo,e_hi)

def local_box(cell,e):
    (r0,r1),(s0,s1),(t0,t1),_=cell
    return ((r0-RSTAR.hi,r1-RSTAR.lo),(1-s1,1-s0),(1-t1,1-t0),e)

def contained(B,A): return all(a0<=b0 and b1<=a1 for (b0,b1),(a0,a1) in zip(B,A))
def mask_hit(B):
    for name,A in MASKS:
      if contained(B,A): return name
    return None

def split_point(lo,hi):
    # fixed navigation split: rationalized geometric midpoint, exact partition after proposal.
    x=math.sqrt(float(lo)*float(hi)); m=F(format(x,'.16g'))
    if not lo<m<hi: m=(lo+hi)/2
    return m

def choose_axis(cell):
    return max(range(4),key=lambda i: math.log(float(cell[i][1]/cell[i][0])))
def split(cell,axis):
    lo,hi=cell[axis];m=split_point(lo,hi); Lc=list(cell);Rc=list(cell);Lc[axis]=(lo,m);Rc[axis]=(m,hi);return tuple(Lc),tuple(Rc),m

def run(max_leaves=10000):
    stack=[(HULL,'')]; leaves=[]; stats={'SAFE_BY_MASK':0,'SAFE_BY_TRACE':0,'UNRESOLVED':0}; max_tr=F(0); start=time.monotonic()
    while stack:
      cell,path=stack.pop()
      tr,e=trace_upper_and_e(cell); max_tr=max(max_tr,tr); B=local_box(cell,e); hit=mask_hit(B)
      if hit:
        leaves.append({'path':path,'cell':[[str(a),str(b)] for a,b in cell],'reason':'SAFE_BY_MASK','dependency':hit,'trace_upper':str(tr)});stats['SAFE_BY_MASK']+=1;continue
      if tr<=2*TAU:
        leaves.append({'path':path,'cell':[[str(a),str(b)] for a,b in cell],'reason':'SAFE_BY_TRACE','trace_upper':str(tr)});stats['SAFE_BY_TRACE']+=1;continue
      # If splitting would push total frontier+terminals past budget, retain unresolved.
      if len(leaves)+len(stack)+2>max_leaves:
        leaves.append({'path':path,'cell':[[str(a),str(b)] for a,b in cell],'reason':'UNRESOLVED','trace_upper':str(tr),'local_box':[[str(a),str(b)] for a,b in B]});stats['UNRESOLVED']+=1;continue
      ax=choose_axis(cell); l,r,m=split(cell,ax); stack.append((r,path+str(ax)+'R'));stack.append((l,path+str(ax)+'L'))
    out={'task':'R7N-019','threshold':'67/250','root_hull':[[str(a),str(b)] for a,b in HULL],
         'max_leaf_budget':max_leaves,'leaf_count':len(leaves),'stats':stats,'complete':stats['UNRESOLVED']==0,
         'elapsed_seconds':time.monotonic()-start,'max_trace_upper_seen':str(max_tr),
         'proof_methods':['canonical sigma-safe mask containment','PSD trace criterion trace(M4)<=2*tau0'],
         'split_policy':'axis with largest log endpoint ratio; exact rationalized geometric-midpoint partition','leaves':leaves}
    (ROOT/'checkpoints/R7N-019_trace_cover.json').write_text(json.dumps(out,indent=2)+'\n')
    return out
if __name__=='__main__':
 out=run(); print(json.dumps({k:v for k,v in out.items() if k!='leaves'},indent=2))

def pair_prob_product_upper(wb,i,j):
    li,ui=wb[i]; lj,uj=wb[j]; R=sum((wb[k][0] for k in range(len(wb)) if k not in (i,j)),F(0))
    def f(x,y): return x*y/(x+y+R)**2
    vals=[f(li,lj),f(li,uj),f(ui,lj),f(ui,uj)]
    # Along x fixed, optimum y=x+R clipped to [lj,uj]; vice versa.
    for x in (li,ui):
        y=min(max(x+R,lj),uj); vals.append(f(x,y))
    for y in (lj,uj):
        x=min(max(y+R,li),ui); vals.append(f(x,y))
    return max(vals)

def trace_upper_tight_and_e(cell):
    wb=weight_bounds(cell); num=F(0)
    for i in range(7):
      for j in range(i+1,7): num += pair_prob_product_upper(wb,i,j)*D2HI[i][j]
    A_lo=sum(wb[i][0] for i in range(4));A_hi=sum(wb[i][1] for i in range(4))
    O_lo=sum(wb[i][0] for i in range(4,7));O_hi=sum(wb[i][1] for i in range(4,7))
    e_lo=O_lo/(A_hi+O_lo) if O_lo else F(0); e_hi=O_hi/(A_lo+O_hi)
    return num,(e_lo,e_hi)
