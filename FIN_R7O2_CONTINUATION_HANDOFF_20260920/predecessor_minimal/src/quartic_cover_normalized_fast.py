from pathlib import Path
from fractions import Fraction as F
import sys,json,heapq,itertools,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR))
import scientific_rechecks as sr
ROOTS=json.load(open(ROOT/'inputs/FR223_20260916/results/R7P-089_quartic_roots.json'))['roots']
RAD=F(1,20); TWOPI=2*sr.iv.pi

def terms():
 rt3=sr.iv.sqrt(sr.I(3));r3,r4,r5,z6=map(sr.I,['0.1131879146','0.1698528641','0.2269339093','-0.3380663037']);a,b,c,d=[x/(2*rt3) for x in [r3,r4,r5,z6]]
 cub=[(6*a*a*d,(2,0,0)),(12*a*b*c,(1,1,1)),(2*b**3,(0,3,0))]
 qua=[(2*a**4,(4,0,0)),(24*a*b*b*c,(1,-2,1)),(48*a*b*c*d,(-1,1,1)),(8*a*c**3,(1,0,-3)),(24*b*c*c*d,(0,1,-2))]
 return [(coef/sr.I(6),n) for coef,n in cub]+[(coef/sr.I(24),n) for coef,n in qua]
TERMS=terms()

def inner_normalized_collars():
 out=[]
 for r in ROOTS:
  axes=[]
  for x in r['phase']:
   c=sr.I(F(str(x))); L=(c-sr.I(RAD))/TWOPI; U=(c+sr.I(RAD))/TWOPI
   Ll,Lh=sr.bounds(L); Ul,Uh=sr.bounds(U)
   axes.append((Lh,Ul)) # conservative inner interval
  out.append((r['id'],axes))
 return out
COLLARS=inner_normalized_collars()

def in_collar_fast(box):
 for rid,axes in COLLARS:
  okall=True
  for (a,b),(lo,hi) in zip(box,axes):
   ok=False
   for n in (-1,0,1):
    if lo+n<=a and b<=hi+n:ok=True;break
   if not ok:okall=False;break
  if okall:return rid
 return None

def grad_box(box):
 z=[sr.iv.mpf([sr.I(a),sr.I(b)]) for a,b in box]
 g=[sr.I(0),sr.I(0),sr.I(0)]
 for coef,n in TERMS:
  ang=TWOPI*sum((sr.I(n[k])*z[k] for k in range(3)),sr.I(0));sn=sr.iv.sin(ang)
  for k in range(3):g[k]+=-coef*sr.I(n[k])*sn
 return g

def classify(box):
 g=grad_box(box)
 for k,v in enumerate(g):
  lo,hi=sr.bounds(v)
  if lo>0 or hi<0:return 'GRADIENT',k,[str(lo),str(hi)]
 rid=in_collar_fast(box)
 if rid is not None:return 'ROOT_COLLAR',rid,None
 return 'UNRESOLVED',None,None

def run(budget=10000):
 heap=[];counter=0
 for idx in itertools.product(range(4),repeat=3):
  box=tuple((F(i,4),F(i+1,4)) for i in idx);heapq.heappush(heap,(-F(1,64),counter,box,''));counter+=1
 safe=[];roots=[];processed=0;st=time.time()
 while heap and processed<budget:
  _,_,box,path=heapq.heappop(heap);processed+=1;reason,data,iv=classify(box)
  if reason=='GRADIENT':safe.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'gradient_component':data,'interval':iv});continue
  if reason=='ROOT_COLLAR':roots.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'root_id':data});continue
  widths=[b-a for a,b in box];ax=max(range(3),key=lambda i: widths[i]);a,b=box[ax];m=(a+b)/2
  for tag,ab in [('L',(a,m)),('R',(m,b))]:
   nb=list(box);nb[ax]=ab;nb=tuple(nb);vol=(nb[0][1]-nb[0][0])*(nb[1][1]-nb[1][0])*(nb[2][1]-nb[2][0]);heapq.heappush(heap,(-vol,counter,nb,path+str(ax)+tag));counter+=1
 unresolved=[{'path':p,'box':[[str(a),str(b)] for a,b in b]} for _,_,b,p in heap]
 out={'task':'R7N-037','method':'exact normalized dyadic torus; exact-decimal resonance intervals; precomputed conservative inner radius-0.05 certified collars','budget':budget,'processed':processed,'safe_leaf_count':len(safe),'root_collar_leaf_count':len(roots),'unresolved_leaf_count':len(unresolved),'complete':not unresolved,'elapsed_seconds':time.time()-st,'safe_leaves':safe,'root_leaves':roots,'unresolved_leaves':unresolved}
 (ROOT/'checkpoints/R7N-037_quartic_cover_normalized.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if not k.endswith('_leaves')},indent=2))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 10000)
