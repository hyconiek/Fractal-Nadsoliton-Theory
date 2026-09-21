from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,heapq,itertools,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR));import scientific_rechecks as sr
sr.iv.dps=30
sys.path.insert(0,str(ROOT/'src'));import k16_surrogate_eval as ks
FULL=json.load(open(ROOT/'inputs/FR223_20260916/certificates/R7P-092_full_phase_roots.json'))['roots']
COL=json.load(open(ROOT/'results/R7N-041_full_uniqueness_collars.json'))['roots']
TWOPI=2*sr.iv.pi

def inner_collars():
 out=[]
 for rec,col in zip(FULL,COL):
  rad=F(col['certified_radius']); axes=[]
  for x in rec['phase']:
   c=sr.I(F(str(x))); L=(c-sr.I(rad))/TWOPI; U=(c+sr.I(rad))/TWOPI
   Ll,Lh=sr.bounds(L); Ul,Uh=sr.bounds(U); axes.append((Lh,Ul))
  out.append((rec['quartic_id'],axes,rad))
 return out
COLLARS=inner_collars()

def in_collar(box):
 for rid,axes,rad in COLLARS:
  okall=True
  for (a,b),(lo,hi) in zip(box,axes):
   ok=False
   for n in (-1,0,1):
    if lo+n<=a and b<=hi+n: ok=True;break
   if not ok:okall=False;break
  if okall:return rid
 return None

def classify(box):
 ok,k,iv=ks.classify(box)
 if ok:return 'FULL_SAFE_BY_K16',k,iv
 rid=in_collar(box)
 if rid is not None:return 'FULL_ROOT_COLLAR',rid,None
 return 'UNRESOLVED',None,None

def fresh_heap():
 heap=[];counter=0
 for idx in itertools.product(range(4),repeat=3):
  box=tuple((F(i,4),F(i+1,4)) for i in idx);heapq.heappush(heap,(-F(1,64),counter,box,''));counter+=1
 return heap,counter,[],[]

def load_checkpoint():
 p=ROOT/'checkpoints/R7N-044_full_cover_k16.json'
 if not p.exists():return fresh_heap()
 d=json.load(open(p));heap=[];counter=0
 for rec in d['unresolved_leaves']:
  box=tuple((F(a),F(b)) for a,b in rec['box']);vol=(box[0][1]-box[0][0])*(box[1][1]-box[1][0])*(box[2][1]-box[2][0]);heapq.heappush(heap,(-vol,counter,box,rec['path']));counter+=1
 return heap,counter,d['safe_leaves'],d['root_leaves']

def run(additional_budget=10000):
 heap,counter,safe,roots=load_checkpoint(); processed0=0
 p=ROOT/'checkpoints/R7N-044_full_cover_k16.json'
 if p.exists(): processed0=json.load(open(p)).get('processed_total',json.load(open(p)).get('processed',0))
 processed=0;st=time.time()
 while heap and processed<additional_budget:
  _,_,box,path=heapq.heappop(heap);processed+=1
  reason,data,iv=classify(box)
  if reason=='FULL_SAFE_BY_K16':safe.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'gradient_component':data,'surrogate_interval':iv});continue
  if reason=='FULL_ROOT_COLLAR':roots.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'root_id':data});continue
  widths=[b-a for a,b in box];ax=max(range(3),key=lambda i:widths[i]);a,b=box[ax];m=(a+b)/2
  for tag,ab in [('L',(a,m)),('R',(m,b))]:
   nb=list(box);nb[ax]=ab;nb=tuple(nb);vol=(nb[0][1]-nb[0][0])*(nb[1][1]-nb[1][0])*(nb[2][1]-nb[2][0]);heapq.heappush(heap,(-vol,counter,nb,path+str(ax)+tag));counter+=1
 unresolved=[{'path':path,'box':[[str(a),str(b)] for a,b in box]} for _,_,box,path in heap]
 out={'task':'R7N-044','method':'K16 retained 25-resonance interval gradient plus rigorous uniform full-gradient error bound; full log-mgf injectivity collars','epsilon':str(ks.EPS),'processed_total':processed0+processed,'processed_last_chunk':processed,'safe_leaf_count':len(safe),'root_collar_leaf_count':len(roots),'unresolved_leaf_count':len(unresolved),'complete':not unresolved,'elapsed_seconds_last_chunk':time.time()-st,'safe_leaves':safe,'root_leaves':roots,'unresolved_leaves':unresolved}
 p.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if not k.endswith('_leaves')},indent=2))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 10000)
