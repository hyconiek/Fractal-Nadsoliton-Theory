"""R7P-091 bounded interval cover of the quartic phase-torus complement."""
from __future__ import annotations
import json,math,itertools,heapq
from pathlib import Path
import numpy as np, mpmath as mp
from .phase_cumulants import _phase_term_data
FIXTURE=(0.1131879146,0.1698528641,0.2269339093,-0.3380663037)
ROOT_RADIUS=0.05

def _terms():
 _,_,cubic,quartic=_phase_term_data(*FIXTURE)
 return [(float(c/6),np.asarray(n,float)) for c,n in cubic]+[(float(c/24),np.asarray(n,float)) for c,n in quartic]

def _grad_intervals(box,terms):
 mp.iv.dps=30; out=[]
 vars=[mp.iv.mpf([repr(float(a)),repr(float(b))]) for a,b in box]
 for i in range(3):
  s=mp.iv.mpf('0')
  for coef,n in terms:
   ang=sum(mp.iv.mpf(repr(float(n[k])))*vars[k] for k in range(3))
   s += -mp.iv.mpf(repr(float(coef*n[i])))*mp.iv.sin(ang)
  out.append((float(s.a),float(s.b)))
 return out

def _inside_root_box(box,roots,radius=ROOT_RADIUS):
 if max(b-a for a,b in box)>2*radius:return False
 corners=list(itertools.product(*[(a,b) for a,b in box]))
 for r in roots:
  ok=True
  for c in corners:
   d=np.abs(np.asarray(c)-r); d=np.minimum(d,2*math.pi-d)
   if float(np.max(d))>radius: ok=False; break
  if ok:return True
 return False

def run(root_catalog,out_path,budget=2000,initial_n=4):
 roots=[np.asarray(r['phase'],float) for r in json.load(open(root_catalog))['roots']];terms=_terms()
 edges=np.linspace(0,2*math.pi,initial_n+1)
 heap=[];counter=0
 for idx in itertools.product(range(initial_n),repeat=3):
  box=tuple((float(edges[i]),float(edges[i+1])) for i in idx)
  vol=np.prod([b-a for a,b in box]);heapq.heappush(heap,(-vol,counter,box));counter+=1
 safe=[];inside=[];processed=0
 while heap and processed<budget:
  _,_,box=heapq.heappop(heap);processed+=1
  if _inside_root_box(box,roots):inside.append(box);continue
  gi=_grad_intervals(box,terms)
  reason=next((i for i,(lo,hi) in enumerate(gi) if lo>0 or hi<0),None)
  if reason is not None:
   safe.append({'box':box,'gradient_component':reason,'interval':gi[reason]});continue
  widths=[b-a for a,b in box];k=int(np.argmax(widths));a,b=box[k];m=(a+b)/2
  for ab in ((a,m),(m,b)):
   nb=list(box);nb[k]=ab;nb=tuple(nb);vol=np.prod([y-x for x,y in nb]);heapq.heappush(heap,(-vol,counter,nb));counter+=1
 unresolved=[x[2] for x in heap]
 result={'method':'bounded mpmath.iv componentwise-gradient exclusion','budget':budget,'processed':processed,
         'initial_n':initial_n,'root_neighborhood_radius':ROOT_RADIUS,'safe_leaf_count':len(safe),
         'root_neighborhood_leaf_count':len(inside),'unresolved_leaf_count':len(unresolved),
         'complete':len(unresolved)==0,'safe_leaves':safe,'root_leaves':inside,'unresolved_leaves':unresolved}
 Path(out_path).write_text(json.dumps(result,indent=2)+'\n');return result
if __name__=='__main__':
 import sys
 o=run(sys.argv[1],sys.argv[2],int(sys.argv[3]) if len(sys.argv)>3 else 2000)
 print(json.dumps({k:v for k,v in o.items() if not k.endswith('leaves')},indent=2))
