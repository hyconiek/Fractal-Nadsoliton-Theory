from pathlib import Path
from fractions import Fraction as F
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR))
import scientific_rechecks as sr
OLD=ROOT/'inputs/FR223_20260916/results/R7P-091_quartic_complement_partial.json'

def exact_terms():
 # a,b,c,d = respective exact-decimal amplitudes /(2 sqrt(3)); d has z6 sign.
 rt3=sr.iv.sqrt(sr.I(3)); two=sr.I(2)
 r3,r4,r5,z6=map(sr.I,['0.1131879146','0.1698528641','0.2269339093','-0.3380663037'])
 a,b,c,d=[x/(two*rt3) for x in [r3,r4,r5,z6]]
 cubic=[(6*a*a*d,(2,0,0)),(12*a*b*c,(1,1,1)),(2*b**3,(0,3,0))]
 quartic=[(2*a**4,(4,0,0)),(24*a*b*b*c,(1,-2,1)),(48*a*b*c*d,(-1,1,1)),(8*a*c**3,(1,0,-3)),(24*b*c*c*d,(0,1,-2))]
 return [(coef/sr.I(6),n) for coef,n in cubic]+[(coef/sr.I(24),n) for coef,n in quartic]
TERMS=exact_terms()
def boxvar(a,b):
 a,b=F.from_float(float(a)),F.from_float(float(b));return sr.iv.mpf([sr.I(a),sr.I(b)])
def grad_box(box):
 vars=[boxvar(a,b) for a,b in box];g=[sr.I(0),sr.I(0),sr.I(0)]
 for coef,n in TERMS:
  ang=sum((sr.I(n[k])*vars[k] for k in range(3)),sr.I(0));sn=sr.iv.sin(ang)
  for k in range(3):g[k]+=-coef*sr.I(n[k])*sn
 return g

def run():
 d=json.load(open(OLD));rows=[];lost=[];st=time.time()
 for i,rec in enumerate(d['safe_leaves']):
  g=grad_box(rec['box']);ints=[];reason=None
  for k,v in enumerate(g):
   lo,hi=sr.bounds(v);ints.append([str(lo),str(hi)])
   if reason is None and (lo>0 or hi<0):reason=k
  rows.append({'old_index':i,'box':rec['box'],'old_reason':rec['gradient_component'],'new_reason':reason,'gradient_intervals':ints,'safe':reason is not None})
  if reason is None:lost.append(i)
 out={'task':'R7N-036','method':'exact-decimal resonance coefficients + true sqrt(3), exact archived binary-float boxes','source_safe_count':len(rows),'revalidated_safe_count':sum(x['safe'] for x in rows),'lost_safe_count':len(lost),'lost_indices':lost,'elapsed_seconds':time.time()-st,'cells':rows}
 (ROOT/'results/R7N-036_revalidated_old_safe_terms.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='cells'},indent=2))
if __name__=='__main__':run()
