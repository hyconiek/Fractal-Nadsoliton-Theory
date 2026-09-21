from pathlib import Path
from fractions import Fraction as F
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign')
IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR))
import scientific_rechecks as sr
OLD=ROOT/'inputs/FR223_20260916/results/R7P-091_quartic_complement_partial.json'

def exact_float(x): return F.from_float(float(x))

def run():
 d=json.load(open(OLD)); out=[];lost=[];st=time.time()
 for idx,rec in enumerate(d['safe_leaves']):
  box=[[exact_float(a),exact_float(b)] for a,b in rec['box']]
  phi=[sr.I(a)+sr.iv.mpf([sr.I(b-a),sr.I(b-a)])*0 for a,b in []] # unused
  vars=[]
  for a,b in box:
   vars.append(sr.iv.mpf(a.numerator)/a.denominator + sr.iv.mpf([0,1])*(sr.iv.mpf((b-a).numerator)/(b-a).denominator))
  g,_=sr.phase_FH(vars,'quartic')
  intervals=[];reason=None
  for j,v in enumerate(g):
   lo,hi=sr.bounds(v);intervals.append([str(lo),str(hi)])
   if reason is None and (lo>0 or hi<0):reason=j
  row={'old_index':idx,'box_binary_exact':[[str(a),str(b)] for a,b in box],'old_reason':rec['gradient_component'],'new_reason':reason,'gradient_intervals':intervals,'safe':reason is not None}
  out.append(row)
  if reason is None:lost.append(idx)
 res={'task':'R7N-036','source_safe_count':len(d['safe_leaves']),'revalidated_safe_count':sum(x['safe'] for x in out),'lost_safe_count':len(lost),'lost_indices':lost,'elapsed_seconds':time.time()-st,'cells':out,
      'endpoint_note':'Boxes use exact binary values of the archived float endpoints. Exact torus seam beyond archived 2*pi float is handled separately.'}
 (ROOT/'results/R7N-036_revalidated_old_safe.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({k:v for k,v in res.items() if k!='cells'},indent=2))
if __name__=='__main__':run()
