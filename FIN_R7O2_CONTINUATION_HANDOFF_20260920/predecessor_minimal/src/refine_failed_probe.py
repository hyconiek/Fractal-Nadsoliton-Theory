from pathlib import Path
from fractions import Fraction as F
from functools import reduce
import sys,json,math,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import target_p_trace_cover as tc
import target_p_trace_tight_rounded as tr
import target_p_e2_bound as e2
import compression_interval_probe as cp
TAU=F(67,250);TAU2=TAU*TAU
SRC=ROOT/'checkpoints/R7N-020_trace_e2_compression_v1.json'
def parse(x):return tuple(tuple(map(F,p)) for p in x['cell'])
def vol(c):return reduce(lambda a,b:a*b,[hi-lo for lo,hi in c],F(1))
def classify(c):
 t,_=tr.trace_e(c)
 if t<=2*TAU:return 'TRACE',None
 u=e2.e2_upper(c)
 if u<=TAU2:return 'E2',None
 r=cp.certify(c)
 if r['ok']:return 'COMP',r['reason']
 return 'FAIL',None
D=json.load(open(SRC));arr=sorted(D['compression_failed'],key=lambda x:vol(parse(x)),reverse=True)[:100]
out=[];st=time.time()
for x in arr:
 c=parse(x);ax=tc.choose_axis(c);L,R,m=tc.split(c,ax)
 for tag,ch in [('L',L),('R',R)]:
  reason,sub=classify(ch);out.append({'parent':x['path'],'child':tag,'axis':ax,'reason':reason,'subreason':sub,'cell':[[str(a),str(b)] for a,b in ch]})
res={'task':'R7N-021-one-split-probe','parents':len(arr),'children':len(out),'counts':{k:sum(z['reason']==k for z in out) for k in ['TRACE','E2','COMP','FAIL']},'elapsed_seconds':time.time()-st,'results':out}
json.dump(res,open(ROOT/'results/R7N-021_one_split_probe100.json','w'),indent=2);print(json.dumps({k:v for k,v in res.items() if k!='results'},indent=2))
