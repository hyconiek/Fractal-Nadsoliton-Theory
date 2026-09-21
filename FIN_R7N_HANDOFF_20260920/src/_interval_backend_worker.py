from fractions import Fraction as F
from pathlib import Path
import sys,json,time
sys.set_int_max_str_digits(0)
ROOT=Path('/mnt/data/r7n_repo_root'); SRC=ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223'; sys.path.insert(0,str(SRC/'src'))
mode=sys.argv[1]
if mode=='bounded':
 sys.path.insert(0,str(ROOT)); from fin_rank7_intake_review.verify_continuation import bounded_rationals; bounded_rationals(60)
import frontier_shifted_boxes as chk
box=[[F(x) for x in p] for p in json.loads(sys.stdin.read())]
t=time.perf_counter(); o=chk.raw_shifted_box(box); dt=time.perf_counter()-t
def pack(I): return [str(I.lo),str(I.hi)]
def width(I): return str(I.hi-I.lo)
print(json.dumps({'status':o['status'],'reason':o['reason'],'seconds':dt,'P':pack(o['P']),'P1':pack(o['P1']),'c2':pack(o['c2']),'widths':{k:width(o[k]) for k in ['P','P1','c2']}}))
