from pathlib import Path
from fractions import Fraction as F
import json,sys,time,subprocess,os,tempfile
sys.set_int_max_str_digits(0)
CAMP=Path('/mnt/data/fin_rank7_next_campaign'); AUD=Path('/mnt/data/fin_rank7_intake_review_unpacked/fin_rank7_intake_review')
union=json.load(open(AUD/'FR223_union_replay.json')); reps=json.load(open(AUD/'FR223_subdivision_repairs.json'))
rows={r['name']:r for r in union['rows']}
corpus=[('FR26_direct',rows['FR26']['exact_decimal_box']),('FR32_leaf0',reps[0]['leaves'][0]['box']),('FR54_leaf0',next(r for r in reps if r['name']=='FR54')['leaves'][0]['box'])]
worker=CAMP/'src/_interval_backend_worker.py'
worker.write_text(r'''from fractions import Fraction as F
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
''')
out={'task':'R7N-006','corpus':[],'policy':'Use exact Fraction-QI by default for small proof objects; use one isolated process with outward 10^-60 constructor rounding for large replay/search. Never initialize the wrapper twice in one process.'}
for name,box in corpus:
 rec={'name':name,'box':box}
 for mode in ['exact','bounded']:
  cp=subprocess.run([sys.executable,str(worker),mode],input=json.dumps(box),text=True,capture_output=True,timeout=120)
  if cp.returncode: raise RuntimeError((name,mode,cp.stderr))
  rec[mode]=json.loads(cp.stdout)
 # containment of bounded result around exact for each scalar output
 rec['bounded_contains_exact']={}
 for k in ['P','P1','c2']:
  elo,ehi=map(F,rec['exact'][k]); blo,bhi=map(F,rec['bounded'][k]); rec['bounded_contains_exact'][k]=(blo<=elo<=ehi<=bhi)
 out['corpus'].append(rec)
(CAMP/'results/R7N-006_arithmetic_benchmark.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
