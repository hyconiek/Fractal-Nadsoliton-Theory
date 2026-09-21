from pathlib import Path
from fractions import Fraction as F
import sys,json,time,hashlib
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916';sys.path.insert(0,str(H/'src'))
import frontier_shifted_boxes as fsb
D=json.load(open(H/'results/FR223_ACTIVE_MASK_LEDGER.json'));by={m['name']:m for m in D['shifted_masks']}
checker=H/'src/frontier_shifted_boxes.py'; checker_sha=hashlib.sha256(checker.read_bytes()).hexdigest()
def parent(name):
 m=by[name];return [(F(str(m['x'][0])),F(str(m['x'][1]))),(F(str(m['u'][0])),F(str(m['u'][1]))),(F(str(m['v'][0])),F(str(m['v'][1]))),(F(0),F(str(m['e'])))]
def endpoint(b):return [[str(a),str(c)] for a,c in b]
def make(name):
 b=parent(name)
 specs={
 'FR32':(0,[(F(0),F(1,2)),(F(1,2),F(1))]),
 'FR48':(1,[(F(0),F(1,2)),(F(1,2),F(1))]),
 'FR52':(0,[(F(0),F(1,2)),(F(1,2),F(1))]),
 'FR54':(1,[(F(0),F(1,8)),(F(1,8),F(1,4)),(F(1,4),F(1,2)),(F(1,2),F(3,4)),(F(3,4),F(1))]),
 'FR56':(2,[(F(0),F(1,2)),(F(1,2),F(1))]),
 'FR58':(1,[(F(0),F(1,2)),(F(1,2),F(1))]),
 'FR60':(1,[(F(0),F(1,2)),(F(1,2),F(1))]),}
 ax,fracs=specs[name];a,c=b[ax];ls=[]
 for lo,hi in fracs:
  z=list(b);z[ax]=(a+(c-a)*lo,a+(c-a)*hi);t=time.time();r=fsb.raw_shifted_box(z);dt=time.time()-t
  rec={'fraction':[str(lo),str(hi)],'box':endpoint(z),'status':r['status'],'reason':r['reason'],'runtime_s':dt,
       'P_hi_float':float(r['P'].hi),'P1_lo_float':float(r['P1'].lo),'c2_lo_float':float(r['c2'].lo)}
  print(name,rec['fraction'],rec['status'],rec['reason'],flush=True);ls.append(rec)
 out={'parent':name,'parent_box':endpoint(b),'split_axis':ax,'leaves':ls,'all_certified':all(x['status']=='INTERVAL_CERTIFIED' for x in ls),
      'checker':'inputs/FR223_20260916/src/frontier_shifted_boxes.py::raw_shifted_box','checker_sha256':checker_sha,
      'threshold':'sigma from frozen strict interval provider','analytic_lemma':'R7P-044 / 3x3 Schur characteristic criterion: c2>0 and (P<=0 or P1>=0)'}
 json.dump(out,open(ROOT/'certificates'/f'{name}_fresh_partition.json','w'),indent=2);return out
if __name__=='__main__':make(sys.argv[1])
