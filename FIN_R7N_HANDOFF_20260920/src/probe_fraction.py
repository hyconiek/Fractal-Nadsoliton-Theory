from pathlib import Path
from fractions import Fraction as F
import sys,json,time
sys.set_int_max_str_digits(0)
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916';sys.path.insert(0,str(H/'src'))
import frontier_shifted_boxes as fsb
D=json.load(open(H/'results/FR223_ACTIVE_MASK_LEDGER.json'));by={m['name']:m for m in D['shifted_masks']}
def box(name):
 m=by[name];return [(F(str(m['x'][0])),F(str(m['x'][1]))),(F(str(m['u'][0])),F(str(m['u'][1]))),(F(str(m['v'][0])),F(str(m['v'][1]))),(F(0),F(str(m['e'])))]
name=sys.argv[1];ax=int(sys.argv[2]);fractions=[(F(sys.argv[i]),F(sys.argv[i+1])) for i in range(3,len(sys.argv),2)];b=box(name);a,c=b[ax]
for lo,hi in fractions:
 z=list(b);z[ax]=(a+(c-a)*lo,a+(c-a)*hi);t=time.time();r=fsb.raw_shifted_box(z);print(str(lo),str(hi),r['status'],r['reason'],'P_hi',float(r['P'].hi),'P1_lo',float(r['P1'].lo),'secs',round(time.time()-t,3),flush=True)
