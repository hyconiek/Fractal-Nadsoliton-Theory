from pathlib import Path
import sys
sys.set_int_max_str_digits(0)
from fractions import Fraction as F
import json,time
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916'
sys.path.insert(0,str(H/'src'));sys.path.insert(0,str(ROOT/'src'))
import frontier_shifted_boxes as old
import generic_threshold_shifted as new
ledger=json.load(open(H/'results/FR223_ACTIVE_MASK_LEDGER.json'))
by={m['name']:m for m in ledger['shifted_masks']}
def box(name):
 m=by[name];return [(F(str(m['x'][0])),F(str(m['x'][1]))),(F(str(m['u'][0])),F(str(m['u'][1]))),(F(str(m['v'][0])),F(str(m['v'][1]))),(F(0),F(str(m['e'])))]
def tup(r): return (r['status'],r['reason'],str(r['P'].lo),str(r['P'].hi),str(r['P1'].lo),str(r['P1'].hi),str(r['c2'].lo),str(r['c2'].hi))
selected=sys.argv[1:]
check_default=('--no-default' not in selected)
selected=[x for x in selected if not x.startswith('--')] or ['FR32','FR48','FR52','FR54','FR56','FR58','FR60']
if check_default:
 a=old.raw_shifted_box(box('FR26')); b=new.raw_shifted_box(box('FR26')); assert tup(a)==tup(b)
outpath=ROOT/'results/R7N-012_014_generic_threshold_tests.json'
if outpath.exists(): out=json.load(open(outpath))
else: out={'default_sigma_matches_historical_FR26':False,'tau0':'67/250','whole_mask_tau0':{}}
if check_default: out['default_sigma_matches_historical_FR26']=True
for name in selected:
 t=time.time();r=new.raw_shifted_box(box(name),F(67,250));out['whole_mask_tau0'][name]={'status':r['status'],'reason':r['reason'],'runtime_s':time.time()-t,'P_hi':str(r['P'].hi),'P1_lo':str(r['P1'].lo),'c2_lo':str(r['c2'].lo)}
 print(name,out['whole_mask_tau0'][name],flush=True)
json.dump(out,open(ROOT/'results/R7N-012_014_generic_threshold_tests.json','w'),indent=2)
