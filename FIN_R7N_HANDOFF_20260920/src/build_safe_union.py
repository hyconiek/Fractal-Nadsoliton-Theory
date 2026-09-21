from pathlib import Path
from fractions import Fraction as F
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1];H=ROOT/'inputs/FR223_20260916'
ledger_path=H/'results/FR223_ACTIVE_MASK_LEDGER.json';D=json.load(open(ledger_path))
failed={'FR32','FR48','FR52','FR54','FR56','FR58','FR60'}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def centered_box(m):
 return {'x':[str(F(str(-m['rx']))),str(F(str(m['rx'])))], 'u':['0',str(F(str(m['ru'])))], 'v':['0',str(F(str(m['rv'])))], 'e':['0',str(F(str(m['e'])))]}
def shifted_box(m):
 return {'x':[str(F(str(m['x'][0]))),str(F(str(m['x'][1])))], 'u':[str(F(str(m['u'][0]))),str(F(str(m['u'][1])))], 'v':[str(F(str(m['v'][0]))),str(F(str(m['v'][1])))], 'e':['0',str(F(str(m['e'])))]}
direct=[]
for m in D['centered_masks']:
 direct.append({'name':m['name'],'kind':'centered','box':centered_box(m),'threshold':'sigma','checker':'inputs/FR223_20260916/src/frontier_local_boxes.py', 'checker_sha256':sha(H/'src/frontier_local_boxes.py'),'proof_status':'HISTORICAL_FR223_DIRECT_DOMAIN; not fully fresh-replayed in R7N','analytic_lemma':'R7P-044 local characteristic criterion'})
for m in D['shifted_masks']:
 if m['name'] in failed: continue
 direct.append({'name':m['name'],'kind':'shifted','box':shifted_box(m),'threshold':'sigma','checker':'inputs/FR223_20260916/src/frontier_shifted_boxes.py','checker_sha256':sha(H/'src/frontier_shifted_boxes.py'),'proof_status':'HISTORICAL_FR223_DIRECT_DOMAIN; not fully fresh-replayed in R7N','analytic_lemma':'R7P-044 shifted characteristic criterion'})
repairs=[]
for name in sorted(failed,key=lambda x:int(x[2:])):
 p=ROOT/'certificates'/f'{name}_fresh_partition.json';r=json.load(open(p));assert r['all_certified']
 for i,L in enumerate(r['leaves']):
  repairs.append({'parent':name,'leaf_index':i,'box':{'x':L['box'][0],'u':L['box'][1],'v':L['box'][2],'e':L['box'][3]},'threshold':'sigma','checker':r['checker'],'checker_sha256':r['checker_sha256'],'proof_status':'FRESH_R7N_INTERVAL_CERTIFIED','reason':L['reason'],'analytic_lemma':r['analytic_lemma'],'certificate':str(p.relative_to(ROOT))})
out={'id':'R7N-safe-union-v1-alternative','threshold':'sigma','navigation_buffer_used':False,
 'source_ledger':'inputs/FR223_20260916/results/FR223_ACTIVE_MASK_LEDGER.json','source_ledger_sha256':sha(ledger_path),
 'historical_expected':{'geometric_domains':106,'direct_domains':99,'repair_leaves':18,'source':'2026-09-19 plan; original consolidated repair trees unavailable'},
 'this_reconstruction':{'geometric_domains':106,'direct_domains':len(direct),'fresh_repair_leaves':len(repairs)},
 'direct_domains':direct,'repair_leaves':repairs,
 'scope_note':'This is a fresh alternative reconstruction of the seven failed whole boxes using 17 certified leaves. It must not be called the missing historical 18-leaf consolidated registry.'}
assert len(direct)==99 and len(repairs)==17
json.dump(out,open(ROOT/'results/safe_union_v1.json','w'),indent=2)
print(json.dumps({'direct':len(direct),'repair_leaves':len(repairs),'navigation_buffer_used':False},indent=2))
