from pathlib import Path
import json,hashlib,sys
ROOT=Path(__file__).resolve().parent
IMM=[ROOT/'inputs/FR223_20260916/src/frontier_shifted_boxes.py',ROOT/'inputs/FR223_20260916/results/FR223_ACTIVE_MASK_LEDGER.json',ROOT/'inputs/FIN_Rank7_Next_Campaign_Plan_EN_20260919.md']
def hashes():return {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in IMM}
def check_registry():
 r=json.load(open(ROOT/'results/safe_union_v1.json'));assert r['navigation_buffer_used'] is False;assert len(r['direct_domains'])==99;assert len(r['repair_leaves'])==17
 for p in sorted({x['parent'] for x in r['repair_leaves']}):
  c=json.load(open(ROOT/'certificates'/f'{p}_fresh_partition.json'));assert c['all_certified'];assert all(x['status']=='INTERVAL_CERTIFIED' for x in c['leaves'])
 return 'safe_union_structural_and_fresh_repairs_PASS'
def check_target():
 d=json.load(open(ROOT/'results/R7N-017_target_p_implication.json'));assert d['sigma_lt_tau0'] and d['local_event_below_gain_endpoint'];return 'target_P_logic_PASS'
def check_phase():
 d=json.load(open(ROOT/'results/R7N-033_phase_cover_audit.json'));assert d['old_partition']['safe_leaf_count']==396 and d['old_partition']['unresolved_leaf_count']==1272 and d['old_partition']['root_neighborhood_leaf_count']==0;return 'phase_provenance_PASS'
checks=[check_registry,check_target,check_phase]
before=hashes();a=[f() for f in checks];mid=hashes();b=[f() for f in reversed(checks)];after=hashes();assert before==mid==after
out={'status':'PASS','forward':a,'reverse':b,'immutable_hashes_unchanged':True,'layers':{'schema_saved_evidence':'PASS','fresh_mathematical_recompute':'limited to individual repair certificates run separately; this verifier does not regenerate them','historical_FR223_direct_99':'not fully fresh-replayed'}}
print(json.dumps(out,indent=2))
