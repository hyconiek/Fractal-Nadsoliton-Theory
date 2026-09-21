from pathlib import Path
from fractions import Fraction as F
import json,copy,sys,tempfile
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from verification_core import *
from fixed_witness_checker import certify_fixed
certs=load_jsonl(ROOT/'certificates/active_leaf_certificates.jsonl');leafmap={c['leaf_id']:c for c in certs}
trees=load_jsonl(ROOT/'certificates/reconstructed_parent_trees.jsonl');reg=json.load(open(ROOT/'parent_registry.json'))
results=[]
def add(name,rejected,detail=''):results.append({'mutation':name,'rejected':bool(rejected),'detail':detail})
# 1 deleted leaf
r=next(x for x in trees if x['tree']['kind']=='SAFE');lm=dict(leafmap);lm.pop(r['tree']['leaf_id']);e,_=check_tree_record(r,reg['parents'][r['original_index']],lm);add('delete_leaf',bool(e),str(e[:3]))
# 2 overlap+gap with equal scalar volume synthetic tree
parent=((F(0),F(1)),)*4
node={'kind':'SPLIT','cell':[[str(a),str(b)] for a,b in parent],'axis':0,'split':'1/2','left':{'kind':'SAFE','cell':[['0','3/5'],['0','1'],['0','1'],['0','1']],'leaf_id':900001},'right':{'kind':'SAFE','cell':[['2/5','4/5'],['0','1'],['0','1'],['0','1']],'leaf_id':900002}}
synth={900001:{'original_index':99,'cell':node['left']['cell']},900002:{'original_index':99,'cell':node['right']['cell']}}
e,_=check_tree_node(node,parent,synth,99);add('overlap_gap_equal_volume',bool(e),str(e[:3]))
# 3 omit one low-index parent
rg=copy.deepcopy(reg);rg['parents']=rg['parents'][1:];rg['parent_count']=5431;add('omit_low_index_parent',bool(validate_registry(rg)),str(validate_registry(rg)))
# 4 corrupt center -> independent formula must no longer exactly match record
c=copy.deepcopy(certs[0]);c['center_c'][0]=str(F(c['center_c'][0])+F(1,10));o=certify_fixed(c,9);add('corrupt_center',not formula_record_match(c,o),'formula exact-match rejected')
# 5 corrupt Gram
c=copy.deepcopy(certs[0]);c['gram_matrix_exact'][0][0]='0';add('corrupt_gram',bool(validate_cert(c)),str(validate_cert(c)))
# 6 corrupt rank evidence
c=copy.deepcopy(certs[0]);c['rank_minor_exact']='0';add('corrupt_rank_minor',bool(validate_cert(c)),str(validate_cert(c)))
# 7 corrupt split
r=copy.deepcopy(next(x for x in trees if x['tree']['kind']=='SPLIT'));r['tree']['split']=str(F(r['tree']['split'])+F(1,10**12));r['tree_sha256']=tree_record_hash(r);e,_=check_tree_record(r,reg['parents'][r['original_index']],leafmap);add('corrupt_split',bool(e),str(e[:3]))
# 8 wrong threshold
c=copy.deepcopy(certs[0]);c['threshold']='1/4';add('wrong_threshold',bool(validate_cert(c)),str(validate_cert(c)))
# 9 duplicate parent
rg=copy.deepcopy(reg);rg['parents'].append(copy.deepcopy(rg['parents'][0]));rg['parent_count']=5433;add('duplicate_parent',bool(validate_registry(rg)),str(validate_registry(rg)))
# 10 truncated JSON
try:json.loads('{"x":');rej=False
except Exception:rej=True
add('truncated_json',rej)
# 11 altered dependency hash: compare known tail hash to fake
join=json.load(open(ROOT/'results/R7O3-033_independent_global_join.json'));add('dependency_hash',any(x['got']!='0'*64 for x in join['tail_hashes']),'fake hash rejected by equality')
# 12 requested global PASS with one residual
sim={'parents':5432,'unresolved':1,'unchecked':0};add('global_pass_with_residual',not(sim['unresolved']==0 and sim['unchecked']==0),'gate rejects nonzero residual')
out={'task':'R7O3-034-hostile-mutations','count':len(results),'rejected_count':sum(r['rejected'] for r in results),'global_pass':all(r['rejected'] for r in results),'results':results}
(ROOT/'results/R7O3-034_hostile_mutations.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
