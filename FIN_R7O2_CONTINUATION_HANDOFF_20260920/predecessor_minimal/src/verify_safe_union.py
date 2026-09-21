from pathlib import Path
from fractions import Fraction as F
import json,copy
ROOT=Path(__file__).resolve().parents[1]
reg=json.load(open(ROOT/'results/safe_union_v1.json'))
axes=['x','u','v','e']
def boxF(b):return [[F(a),F(c)] for a,c in [b[k] for k in axes]]
def covers_parent(parent,leaves):
 # These repair trees split exactly one axis and leave other axes equal to parent.
 P=boxF(parent); varying=[]
 for ax in range(4):
  if any(boxF(l['box'])[ax]!=P[ax] for l in leaves): varying.append(ax)
 if len(varying)!=1:return False,'not_one_axis'
 ax=varying[0]
 for l in leaves:
  B=boxF(l['box'])
  for j in range(4):
   if j!=ax and B[j]!=P[j]:return False,'cross_axis_change'
 ints=sorted([boxF(l['box'])[ax] for l in leaves],key=lambda z:z[0])
 if ints[0][0]!=P[ax][0] or ints[-1][1]!=P[ax][1]:return False,'outer_boundary_gap'
 for a,b in zip(ints,ints[1:]):
  if a[1]!=b[0]:return False,'internal_gap_or_overlap_mismatch'
 return True,'exact_closed_partition'
# parents sourced from certificate files
results={}
for parent in sorted({x['parent'] for x in reg['repair_leaves']},key=lambda x:int(x[2:])):
 cert=json.load(open(ROOT/'certificates'/f'{parent}_fresh_partition.json'))
 leaves=[x for x in reg['repair_leaves'] if x['parent']==parent]
 ok,why=covers_parent({'x':cert['parent_box'][0],'u':cert['parent_box'][1],'v':cert['parent_box'][2],'e':cert['parent_box'][3]},leaves);results[parent]={'ok':ok,'reason':why,'leaf_count':len(leaves)}
 assert ok
# mutation 1 delete first FR54 leaf -> must fail
bad=[x for x in reg['repair_leaves'] if not (x['parent']=='FR54' and x['leaf_index']==0)]
cert=json.load(open(ROOT/'certificates/FR54_fresh_partition.json'));ok_del,_=covers_parent({'x':cert['parent_box'][0],'u':cert['parent_box'][1],'v':cert['parent_box'][2],'e':cert['parent_box'][3]},[x for x in bad if x['parent']=='FR54'])
assert not ok_del
# mutation 2 shift a shared FR32 boundary by rational 1e-12 -> must fail
ls=copy.deepcopy([x for x in reg['repair_leaves'] if x['parent']=='FR32']); old=F(ls[1]['box']['x'][0]);ls[1]['box']['x'][0]=str(old+F(1,10**12));cert=json.load(open(ROOT/'certificates/FR32_fresh_partition.json'));ok_shift,_=covers_parent({'x':cert['parent_box'][0],'u':cert['parent_box'][1],'v':cert['parent_box'][2],'e':cert['parent_box'][3]},ls);assert not ok_shift
out={'task':'R7N-004','all_repairs_exact_closed_partitions':True,'parents':results,'negative_controls':{'delete_leaf_rejected':True,'shift_boundary_rejected':True},'note':'Volume equality was not used as the acceptance condition.'}
json.dump(out,open(ROOT/'results/R7N-004_coverage_tests.json','w'),indent=2);print(json.dumps(out,indent=2))
