from pathlib import Path
from fractions import Fraction as F
import json,hashlib
ROOT=Path(__file__).resolve().parents[1]
HULL=((F(1,900),F(1)),(F(1,128),F(1)),(F(1,9),F(1)),(F(1,1000000),F(1)))
def C(cell): return tuple((F(a),F(b)) for a,b in cell)
def toks(path):
 if len(path)%2: raise ValueError('odd path')
 return [(int(path[i]),path[i+1]) for i in range(0,len(path),2)]
def build(rows):
 root={}
 for row in rows:
  node=root
  for tok in toks(row['path']): node=node.setdefault(tok,{})
  if '_leaf' in node or any(k!='_leaf' for k in node):
   if '_leaf' in node: raise ValueError('duplicate path')
  node['_leaf']=row
 return root
def recur(node,prefix=''):
 leaf=node.get('_leaf')
 kids=[k for k in node if k!='_leaf']
 if leaf is not None:
  if kids: raise ValueError(f'leaf has descendants {prefix}')
  return C(leaf['cell']),[leaf]
 if not kids: raise ValueError(f'empty node {prefix}')
 axes={a for a,s in kids}; sides={s for a,s in kids}
 if len(axes)!=1 or sides!={'L','R'} or len(kids)!=2: raise ValueError(f'bad children {prefix} {kids}')
 ax=next(iter(axes)); L,lr=recur(node[(ax,'L')],prefix+str(ax)+'L'); R,rr=recur(node[(ax,'R')],prefix+str(ax)+'R')
 for k in range(4):
  if k==ax:
   if L[k][1]!=R[k][0]: raise ValueError(f'boundary mismatch {prefix} axis{ax}')
  elif L[k]!=R[k]: raise ValueError(f'cross-axis mismatch {prefix} dim{k}')
 P=list(L); P[ax]=(L[ax][0],R[ax][1]); return tuple(P),lr+rr
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check():
 rows=[json.loads(x) for x in (ROOT/'inherited/r7n_final_compact_partition.jsonl').read_text().splitlines() if x.strip()]
 box,leaves=recur(build(rows))
 reg=json.load(open(ROOT/'parent_registry.json'))
 residual=sorted((r for r in rows if not r['safe']),key=lambda x:x['path'])
 rp=sorted(reg['parents'],key=lambda x:x['original_path'])
 iderrs=[]
 if len(residual)!=len(rp): iderrs.append(['count',len(residual),len(rp)])
 else:
  for a,b in zip(residual,rp):
   if a['path']!=b['original_path'] or C(a['cell'])!=C(b['original_cell']): iderrs.append(['parent',a['path'],b['original_path']]); break
 expected={
 'proofs/R7N-010_residual_hull.md':'1fb5382da21501b33066edf5cb711cc0814fd5bfd46ac90872866b3b7075fbfe',
 'inputs/FR223_20260916/results/FR1_residual_tail_upgrade.json':'3ee61056126b3fc269dea43c902e086597be49cabea5ab799834bd8673814789',
 'inputs/FR223_20260916/results/FR42_global_large_J6_tail_5x.json':'c92ff8ac458b35f039ca41bf680fa3ba4625f4b7ce96b10adf4fa9ed2f916010'}
 tails=[]
 for rel,want in expected.items():
  got=sha(ROOT/'inherited'/rel);tails.append({'path':rel,'want':want,'got':got,'ok':got==want})
 counts={}
 for r in rows:counts[r['source']]=counts.get(r['source'],0)+1
 return {'leaf_count':len(rows),'safe_count':sum(r['safe'] for r in rows),'residual_count':sum(not r['safe'] for r in rows),'counts':counts,'root_box':[[str(a),str(b)] for a,b in box],'exact_tree_pass':box==HULL,'r7o3_parent_identity_pass':not iderrs,'identity_errors':iderrs,'tail_hashes':tails,'tails_pass':all(x['ok'] for x in tails),'global_join_pass':box==HULL and not iderrs and all(x['ok'] for x in tails)}
if __name__=='__main__': print(json.dumps(check(),indent=2))
