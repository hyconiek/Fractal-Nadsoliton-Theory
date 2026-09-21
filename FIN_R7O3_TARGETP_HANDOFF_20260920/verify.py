from pathlib import Path
import json,sys,hashlib
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT/'src'))
from verification_core import *
from global_join_checker import check as check_join

def main():
 errs=[]
 try:reg=json.load(open(ROOT/'parent_registry.json'))
 except Exception as ex:return 2
 errs += [('registry',x) for x in validate_registry(reg)]
 certs=load_jsonl(ROOT/'certificates/active_leaf_certificates.jsonl');leafmap={c['leaf_id']:c for c in certs}
 if len(certs)!=12425 or len(leafmap)!=12425:errs.append(('cert_count',len(certs),len(leafmap)))
 for c in certs:
  e=validate_cert(c)
  if e:errs.append(('cert_schema',c.get('leaf_id'),e));break
 trees=load_jsonl(ROOT/'certificates/reconstructed_parent_trees.jsonl')
 if len(trees)!=5432:errs.append(('tree_count',len(trees)))
 trmap={r['original_index']:r for r in trees}
 for p in reg['parents']:
  i=p['original_index'];r=trmap.get(i)
  if r is None:errs.append(('missing_tree',i));break
  e,lids=check_tree_record(r,p,leafmap)
  if e:errs.append(('tree',i,e[:10]));break
  got={leafmap[x]['certificate_sha256'] for x in lids};want=set(p['leaf_certificate_sha256s'])
  if got!=want:errs.append(('tree_leaf_set',i));break
 try:
  replay=json.load(open(ROOT/'results/R7O3-037_clean_directory_math_replay.json'))
  if not replay.get('global_pass') or replay.get('pass_count')!=12425:errs.append(('clean_math_replay',replay.get('pass_count'),replay.get('fail_count')))
 except Exception as ex:errs.append(('clean_math_replay_missing',str(ex)))
 try:
  join=check_join()
  if not join['global_join_pass']:errs.append(('global_join',join))
 except Exception as ex:errs.append(('global_join_exception',str(ex)))
 try:
  mut=json.load(open(ROOT/'results/R7O3-034_hostile_mutations.json'))
  if not mut.get('global_pass') or mut.get('rejected_count')!=mut.get('count'):errs.append(('mutations',mut))
 except Exception as ex:errs.append(('mutations_missing',str(ex)))
 out={'task':'R7O3-final-readonly-verifier','global_pass':not errs,'error_count':len(errs),'errors':errs[:50],'counts':{'parents':len(reg.get('parents',[])),'active_safe_leaves':len(certs),'clean_formula_replay_pass':12425 if not any(x[0]=='clean_math_replay' for x in errs) else None,'r7n_prior_safe_leaves':13231}}
 print(json.dumps(out,indent=2));return 0 if not errs else 1
if __name__=='__main__':raise SystemExit(main())
