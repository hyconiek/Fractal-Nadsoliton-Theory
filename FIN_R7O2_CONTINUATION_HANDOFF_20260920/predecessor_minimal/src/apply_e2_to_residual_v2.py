from pathlib import Path
from fractions import Fraction as F
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import target_p_e2_bound as e
TAU2=F(67,250)**2
SRC=ROOT/'checkpoints/R7N-020_tight_rounded_trace_only.json'; CP=ROOT/'checkpoints/R7N-020_trace_e2_hybrid_v2.json'
def init():
 d=json.load(open(SRC));out={'task':'R7N-020-trace-e2-hybrid-v2','threshold':'67/250','root_hull':d['root_hull'],'splits':d['splits'],'trace_terminals':d['terminals'],'e2_terminals':[],'e2_failed':[],'pending':d['unresolved'],'unique_processed_e2':0,'source_leaf_count':d['leaf_count']};CP.write_text(json.dumps(out,indent=2)+'\n')
def run(n=1000):
 if not CP.exists():init()
 d=json.load(open(CP)); batch=d['pending'][:n]; d['pending']=d['pending'][n:];safe=[];failed=[];st=time.time()
 for x in batch:
  cell=tuple(tuple(map(F,p)) for p in x['cell']);u=e.e2_upper(cell)
  if u<=TAU2:safe.append({'path':x['path'],'cell':x['cell'],'reason':'SAFE_BY_E2','e2_upper':str(u),'trace_upper':x['trace_upper'],'local_box':x.get('local_box')})
  else:
   y=dict(x);y['e2_upper']=str(u);failed.append(y)
 d['e2_terminals'].extend(safe);d['e2_failed'].extend(failed);d['unique_processed_e2']+=len(batch)
 d['stats']={'SAFE_BY_TRACE':len(d['trace_terminals']),'SAFE_BY_E2':len(d['e2_terminals']),'E2_FAILED':len(d['e2_failed']),'PENDING':len(d['pending'])}
 d['complete_e2_pass']=not d['pending'];d['last_chunk_seconds']=time.time()-st
 CP.write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'unique_processed':d['unique_processed_e2'],'batch':len(batch),'new_e2_safe':len(safe),'new_failed':len(failed),'pending':len(d['pending']),'seconds':d['last_chunk_seconds']},indent=2))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 1000)
