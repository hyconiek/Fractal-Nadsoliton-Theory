from pathlib import Path
from fractions import Fraction as F
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import compression_interval_probe as c
SRC=ROOT/'checkpoints/R7N-020_trace_e2_hybrid_v2.json'; CP=ROOT/'checkpoints/R7N-020_trace_e2_compression_v1.json'
def init():
 d=json.load(open(SRC)); assert d.get('complete_e2_pass') and not d['pending']
 out={'task':'R7N-020-trace-e2-compression-v1','threshold':'67/250','root_hull':d['root_hull'],
      'splits':d['splits'],'trace_terminals':d['trace_terminals'],'e2_terminals':d['e2_terminals'],
      'compression_terminals':[],'compression_failed':[],'pending':d['e2_failed'],
      'unique_processed_compression':0,'source_leaf_count':d['source_leaf_count'],
      'proof_method':'cell-adaptive rational rank-3 B; exact-rational rank minor; interval pairwise covariance; Sylvester PD of B^T(tau I-M4)B'}
 CP.write_text(json.dumps(out,indent=2)+'\n')
def run(n=250):
 if not CP.exists():init()
 d=json.load(open(CP)); batch=d['pending'][:n];d['pending']=d['pending'][n:];safe=[];failed=[];st=time.time()
 for x in batch:
  cell=tuple(tuple(map(F,p)) for p in x['cell']);r=c.certify(cell)
  rec={'path':x['path'],'cell':x['cell'],'local_box':x.get('local_box'),'trace_upper':x.get('trace_upper'),'e2_upper':x.get('e2_upper'),**r}
  (safe if r['ok'] else failed).append(rec)
 d['compression_terminals'].extend(safe);d['compression_failed'].extend(failed);d['unique_processed_compression']+=len(batch)
 d['stats']={'SAFE_BY_TRACE':len(d['trace_terminals']),'SAFE_BY_E2':len(d['e2_terminals']),'SAFE_BY_COMPRESSION':len(d['compression_terminals']),'COMPRESSION_FAILED':len(d['compression_failed']),'PENDING':len(d['pending'])}
 d['complete_compression_pass']=not d['pending'];d['last_chunk_seconds']=time.time()-st
 CP.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps({'unique_processed':d['unique_processed_compression'],'batch':len(batch),'new_safe':len(safe),'new_failed':len(failed),'pending':len(d['pending']),'seconds':d['last_chunk_seconds']},indent=2))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 250)
