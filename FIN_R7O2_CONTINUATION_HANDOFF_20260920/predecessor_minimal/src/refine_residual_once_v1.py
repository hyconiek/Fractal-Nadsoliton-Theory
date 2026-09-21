from pathlib import Path
from fractions import Fraction as F
import sys,json,time,math
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); sys.path.insert(0,str(ROOT/'src'))
import target_p_trace_cover as tc
import target_p_trace_tight_rounded as tr
import target_p_e2_bound as e2
import compression_interval_probe as cp
TAU=F(67,250); TAU2=TAU*TAU
SRC=ROOT/'checkpoints/R7N-020_trace_e2_compression_v1.json'
OUT=ROOT/'checkpoints/R7N-021_refine_once_v1.json'

def parse(cell): return tuple(tuple(map(F,p)) for p in cell)
def enc(c): return [[str(a),str(b)] for a,b in c]
def classify(c):
    t,_=tr.trace_e(c)
    if t<=2*TAU: return {'reason':'SAFE_BY_TRACE','trace_upper':str(t)}
    u=e2.e2_upper(c)
    if u<=TAU2: return {'reason':'SAFE_BY_E2','trace_upper':str(t),'e2_upper':str(u)}
    r=cp.certify(c)
    if r['ok']:
        return {'reason':'SAFE_BY_COMPRESSION','trace_upper':str(t),'e2_upper':str(u),
                'compression_reason':r['reason'],'basis_den':r['basis_den'],'basis_num':r['basis_num'],
                'rank_rows':r['rank_rows'],'rank_det':r['rank_det'],'d1':r['d1'],'d2':r['d2'],'d3':r['d3'],
                'gersh_lower':r['gersh_lower']}
    return {'reason':'UNRESOLVED','trace_upper':str(t),'e2_upper':str(u),
            'compression_reason':r['reason'],'center_eigs':r['center_eigs'],'rank_det':r['rank_det'],
            'd1':r['d1'],'d2':r['d2'],'d3':r['d3'],'gersh_lower':r['gersh_lower']}

def init():
    d=json.load(open(SRC)); assert d['complete_compression_pass'] and not d['pending']
    out={'task':'R7N-021-refine-once-v1','threshold':'67/250','root_hull':d['root_hull'],'source_leaf_count':d['source_leaf_count'],
         'inherited_trace_count':len(d['trace_terminals']),'inherited_e2_count':len(d['e2_terminals']),
         'inherited_compression_count':len(d['compression_terminals']),
         'parent_residual_count':len(d['compression_failed']),'pending_parents':d['compression_failed'],
         'refined_safe':[],'refined_failed':[],'processed_parents':0,
         'policy':'one exact binary split per residual parent; axis=max endpoint log-ratio; classify children trace->e2->cell-adaptive compression'}
    OUT.write_text(json.dumps(out,indent=2)+'\n')

def run(n=100):
    if not OUT.exists(): init()
    d=json.load(open(OUT)); batch=d['pending_parents'][:n]; d['pending_parents']=d['pending_parents'][n:]
    safe=[]; failed=[]; st=time.time()
    for x in batch:
        c=parse(x['cell']); ax=tc.choose_axis(c); L,R,m=tc.split(c,ax)
        for side,ch in [('L',L),('R',R)]:
            r=classify(ch); rec={'parent_path':x['path'],'path':x['path']+str(ax)+side,'axis':ax,'split':str(m),'cell':enc(ch),**r}
            (safe if r['reason']!='UNRESOLVED' else failed).append(rec)
    d['refined_safe'].extend(safe); d['refined_failed'].extend(failed); d['processed_parents']+=len(batch)
    d['stats']={'processed_parents':d['processed_parents'],'refined_children_total':len(d['refined_safe'])+len(d['refined_failed']),
                'SAFE_BY_TRACE':sum(x['reason']=='SAFE_BY_TRACE' for x in d['refined_safe']),
                'SAFE_BY_E2':sum(x['reason']=='SAFE_BY_E2' for x in d['refined_safe']),
                'SAFE_BY_COMPRESSION':sum(x['reason']=='SAFE_BY_COMPRESSION' for x in d['refined_safe']),
                'UNRESOLVED':len(d['refined_failed']),'pending_parents':len(d['pending_parents'])}
    d['complete']=not d['pending_parents']; d['last_chunk_seconds']=time.time()-st
    OUT.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'batch_parents':len(batch),'new_safe_children':len(safe),'new_failed_children':len(failed),'seconds':d['last_chunk_seconds'],**d['stats']},indent=2))
if __name__=='__main__': run(int(sys.argv[1]) if len(sys.argv)>1 else 100)
