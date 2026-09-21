from pathlib import Path
from fractions import Fraction as F
import sys,json,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); sys.path.insert(0,str(ROOT/'src'))
import target_p_trace_cover as tc
import target_p_trace_tight_rounded as tr
import target_p_e2_bound as e2
import compression_centered_moment as cm
import compression_interval_probe as cp
TAU=F(67,250); TAU2=TAU*TAU
OUT=ROOT/'checkpoints/R7N-021_t_refine_once_v2.json'; AXIS=2

def parse(cell): return tuple(tuple(map(F,p)) for p in cell)
def enc(c): return [[str(a),str(b)] for a,b in c]
def classify(c):
    t,_=tr.trace_e(c)
    if t<=2*TAU: return {'reason':'SAFE_BY_TRACE','trace_upper':str(t)}
    u=e2.e2_upper(c)
    if u<=TAU2: return {'reason':'SAFE_BY_E2','trace_upper':str(t),'e2_upper':str(u)}
    q=cm.certify(c)
    if q['ok']:
        return {'reason':'SAFE_BY_CENTERED_MOMENT','trace_upper':str(t),'e2_upper':str(u),'centered_reason':q['reason'],
                'basis_den':q['basis_den'],'basis_num':q['basis_num'],'rank_rows':q['rank_rows'],'rank_det':q['rank_det'],
                'center_c':q['center_c'],'d1':q['d1'],'d2':q['d2'],'d3':q['d3'],'gersh_lower':q['gersh_lower']}
    r=cp.certify(c)
    if r['ok']:
        return {'reason':'SAFE_BY_COMPRESSION','trace_upper':str(t),'e2_upper':str(u),
                'compression_reason':r['reason'],'basis_den':r['basis_den'],'basis_num':r['basis_num'],
                'rank_rows':r['rank_rows'],'rank_det':r['rank_det'],'d1':r['d1'],'d2':r['d2'],'d3':r['d3'],
                'gersh_lower':r['gersh_lower']}
    return {'reason':'UNRESOLVED','trace_upper':str(t),'e2_upper':str(u),
            'centered_reason':q['reason'],'compression_reason':r['reason'],'center_eigs':r['center_eigs'],'rank_det':r['rank_det'],
            'd1':r['d1'],'d2':r['d2'],'d3':r['d3'],'gersh_lower':r['gersh_lower']}

def run(n=150):
    d=json.load(open(OUT)); batch=d['pending_parents'][:n]; d['pending_parents']=d['pending_parents'][n:]
    safe=[]; failed=[]; st=time.time()
    for x in batch:
        c=parse(x['cell']); L,R,m=tc.split(c,AXIS)
        for side,ch in [('L',L),('R',R)]:
            r=classify(ch); rec={'parent_path':x['path'],'path':x['path']+str(AXIS)+side,'axis':AXIS,'split':str(m),'cell':enc(ch),**r}
            (safe if r['reason']!='UNRESOLVED' else failed).append(rec)
    d['refined_safe'].extend(safe); d['refined_failed'].extend(failed); d['processed_parents']+=len(batch)
    reasons=['SAFE_BY_TRACE','SAFE_BY_E2','SAFE_BY_CENTERED_MOMENT','SAFE_BY_COMPRESSION']
    d['stats']={'processed_parents':d['processed_parents'],'refined_children_total':len(d['refined_safe'])+len(d['refined_failed']),
                **{r:sum(x['reason']==r for x in d['refined_safe']) for r in reasons},
                'UNRESOLVED':len(d['refined_failed']),'pending_parents':len(d['pending_parents'])}
    d['complete']=not d['pending_parents']; d['last_chunk_seconds']=time.time()-st
    tmp=OUT.with_suffix('.tmp'); tmp.write_text(json.dumps(d,indent=2)+'\n'); tmp.replace(OUT)
    print(json.dumps({'batch_parents':len(batch),'new_safe_children':len(safe),'new_failed_children':len(failed),'seconds':d['last_chunk_seconds'],**d['stats']},indent=2))
if __name__=='__main__': run(int(sys.argv[1]) if len(sys.argv)>1 else 150)
