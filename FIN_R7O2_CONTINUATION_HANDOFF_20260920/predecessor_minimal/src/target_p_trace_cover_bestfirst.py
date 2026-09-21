from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,math,time,heapq,itertools
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import target_p_trace_cover as base
TAU=base.TAU

def log_volume(cell):
    return sum(math.log(float(hi/lo)) for lo,hi in cell)

def classify(cell):
    tr,e=base.trace_upper_and_e(cell); B=base.local_box(cell,e); hit=base.mask_hit(B)
    if hit:return 'SAFE_BY_MASK',{'dependency':hit,'trace_upper':str(tr),'local_box':[[str(a),str(b)] for a,b in B]}
    if tr<=2*TAU:return 'SAFE_BY_TRACE',{'trace_upper':str(tr),'local_box':[[str(a),str(b)] for a,b in B]}
    return 'UNRESOLVED',{'trace_upper':str(tr),'local_box':[[str(a),str(b)] for a,b in B]}

def run(max_leaves=10000):
    start=time.monotonic(); counter=itertools.count(); root=base.HULL
    # frontier max-heap by log-volume, only unresolved cells live here.
    heap=[]; terminals=[]; splits=[]; checks=0
    reason,meta=classify(root);checks+=1
    if reason!='UNRESOLVED':terminals.append({'path':'','cell':[[str(a),str(b)] for a,b in root],'reason':reason,**meta})
    else: heapq.heappush(heap,(-log_volume(root),next(counter),'',root,meta))
    leaf_count=len(terminals)+len(heap)
    while heap and leaf_count<max_leaves:
        _,_,path,cell,meta=heapq.heappop(heap)
        ax=base.choose_axis(cell); left,right,m=base.split(cell,ax)
        splits.append({'path':path,'axis':ax,'split':str(m)})
        # one unresolved leaf replaced by two => leaf_count +1, regardless child statuses
        leaf_count+=1
        for side,ch in [('L',left),('R',right)]:
            p=path+str(ax)+side; reason,mta=classify(ch);checks+=1
            if reason=='UNRESOLVED': heapq.heappush(heap,(-log_volume(ch),next(counter),p,ch,mta))
            else: terminals.append({'path':p,'cell':[[str(a),str(b)] for a,b in ch],'reason':reason,**mta})
    unresolved=[]
    while heap:
        _,_,path,cell,meta=heapq.heappop(heap);unresolved.append({'path':path,'cell':[[str(a),str(b)] for a,b in cell],'reason':'UNRESOLVED',**meta})
    stats={'SAFE_BY_MASK':sum(x['reason']=='SAFE_BY_MASK' for x in terminals),'SAFE_BY_TRACE':sum(x['reason']=='SAFE_BY_TRACE' for x in terminals),'UNRESOLVED':len(unresolved)}
    out={'task':'R7N-020-best-first-trace-cover','threshold':'67/250','root_hull':[[str(a),str(b)] for a,b in root],
         'max_leaf_budget':max_leaves,'leaf_count':len(terminals)+len(unresolved),'classifier_calls':checks,'stats':stats,'complete':not unresolved,
         'elapsed_seconds':time.monotonic()-start,'split_policy':'best-first largest log-volume unresolved cell; within cell split axis with largest log ratio at exact rationalized geometric midpoint',
         'splits':splits,'terminals':terminals,'unresolved':unresolved}
    (ROOT/'checkpoints/R7N-020_bestfirst_trace_cover.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ['splits','terminals','unresolved']},indent=2))
    if unresolved:
      print('residual log-volume range',min(log_volume(tuple(tuple(map(F,p)) for p in x['cell'])) for x in unresolved),max(log_volume(tuple(tuple(map(F,p)) for p in x['cell'])) for x in unresolved))
      print('residual trace upper min/max',min(float(F(x['trace_upper'])) for x in unresolved),max(float(F(x['trace_upper'])) for x in unresolved))
    return out
if __name__=='__main__':run()
