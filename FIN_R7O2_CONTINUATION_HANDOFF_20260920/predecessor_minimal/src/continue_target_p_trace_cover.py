from pathlib import Path
from fractions import Fraction as F
import sys,json,math,time,heapq,itertools
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import target_p_trace_cover as base
import target_p_trace_tight_rounded as tight
TAU=base.TAU
CP=ROOT/'checkpoints/R7N-020_tight_rounded_trace_only.json'
def cellof(x):return tuple(tuple(map(F,p)) for p in x['cell'])
def log_volume(cell):return sum(math.log(float(hi/lo)) for lo,hi in cell)
def classify(cell):
 tr,e=tight.trace_e(cell);B=base.local_box(cell,e)
 if tr<=2*TAU:return 'SAFE_BY_TRACE',{'trace_upper':str(tr),'local_box':[[str(a),str(b)] for a,b in B]}
 return 'UNRESOLVED',{'trace_upper':str(tr),'local_box':[[str(a),str(b)] for a,b in B]}
def resume(target):
 d=json.load(open(CP)); terminals=d['terminals'];splits=d['splits'];unres=d['unresolved'];counter=itertools.count();heap=[]
 for x in unres:heapq.heappush(heap,(-log_volume(cellof(x)),next(counter),x['path'],cellof(x),{k:v for k,v in x.items() if k not in ['path','cell','reason']}))
 leaf_count=len(terminals)+len(heap);checks=0;start=time.monotonic()
 while heap and leaf_count<target:
  _,_,path,cell,meta=heapq.heappop(heap);ax=base.choose_axis(cell);left,right,m=base.split(cell,ax);splits.append({'path':path,'axis':ax,'split':str(m)});leaf_count+=1
  for side,ch in [('L',left),('R',right)]:
   p=path+str(ax)+side;reason,mta=classify(ch);checks+=1
   if reason=='UNRESOLVED':heapq.heappush(heap,(-log_volume(ch),next(counter),p,ch,mta))
   else:terminals.append({'path':p,'cell':[[str(a),str(b)] for a,b in ch],'reason':reason,**mta})
 unresolved=[]
 while heap:
  _,_,path,cell,meta=heapq.heappop(heap);unresolved.append({'path':path,'cell':[[str(a),str(b)] for a,b in cell],'reason':'UNRESOLVED',**meta})
 stats={'SAFE_BY_MASK':0,'SAFE_BY_TRACE':len(terminals),'UNRESOLVED':len(unresolved)}
 d.update({'max_leaf_budget':target,'leaf_count':len(terminals)+len(unresolved),'classifier_calls':d.get('classifier_calls',0)+checks,'stats':stats,'complete':not unresolved,'elapsed_seconds_total_chunks':d.get('elapsed_seconds_total_chunks',d.get('elapsed_seconds',0))+time.monotonic()-start,'splits':splits,'terminals':terminals,'unresolved':unresolved})
 CP.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps({'target':target,'leaf_count':d['leaf_count'],'new_classifier_calls':checks,'stats':stats,'chunk_seconds':time.monotonic()-start},indent=2))
 return d
if __name__=='__main__':resume(int(sys.argv[1]))
