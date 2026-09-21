from pathlib import Path
from fractions import Fraction as F
import sys,json,heapq,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign');sys.path.insert(0,str(ROOT/'src'))
import quartic_cover_normalized_fast as q
CP=ROOT/'checkpoints/R7N-037_quartic_cover_normalized.json'
def run(extra=10000):
 d=json.load(open(CP)); heap=[];counter=0
 for x in d['unresolved_leaves']:
  box=tuple(tuple(map(F,p)) for p in x['box']);vol=(box[0][1]-box[0][0])*(box[1][1]-box[1][0])*(box[2][1]-box[2][0]);heapq.heappush(heap,(-vol,counter,box,x['path']));counter+=1
 safe=d['safe_leaves'];roots=d['root_leaves'];done=0;st=time.time()
 while heap and done<extra:
  _,_,box,path=heapq.heappop(heap);done+=1;reason,data,iv=q.classify(box)
  if reason=='GRADIENT':safe.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'gradient_component':data,'interval':iv});continue
  if reason=='ROOT_COLLAR':roots.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'root_id':data});continue
  widths=[b-a for a,b in box];ax=max(range(3),key=lambda i: widths[i]);a,b=box[ax];m=(a+b)/2
  for tag,ab in [('L',(a,m)),('R',(m,b))]:
   nb=list(box);nb[ax]=ab;nb=tuple(nb);vol=(nb[0][1]-nb[0][0])*(nb[1][1]-nb[1][0])*(nb[2][1]-nb[2][0]);heapq.heappush(heap,(-vol,counter,nb,path+str(ax)+tag));counter+=1
 unresolved=[{'path':p,'box':[[str(a),str(b)] for a,b in b]} for _,_,b,p in heap]
 d.update({'budget':d.get('budget',d.get('processed',0))+extra,'processed':d.get('processed',0)+done,'safe_leaf_count':len(safe),'root_collar_leaf_count':len(roots),'unresolved_leaf_count':len(unresolved),'complete':not unresolved,'elapsed_seconds_last_chunk':time.time()-st,'safe_leaves':safe,'root_leaves':roots,'unresolved_leaves':unresolved})
 d['elapsed_seconds_total']=d.get('elapsed_seconds_total',d.get('elapsed_seconds',0))+d['elapsed_seconds_last_chunk'];CP.write_text(json.dumps(d,indent=2)+'\n')
 print(json.dumps({'processed_total':d['processed'],'chunk':done,'safe':len(safe),'root_collars':len(roots),'unresolved':len(unresolved),'seconds':d['elapsed_seconds_last_chunk']},indent=2))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 10000)
