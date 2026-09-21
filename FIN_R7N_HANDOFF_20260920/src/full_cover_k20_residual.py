from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import sys,json,heapq,time
ROOT=Path('/mnt/data/fin_rank7_next_campaign'); IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR));import scientific_rechecks as sr
sr.iv.dps=30
sys.path.insert(0,str(ROOT/'src'));import k20_surrogate_eval as ks
FULL=json.load(open(ROOT/'inputs/FR223_20260916/certificates/R7P-092_full_phase_roots.json'))['roots']
COL=json.load(open(ROOT/'results/R7N-041_full_uniqueness_collars.json'))['roots']
TWOPI=2*sr.iv.pi
SRC=ROOT/'checkpoints/R7N-044_K20_residual.json'; OUT=ROOT/'checkpoints/R7N-044_K20_adaptive.json'

def inner_collars():
    out=[]
    for rec,col in zip(FULL,COL):
        rad=F(col['certified_radius']); axes=[]
        for x in rec['phase']:
            c=sr.I(F(str(x))); L=(c-sr.I(rad))/TWOPI; U=(c+sr.I(rad))/TWOPI
            Ll,Lh=sr.bounds(L); Ul,Uh=sr.bounds(U); axes.append((Lh,Ul))
        out.append((rec['quartic_id'],axes,rad))
    return out
COLLARS=inner_collars()

def in_collar(box):
    for rid,axes,rad in COLLARS:
        okall=True
        for (a,b),(lo,hi) in zip(box,axes):
            ok=False
            for n in (-1,0,1):
                if lo+n<=a and b<=hi+n: ok=True; break
            if not ok: okall=False; break
        if okall:return rid
    return None

def classify(box):
    ok,k,iv=ks.classify(box)
    if ok:return 'FULL_SAFE_BY_K20',k,iv
    rid=in_collar(box)
    if rid is not None:return 'FULL_ROOT_COLLAR',rid,None
    return 'UNRESOLVED',None,None

def initial_heap():
    d=json.load(open(SRC)); assert d['complete']
    heap=[];counter=0
    for rec in d['failed_leaves']:
        box=tuple((F(a),F(b)) for a,b in rec['box']);vol=1
        for a,b in box:vol*=b-a
        heapq.heappush(heap,(-vol,counter,box,rec['path']));counter+=1
    return heap,counter,[],[]

def load():
    if not OUT.exists(): return initial_heap(),0
    d=json.load(open(OUT));heap=[];counter=0
    for rec in d['unresolved_leaves']:
        box=tuple((F(a),F(b)) for a,b in rec['box']);vol=1
        for a,b in box:vol*=b-a
        heapq.heappush(heap,(-vol,counter,box,rec['path']));counter+=1
    return (heap,counter,d['safe_leaves'],d['root_leaves']),d['processed_total']

def run(n=2000):
    (heap,counter,safe,roots),processed0=load();st=time.time();processed=0
    while heap and processed<n:
        _,_,box,path=heapq.heappop(heap);processed+=1
        reason,data,iv=classify(box)
        if reason=='FULL_SAFE_BY_K20':
            safe.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'gradient_component':data,'surrogate_interval':iv});continue
        if reason=='FULL_ROOT_COLLAR':
            roots.append({'path':path,'box':[[str(a),str(b)] for a,b in box],'root_id':data});continue
        widths=[b-a for a,b in box];ax=max(range(3),key=lambda i:widths[i]);a,b=box[ax];m=(a+b)/2
        for tag,ab in [('L',(a,m)),('R',(m,b))]:
            nb=list(box);nb[ax]=ab;nb=tuple(nb);vol=1
            for x,y in nb:vol*=y-x
            heapq.heappush(heap,(-vol,counter,nb,path+str(ax)+tag));counter+=1
    unresolved=[{'path':path,'box':[[str(a),str(b)] for a,b in box]} for _,_,box,path in heap]
    out={'task':'R7N-044-K20-adaptive','method':'K20 45-resonance interval gradient + rigorous uniform full-gradient error; certified full-root collars',
         'epsilon':str(ks.EPS),'source_failed_count':json.load(open(SRC))['failed_count'],'processed_total':processed0+processed,'processed_last_chunk':processed,
         'safe_leaf_count':len(safe),'root_collar_leaf_count':len(roots),'unresolved_leaf_count':len(unresolved),'complete':not unresolved,
         'elapsed_seconds_last_chunk':time.time()-st,'safe_leaves':safe,'root_leaves':roots,'unresolved_leaves':unresolved}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if not k.endswith('_leaves')},indent=2))
if __name__=='__main__':run(int(sys.argv[1]) if len(sys.argv)>1 else 2000)
