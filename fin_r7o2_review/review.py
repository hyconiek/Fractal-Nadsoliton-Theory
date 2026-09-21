"""R7O2 integration audit. Archives are immutable; outputs are written here."""
from pathlib import Path
from fractions import Fraction as F
from collections import defaultdict
from itertools import combinations
import argparse, hashlib, importlib.util, json, math, sys, time, types
import multiprocessing as mp

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'FIN_R7O2_CONTINUATION_HANDOFF_20260920'
PACK=SOURCE/'fin_rank7_targetp_physical_full_20260920'
PRE=ROOT/'FIN_R7N_HANDOFF_20260920'

def save(name,d): (HERE/name).write_text(json.dumps(d,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def box(row):return tuple(tuple(map(F,p)) for p in row['cell'])
def volume(b):return math.prod(hi-lo for lo,hi in b)

def partition(parent,rows):
    trie={}
    for row in rows:
        node=trie
        for tag in row['path']:
            assert tag in 'LR' and 'leaf' not in node
            node=node.setdefault(tag,{})
        assert not node;node['leaf']=box(row)
    def visit(node):
        if 'leaf' in node:
            assert len(node)==1;return node['leaf']
        assert set(node)=={'L','R'}
        L,R=visit(node['L']),visit(node['R'])
        axes=[k for k in range(4) if L[k]!=R[k]];assert len(axes)==1
        k=axes[0];assert L[k][0]<L[k][1]==R[k][0]<R[k][1]
        result=list(L);result[k]=(L[k][0],R[k][1]);return tuple(result)
    assert visit(trie)==parent

def registry():
    failures=[];entries=0
    for line in (SOURCE/'MANIFEST.sha256').read_text().splitlines():
        if not line.strip():continue
        want,rel=line.split(maxsplit=1);p=SOURCE/rel.lstrip('*');entries+=1
        if not p.exists() or sha(p)!=want:failures.append(rel)
    assert not failures,failures
    data=load(PACK/'inputs/target_p_residual_5432.json')
    prior=load(PRE/'checkpoints/R7N-021_t_refine_second_cheap_v1.json')
    assert data['refined_failed']==prior['refined_failed'] and data['root_hull']==prior['root_hull']
    original=data['refined_failed'];assert len(original)==5432
    parents={};hashes={}
    idx=load(PACK/'checkpoints/full_cover_chunks/index.json')
    for rec in idx['chunks']:
        p=PACK/'checkpoints/full_cover_chunks'/rec['file'];assert sha(p)==rec['sha256'];hashes[str(p.relative_to(SOURCE))]=sha(p)
        d=load(p);s=defaultdict(list);u=defaultdict(list)
        for row in d['safe_terminal_leaves']:s[row['original_index']].append(row)
        for row in d['unresolved_terminal_leaves']:u[row['original_index']].append(row)
        for meta in d['processed_parents']:
            i=meta['index'];assert i not in parents
            parents[i]=dict(safe=s[i],unresolved=u[i],source=str(p.relative_to(SOURCE)))
    for name,key in [('stream_v2','stream_rows'),('depth5_stream','rows')]:
        index=load(PACK/'checkpoints'/name/'index.json')
        for rec in index[key]:
            p=PACK/'checkpoints'/name/'parents'/rec['file'];assert sha(p)==rec['sha256'];hashes[str(p.relative_to(SOURCE))]=sha(p)
            d=load(p);i=d['index'];assert i==rec['index'] and i not in parents
            assert d['parent_cell']==original[i]['cell']
            parents[i]=dict(safe=d['safe_terminal_leaves'],unresolved=d['unresolved_terminal_leaves'],source=str(p.relative_to(SOURCE)))
    for repair in idx['repairs']:
        p=(PACK/'checkpoints/full_cover_chunks'/repair['file']).resolve();assert sha(p)==repair['sha256'];hashes[str(p.relative_to(SOURCE))]=sha(p)
        d=load(p);assert not d['unresolved_terminal_leaves']
        for i in repair['parent_indices']:
            assert i in parents
            parents[i]=dict(safe=[x for x in d['safe_terminal_leaves'] if x['original_index']==i],unresolved=[],source=str(p.relative_to(SOURCE)))
    safe=[];residual=[];closed=[];repair=[]
    for i,parent in sorted(parents.items()):
        for row in parent['safe']+parent['unresolved']:assert row['original_index']==i
        partition(box(original[i]),parent['safe']+parent['unresolved'])
        assert all(x['ok'] for x in parent['safe']) and all(not x['ok'] for x in parent['unresolved'])
        safe.extend(parent['safe']);residual.extend(parent['unresolved'])
        (repair if parent['unresolved'] else closed).append(i)
    missing=sorted(set(range(5432))-set(parents));hull=tuple(tuple(map(F,p)) for p in data['root_hull']);hv=volume(hull)
    pv=lambda ids:sum((volume(box(original[i])) for i in ids),F(0))/hv
    rv=sum((volume(box(r)) for r in residual),F(0))/hv
    out=dict(manifest_entries=entries,manifest_pass=True,source_parent_list_matches=True,
             recorded_processed=len(parents),claimed_closed=len(closed),repair_parents=repair,
             missing=missing,safe_leaf_count=len(safe),unresolved_leaf_count=len(residual),
             closed_parent_fraction=str(pv(closed)),repair_parent_fraction=str(pv(repair)),
             unprocessed_fraction=str(pv(missing)),unresolved_terminal_fraction=str(rv),
             geometry_checked=True,source_sha256=hashes)
    save('registry.json',out);save('safe_leaves.json',safe)
    print(json.dumps({k:v for k,v in out.items() if k not in ['source_sha256','missing','repair_parents']},indent=2))

def backend(fast=True):
    h=PRE/'inputs/FR223_20260916/src'
    sys.path[:0]=[str(PRE/'src'),str(h),str(ROOT)]
    if fast:
        import intervals
        from fin_r7o2_review.intervals_fast import Interval
        intervals.QI=Interval
    else:
        from fin_rank7_intake_review.verify_continuation import bounded_rationals
        bounded_rationals(12)
    def relocated(name,path):
        code=path.read_text().replace('/mnt/data/fin_rank7_next_campaign',str(PRE)).replace('/mnt/data/r7n_full_restore/FIN_R7N_HANDOFF_20260920',str(PRE))
        m=types.ModuleType(name);m.__file__=str(path);sys.modules[name]=m
        exec(compile(code,str(path),'exec'),m.__dict__);return m
    relocated('target_p_trace_cover',PRE/'src/target_p_trace_cover.py')
    old=relocated('compression_interval_probe',PRE/'src/compression_interval_probe.py')
    pc=relocated('physical_centered_moment',PACK/'src/physical_centered_moment.py')
    return pc,old

def determinant(A):return A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])

def certify_fixed(row,pc,old,center=None):
    import numpy as np
    from intervals import QI
    cbox=box(row);den=int(row['basis_den']);B=[[F(int(x),den) for x in r] for r in row['basis_num']]
    minors=[determinant([B[i] for i in ids]) for ids in combinations(range(4),3)];assert any(minors)
    # Missing original c is replaced by a newly proposed, frozen rational c.
    if center is None:
        mid=np.array([float((a+b)/2) for a,b in cbox]);p=pc.off_face.p_from_aligned_compact(*mid)
        pa=[p[0],p[4]+p[8],p[6],p[2]+p[10],p[3]+p[9],p[5]+p[7],p[1]+p[11]]
        c=[F(format(sum(pa[i]*sum(float(B[r][k])*float((old.OBS[i][r].lo+old.OBS[i][r].hi)/2) for r in range(4)) for i in range(7)),'.16g')) for k in range(3)]
    else:
        c=list(map(F,center));assert len(c)==3
    bounds=pc._internal_bounds(cbox);center,rads=pc._midpoint_bounds(bounds)
    full=pc._moment_jets(bounds,B,c);cen=pc._moment_jets(center,B,c)
    E=[[pc._range_from_jets(full[i][j],cen[i][j],rads) for j in range(3)] for i in range(3)]
    K=[[QI(F(67,250)*sum(B[r][i]*B[r][j] for r in range(4)))-E[i][j] for j in range(3)] for i in range(3)]
    sylv,gersh,d1,d2,d3,g=pc._pd_result(K)
    return dict(original_index=row['original_index'],path=row['path'],ok=bool(sylv or gersh),
                center_c=list(map(str,c)),basis_num=row['basis_num'],basis_den=den,
                rank_minor=str(next(x for x in minors if x)),
                pd_bounds=[[str(v.lo),str(v.hi)] for v in [d1,d2,d3]],gersh_lower=str(min(g)),
                proof='Replacement fixed-B/new-rational-c certificate; outward interval arithmetic; no recorded float sign used.')

def worker_init():
    global WORKER_PC,WORKER_OLD
    WORKER_PC,WORKER_OLD=backend(True)

def worker(row):return certify_fixed(row,WORKER_PC,WORKER_OLD,row.get('_audit_center'))

def replay(limit=None,fast=True,workers=1):
    rows=load(HERE/'safe_leaves.json');out=[];start=time.monotonic();pool=None
    assert limit is None or 0<limit<=len(rows)
    output='leaf_replay.json' if limit is None else ('sample_binary_replay.json' if fast else 'sample_rational_replay.json')
    prior=HERE/'leaf_replay.json'
    if prior.exists():
        saved=load(prior)
        if saved.get('complete'):
            assert len(saved['certificates'])==len(rows)
            for row,c in zip(rows,saved['certificates']):
                assert (row['original_index'],row['path'])==(c['original_index'],c['path'])
                row['_audit_center']=c['center_c']
    if workers>1:
        assert fast and workers<=4
        pool=mp.Pool(workers,initializer=worker_init)
        iterator=pool.imap(worker,rows[:limit],chunksize=5)
    else:
        pc,old=backend(fast);iterator=(certify_fixed(row,pc,old,row.get('_audit_center')) for row in rows[:limit])
    for i,result in enumerate(iterator):
        out.append(result)
        if (i+1)%100==0 or i+1==len(rows[:limit]):
            save(output,dict(processed=i+1,total=len(rows),complete=i+1==len(rows),
                 arithmetic='Outward binary64 intervals, exact rational exposed endpoints' if fast else 'Outward rational 10^-12 intervals',
                 failed=[j for j,x in enumerate(out) if not x['ok']],seconds=time.monotonic()-start,certificates=out))
            print(i+1,'/',len(rows),'failures',sum(not x['ok'] for x in out),'seconds',round(time.monotonic()-start,1),flush=True)
    if pool is not None:pool.close();pool.join()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['registry','replay']);p.add_argument('--limit',type=int);p.add_argument('--rational',action='store_true');p.add_argument('--workers',type=int,default=1);a=p.parse_args()
    registry() if a.stage=='registry' else replay(a.limit,not a.rational,a.workers)
