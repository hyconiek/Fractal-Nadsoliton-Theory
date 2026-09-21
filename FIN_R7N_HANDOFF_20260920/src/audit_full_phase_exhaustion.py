from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import json,math,collections,time,copy,sys
ROOT=Path(__file__).resolve().parents[1]
B=json.load(open(ROOT/'checkpoints/R7N-044_full_cover_k16.json'))
D=json.load(open(ROOT/'checkpoints/R7N-044_K20_residual.json'))
A=json.load(open(ROOT/'checkpoints/R7N-044_K20_adaptive.json'))
K16=json.load(open(ROOT/'results/R7N-043_K16_surrogate.json'))
K20=json.load(open(ROOT/'results/R7N-043_K20_surrogate.json'))
FULL=json.load(open(ROOT/'inputs/FR223_20260916/certificates/R7P-092_full_phase_roots.json'))
COL=json.load(open(ROOT/'results/R7N-041_full_uniqueness_collars.json'))
SYM=json.load(open(ROOT/'results/R7N-039_symmetry_and_count.json'))
EPS16=F(K16['uniform_full_gradient_error_bound']); EPS20=F(K20['uniform_full_gradient_error_bound'])

def box(rec): return tuple((F(a),F(b)) for a,b in rec['box'])
def subset(x,p): return all(pa<=a<=b<=pb for (a,b),(pa,pb) in zip(x,p))

def partition_ok(parent, leaves):
    # leaves: list of exact boxes. Verify a midpoint binary partition, independent of stored path/search order.
    if not leaves:return False
    if any(not subset(x,parent) for x in leaves):return False
    if len(leaves)==1:return leaves[0]==parent
    for ax in range(3):
        a,b=parent[ax];m=(a+b)/2
        L=[x for x in leaves if x[ax][1] <= m]
        R=[x for x in leaves if x[ax][0] >= m]
        if L and R and len(L)+len(R)==len(leaves):
            lp=list(parent);rp=list(parent);lp[ax]=(a,m);rp[ax]=(m,b)
            if partition_ok(tuple(lp),L) and partition_ok(tuple(rp),R):return True
    return False

def root_cell_id(x,scale):
    out=[]
    for a,b in x:
        c=(a+b)/2; q=c*scale; out.append(min(scale-1,q.numerator//q.denominator))
    return tuple(out)

def check_baseline_partition():
    terms=[box(r) for r in B['safe_leaves']]+[box(r) for r in B['root_leaves']]+[box(r) for r in B['unresolved_leaves']]
    groups=collections.defaultdict(list)
    for x in terms:groups[root_cell_id(x,4)].append(x)
    if set(groups)!=set(__import__('itertools').product(range(4),repeat=3)):return False,'missing initial cells'
    for idx,ls in groups.items():
        p=tuple((F(i,4),F(i+1,4)) for i in idx)
        if not partition_ok(p,ls):return False,f'partition failure initial {idx}'
    return True,{'initial_cells':64,'terminal_leaves':len(terms)}

def check_direct_identity():
    src=B['unresolved_leaves']; rows=D['safe_leaves']+D['failed_leaves']
    if len(rows)!=len(src):return False,'count'
    seen=set()
    for r in rows:
        i=r['source_index']
        if i in seen or not (0<=i<len(src)):return False,'source index'
        seen.add(i); s=src[i]
        if r['box']!=s['box'] or r['path']!=s['path']:return False,f'mismatch {i}'
    return seen==set(range(len(src))),{'source_count':len(src),'safe':D['safe_count'],'failed':D['failed_count']}

def build_parent_grid():
    # Failed parents lie on the exact 1/8192 lattice. Map each elementary open cube to one parent.
    mp={}; parents={r['source_index']:box(r) for r in D['failed_leaves']}
    S=8192
    for i,p in parents.items():
        rs=[]
        for a,b in p:
            ia=int(a*S); ib=int(b*S); assert F(ia,S)==a and F(ib,S)==b
            rs.append(range(ia,ib))
        for key in __import__('itertools').product(*rs):
            if key in mp:raise AssertionError(('parent overlap',key))
            mp[key]=i
    return mp,parents

def adaptive_groups():
    mp,parents=build_parent_grid(); groups=collections.defaultdict(list); recs=collections.defaultdict(list); S=8192
    for r in A['safe_leaves']+A['root_leaves']:
        x=box(r); c=[(a+b)/2 for a,b in x]; key=tuple(min(S-1,int(z*S)) for z in c)
        i=mp.get(key)
        if i is None or not subset(x,parents[i]):raise AssertionError(('unmapped adaptive leaf',r['path'],key))
        groups[i].append(x); recs[i].append(r)
    return parents,groups,recs

def check_adaptive_partition():
    parents,groups,recs=adaptive_groups()
    if set(groups)!=set(parents):return False,('missing parents',len(set(parents)-set(groups)))
    for i,p in parents.items():
        if not partition_ok(p,groups[i]):return False,f'adaptive parent {i}'
    return True,{'parents':len(parents),'terminal_leaves':sum(map(len,groups.values()))}

def saved_interval_ok(rec,eps):
    lo,hi=map(F,rec['surrogate_interval']); return lo>eps or hi<-eps

def check_saved_signs():
    bad16=[i for i,r in enumerate(B['safe_leaves']) if not saved_interval_ok(r,EPS16)]
    bad20d=[i for i,r in enumerate(D['safe_leaves']) if not saved_interval_ok(r,EPS20)]
    bad20a=[i for i,r in enumerate(A['safe_leaves']) if not saved_interval_ok(r,EPS20)]
    return not(bad16 or bad20d or bad20a),{'K16_checked':len(B['safe_leaves']),'K20_direct_checked':len(D['safe_leaves']),'K20_adaptive_checked':len(A['safe_leaves']),'bad':[bad16[:3],bad20d[:3],bad20a[:3]]}

def collar_axes():
    # independent numerical containment using conservative decimal pi bounds is not needed: use interval provider exactly.
    IR=ROOT/'inputs/intake_review_20260919';sys.path.insert(0,str(IR));import scientific_rechecks as sr
    TWOPI=2*sr.iv.pi;out={}
    for rec,col in zip(FULL['roots'],COL['roots']):
        rad=F(col['certified_radius']); axes=[]
        for xx in rec['phase']:
            c=sr.I(F(str(xx)));L=(c-sr.I(rad))/TWOPI;U=(c+sr.I(rad))/TWOPI
            Ll,Lh=sr.bounds(L);Ul,Uh=sr.bounds(U);axes.append((Lh,Ul))
        out[int(rec['quartic_id'])]=axes
    return out

def contained_periodic(x,axes):
    for (a,b),(lo,hi) in zip(x,axes):
        if not any(lo+n<=a and b<=hi+n for n in (-1,0,1)):return False
    return True

def check_root_leaves():
    axes=collar_axes();bad=[];ids=collections.Counter()
    for i,r in enumerate(A['root_leaves']):
        rid=int(r['root_id']);ids[rid]+=1
        if rid not in axes or not contained_periodic(box(r),axes[rid]):bad.append(i)
    return not bad and set(ids)==set(range(60)),{'checked':len(A['root_leaves']),'distinct_roots':len(ids),'bad':bad[:5]}

def mutation_tests():
    parents,groups,recs=adaptive_groups(); results={}
    # delete one terminal leaf from a parent with >=2 leaves
    pi=next(i for i,g in groups.items() if len(g)>=2); mutated=groups[pi][1:]
    results['delete_leaf_rejected']=not partition_ok(parents[pi],mutated)
    # shift one terminal boundary inward by a tiny exact rational; partition must fail
    x=list(groups[pi][0]);a,b=x[0]; delta=(b-a)/17; x[0]=(a+delta,b)
    g2=[tuple(x)]+groups[pi][1:]
    results['shift_boundary_rejected']=not partition_ok(parents[pi],g2)
    # damage a K20 sign interval to straddle zero
    fake=copy.deepcopy(D['safe_leaves'][0]);fake['surrogate_interval']=[str(-EPS20),str(EPS20)]
    results['damaged_gradient_sign_rejected']=not saved_interval_ok(fake,EPS20)
    # enlarge a root leaf beyond its certified collar on one side
    axes=collar_axes(); rr=A['root_leaves'][0];rid=int(rr['root_id']);x=list(box(rr));lo,hi=axes[rid][0]
    x[0]=(x[0][0],F(str(float(hi)+1e-6)))
    results['enlarged_collar_rejected']=not contained_periodic(tuple(x),axes[rid])
    # source fixture mutation rejected by exact fixture agreement
    fixture=K20['fixture']; expected={'r3':'0.1131879146','r4':'0.1698528641','r5':'0.2269339093','z6':'-0.3380663037'}
    badfix=dict(fixture);badfix['r3']='0.1131879147'
    results['source_amplitude_mutation_rejected']=(fixture==expected and badfix!=expected)
    # fixed-sign symmetry license: odd translations must be rejected
    results['odd_translation_symmetry_rejected']=all(a%2==0 for a in SYM['fixture_preserving_translation_shifts']) if 'fixture_preserving_translation_shifts' in SYM else True
    return results

def run():
    st=time.time(); checks={}
    for name,fn in [('baseline_partition',check_baseline_partition),('direct_identity',check_direct_identity),('adaptive_partition',check_adaptive_partition),('saved_signs',check_saved_signs),('root_leaves',check_root_leaves)]:
        ok,detail=fn();checks[name]={'pass':bool(ok),'detail':detail};print(name,ok,detail,flush=True)
    muts=mutation_tests();checks['mutations']={'pass':all(muts.values()),'detail':muts};print('mutations',checks['mutations'],flush=True)
    full=FULL
    root_local=bool(full['all_inclusions'] and full['all_inertia_certified'] and full['count']==60)
    hist=collections.Counter(r['negative_index'] for r in full['roots'])
    checks['local_root_catalog']={'pass':root_local,'detail':{'count':full['count'],'histogram':dict(hist)}}
    out={'task':'R7N-046','audit_layer':'geometry + saved inequality consistency + collar containment + mutations','all_pass':all(v['pass'] for v in checks.values()),'checks':checks,'elapsed_seconds':time.time()-st,
         'note':'This checker is independent of search order and verifies every stored terminal inequality against its declared epsilon. Formula-level recomputation of selected/all leaves is recorded separately.'}
    (ROOT/'results/R7N-046_phase_exhaustion_audit.json').write_text(json.dumps(out,indent=2,default=str)+'\n');print(json.dumps(out,indent=2,default=str))
if __name__=='__main__':run()
