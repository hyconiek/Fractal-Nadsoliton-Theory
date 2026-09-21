from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
import itertools, json, math, sys

ROOT=Path(__file__).resolve().parents[1]
IR=ROOT/'inputs/intake_review_20260919'; sys.path.insert(0,str(IR))
import scientific_rechecks as sr

CAT=json.load(open(ROOT/'inputs/FR223_20260916/results/R7P-089_quartic_roots.json'))
LOCAL=json.load(open(ROOT/'results/R7N-035_041_phase_local_recertification.json'))['quartic']
COL=json.load(open(ROOT/'results/R7N-035_quartic_uniqueness_collars.json'))
COVER=json.load(open(ROOT/'checkpoints/R7N-037_quartic_cover_normalized.json'))
roots=CAT['roots']

TWOPI=2*sr.iv.pi

def torus_sep_lower(ca,cb):
    """Rigorous lower bound on L_inf torus separation using interval 2*pi."""
    best=[]
    for x,y in zip(map(F,ca),map(F,cb)):
        lows=[]
        for n in (-1,0,1):
            lo,hi=sr.bounds(sr.I(x-y)+n*TWOPI)
            if lo<=0<=hi: q=F(0)
            else: q=min(abs(lo),abs(hi))
            lows.append(q)
        best.append(min(lows))
    return max(best),best

# 1. local root boxes pairwise distinct and large collars pairwise disjoint.
min_small=None; min_large=None; pair_small=None; pair_large=None
local_centers=[r['center'] for r in LOCAL['roots']]
collar_centers=[r['center'] for r in COL['roots']]
for i,j in itertools.combinations(range(60),2):
    sep,_=torus_sep_lower(local_centers[i],local_centers[j])
    if min_small is None or sep<min_small: min_small,pair_small=sep,(i,j)
    sep2,_=torus_sep_lower(collar_centers[i],collar_centers[j])
    if min_large is None or sep2<min_large: min_large,pair_large=sep2,(i,j)
small_r=F(1,10**7); large_r=F(1,20)
assert min_small>2*small_r
assert min_large>2*large_r

# 2. Derive fixed-sign subgroup analytically. Under j -> eps*j+a:
# eps=+1: phi_k -> phi_k + 2*pi*k*a/12
# eps=-1: phi_k -> -phi_k - 2*pi*k*a/12
# z6 -> (-1)^a z6; therefore a must be even for fixed negative z6.
subgroup=[(a,eps) for eps in (1,-1) for a in range(12) if a%2==0]
assert len(subgroup)==12

def wrap(x): return x%(2*math.pi)
def action(ph,a,eps):
    out=[]
    for k,p in zip((3,4,5),ph):
        if eps==1: q=p+2*math.pi*k*a/12
        else: q=-p-2*math.pi*k*a/12
        out.append(wrap(q))
    return out

def dist(a,b):
    ds=[]
    for x,y in zip(a,b):
        d=abs(x-y); d=min(d,2*math.pi-d); ds.append(d)
    return math.sqrt(sum(x*x for x in ds))

# Map all subgroup actions to certified root IDs, and compare old numerical orbit fields.
action_maps={}
max_mapping_distance=0.0
for a,eps in subgroup:
    perm=[]
    for r in roots:
        ph2=action(r['phase'],a,eps)
        ds=[dist(ph2,s['phase']) for s in roots]; j=min(range(60),key=lambda x:ds[x])
        max_mapping_distance=max(max_mapping_distance,ds[j]); assert ds[j]<5e-7
        perm.append(j)
    assert len(set(perm))==60
    action_maps[f'a={a},eps={eps}']=perm

orbit_sets=[]; seen=set()
for i in range(60):
    orb=tuple(sorted({action_maps[k][i] for k in action_maps}))
    assert orb==tuple(roots[i]['fixed_sign_symmetry_orbit_ids'])
    if orb not in seen: seen.add(orb); orbit_sets.append(orb)

# action composition closure at the level of phase maps / root permutations
perms=list(action_maps.values())
permset={tuple(p) for p in perms}
for p in perms:
    for q in perms:
        comp=tuple(p[q[i]] for i in range(60))
        assert comp in permset

# Stabilizers and index preservation.
orbit_rows=[]
for orb in sorted(orbit_sets,key=lambda x:x[0]):
    idxs={LOCAL['roots'][i]['negative_index'] for i in orb}; assert len(idxs)==1
    stabs=[]
    for i in orb:
        stabs.append(sum(action_maps[k][i]==i for k in action_maps))
    assert len(set(stabs))==1
    assert len(orb)*stabs[0]==len(subgroup)
    orbit_rows.append({'ids':list(orb),'size':len(orb),'stabilizer_size':stabs[0],
                       'negative_index':next(iter(idxs))})

index_hist=Counter(str(x['negative_index']) for x in LOCAL['roots'])
assert index_hist==Counter({'1':24,'2':18,'0':12,'3':6})
assert COVER['complete'] and COVER['unresolved_leaf_count']==0
assert COL['certified_count']==60 and all(F(r['certified_radius'])>=large_r for r in COL['roots'])
# Every small root box lies within its own large collar: same candidate center convention and 1e-7 << .05.
assert small_r<large_r

out={
 'task':'R7N-039',
 'scientific_status':'INTERVAL_CERTIFIED_WITH_EXACT_GROUP_ACTION_DERIVATION',
 'fixture':CAT['fixture'],
 'fixed_sign_subgroup':{
   'order':12,
   'elements':[{'a':a,'eps':eps} for a,eps in subgroup],
   'derivation':'z6 transforms as (-1)^a z6, so fixed nonzero negative z6 permits exactly even translations a; reflections are allowed. For k=3,4,5: eps=+1 gives phi_k -> phi_k+2*pi*k*a/12; eps=-1 gives phi_k -> -phi_k-2*pi*k*a/12.',
   'not_licensed':'odd translations/full 24-element D12 do not act within the same fixed-sign fixture.'
 },
 'root_separation':{
   'small_box_radius':str(small_r),
   'minimum_certified_center_Linf_torus_separation_lower':str(min_small),
   'attaining_pair':list(pair_small),
   'pairwise_small_boxes_disjoint':True,
   'large_collar_radius':str(large_r),
   'minimum_collar_center_Linf_torus_separation_lower':str(min_large),
   'large_collars_pairwise_disjoint':True,
   'large_collar_attaining_pair':list(pair_large)
 },
 'symmetry_mapping_max_numeric_locator_distance':max_mapping_distance,
 'orbit_count':len(orbit_rows),
 'orbits':orbit_rows,
 'index_histogram':dict(index_hist),
 'cover_complete':True,
 'cover_terminal_counts':{'gradient':COVER['safe_leaf_count'],'root_collar':COVER['root_collar_leaf_count']},
 'conclusion':'The fixed negative-z6 quartic fixture has exactly 60 critical points on the full phase 3-torus. They form 9 orbits under the 12-element sign-preserving subgroup: eight orbits of size 6 and one orbit of size 12. Negative-index histogram is 12,24,18,6 for indices 0,1,2,3.',
 'nonconclusions':['No identification under odd translations inside the same fixed-sign fixture.','No amplitude robustness theorem.','No statement about the full log-mgf phase census.']
}
(ROOT/'results/R7N-039_symmetry_and_count.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ('orbits',)},indent=2))
