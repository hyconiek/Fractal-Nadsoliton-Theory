"""Read-only scientific replay of the supplied FR223 sources.

Only this audit directory receives generated results; source archives are immutable.
"""
from pathlib import Path
from fractions import Fraction as F
import argparse
import hashlib
import json
import sys
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223'
sys.path.insert(0,str(SOURCE/'src'))

def save(name,data):
    (HERE/name).write_text(json.dumps(data,indent=2)+'\n')

def hashes(p):
    return {str(f.relative_to(p)):hashlib.sha256(f.read_bytes()).hexdigest()
            for f in sorted(p.rglob('*')) if f.is_file()
            and '__pycache__' not in f.parts and '.pytest_cache' not in f.parts}

def inventory():
    paths=[ROOT/'FIN_rank7_followup_handoff_bundle'/'fin_rank7_followup',ROOT/'fin_rank7_followup',
           ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260915',SOURCE]
    hs=[hashes(p) for p in paths]; pairs=[]
    for i in range(3):
        a,b=hs[i:i+2];shared=set(a)&set(b)
        pairs.append(dict(older=str(paths[i].relative_to(ROOT)),newer=str(paths[i+1].relative_to(ROOT)),
            older_count=len(a),newer_count=len(b),identical=sum(a[k]==b[k] for k in shared),
            changed=sorted(k for k in shared if a[k]!=b[k]),missing=sorted(set(a)-set(b)),added=sorted(set(b)-set(a))))
    manifests=[]
    for p,name in [(paths[2],'CONTINUATION_MANIFEST.sha256'),(SOURCE,'MANIFEST_FR223.sha256')]:
        missing=[];changed=[];count=0
        for line in (p/name).read_text().splitlines():
            want,rel=line.split(maxsplit=1);rel=rel.lstrip('*');f=p/rel;count+=1
            if not f.exists():missing.append(rel)
            elif hashlib.sha256(f.read_bytes()).hexdigest()!=want:changed.append(rel)
        manifests.append(dict(package=str(p.relative_to(ROOT)),manifest=name,count=count,missing=missing,changed=changed))
    out=dict(comparisons=pairs,manifests=manifests,source_sha256=hs[-1])
    save('continuation_inventory.json',out)
    print(json.dumps(dict(comparisons=pairs,manifests=manifests),indent=2))

def pair(I,digits=24):
    q=10**digits
    return [str(F((I.lo*q).__floor__(),q)),str(F((I.hi*q).__ceil__(),q))]

def bounded_rationals(digits=60):
    """Outward rational rounding at every QI construction, never to nearest.

    This prevents enormous exact denominators. It only widens enclosures,
    hence can cause a false rejection but cannot create a false sign proof.
    """
    import intervals
    q=10**digits
    original=intervals.QI.__init__
    def initialize(self,lo,hi=None):
        original(self,lo,hi)
        self.lo=F((self.lo*q).__floor__(),q)
        self.hi=F((self.hi*q).__ceil__(),q)
    intervals.QI.__init__=initialize
    for a,b in [(F(1,3),F(2,3)),(F(-2,7),F(-1,7)),(F(0),F(0))]:
        v=intervals.QI(a,b);assert v.lo<=a<=b<=v.hi

def boxes():
    prior=HERE/'FR223_union_replay.json'
    if prior.exists() and not (HERE/'FR223_union_unrounded_partial.json').exists():
        save('FR223_union_unrounded_partial.json',json.loads(prior.read_text()))
    bounded_rationals()
    import frontier_local_boxes as centered
    import frontier_shifted_boxes as shifted
    ledger=json.loads((SOURCE/'results/FR223_ACTIVE_MASK_LEDGER.json').read_text())
    source=json.loads((SOURCE/'results/FR223_post_FR222_residual_search.json').read_text())
    assert ledger['centered_masks']==source['centered_masks']
    assert ledger['shifted_masks']==source['shifted_masks']
    rows=[]
    for kind in ['centered','shifted']:
        for mask in ledger[kind+'_masks']:
            start=time.monotonic();toF=lambda x:F(str(x))
            if kind=='centered':
                radii=[toF(mask[k]) for k in ['rx','ru','rv','e']]
                result=centered.raw_box(*radii)
                bounds=[[-radii[0],radii[0]]]+[[F(0),r] for r in radii[1:]]
                signs=result['signs']
                proof={'signs':signs,'boundary_schur':[pair(v) for v in result['boundary_schur']],
                       'endpoint_schur':[pair(v) for v in result['endpoint_schur']]}
            else:
                bounds=[list(map(toF,mask[k])) for k in ['x','u','v']]+[[F(0),toF(mask['e'])]]
                result=shifted.raw_shifted_box(bounds)
                proof={k:pair(result[k]) for k in ['P','P1','c2']}
                proof['reason']=result['reason']
            rstar=result['rstar']
            assert rstar.lo+bounds[0][0]>0 and rstar.hi+bounds[0][1]<1
            assert all(F(0)<=lo<=hi<F(1) for lo,hi in bounds[1:3])
            assert F(0)<=bounds[3][0]<=bounds[3][1]<F(1,2)
            row=dict(name=mask['name'],kind=kind,status=result['status'],
                     exact_decimal_box=[[str(x) for x in v] for v in bounds],proof=proof,
                     seconds=time.monotonic()-start)
            rows.append(row)
            out=dict(source='FR223_ACTIVE_MASK_LEDGER.json',navigation_buffer_used=False,
                     arithmetic='Fraction intervals, outward rounded at every operation to denominator 10^60',
                     finished=len(rows)==106,count=len(rows),passed=sum(x['status']=='INTERVAL_CERTIFIED' for x in rows),rows=rows)
            save('FR223_union_replay.json',out)
            print(mask['name'],row['status'],round(row['seconds'],2),flush=True)
    assert len(rows)==106

def mask_coverage():
    data=json.loads((HERE/'FR223_union_replay.json').read_text());assert data['finished']
    good=[x for x in data['rows'] if x['status']=='INTERVAL_CERTIFIED']
    def box(row):return [tuple(map(F,p)) for p in row['exact_decimal_box']]
    def subtract(B,A):
        inter=[(max(b[0],a[0]),min(b[1],a[1])) for a,b in zip(A,B)]
        if any(lo>=hi for lo,hi in inter):return [B]
        out=[];core=list(B)
        for k,(lo,hi) in enumerate(inter):
            if core[k][0]<lo:
                part=list(core);part[k]=(core[k][0],lo);out.append(part)
                core[k]=(lo,core[k][1])
            if core[k][1]>hi:
                part=list(core);part[k]=(hi,core[k][1]);out.append(part)
                core[k]=(core[k][0],hi)
        return out
    rows=[]
    for r in data['rows']:
        if r['status']=='INTERVAL_CERTIFIED':continue
        remaining=[box(r)]
        for g in good:
            remaining=[piece for B in remaining for piece in subtract(B,box(g))]
            if not remaining:break
        rows.append(dict(name=r['name'],covered_by_other_accepted_masks=not remaining,
                         remaining_rectangles=[[[str(x) for x in v] for v in b] for b in remaining]))
    save('FR223_failed_mask_coverage.json',rows)
    print(json.dumps([dict(name=r['name'],covered=r['covered_by_other_accepted_masks'],remaining=len(r['remaining_rectangles'])) for r in rows],indent=2))

def repair_boxes():
    bounded_rationals()
    import frontier_shifted_boxes as checker
    data=json.loads((HERE/'FR223_union_replay.json').read_text());assert data['finished']
    records=[]
    for row in data['rows']:
        if row['status']=='INTERVAL_CERTIFIED':continue
        root=[tuple(map(F,p)) for p in row['exact_decimal_box']]
        stack=[(root,0,'')];leaves=[]
        while stack:
            box,depth,path=stack.pop();out=checker.raw_shifted_box(box)
            if out['status']=='INTERVAL_CERTIFIED' or depth==4:
                leaves.append(dict(path=path,box=[[str(x) for x in v] for v in box],
                                   status=out['status'],reason=out['reason'],
                                   bounds={k:pair(out[k]) for k in ['P','P1','c2']}))
            else:
                axis=depth%4;lo,hi=box[axis];mid=(lo+hi)/2
                left=list(box);right=list(box);left[axis]=(lo,mid);right[axis]=(mid,hi)
                stack.extend([(right,depth+1,path+str(axis)+'R'),(left,depth+1,path+str(axis)+'L')])
        records.append(dict(name=row['name'],status='PASS' if all(x['status']=='INTERVAL_CERTIFIED' for x in leaves) else 'UNRESOLVED',
                            original_box=row['exact_decimal_box'],leaves=leaves))
        save('FR223_subdivision_repairs.json',records)
        print(row['name'],records[-1]['status'],len(leaves),flush=True)

def tails():
    import boundary_cover as bc
    import boundary_checker as chk
    import boundary_ising as bi
    import frontier_j6_tail as j6
    import frontier_residual_tails as fr1
    import frontier_local_boxes as local
    from intervals import QI,sqrt_interval
    # FR1 coefficients are sums of squared projected distances; certify margins.
    one=fr1.build()
    for name in ['large_J3_tail','large_J4_tail','large_J5_tail']:
        assert F(one[name]['strict_gap_interval'][0])>0
    save('FR1_replay.json',one)
    # Correct the source's unjustified cyclic-distance claim for C4: enumerate
    # all 66 pairs, since the cosine-only restriction is not rotation invariant.
    L=bi.strict_intervals();rt=sqrt_interval(QI(3),40)
    c3=[1,0,-1,0]*3;c4=[QI(1),QI(F(-1,2)),QI(F(-1,2))]*4
    c5=[QI(1),-rt/2,QI(F(1,2)),QI(0),QI(F(-1,2)),rt/2,
        QI(-1),rt/2,QI(F(-1,2)),QI(0),QI(F(1,2)),-rt/2]
    vectors=[[QI(c3[j]),c4[j],c5[j],QI((-1)**j)] for j in range(12)]
    distances=[]
    for i in range(12):
        for j in range(i):
            d=sum((L[k]/den*(vectors[i][a]-vectors[j][a])**2
                   for a,(k,den) in enumerate([(3,6),(4,6),(5,6),(6,12)])),QI(0))
            distances.append(d.hi)
    assert len(distances)==66 and max(distances)<5
    saved=json.loads((SOURCE/'results/FR42_global_large_J6_tail_5x.json').read_text())
    tree=saved['boundary_shift']['cover'];delta=F(1,200000)
    root=j6.shifted_root(delta);trie=chk.build_trie(tree['leaves'])
    tau=bi.interval_signs()['t'];R=(QI(1)-QI(*tau))/(QI(1)+QI(*tau))
    radii={'LOCAL_FR9':[F(1,6500),F(1,2432),F(1,8192),F(1,100000)],
           'LOCAL_FR16':[F(1,5000),F(1,8192),F(1,5600),F(1,100000)]}
    for key,v in radii.items():
        assert local.raw_box(*v)['status']=='INTERVAL_CERTIFIED',key
    rows=[]
    def walk(node,box,A,B,path=''):
        leaf=node['leaf']
        if leaf is not None:
            assert leaf['path']==path and leaf['box']==bc._box_json(box)
            a=bc.bernstein_bounds(A);b=bc.bernstein_bounds(B);reason=leaf['reason']
            assert leaf['A_bounds']==bc.qi_json(a) and leaf['B_bounds']==bc.qi_json(b)
            if reason=='SAFE_A_NONPOS':assert a.hi<=0
            elif reason=='SAFE_B_NONNEG':assert b.lo>=0
            else:
                rx,ru,rv,re=radii[reason]
                assert box[0][0]>=R.hi-rx and box[0][1]<=R.lo+rx
                assert 1-box[1][0]<=ru and 1-box[2][0]<=rv
                assert F(1,1000001)<=re
            rows.append(dict(path=path,reason=reason))
            return
        children=node['children'];assert len(children)==2
        axes={a for a,s in children};sides={s for a,s in children}
        assert len(axes)==1 and sides=={'L','R'}
        axis=next(iter(axes));al,ar=bc.split_bernstein(A,root['degrees_A'],axis)
        bl,br=bc.split_bernstein(B,root['degrees_B'],axis);xl,xr=bc.split_box(box,axis)
        walk(children[(axis,'L')],xl,al,bl,path+str(axis)+'L')
        walk(children[(axis,'R')],xr,ar,br,path+str(axis)+'R')
    walk(trie,((F(0),F(1)),)*3,root['A'],root['B'])
    assert len(rows)==637 and tree['unresolved_count']==0
    assert F(bi.interval_signs()['trace_margin'][0])-3*delta>0
    assert F(5,1000001)<delta
    out=dict(status='PASS',leaves=len(rows),all_66_pair_diameter_squared_upper=str(max(distances)),
             perturbation_upper='5/1000001',boundary_reserve=str(delta),
             y_upper='1/1000000',local_boxes={k:list(map(str,v)) for k,v in radii.items()},
             source_correction='C4 pair distances are not assumed cyclic; all 66 pairs checked.',rows=rows)
    save('FR42_replay.json',out);print('FR1 + FR42 PASS: 637 leaves, 66 diameter pairs',flush=True)

def geometry():
    """Independently derive FR1 bounds by explicit projected feature distances."""
    bounded_rationals()
    import boundary_ising as bi
    from intervals import QI,sqrt_interval
    L=bi.strict_intervals();rt=sqrt_interval(QI(3),40)
    c3=[1,0,-1,0]*3;c4=[QI(1),QI(F(-1,2)),QI(F(-1,2))]*4
    c5=[QI(1),-rt/2,QI(F(1,2)),QI(0),QI(F(-1,2)),rt/2,
        QI(-1),rt/2,QI(F(-1,2)),QI(0),QI(F(1,2)),-rt/2]
    scales=[sqrt_interval(L[k]/d,40) for k,d in [(3,6),(4,6),(5,6),(6,12)]]
    features=[[scales[a]*v for a,v in enumerate([QI(c3[j]),c4[j],c5[j],QI((-1)**j)])] for j in range(12)]
    def dot(v,w):return sum((a*b for a,b in zip(v,w)),QI(0))
    diffs=[[a-b for a,b in zip(v,features[0])] for v in features]
    def projected(j,axis,remove_alt=False):
        v=diffs[axis];d=diffs[j]
        x=dot(d,d)-dot(d,v)**2/dot(v,v)
        return x-d[3]**2 if remove_alt else x
    a0=F(1,30);s0=F(1,128);t0=F(1,9)
    a_bound=sum((projected(j,4)*a0 for j in [1,3,5,7,9,11]),QI(0))
    a_bound+=sum((projected(j,4)*a0*a0 for j in [2,6,10]),QI(0))
    s_bound=L[6]/12+sum((projected(j,6,True)*s0 for j in range(12) if j%3),QI(0))
    groups={1:[2,10],2:[3,9],3:[1,11,4,8],4:[6]}
    t_bound=sum((projected(j,5)*t0**power for power,js in groups.items() for j in js),QI(0))
    sigma=(2*L[3]*(L[4]+L[5])-L[4]*L[5])/(24*L[3])
    rows={}
    for name,b in [('large_J3',a_bound),('large_J4',s_bound),('large_J5',t_bound)]:
        gap=sigma-b;assert gap.lo>0
        rows[name]=dict(bound=pair(b),strict_gap=pair(gap))
    save('FR1_independent_geometry.json',dict(status='PASS',rows=rows,
         method='Explicit normalized C4 projections; no imported closed-form coefficient used.'))
    print('FR1 independent geometric bounds PASS',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['inventory','boxes','tails','mask_coverage','repair_boxes','geometry']);args=p.parse_args()
    globals()[args.stage]()
