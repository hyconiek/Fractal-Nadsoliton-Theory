"""Independent intake of the R7N handoff. Supplied archives are read-only."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile
import shutil
import time
from functools import lru_cache

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'FIN_R7N_HANDOFF_20260920'
sys.path.insert(0,str(ROOT))
from fin_rank7_intake_review import scientific_rechecks as sr
sr.iv.dps=55

def load(rel):return json.loads((SOURCE/rel).read_text())
def save(name,data): (HERE/name).write_text(json.dumps(data,indent=2)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inventory():
    out={}
    for name in ['MANIFEST.sha256','INPUT_MANIFEST.sha256']:
        wrong=[];missing=[];n=0
        for line in (SOURCE/name).read_text().splitlines():
            if not line.strip():continue
            want,rel=line.split(maxsplit=1);rel=rel.lstrip('*');n+=1
            p=SOURCE/rel
            if not p.exists():missing.append(rel)
            elif sha(p)!=want:wrong.append(rel)
        out[name]=dict(entries=n,missing=missing,changed=wrong)
    copied=SOURCE/'inputs/intake_review_20260919'
    comparisons=[]
    for p in copied.iterdir():
        if p.is_file() and (ROOT/'fin_rank7_intake_review'/p.name).exists():
            comparisons.append(dict(path=p.name,identical=sha(p)==sha(ROOT/'fin_rank7_intake_review'/p.name)))
    out['intake_dependency_comparison']=comparisons
    aliases={
        'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223.zip':ROOT/'FIN_rank7_CONTINUATION_HANDOFF_20260916_FR223.zip',
        'FIN_Rank7_Next_Campaign_Plan_EN_20260919.md':ROOT/'FIN_Rank7_Next_Campaign_Plan_EN_20260919.md',
        'AGENTS(1).md':SOURCE/'inputs/AGENTS_preconsolidation.md'}
    resolved=[]
    for line in (SOURCE/'INPUT_MANIFEST.sha256').read_text().splitlines():
        want,old=line.split(maxsplit=1);p=aliases[Path(old).name]
        resolved.append(dict(original_path=old,local_path=str(p.relative_to(ROOT)),matches=p.exists() and sha(p)==want))
    assert all(r['matches'] for r in resolved)
    out['resolved_original_inputs']=resolved
    out['archive_sha256']={str(p.relative_to(SOURCE)):sha(p) for p in sorted(SOURCE.rglob('*')) if p.is_file()}
    save('inventory.json',out);print(json.dumps({k:v for k,v in out.items() if k!='archive_sha256'},indent=2))

def portable():
    if (HERE/'portable_replay.json').exists():
        previous=json.loads((HERE/'portable_replay.json').read_text())
        if previous['returncode']!=0:save('portable_sanitized_copy_attempt.json',previous)
    parent=Path(tempfile.mkdtemp(prefix='fin_r7n_review_'));copy=parent/'fin_rank7_next_campaign'
    # The supplied manifest includes four pytest-cache files; an exact copy is
    # required for its own manifest check. Bytecode writing remains disabled.
    shutil.copytree(SOURCE,copy)
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    start=time.monotonic()
    proc=subprocess.run([sys.executable,str(copy/'portable_verify.py')],cwd=copy,env=env,capture_output=True,text=True,timeout=180)
    out=dict(command='python portable_verify.py in isolated copy',returncode=proc.returncode,
             stdout=proc.stdout,stderr=proc.stderr,seconds=time.monotonic()-start,copy=str(copy))
    save('portable_replay.json',out);print(json.dumps(out,indent=2))

def polynomial_mul(A,B):
    C={}
    for e,a in A.items():
        for f,b in B.items():
            g=tuple(x+y for x,y in zip(e,f));C[g]=C.get(g,sr.I(0))+a*b
    return C

def reconstruct_surrogates():
    """Reconstruct log(E exp(alpha h)) coefficients by the E'=K'E recurrence."""
    amp=list(map(sr.I,['0.1131879146','0.1698528641','0.2269339093','-0.3380663037']))
    norm=2*sr.iv.sqrt(3)
    base=[]
    for i,k in enumerate([3,4,5]):
        e=[0,0,0];e[i]=1;base.append((k,tuple(e),amp[i]/norm))
        e[i]=-1;base.append((12-k,tuple(e),amp[i]/norm))
    base.append((6,(0,0,0),amp[3]/norm))
    state={(0,(0,0,0)):sr.I(1)};E={0:{(0,0,0):sr.I(1)}};K={};aggregates={};aggregate={}
    for n in range(1,21):
        new={}
        for (lab,e),v in state.items():
            for k,f,a in base:
                key=((lab+k)%12,tuple(x+y for x,y in zip(e,f)))
                new[key]=new.get(key,sr.I(0))+v*a/n
        state=new;E[n]={e:v for (lab,e),v in state.items() if lab==0}
        out=dict(E[n])
        for k in range(1,n):
            for e,v in polynomial_mul(K[k],E[n-k]).items():out[e]=out.get(e,sr.I(0))-sr.I(F(k,n))*v
        K[n]=out
        for e,v in out.items():aggregate[e]=aggregate.get(e,sr.I(0))+v
        if n in [4,16,20]:aggregates[n]=dict(aggregate)
        if n%4==0:print('independent log series',n,flush=True)
    reports={}
    for N in [16,20]:
        source=load(f'results/R7N-043_K{N}_surrogate.json');actual=aggregates[N]
        assert source['fixture']=={'r3':'0.1131879146','r4':'0.1698528641','r5':'0.2269339093','z6':'-0.3380663037'}
        retained={tuple(r['frequency']):r for r in source['terms']};drop=[F(0)]*3
        for e,v in actual.items():
            if e==(0,0,0) or next(x for x in e if x)!=abs(next(x for x in e if x)):continue
            opp=tuple(-x for x in e);va,vb=sr.bounds(v);wa,wb=sr.bounds(actual[opp])
            lo,hi=max(va,wa),min(vb,wb);assert lo<=hi
            if e in retained:
                sl,sh=map(F,retained[e]['coefficient_interval']);assert sl<=lo<=hi<=sh,(N,e)
            else:
                for k in range(3):drop[k]+=2*max(abs(lo),abs(hi))*abs(e[k])
        assert set(retained)<=set(actual)
        R=sr.I(F(source['complex_analytic_tail']['R']))
        assert sr.bounds(R)[0]>1
        H=sum(amp[:3],sr.I(0))/sr.iv.sqrt(3)+abs(amp[3])/sr.iv.sqrt(12)
        RH=R*H;assert sr.bounds(RH-sr.iv.pi/2)[1]<0
        M=R*(amp[2]/sr.iv.sqrt(3))*sr.iv.exp(2*RH)/sr.iv.cos(RH)
        tail=M*R**(-N-1)/(1-1/R)
        eps=max(drop)+sr.bounds(tail)[1]
        # Use our recomputed bound, never a recorded PASS flag or unverified epsilon.
        reports[str(N)]=dict(retained=len(retained),drop_upper=[str(v) for v in drop],
                            tail=sr.serialize(tail),epsilon=str(eps),source_epsilon=source['uniform_full_gradient_error_bound'],
                            source_coefficients_enclose_independent_reconstruction=True)
    # Exact quartic coefficients are independently available from the same construction.
    qterms=[]
    for e,v in aggregates[4].items():
        if e==(0,0,0) or next(x for x in e if x)<0:continue
        lo,hi=sr.bounds(v)
        if lo==hi==0:continue
        qterms.append(dict(frequency=list(e),coefficient_interval=sr.serialize(2*v)))
    reports['quartic_terms']=qterms
    save('surrogate_reconstruction.json',reports)
    print(json.dumps({k:{kk:vv for kk,vv in v.items() if kk in ['retained','source_coefficients_enclose_independent_reconstruction']} for k,v in reports.items() if k!='quartic_terms'},indent=2))

def interval_pair(x):
    lo,hi=map(F,x);assert lo<=hi
    return sr.I((lo+hi)/2)+sr.I((hi-lo)/2)*sr.iv.mpf([-1,1])

def terms_for(kind):
    if kind=='quartic':
        data=json.loads((HERE/'surrogate_reconstruction.json').read_text())['quartic_terms']
        return [(interval_pair(r['coefficient_interval']),tuple(r['frequency'])) for r in data]
    data=load(f'results/R7N-043_K{kind}_surrogate.json')['terms']
    return [(2*interval_pair(r['coefficient_interval']),tuple(r['frequency'])) for r in data]

def quartet_hessian(phi):
    H=[[sr.I(0) for _ in range(3)] for __ in range(3)]
    for coef,n in terms_for('quartic'):
        v=sr.iv.cos(sum((n[k]*phi[k] for k in range(3)),sr.I(0)))
        for i in range(3):
            for j in range(3):H[i][j]-=coef*n[i]*n[j]*v
    return H

def collars():
    reports={}
    for kind,catfile,colfile in [
        ('quartic','inputs/FR223_20260916/results/R7P-089_quartic_roots.json','results/R7N-035_quartic_uniqueness_collars.json'),
        ('full','inputs/FR223_20260916/certificates/R7P-092_full_phase_roots.json','results/R7N-041_full_uniqueness_collars.json')]:
        cat=load(catfile)['roots'];col=load(colfile)['roots'];rows=[];assert len(cat)==len(col)==60
        for rec,c in zip(cat,col):
            local=sr.phase_root(rec,kind)
            assert list(map(F,local['center']))==list(map(F,c['center']))
            center=list(map(F,c['center']));radius=F(c['certified_radius']);assert radius>F(local['radius'])
            X=[sr.I(x)+sr.I(radius)*sr.iv.mpf([-1,1]) for x in center]
            H=quartet_hessian(X) if kind=='quartic' else sr.phase_FH(X,'full')[1]
            A=[[F(v) for v in row] for row in c['preconditioner']]
            determinant=A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])-A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])+A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0])
            assert determinant
            qs=[]
            for i in range(3):
                q=F(0)
                for j in range(3):
                    e=sr.I(int(i==j))-sum((sr.I(A[i][k])*H[k][j] for k in range(3)),sr.I(0))
                    lo,hi=sr.bounds(e);q+=max(abs(lo),abs(hi))
                qs.append(q)
            assert max(qs)<1,(kind,rec.get('id'),float(max(qs)))
            axes=[]
            for x in center:
                # Inner normalized collar; every contained cell lies in actual collar.
                lower=sr.bounds(sr.I(x-radius)/(2*sr.iv.pi))[1]
                upper=sr.bounds(sr.I(x+radius)/(2*sr.iv.pi))[0]
                axes.append([str(lower),str(upper)])
            rows.append(dict(id=rec.get('id',rec.get('quartic_id')),radius=str(radius),center=c['center'],
                             q_upper=str(max(qs)),preconditioner=c['preconditioner'],normalized_inner_axes=axes,
                             local_root=local))
        min_separation=None
        for a,b in combinations(rows,2):
            diffs=[]
            for x,y in zip(a['center'],b['center']):
                shifts=[]
                for n in [-1,0,1]:
                    lo,hi=sr.bounds(sr.I(F(x)-F(y))+n*2*sr.iv.pi)
                    shifts.append(F(0) if lo<=0<=hi else min(abs(lo),abs(hi)))
                diffs.append(min(shifts))
            sep=max(diffs)-F(a['radius'])-F(b['radius']);assert sep>0
            min_separation=sep if min_separation is None else min(sep,min_separation)
        reports[kind]=dict(roots=rows,count=60,minimum_collar_separation=str(min_separation),
                           indices=dict(Counter(str(r['local_root']['negative_index']) for r in rows)))
        save('collar_replay.json',reports);print(kind,'60 local roots and injectivity collars PASS',flush=True)

def parse_box(r):return tuple(tuple(map(F,p)) for p in r['box'])
def partition(parent,leaves):
    """Check a complete midpoint rectangle partition from geometry, not labels."""
    assert leaves and all(all(p<=a<b<=q for (p,q),(a,b) in zip(parent,x)) for x in leaves)
    if len(leaves)==1:
        assert leaves[0]==parent;return
    for axis in range(3):
        lo,hi=parent[axis];m=(lo+hi)/2
        left=[b for b in leaves if b[axis][1]<=m];right=[b for b in leaves if b[axis][0]>=m]
        if left and right and len(left)+len(right)==len(leaves):
            lp=list(parent);rp=list(parent);lp[axis]=(lo,m);rp[axis]=(m,hi)
            partition(tuple(lp),left);partition(tuple(rp),right);return
    raise AssertionError('Not a complete midpoint partition')

def phase_geometry():
    quart=load('checkpoints/R7N-037_quartic_cover_normalized.json')
    b=load('checkpoints/R7N-044_full_cover_k16.json');d=load('checkpoints/R7N-044_K20_residual.json');a=load('checkpoints/R7N-044_K20_adaptive.json')
    assert not quart['unresolved_leaves'] and not a['unresolved_leaves']
    # The K20 direct list replaces the entire saved K16 unresolved list once.
    refs=d['safe_leaves']+d['failed_leaves'];seen=set()
    for row in refs:
        i=row['source_index'];assert i not in seen;seen.add(i)
        assert row['box']==b['unresolved_leaves'][i]['box']
    assert seen==set(range(len(b['unresolved_leaves'])))
    coll=json.loads((HERE/'collar_replay.json').read_text());out={}
    for kind,safes,rootleaves in [
        ('quartic',quart['safe_leaves'],quart['root_leaves']),
        ('full',b['safe_leaves']+d['safe_leaves']+a['safe_leaves'],b['root_leaves']+a['root_leaves'])]:
        leaves=safes+rootleaves;groups=defaultdict(list)
        for row in leaves:
            box=parse_box(row);cell=tuple(int(4*(lo+hi)/2) for lo,hi in box)
            groups[cell].append(box)
        assert set(groups)==set(product(range(4),repeat=3))
        for key,boxes in groups.items():partition(tuple((F(k,4),F(k+1,4)) for k in key),boxes)
        collars_by_id={int(r['id']):r for r in coll[kind]['roots']};ids=set()
        for row in rootleaves:
            rid=int(row['root_id']);c=collars_by_id[rid];box=parse_box(row);ids.add(rid)
            for (lo,hi),ab in zip(box,c['normalized_inner_axes']):
                x,y=map(F,ab);assert any(x+n<=lo and hi<=y+n for n in [-1,0,1])
        assert ids==set(range(60))
        out[kind]=dict(safe=len(safes),root_leaves=len(rootleaves),complete_partition=True,root_ids=sorted(ids))
        print(kind,'complete partition and collar containments PASS',flush=True)
    save('phase_geometry.json',out)

def symmetry():
    data=json.loads((HERE/'collar_replay.json').read_text())['quartic']['roots'];perms=[];labels=[]
    for eps in [1,-1]:
        for a in range(0,12,2):
            perm=[]
            for row in data:
                image=[];locator=[]
                for k,v in zip([3,4,5],row['center']):
                    phi=sr.I(F(v))+sr.I(F(1,10**7))*sr.iv.mpf([-1,1])
                    im=eps*(phi+2*sr.iv.pi*k*a/12);image.append(im)
                    locator.append(float(sr.mid(im))%(2*math.pi))
                def distance(target):
                    ds=[abs(x-float(F(y))) for x,y in zip(locator,target['center'])]
                    return sum(min(d,2*math.pi-d)**2 for d in ds)
                target=min(data,key=distance);radius=F(target['radius'])
                for im,v in zip(image,target['center']):
                    c=F(v);n=round((float(c)-sr.mid(im))/(2*math.pi))
                    lo,hi=sr.bounds(im+n*2*sr.iv.pi)
                    assert c-radius<=lo<=hi<=c+radius
                perm.append(int(target['id']))
            assert set(perm)==set(range(60));perms.append(perm);labels.append([a,eps])
    ps={tuple(p) for p in perms}
    for p in perms:
        for q in perms:assert tuple(p[q[i]] for i in range(60)) in ps
    orbits={tuple(sorted({p[i] for p in perms})) for i in range(60)}
    assert Counter(map(len,orbits))=={6:8,12:1}
    save('symmetry_replay.json',dict(status='PASS',group_order=12,group_elements=labels,
         orbits=[list(x) for x in sorted(orbits)],method='Exact action, numerical target proposals, interval inclusion in unique-root collars.'))
    print('Quartic exact subgroup orbit classification PASS',flush=True)

def replay_leaves(kind):
    previous=HERE/('leaf_replay_'+kind+'.json')
    if previous.exists() and not (HERE/('leaf_replay_'+kind+'_initial_partial.json')).exists():
        save('leaf_replay_'+kind+'_initial_partial.json',json.loads(previous.read_text()))
    if kind=='quartic':leaves=load('checkpoints/R7N-037_quartic_cover_normalized.json')['safe_leaves'];epsilon=F(0)
    else:
        epsilon=F(json.loads((HERE/'surrogate_reconstruction.json').read_text())[kind]['epsilon'])
        if kind=='16':leaves=load('checkpoints/R7N-044_full_cover_k16.json')['safe_leaves']
        else:leaves=load('checkpoints/R7N-044_K20_residual.json')['safe_leaves']+load('checkpoints/R7N-044_K20_adaptive.json')['safe_leaves']
    sr.iv.dps=35;terms=terms_for(kind);twopi=2*sr.iv.pi;failed=[];minimum=None;start=time.monotonic()
    @lru_cache(maxsize=150000)
    def sine(lo,hi):
        if hi-lo>=1:return sr.iv.mpf([-1,1])
        n=lo.__floor__();return sr.iv.sin(twopi*interval_pair((lo-n,hi-n)))
    for idx,row in enumerate(leaves):
        box=parse_box(row);k=row['gradient_component'];assert k in [0,1,2];g=sr.I(0)
        for coef,n in terms:
            if not n[k]:continue
            lo=sum((a*e if e>=0 else b*e for (a,b),e in zip(box,n)),F(0))
            hi=sum((b*e if e>=0 else a*e for (a,b),e in zip(box,n)),F(0))
            g-=coef*n[k]*sine(lo,hi)
        lo,hi=sr.bounds(g);margin=max(lo-epsilon,-hi-epsilon)
        if margin<=0:failed.append(idx)
        else:minimum=margin if minimum is None else min(minimum,margin)
        if (idx+1)%2000==0 or idx==len(leaves)-1:
            out=dict(kind=kind,processed=idx+1,total=len(leaves),failed=failed,
                     complete=idx+1==len(leaves),minimum_margin=str(minimum),seconds=time.monotonic()-start)
            save('leaf_replay_'+kind+'.json',out);print(kind,idx+1,'/',len(leaves),'failed',len(failed),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('stage',choices=['inventory','portable','reconstruct_surrogates','collars','phase_geometry','symmetry','leaves']);parser.add_argument('--kind',choices=['quartic','16','20']);a=parser.parse_args()
    if a.stage=='leaves':replay_leaves(a.kind)
    else:globals()[a.stage]()
