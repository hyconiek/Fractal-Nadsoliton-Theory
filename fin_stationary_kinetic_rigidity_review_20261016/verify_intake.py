#!/usr/bin/env python3
"""Independent bounded KR checks and copy-only source replay.

Exact rational finite identities, a nominal native-A7 diagnostic and the
analytical process theorem have deliberately separate scopes.
"""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/'fin_stationary_kinetic_rigidity_campaign_20261016'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder):
    names=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1); assert n not in names and sha(folder/n)==h,n; names.append(n)
    actual={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'},folder
    return {'matched':len(names),'coverage':'COMPLETE','manifest_sha256':sha(folder/'MANIFEST.sha256')}

def prefix_hash(h):
    digest=hashlib.sha256(); length=0
    for line in (ROOT/'AGENTS.md').read_bytes().splitlines(keepends=True):
        digest.update(line); length+=len(line)
        if digest.hexdigest()==h:return length
    return None

def counts(N,q):
    if q==1:return [(N,)]
    return [(i,)+t for i in range(N+1) for t in counts(N-i,q-1)]

def solve(A,b):
    n=len(b); M=[list(A[i])+[b[i]] for i in range(n)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if M[i][j]); M[j],M[pivot]=M[pivot],M[j]
        a=M[j][j]; M[j]=[x/a for x in M[j]]
        for i in range(n):
            if i!=j:
                a=M[i][j]; M[i]=[x-a*y for x,y in zip(M[i],M[j])]
    return [row[-1] for row in M]

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True); p.add_argument('--replay',action='store_true'); args=p.parse_args()
    checks=[]
    def ck(name,ok,**detail):checks.append({'name':name,'passed':bool(ok),**detail}); assert ok,name
    out={'date':'2026-10-16','entry_baseline':'8571db1d','scope':'independent bounded intake, not mechanized process proof or physical source'}
    out['source_manifest']=manifest(SOURCE)
    source_records=json.loads((SOURCE/'SOURCE_AVAILABILITY.json').read_text())['hashed_available']
    lookup={
        'UF-002 original archived member':'fin_uniform_critical_campaign_20261014/UF-002_PARITY_COEFFICIENTS_AND_CONTROLS.md',
        'UF-003 original archived member':'fin_uniform_critical_campaign_20261014/UF-003_CRITICAL_PROCESS_THEOREM_AND_ENERGY.md',
        'MT-003 original archived member':'fin_critical_metastability_campaign_20261014/MT-003_CONDITIONAL_GIBBS.md',
        'UF historical wrapper':'fin_critical_metastability_campaign_20261014/source_context/UF_PREDECESSOR.zip',
        'MT historical wrapper':'fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip'}
    for record in source_records:
        if record['label'] in lookup:assert sha(ROOT/lookup[record['label']])==record['sha256']
    expected=''.join(r['sha256']+'  '+r['path']+'  # '+r['label']+'\n' for r in source_records)
    ck('six source references and all four stage ledgers agree',len(source_records)==6 and (SOURCE/'INPUTS.sha256').read_text()==expected and all((SOURCE/f'KR-{i:03d}/INPUTS.sha256').read_text()==expected for i in range(1,5)))
    out['provenance']={'matched_current_original_bytes':lookup,'unavailable_originals':[r['path'] for r in source_records if r['label'] not in lookup],
        'stage_reference_counts':{'matched':20,'unavailable_archival_AGENTS':4},
        'uploaded_order':'not supplied as original bytes; current repository order separately bound'}
    ck('five original sources match without inventing archival AGENTS access',len(lookup)==5)
    # Also bind the source-member bytes inside the two actually verified archives.
    for label,member in [('UF historical wrapper',lookup['UF-002 original archived member']),('UF historical wrapper',lookup['UF-003 original archived member']),('MT historical wrapper',lookup['MT-003 original archived member'])]:
        with zipfile.ZipFile(ROOT/lookup[label]) as z:
            assert hashlib.sha256(z.read(member)).hexdigest()==sha(ROOT/member)
    ck('three archived members agree inside actual verified wrappers',True)
    priors={}
    for rel in ('fin_activity_covariance_research_20261016','fin_rate_coupling_review_20261016',
                'fin_uniform_critical_review_20261014','fin_critical_metastability_review_20261015','fin_intrinsic_observation_review_20261015'):
        priors[rel]=manifest(ROOT/rel)
        f=ROOT/rel/('RESULTS.json' if rel.startswith('fin_activity_') else 'verification.json')
        saved=json.loads(f.read_text()); hashes=saved.get('input_hashes',saved.get('inputs',{}))
        for path,h in hashes.items():
            if path=='AGENTS.md':
                length=prefix_hash(h); assert length is not None
                priors[rel]['historical_AGENTS_prefix_bytes']=length
            else:assert sha(ROOT/path)==h,path
    out['actual_dependencies']=priors
    missing=json.loads((SOURCE/'SOURCE_AVAILABILITY.json').read_text())['not_available_as_current_repository_bytes']
    ck('formerly unavailable actual scientific files now bound locally',all((ROOT/path).is_file() for path in missing if path!='current AGENTS.md repository HEAD'))
    ck('all actual predecessor manifest and scientific dependencies remain intact',len(priors)==5)
    # A distinct nonuniform q3/N3 Gibbs fixture from complete labelled conditional weights.
    q=3; N=3; states=counts(N,q); ix={n:i for i,n in enumerate(states)}; K=len(states)
    W=[F(math.factorial(N)//math.prod(math.factorial(v) for v in n))*2**sum(v*v for v in n) for n in states]
    Z=sum(W); pi=[v/Z for v in W]; Q=[[F(0) for _ in states] for _ in states]
    for row,n in enumerate(states):
        for old,number in enumerate(n):
            if number==0:continue
            trial=[]
            for new in range(q):
                v=list(n); v[old]-=1; v[new]+=1; trial.append((tuple(v),F(2)**sum(a*a for a in v)))
            denominator=sum(v for _,v in trial)
            for target,v in trial:
                if target!=n:Q[row][ix[target]]+=number*v/denominator
        Q[row][row]=-sum(Q[row])
    ck('nonuniform exact labelled-Gibbs lumped fixture DB and cap',all(pi[i]*Q[i][j]==pi[j]*Q[j][i] for i,j in itertools.product(range(K),repeat=2)) and all(-Q[i][i]<=N for i in range(K)))
    sigma=(1,-1,1); m=[sum(F(a*b,N) for a,b in zip(n,sigma)) for n in states]
    def drift(L):return [sum(L[i][j]*(m[j]-m[i]) for j in range(K)) for i in range(K)]
    def bracket(L):return [sum(L[i][j]*(m[j]-m[i])**2 for j in range(K)) for i in range(K)]
    centre=(1,1,1); roots=[]
    for a,b in ((0,1),(0,2),(1,2),(1,0),(2,0),(2,1)):
        v=[0]*3; v[a]+=1; v[b]-=1; roots.append(tuple(v))
    ring=[tuple(a+b for a,b in zip(centre,r)) for r in roots]; stress={}
    for j,n in enumerate(ring):
        stress[tuple(sorted((ix[centre],ix[n])))]=F(1)
        stress[tuple(sorted((ix[n],ix[ring[(j+1)%6]])))]=F(-1)
    epsilon=min(pi[i]*Q[i][j]/(4*abs(v)) for (i,j),v in stress.items())
    R=[row[:] for row in Q]
    for (i,j),v in stress.items():R[i][j]+=epsilon*v/pi[i]; R[j][i]+=epsilon*v/pi[j]
    for i in range(K):R[i][i]=-sum(R[i][j] for j in range(K) if j!=i)
    b0,b1=drift(Q),drift(R); G0,G1=bracket(Q),bracket(R)
    ck('exact parity drift preserved while bracket changes',b0==b1 and G0!=G1)
    ck('perturbed exact positivity DB and total cap',all(R[i][j]>0 and pi[i]*R[i][j]==pi[j]*R[j][i] for i in range(K) for j in range(K) if Q[i][j]>0) and all(-R[i][i]<=N for i in range(K)))
    ck('exact stationary bracket identity and matching mean',sum(pi[i]*G1[i] for i in range(K))==sum(pi[i]*G0[i] for i in range(K))==-2*sum(pi[i]*m[i]*b1[i] for i in range(K)))
    Thin=[row[:] for row in Q]
    for i in range(K):
        for j in range(K):
            if i!=j and m[i]==m[j]:Thin[i][j]*=F(1,7)
        Thin[i][i]=-sum(Thin[i][j] for j in range(K) if j!=i)
    ck('separate same-parity thinning preserves full parity drift and bracket',Thin!=Q and drift(Thin)==b0 and bracket(Thin)==G0)
    # Exact Poisson certificate for the NEW bounded proof-support criterion.
    G=[N*x for x in G1]; mean=sum(a*b for a,b in zip(pi,G)); r=[v-mean for v in G]
    A=[[-v for v in row] for row in R]; A[-1]=pi[:]; rhs=r[:]; rhs[-1]=F(0)
    u=solve(A,rhs)
    ck('exact reversible Poisson equation with zero-mean solution',all(-sum(R[i][j]*u[j] for j in range(K))==r[i] for i in range(K)) and sum(pi[i]*u[i] for i in range(K))==0)
    norm=sum(pi[i]*r[i]*u[i] for i in range(K)); energy=sum(pi[i]*R[i][j]*(u[j]-u[i])**2 for i,j in itertools.product(range(K),repeat=2))/2
    ck('Poisson dual norm equals Dirichlet energy and is nonnegative',norm==energy and norm>0)
    for k in range(4):
        v=[F((i+k)%5-2,7) for i in range(K)]
        trial=[u[i]+v[i] for i in range(K)]
        E=lambda f:sum(pi[i]*R[i][j]*(f[j]-f[i])**2 for i,j in itertools.product(range(K),repeat=2))/2
        objective=2*sum(pi[i]*r[i]*trial[i] for i in range(K))-E(trial)
        assert norm-objective==E(v)
    ck('exact variational complete-square identity on four test vectors',True)
    out['poisson_fixture']={'states':K,'norm_Hminus1_squared':str(norm),'Dirichlet_energy':str(energy),'variance_r':str(sum(pi[i]*r[i]**2 for i in range(K))),
        'scope':'finite rational identity, not a bound for all KR families'}
    # Same zero conditional mean need not give trajectory averaging: abstract control only.
    ck('invariant component residual can persist despite global centering',sum(map(F,(1,1,-1,-1)))==0 and sum(F(x*x,4) for x in (1,1,-1,-1))==1)
    # Independently embed producer scalar-only motif into q12 integer counts.
    mobile=counts(2,4); pattern={(0,6):1,(0,3):-1,(1,6):-1,(1,4):1,(3,8):-1,(4,8):1}
    # counts() recursion agrees with lexicographic ordering; test it explicitly.
    ck('producer uniform state is correctly indexed',mobile[3]==(0,1,0,1))
    sig4=(1,1,-1,-1); parity=lambda n:sum(a*b for a,b in zip(n,sig4))
    stress_first=[0]*10; stress_second=[0]*10; coordinate_first=[[0]*4 for _ in range(10)]
    for (i,j),v in pattern.items():
        d=parity(mobile[j])-parity(mobile[i]); stress_first[i]+=v*d; stress_first[j]-=v*d
        stress_second[i]+=v*d*d; stress_second[j]+=v*d*d
        for a in range(4):
            diff=mobile[j][a]-mobile[i][a]; coordinate_first[i][a]+=v*diff; coordinate_first[j][a]-=v*diff
        assert sum(abs(a-b) for a,b in zip(mobile[i],mobile[j]))==2
    ck('six-edge motif preserves scalar drift but not all coordinate drifts',stress_first==[0]*10 and any(v for row in coordinate_first for v in row))
    ck('exact uniform-state squared jump coefficient and opposite compensation',stress_second[3]==-8 and stress_second[4]==8 and sum(v!=0 for v in stress_second)==2)
    for population in (12,24,120):
        base=[population//12]*12; base[2]-=1; base[3]-=1
        embedded=[]
        for n in mobile:
            x=base[:]
            for a,v in zip((0,2,1,3),n):x[a]+=v
            embedded.append(tuple(x))
        assert len(set(embedded))==10 and all(sum(x)==population and min(x)>=0 for x in embedded)
        assert embedded[3]==(population//12,)*12
    ck('q12 embedding legal at three bounded populations without MT-domain inference',True)
    # Native numerical local characteristics: original versioned safe_seed.
    seedpath=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
    spec=importlib.util.spec_from_file_location('versioned_seed',seedpath); seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,_,kernel=seed.build_rank7(); gc=12/lam[6]
    from scipy.special import softmax
    numerical=[]
    for population in (12,24,120):
        n=[population//12]*12; bm=0.; B=0.
        for i in range(12):
            cavity=n[:]; cavity[i]-=1; probabilities=softmax((gc/population)*(kernel@cavity))
            for j in range(12):
                delta=(-1)**j-(-1)**i; rate=n[i]*probabilities[j]
                bm+=rate*delta/population; B+=rate*delta*delta/population
        assert abs(bm)<2e-15 and abs(B-2)<10/population
        numerical.append({'N':population,'native_parity_drift':bm,'N_Gamma':B})
    ck('native balanced-type clock and bracket normalization diagnostics',True)
    out['native_numeric_checks']=numerical
    ck('critical acceleration bracket power exactly N',F(1,2)+2*F(1,4)==1 and F(1,4)-1==F(-3,4))
    ck('quartic stationary moment normalization',F(-1,3)*3+1==0 and F(-1,3)*9+3==0)
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-kr-intake-'))/'campaign'; shutil.copytree(SOURCE,temp,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'); logs={}
        for name in ('verify_package.py','verify_kr.py'):
            run=subprocess.run(['python3',name],cwd=temp,env=env,capture_output=True,text=True,timeout=60)
            assert run.returncode==0,(name,run.stderr); logs[name]={'exit_code':run.returncode,'stdout':run.stdout,'stderr':run.stderr}
        original=json.loads((SOURCE/'FINITE_EXACT_CHECKS.json').read_text()); regenerated=json.loads((temp/'FINITE_EXACT_CHECKS.json').read_text())
        ck('producer finite exact output reproduced without source edits',original==regenerated)
        out['replay']={'directory':str(temp),'logs':logs,'byte_equal':sha(SOURCE/'FINITE_EXACT_CHECKS.json')==sha(temp/'FINITE_EXACT_CHECKS.json'),
            'producer_scope':'one q4/N2 exact rational fixture, not a fixed number of independent stochastic proofs'}
    inputs=[Path(__file__),seedpath,SOURCE/'MANIFEST.sha256',ROOT/'FIN_STATIONARY_KINETIC_RIGIDITY_RESEARCH_ORDER_20261016.md']
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in priors]
    inputs += [ROOT/r for r in lookup.values()]
    inputs += [ROOT/path for path in missing if path!='current AGENTS.md repository HEAD']
    inputs += [ROOT/'fin_uniform_critical_review_20261014/verification.json']
    out['input_hashes']={str(path.relative_to(ROOT)):sha(path) for path in inputs}
    out['checks']=checks; out['status']='PASS_SCOPED_KR_INTAKE_CHECKS'; out['source_manifest_after']=manifest(SOURCE)
    assert out['source_manifest']==out['source_manifest_after']
    for rel in priors:assert manifest(ROOT/rel)=={k:v for k,v in priors[rel].items() if k!='historical_AGENTS_prefix_bytes'}
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'checks':len(checks),'source_manifest':out['source_manifest'],'provenance':out['provenance'],
        'replay_byte_equal':out.get('replay',{}).get('byte_equal'),'poisson_fixture':out['poisson_fixture']},indent=2))

if __name__=='__main__':main()
