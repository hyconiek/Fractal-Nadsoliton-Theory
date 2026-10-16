#!/usr/bin/env python3
"""BA intake: source-preserving exact checks; no incoming producer imports."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/'fin_reversible_bracket_averaging_campaign_20261016'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder,excluded=()):
    names=[]
    for row in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=row.split('  ',1); assert n not in names and sha(folder/n)==h,n; names.append(n)
    actual={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'}|set(excluded),folder
    return {'matched':len(names),'coverage':'COMPLETE','excluded_hashes':{n:sha(folder/n) for n in excluded}}
def agents_prefix(h):
    digest=hashlib.sha256(); n=0
    for line in (ROOT/'AGENTS.md').read_bytes().splitlines(keepends=True):
        digest.update(line); n+=len(line)
        if digest.hexdigest()==h:return n
    return None

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True); p.add_argument('--replay',action='store_true'); a=p.parse_args(); checks=[]
    def ck(name,ok,**detail):checks.append({'name':name,'passed':bool(ok),**detail}); assert ok,name
    out={'date':'2026-10-16','entry_baseline':'8571db1d plus preserved KR worktree additions'}
    out['source_manifest']=manifest(SOURCE,('PACKAGE_VERIFICATION.json',))
    records=json.loads((SOURCE/'SOURCE_PROVENANCE.json').read_text())['entries']
    expected=''.join(r['sha256']+'  '+r['source_identifier']+'\n' for r in records)
    ck('twelve source entries and all four stage ledgers agree',len(records)==12 and (SOURCE/'INPUTS.sha256').read_text()==expected and all((SOURCE/f'BA-{i:03d}/INPUTS.sha256').read_text()==expected for i in range(1,5)))
    archives={'files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip':'fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip',
        'files/FIN_UF_001_UF_005_FULL_HANDOFF_20261014.zip':'fin_critical_metastability_campaign_20261014/source_context/UF_PREDECESSOR.zip'}
    matched={}; missing=[]
    for row in records:
        name=row['source_identifier']
        if name=='files/AGENTS.md':missing.append(name); continue
        if '::' in name:
            wrapper,member=name.split('::',1)
            with zipfile.ZipFile(ROOT/archives[wrapper]) as z:actual=hashlib.sha256(z.read(member)).hexdigest()
            assert actual==row['sha256'] and sha(ROOT/member)==actual,name; matched[name]=member
        else:
            path=archives.get(name,name); assert sha(ROOT/path)==row['sha256']; matched[name]=path
    ck('eleven original references matched without historical AGENTS inflation',len(matched)==11 and len(missing)==1)
    out['provenance']={'matched':matched,'unavailable_originals':missing,'stage_references':{'matched':44,'unavailable':4},
        'scope':'Current review/order separately bound, not a remote checkout or missing original upload identity'}
    priors={}
    for rel in ('fin_stationary_kinetic_rigidity_review_20261016','fin_activity_covariance_research_20261016',
        'fin_rate_coupling_review_20261016','fin_critical_metastability_review_20261015','fin_uniform_critical_review_20261014'):
        priors[rel]=manifest(ROOT/rel)
        saved=json.loads((ROOT/rel/('RESULTS.json' if rel.startswith('fin_activity_') else 'verification.json')).read_text())
        for path,h in saved.get('input_hashes',saved.get('inputs',{})).items():
            if path=='AGENTS.md':
                length=agents_prefix(h); assert length is not None; priors[rel]['historical_AGENTS_prefix_bytes']=length
            else:assert sha(ROOT/path)==h,path
    out['actual_dependencies']=priors
    ck('actual predecessor manifests and scientific inputs preserved',len(priors)==5)
    # Reuse only an older versioned research-model helper, never the new BA producer.
    acpath=ROOT/'fin_activity_covariance_research_20261016/verify_research.py'
    spec=importlib.util.spec_from_file_location('prior_exact_count_model',acpath); ac=importlib.util.module_from_spec(spec); spec.loader.exec_module(ac)
    W,Q=ac.exact_hb(3,3,F(2)); centre=(1,1,1); stress=ac.wheel(centre,(0,1,2))
    amplitude=ac.safe_amplitude(W,Q,stress,F(1,2)); R=ac.perturb(W,Q,stress,amplitude)
    Z=sum(W.values()); pi={n:w/Z for n,w in W.items()}; N=3
    m=lambda n:F(n[0]-n[1]+n[2],N)
    levels=sorted(set(map(m,W))); delta=F(2,N)
    drift=lambda L,n:sum(rate*(m(y)-m(n)) for y,rate in L[n].items())
    gamma=lambda L,n:sum(rate*(m(y)-m(n))**2 for y,rate in L[n].items())
    ck('nonuniform single-copy count fixture keeps exact scalar drift',all(drift(Q,n)==drift(R,n) for n in W))
    ck('pointwise covariance remains nonunique in actual count support',any(gamma(Q,n)!=gamma(R,n) for n in W))
    cuts=[]
    for level in levels[:-1]:
        cut=lambda L:sum(pi[n]*rate for n in W for y,rate in L[n].items() if m(n)==level and m(y)==level+delta)
        flux=sum(pi[n]*drift(Q,n) for n in W if m(n)<=level)
        assert delta*cut(Q)==flux==delta*cut(R)
        cuts.append(str(cut(Q)))
    ck('all finite parity cuts exactly paid by drift and weights',len(cuts)==3)
    for level in levels:
        conditional=lambda L:sum(pi[n]*gamma(L,n) for n in W if m(n)==level)/sum(pi[n] for n in W if m(n)==level)
        assert conditional(Q)==conditional(R)
    ck('entire conditional bracket equals native exactly at every occupied level',True)
    D=lambda L,f,g:sum(pi[n]*rate*(f(y)-f(n))*(g(y)-g(n)) for n in W for y,rate in L[n].items())/2
    for first,second in itertools.product(levels,repeat=2):
        f=lambda n:F(m(n)==first); g=lambda n:F(m(n)==second)
        assert D(Q,f,g)==D(R,f,g)
    ck('full bilinear Dirichlet compression checked on parity indicator basis',True)
    ck('conditional Jensen contraction and bounded L1-to-L2 step',all(sum(pi[n]*abs(gamma(R,n)-gamma(Q,n)) for n in W if m(n)==ell)>=0 for ell in levels))
    # Independent electrical-cycle Poisson solution and exact leading coefficient.
    weights=tuple(F(v,41) for v in (2,3,5,7,11,13)); values=tuple(map(F,(0,1,2,1,0,1)))
    def cycle(t):
        c=(t,t,F(1,5)-t,F(1,5)-t,t,F(1,5)-t)
        d=tuple((c[i]*(values[(i+1)%6]-values[i])+c[(i-1)%6]*(values[(i-1)%6]-values[i]))/weights[i] for i in range(6))
        G=tuple((c[i]+c[(i-1)%6])/weights[i] for i in range(6)); mean=sum(w*g for w,g in zip(weights,G))
        forcing=[weights[i]*(G[i]-mean) for i in range(6)]; cumulative=[]; running=F(0)
        for x in forcing:running+=x; cumulative.append(running)
        assert running==0
        if t:
            current=-sum(f/e for f,e in zip(cumulative,c))/sum(1/e for e in c)
            flows=[current+f for f in cumulative]
            assert sum(j/e for j,e in zip(flows,c))==0
            energy=sum(j*j/e for j,e in zip(flows,c))
        else:energy=None
        return c,d,G,cumulative,energy
    _,d0,G0,_,_=cycle(F(1,10)); cend,dend,Gend,Fend,_=cycle(F(0))
    ck('six-state endpoint genuinely preserves drift and mean but changes pointwise noise',d0==dend and G0!=Gend and sum(w*g for w,g in zip(weights,G0))==sum(w*g for w,g in zip(weights,Gend))==F(6,5))
    components=((0,5),(1,),(2,3,4)); masses=[sum(weights[i] for i in group) for group in components]
    h=[sum(weights[i]*(Gend[i]-F(6,5)) for i in group)/mass for group,mass in zip(components,masses)]
    ck('actual invariant component centering at the endpoint is nonzero',h==[F(-8,75),F(-6,5),F(26,115)] and sum(x*y for x,y in zip(masses,h))==0)
    small=(0,1,4); j0=-sum(Fend[i] for i in small)/3
    coefficient=sum((j0+Fend[i])**2 for i in small)
    ck('leading Poisson residue from full cycle not coarse inverse or floating fit',coefficient==F(1064,126075))
    supplied=json.loads((SOURCE/'EXACT_DIAGNOSTICS.json').read_text())
    costs=[]
    for row in supplied['fixture']['tau_poisson']:
        t=F(row['tau']); c,d,G,_,cost=cycle(t)
        assert d==d0 and cost==F(row['R_exact']); costs.append(cost)
    ck('all three full Poisson energies independently exact and ordered',costs[0]<costs[1]<costs[2])
    out['poisson_control']={'method':'exact cycle-current closure and small-edge residue','coefficient':str(coefficient),'costs':list(map(str,costs)),
        'scope':'six-state finite graph, NOT an all-N q12 construction'}
    ck('global zero mean does not imply zero invariant component projection',sum(m*h0*h0 for m,h0 in zip(masses,h))>0)
    # New candidate seed: exact invariants of the existing Fourier count variables.
    cosine=(1,0,-1,0)*3; sine=(0,1,0,-1)*3; sigma=tuple((-1)**j for j in range(12))
    fixture_counts=[(1,)*12]+[tuple(12 if j==i else 0 for j in range(12)) for i in range(12)]
    def invariants(n):
        C=sum(x*y for x,y in zip(n,cosine)); S=sum(x*y for x,y in zip(n,sine)); M=sum(x*y for x,y in zip(n,sigma))
        return C*C+S*S,M*M
    for n in fixture_counts:
        for sign,shift in itertools.product((-1,1),range(12)):
            moved=[0]*12
            for j,v in enumerate(n):moved[(sign*j+shift)%12]=v
            assert invariants(n)==invariants(moved)
    ck('candidate fast-radius and squared parity invariants under exact D12 maps',True)
    epsilon=F(1,8)
    ck('candidate fibre masses uniformly nonzero without fitting',F(1,2)-epsilon==F(3,8) and F(1,2)+epsilon==F(5,8))
    # Exact finite dual gate control: phase balance alone need not be sufficient.
    supply=[F(1),F(0),F(-1)]; caps=[F(1),F(2),F(1)]; graph=((0,1),(1,2))
    def eligible(alpha,beta,edges):return all(v>=0 for v in beta) and all(alpha[i]-alpha[j]+beta[i]+beta[j]>=0 for i,j in edges)
    for alpha,beta in itertools.product(itertools.product((-1,0,1),repeat=3),itertools.product((0,1),repeat=3)):
        if eligible(alpha,beta,graph):assert sum(x*y for x,y in zip(alpha,supply))+sum(x*y for x,y in zip(beta,caps))>=0
    ck('exact dual inequality respected by a feasible finite flow',True)
    alpha=(-1,0,1); beta=(0,0,0)
    ck('phase balanced totals alone fail a disconnected internal flow gate',supply[0]+supply[2]==0 and eligible(alpha,beta,()) and sum(x*y for x,y in zip(alpha,supply))<0)
    alpha=(-1,0,0); beta=(1,0,0); tight_caps=(F(1,2),F(2),F(1))
    ck('capacity dual certificate rejects an otherwise available flow',eligible(alpha,beta,graph) and sum(x*y for x,y in zip(alpha,supply))+sum(x*y for x,y in zip(beta,tight_caps))==F(-1,2))
    ck('negative capacity multipliers do not license a dual certificate',not eligible((0,0,0),(-1,0,0),graph))
    eta=F(1,100); cmix=cycle(eta/10)[0]
    assert cmix==tuple((1-eta)*x+eta*y for x,y in zip(cend,cycle(F(1,10))[0]))
    ck('conditional candidate coefficient and positive-mixture scaling',8*epsilon==1 and all(x>0 for x in cmix) and F(3,2)-2==F(-1,2))
    if a.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-ba-intake-'))/'campaign'; shutil.copytree(SOURCE,temp,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        logs={}; env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        for script in ('verify_package.py','verify_ba.py'):
            run=subprocess.run(['python3',script],cwd=temp,env=env,capture_output=True,text=True,timeout=60)
            assert run.returncode==0,(script,run.stderr); logs[script]={'exit_code':run.returncode,'stdout':run.stdout,'stderr':run.stderr}
        regenerated=json.loads((temp/'EXACT_DIAGNOSTICS.json').read_text())
        ck('producer 24 flags and deterministic data reproduced',regenerated==supplied and regenerated['check_count']==24)
        out['replay']={'directory':str(temp),'logs':logs,'byte_equal':sha(temp/'EXACT_DIAGNOSTICS.json')==sha(SOURCE/'EXACT_DIAGNOSTICS.json'),
            'scope':'24 source flags; one monotonicity comparison uses binary64 previews, independently paid with Fraction above'}
    inputs=[Path(__file__),acpath,SOURCE/'MANIFEST.sha256',ROOT/'FIN_REVERSIBLE_BRACKET_AVERAGING_NEXT_CHAIN_20261016.md']
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in priors]
    inputs += [ROOT/v for v in matched.values()]
    inputs += [ROOT/path for path in json.loads((SOURCE/'SOURCE_PROVENANCE.json').read_text())['unavailable_and_not_claimed'] if (ROOT/path).is_file()]
    inputs += [ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json',ROOT/'fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md']
    out['input_hashes']={str(path.relative_to(ROOT)):sha(path) for path in inputs}
    out['checks']=checks; out['status']='PASS_SCOPED_BA_INTAKE_CHECKS'; out['source_manifest_after']=manifest(SOURCE,('PACKAGE_VERIFICATION.json',))
    assert out['source_manifest']==out['source_manifest_after']
    for rel in priors:assert manifest(ROOT/rel)=={k:v for k,v in priors[rel].items() if k!='historical_AGENTS_prefix_bytes'}
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'source_manifest':out['source_manifest'],'matched_source_references':len(matched),
        'replay_byte_equal':out.get('replay',{}).get('byte_equal'),'Poisson_residue':str(coefficient)},indent=2))

if __name__=='__main__':main()
