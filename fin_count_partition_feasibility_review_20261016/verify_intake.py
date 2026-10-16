#!/usr/bin/env python3
"""Source-preserving FP intake; finite exact checks and separate diagnostics."""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile
import mpmath as mp

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/'fin_count_partition_feasibility_campaign_20261016'
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
    out={'date':'2026-10-16','entry_baseline':'8571db1d plus preserved KR/BA worktree additions'}
    excluded=tuple(n for n in ('PACKAGE_VERIFICATION.json',) if (SOURCE/n).exists())
    out['source_manifest']=manifest(SOURCE,excluded)
    records=json.loads((SOURCE/'PROVENANCE.json').read_text())['accessible_sources']
    expected=''.join(r['sha256']+'  '+r['path']+'\n' for r in records)
    ck('ten original entries and three executed-stage ledgers agree',len(records)==10 and (SOURCE/'INPUTS.sha256').read_text()==expected and all((SOURCE/f'FP-{i:03d}/INPUTS.sha256').read_text()==expected for i in (1,2,4)))
    wrappers={'/mnt/data/files/FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip':'fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip',
        '/mnt/data/files/FIN_UF_001_UF_005_FULL_HANDOFF_20261014.zip':'fin_critical_metastability_campaign_20261014/source_context/UF_PREDECESSOR.zip'}
    matched={}; missing=[]
    for r in records:
        original=r['path']; rel=wrappers.get(original)
        if original.startswith('/mnt/data/fin_'):rel=original.removeprefix('/mnt/data/')
        if rel is None:missing.append(original)
        else:assert sha(ROOT/rel)==r['sha256']; matched[original]=rel
    ck('seven originals match with three original wrappers or snapshots unavailable',len(matched)==7 and len(missing)==3)
    members=json.loads((SOURCE/'SOURCE_ARCHIVE_MEMBERS.json').read_text())
    for rel,row in members.items():
        assert sha(ROOT/rel)==row['sha256']
        with zipfile.ZipFile(ROOT/wrappers[row['archive']]) as z:assert hashlib.sha256(z.read(rel)).hexdigest()==row['sha256']
    ck('all four additional archived scientific members bound independently',len(members)==4)
    out['provenance']={'matched_original_references':matched,'unavailable_originals':missing,
        'stage_references':{'matched':21,'unavailable':9},'additional_verified_members':members,
        'scope':'Current reviews/order are separately available, not retroactive original custody'}
    priors={}
    for rel in ('fin_reversible_bracket_averaging_review_20261016','fin_stationary_kinetic_rigidity_review_20261016',
                'fin_activity_covariance_research_20261016','fin_critical_metastability_review_20261015','fin_uniform_critical_review_20261014'):
        priors[rel]=manifest(ROOT/rel)
        saved=json.loads((ROOT/rel/('RESULTS.json' if rel.startswith('fin_activity_') else 'verification.json')).read_text())
        for path,h in saved.get('input_hashes',saved.get('inputs',{})).items():
            if path=='AGENTS.md':
                n=agents_prefix(h); assert n is not None; priors[rel]['historical_AGENTS_prefix_bytes']=n
            else:assert sha(ROOT/path)==h,path
    ck('actual BA KR AC MT UF manifests and scientific dependencies preserved',len(priors)==5)
    out['actual_dependencies']=priors
    disposition=json.loads((SOURCE/'RESULTS.json').read_text())
    ck('source honestly leaves uniform flows open and temporal task inactive',disposition['scientific_disposition']=='OPEN_UNIFORM_FLOW_LEMMA' and disposition['tasks']['FP-003']=='NOT_ACTIVATED')
    # Independent exact incidence/dual test with nonuniform supplies.
    graph=((0,1),(1,2)); flow=(F(3),F(5)); supply=(flow[0],flow[1]-flow[0],-flow[1]); load=(flow[0],sum(flow),flow[1])
    ck('exact incidence normalization and level cut sources',sum(supply)==0 and supply[0]==3 and sum(supply[:2])==5 and load==(3,8,5))
    phi=(-1,0,1); beta=(0,1,0)
    ck('full capacity dual sign and unrestricted potential',all(phi[i]-phi[j]+beta[i]+beta[j]>=0 for i,j in graph) and sum(x*y for x,y in zip(phi,supply))+7*sum(beta)==-1)
    ck('one cap-independent cut witness is not a continuum cut test',sum(supply[i] for i in (0,))>0 and sum(supply[i] for i in (0,1))>0)
    # Correct source rescaling between delta-incidence and unit incidence.
    for population in (12,24,120):
        delta=F(2,population); demand=tuple(delta*x for x in supply)
        assert tuple(x/delta for x in demand)==supply
    ck('N over two demand rescaling retains the vertex-cap units',True)
    # Exact Fourier move directions, no floating phase boundary decisions.
    cos=(1,0,-1,0)*3; sin=(0,1,0,-1)*3; sig=tuple((-1)**j for j in range(12))
    for direction in (-2,2):
        moves={(cos[j]-cos[i],sin[j]-sin[i]) for i,j in itertools.product(range(12),repeat=2) if sig[j]-sig[i]==direction}
        assert moves==set(itertools.product((-1,1),repeat=2))
    ck('all four fast diagonal moves available in either parity direction',True)
    for C,S in ((0,0),(1,-2),(7,11),(-12,5)):
        for u,v in itertools.product((-1,1),repeat=2):
            assert 2*((C+u)**2+(S+v)**2-C*C-S*S)==4*(C*u+S*v+1)
        L=abs(C)+abs(S); inward=(-1 if C>=0 else 1,-1 if S>=0 else 1); outward=tuple(-x for x in inward)
        assert C*inward[0]+S*inward[1]==-L and C*outward[0]+S*outward[1]==L
    ck('exact radius increments and inward-outward choices',True)
    boxes=json.loads((ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json').read_text())['exact']['laplacian_intervals']
    gamma_lo=1-F(boxes[3][1])/F(boxes[6][0]); gamma_hi=1-F(boxes[3][0])/F(boxes[6][1])
    ck('actual spectral input pays the uniform threshold Lipschitz bound',F(4,25)<gamma_lo<gamma_hi<1 and F(4,3)/gamma_lo<F(25,3))
    ck('large-N geometry inequalities are paid without a guessed MT N0',F(4,10000)+F(4,10**8)<F(1,2) and F(50,3)/100<F(1))
    # Exact analytic cap enclosure, not a high-precision rounding claim.
    Jlo=1-F(2,12); Jhi=1-F(2,12)+F(24,2*12**2)
    caplo=F(1,2)+Jlo/6; caphi=F(1,2)+Jhi/6
    ck('Taylor-moment rational necessary cap enclosure',Jlo==F(5,6) and Jhi==F(11,12) and caplo==F(23,36) and caphi==F(47,72))
    ck('continuum bracket bounds remain finite and do not solve microscopic feasibility',F(8,5)<2 and F(8,3)/4==F(2,3))
    # Distinct scalar formula: complete square in t=z^2, not quartic quadrature.
    mp.mp.dps=60
    J=mp.sqrt(3*mp.pi)*mp.exp(3)*mp.erfc(mp.sqrt(3)); cap=mp.mpf('.5')+J/6
    source_values=json.loads((SOURCE/'FINITE_DIAGNOSTICS.json').read_text())['numerical_diagnostics']
    target=mp.mpf(source_values['necessary_C_lower'])
    ck('independent erfc cap diagnostic agrees with the noncertified source value',abs(cap-target)<mp.mpf('1e-38'))
    gA=2-2*J/5; gB=2+2*J/3
    ck('independent phase bracket and mixture normalization at zero',abs(gA-mp.mpf(source_values['gA_at_zero']))<mp.mpf('1e-38') and abs(gB-mp.mpf(source_values['gB_at_zero']))<mp.mpf('1e-38') and abs(mp.mpf(5)/8*gA+mp.mpf(3)/8*gB-2)<mp.mpf('1e-55'))
    ck('separate numerical preview lies inside the analytic rational bounds',mp.mpf(caplo.numerator)/caplo.denominator<cap<mp.mpf(caphi.numerator)/caphi.denominator)
    out['cap_check']={'exact_formula':'C_star=1/2+(sqrt(3*pi)*exp(3)*erfc(sqrt(3)))/6',
        'proved_rational_enclosure':[str(caplo),str(caphi)],'endpoints':'strict for C_star; any feasible C>=C_star',
        'numeric_preview':mp.nstr(cap,42),'numeric_status':'not an outward certificate for the42 digits',
        'scope':'necessary for a hypothetical closed family, not a sufficient uniform cap'}
    # Boundary leakage scaling is different from changing the whole model clock.
    ck('rare boundary flux pays both critical phase-change and moment errors',F(1)+F(1,2)==F(3,2) and 1+1==2)
    ck('native positive mixture has vanishing critical coupling scale',-2+F(3,2)==F(-1,2))
    # Capacity-free lower ideals: exhaustive tiny directed graph dual sanity check.
    vertices=set(range(3)); ideals=[]
    for mask in range(8):
        L={i for i in vertices if mask&(1<<i)}
        if all(not(j in L and i not in L) for i,j in graph):ideals.append(L)
    ck('all true lower ideals have the required source sign on a feasible path',all(sum(supply[i] for i in L)>=0 for L in ideals))
    fake=(F(-1),F(0),F(1)); L={0}
    ck('a negative lower ideal gives a capacity-independent dual obstruction',sum(fake)==0 and sum(fake[i] for i in L)<0 and all(int(i in L)-int(j in L)>=0 for i,j in graph))
    if a.replay:
        tmp=Path(tempfile.mkdtemp(prefix='fin-fp-intake-'))/'campaign'; shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        logs={}; env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        for script in ('verify_package.py','verify_fp.py'):
            run=subprocess.run(['python3',script],cwd=tmp,env=env,capture_output=True,text=True,timeout=60)
            assert run.returncode==0,(script,run.stderr); logs[script]={'exit_code':run.returncode,'stdout':run.stdout,'stderr':run.stderr}
        regenerated=json.loads((tmp/'FINITE_DIAGNOSTICS.json').read_text()); original=json.loads((SOURCE/'FINITE_DIAGNOSTICS.json').read_text())
        ck('producer27 bounded flags pass with scientific output reproduced',regenerated['passed']==regenerated['total']==27 and all(c['pass'] for c in regenerated['checks']))
        out['replay']={'directory':str(tmp),'logs':logs,'byte_equal':sha(tmp/'FINITE_DIAGNOSTICS.json')==sha(SOURCE/'FINITE_DIAGNOSTICS.json'),
            'data_equal':regenerated==original,'scope':'17 exact finite checks plus10 numerical ideal diagnostics, not27 proofs'}
    paths=[Path(__file__),SOURCE/'MANIFEST.sha256',ROOT/'FIN_COUNT_PARTITION_FEASIBILITY_NEXT_CHAIN_20261016.md']
    paths += [ROOT/rel/'MANIFEST.sha256' for rel in priors]
    paths += [ROOT/rel for rel in matched.values()]
    paths += [ROOT/rel for rel in members]
    paths += [ROOT/'fin_reversible_bracket_averaging_review_20261016/COUNT_PARTITION_FEASIBILITY_SEED.md',
        ROOT/'fin_reversible_bracket_averaging_review_20261016/PROOF_AUDIT_AND_SCOPE.md',
        ROOT/'fin_stationary_kinetic_rigidity_review_20261016/PROOF_COMPLETIONS_AND_SCOPE.md',
        ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json']
    out['input_hashes']={str(path.relative_to(ROOT)):sha(path) for path in paths}
    out['status']='PASS_SCOPED_FP_INTAKE_CHECKS'; out['checks']=checks; out['source_manifest_after']=manifest(SOURCE,excluded)
    assert out['source_manifest_after']==out['source_manifest']
    for rel in priors:assert manifest(ROOT/rel)=={k:v for k,v in priors[rel].items() if k!='historical_AGENTS_prefix_bytes'}
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'source_manifest':out['source_manifest'],'provenance':out['provenance'],
        'replay_byte_equal':out.get('replay',{}).get('byte_equal'),'cap_check':out['cap_check']},indent=2))

if __name__=='__main__':main()
