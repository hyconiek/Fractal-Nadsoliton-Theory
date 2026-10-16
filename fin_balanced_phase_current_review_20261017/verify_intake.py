#!/usr/bin/env python3
"""Independent bounded PC intake, preserving all incoming source bytes."""
import argparse
from collections import Counter,defaultdict
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

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent;SOURCE=ROOT/'fin_balanced_phase_current_campaign_20261016'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj
def types(n,q):
    if q==1:yield(n,);return
    for first in range(n+1):
        for tail in types(n-first,q-1):yield(first,)+tail
def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--replay',action='store_true');args=p.parse_args()
    helper=module(ROOT/'fin_rbcut_weighted_degree_review_20261016/verify_intake.py','accepted_hash_helpers');checks=[]
    def ck(name,condition,**detail):checks.append({'name':name,'passed':bool(condition),**detail});assert condition,name
    out={'date':'2026-10-17','entry_HEAD':'4aa415bc','source':'fin_balanced_phase_current_campaign_20261016'}
    out['source_manifest']=helper.manifest(SOURCE,('check_package.py','PACKAGE_VERIFICATION.json'))
    pk=json.loads((SOURCE/'PACKAGE_VERIFICATION.json').read_text());assert pk['check_package_sha256']==sha(SOURCE/'check_package.py')
    ck('all26 scientific entries and separately excluded checker metadata bound',out['source_manifest']['matched']==26)
    aliases={
        'FP_CONDUCTANCE':ROOT/'fin_count_partition_feasibility_campaign_20261016/FP-002/CONDUCTANCE_CONSTRUCTION_OR_DUAL.md',
        'RB_ORBIT_LIFT':helper.STAGES['orbit']/'ORBIT_QUOTIENT_AND_LIFT_THEOREM.md',
        'RB_WEIGHT_RATIO':helper.STAGES['degree']/'UNIFORM_NATIVE_EDGE_ELLIPTICITY.md',
        'SC001_LOCAL_THEOREM':ROOT/'fin_signed_interface_orbit_cut_campaign_20261016/SC-001/INCOMING_COMPENSATION_PROOF_OR_OPEN_LEMMA.md',
        'SC_HANDOFF':ROOT/'fin_signed_interface_orbit_cut_campaign_20261016/HANDOFF.md'}
    matched={};missing={}
    for row in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,name=row.split('  ',1);path=aliases.get(name,ROOT/name)
        if path.is_file():assert sha(path)==h,name;matched[name]=str(path.relative_to(ROOT))
        else:missing[name]={'sha256':h,'matching_current_AGENTS_prefix_bytes':helper.prefix(h) if name=='ARCHIVED_AGENTS_NOT_CURRENT' else None}
    out['original_inputs']={'matched':matched,'unavailable_original_files':missing}
    ck('seven actual input references match while old original snapshot is not invented',len(matched)==7 and set(missing)=={'ARCHIVED_AGENTS_NOT_CURRENT'})
    mtzip=ROOT/'fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip';members=[]
    provenance=json.loads((SOURCE/'SOURCE_PROVENANCE.json').read_text())
    with zipfile.ZipFile(mtzip) as z:
        for name,row in provenance['actual_hashed_inputs'].items():
            if row.get('archive_member'):
                assert hashlib.sha256(z.read(row['archive_member'])).hexdigest()==row['sha256']==sha(ROOT/name);members.append(name)
    out['actual_MT_member_binding']={'wrapper_sha256':sha(mtzip),'members':members,
        'scope':'actual existing wrapper and listed member bytes, not a supplied original PC outer ZIP'}
    ck('both declared MT member contents match actual files and existing archive',len(members)==2)
    previous={}
    for rel in ('fin_signed_interface_orbit_cut_review_20261016','fin_balanced_phase_cut_research_20261016',
                'fin_rbcut_weighted_degree_review_20261016','fin_critical_metastability_review_20261015',
                'fin_stationary_kinetic_rigidity_review_20261016'):
        previous[rel]=helper.manifest(ROOT/rel)
        saved=json.loads((ROOT/rel/('RESULTS.json' if 'research' in rel else 'verification.json')).read_text())
        for path,h in saved.get('input_hashes',{}).items():
            if path=='AGENTS.md':assert helper.prefix(h) is not None
            else:assert sha(ROOT/path)==h,path
    out['actual_dependencies']=previous
    ck('actual PC03 PC04 and current predecessor proof dependencies are bound now',len(previous)==5)
    # New native weight-ratio proof, not the much weaker a0 rate bound.
    pmin=F(103,1250);pmax=F(158,1875);q=F(1,50)+F(18,10000)
    ratio=pmin/(pmax+F(1,10000))*(1-q)
    ck('rational all-local edge weight lower bound independently recomputed',ratio==F(1511319,1581875) and ratio>F(9,10))
    ck('outside-cylinder negative control cannot inherit the bound',F(0)/(pmax+F(1,10000))*(1-q)<F(9,10))
    windows={}
    for P,lo,hi in [('A',F(1,2),F(5,8)),('B',F(3,8),F(1,2))]:
        windows[P]=(lo/(4*hi),hi/(4*lo))
    ck('phase current quarter factor and analytic bounds exact',windows=={'A':(F(1,5),F(5,16)),'B':(F(3,16),F(1,3))})
    ck('finite current window has strict asymptotic margins',F(9,50)<min(x[0] for x in windows.values()) and F(7,20)>max(x[1] for x in windows.values()))
    row=F(357,1000);column=F(9,50);neighbour=F(9,10)
    ck('nondeparture drift and row allowance are conservative',F(1,200)*2+F(24,10000)==F(124,10000) and F(124,20000)<F(7,1000) and row==F(7,20)+F(7,1000))
    right=9*column*neighbour/4
    ck('restricted Hall weighted incidence margin exact',right-row==F(3,400))
    ck('conditional E singleton margin only if nine successors exist',9*column*neighbour-(F(7,20)+1)==F(27,250))
    C=(1,0,-1,0)*3;S=(0,1,0,-1)*3
    for sign in (-2,2):
        cnt=Counter((C[b]-C[a],S[b]-S[a]) for a,b in itertools.product(range(12),repeat=2) if (-1)**b-(-1)**a==sign)
        assert cnt==Counter({d:9 for d in itertools.product((-1,1),repeat=2)})
    ck('q12 diagonal multiplicities both parity directions exact',True)
    # Full exact iid stratum coupling at different N from the producer fixture.
    N=3;states=list(types(N,12));fibres={E:[n for n in states if sum(n[::2])==E] for E in range(N+1)};edge_count=0
    mult=lambda n:F(math.factorial(N),math.prod(math.factorial(v) for v in n))
    for E in range(N):
        incoming=defaultdict(F);O=N-E
        for n in fibres[E]:
            mu=mult(n)/(math.comb(N,E)*6**N);rowmass=F(0)
            for a,b in itertools.product(range(1,12,2),range(0,12,2)):
                if not n[a]:continue
                y=list(n);y[a]-=1;y[b]+=1;v=mu*F(n[a],6*O);incoming[tuple(y)]+=v;rowmass+=v;edge_count+=1
            assert rowmass==mu
        assert all(incoming[y]==mult(y)/(math.comb(N,E+1)*6**N) for y in fibres[E+1])
    ck('all iid q12 N3 stratum coupling rows and columns exact',len(states)==364,
       edges=edge_count,scope='untruncated g0 baseline, not critical MT or phase-H coupling')
    ck('singleton tests cannot substitute for dense Hall',F(1)<=F(3,2) and F(2)>F(3,2))
    ck('nine neighbours and coarse windows still allow a generic dense failure',9*row>9*column,
       scope='complete9by9 logical control, not an actual q12 adverse family; overlap9 violates source hypothesis')
    rr=(F(1),F(1),F(3),F(3));cc=(F(2),F(2),F(4));adj=({0,2},{1,2},{0},{1});deficits=[]
    for mask in range(16):
        picked={i for i in range(4) if mask>>i&1};nbr=set().union(*(adj[i] for i in picked)) if picked else set()
        delta=sum(rr[i] for i in picked)-sum(cc[j] for j in nbr);deficits.append((mask,delta))
    invariant=[v for mask,v in deficits if bool(mask&1)==bool(mask&2) and bool(mask&4)==bool(mask&8)]
    ck('invariant maximum Hall deficit agrees on all finite control cuts',max(v for _,v in deficits)==max(invariant)==2,
       scope='exact small symmetric transport control, not actual dense FIN proof')
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-pc-intake-'));copy=temp/'campaign';shutil.copytree(SOURCE,copy)
        wrapper=temp/'derived_replay_wrapper.zip'
        with zipfile.ZipFile(wrapper,'w',zipfile.ZIP_DEFLATED) as z:
            for path in sorted(copy.rglob('*')):
                if path.is_file():z.write(path,SOURCE.name+'/'+str(path.relative_to(copy)))
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1');logs={}
        for script,extra in [('check_package.py',[str(wrapper)]),('verify_pc.py',[])]:
            run=subprocess.run(['python3',script,*extra],cwd=copy,env=env,capture_output=True,text=True,timeout=60)
            assert run.returncode==0,(script,run.stderr);logs[script]={'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr}
        equal=sha(copy/'verification.json')==sha(SOURCE/'verification.json')
        ck('scientific writer replays byteidentically on copy',equal)
        out['replay']={'directory':str(temp),'logs':logs,'byte_equal':equal,
            'wrapper_scope':'locally derived wrapper to exercise original checker, NOT original uploaded ZIP custody',
            'reported_categories':12,'count_field':'producer hardcoded category label, not twelve independent proofs'}
    inputs=[Path(__file__),SOURCE/'MANIFEST.sha256',ROOT/'AGENTS.md',
        ROOT/'FIN_BALANCED_PHASE_CURRENT_NEXT_CHAIN_20261016.md',
        ROOT/'FIN_PC_PHYSICAL_COMPATIBILITY_NEXT_CHAIN_20261017.md',
        ROOT/'fin_q12_physical_bridge_research_20261017/MANIFEST.sha256']
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in previous]
    out['input_hashes']={str(path.relative_to(ROOT)):sha(path) for path in inputs}
    out['status']='PASS_SCOPED_PC_INTAKE_AND_PHYSICAL_BOUNDARY';out['checks']=checks
    out['source_manifest_after']=helper.manifest(SOURCE,('check_package.py','PACKAGE_VERIFICATION.json'));assert out['source_manifest_after']==out['source_manifest']
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'source_entries':out['source_manifest']['matched'],'original_reference_matches':len(matched)},indent=2))

if __name__=='__main__':main()
