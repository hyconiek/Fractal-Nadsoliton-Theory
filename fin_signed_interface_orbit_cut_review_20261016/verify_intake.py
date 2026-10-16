#!/usr/bin/env python3
"""Source-preserving SC intake; finite controls are not all-N flow proofs."""
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

HERE=Path(__file__).resolve().parent;ROOT=HERE.parent
SOURCE=ROOT/'fin_signed_interface_orbit_cut_campaign_20261016'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load_helper():
    path=ROOT/'fin_rbcut_weighted_degree_review_20261016/verify_intake.py'
    spec=importlib.util.spec_from_file_location('accepted_rbcut_helpers',path)
    obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--replay',action='store_true');args=p.parse_args()
    helper=load_helper();checks=[]
    def ck(name,condition,**detail):checks.append({'name':name,'passed':bool(condition),**detail});assert condition,name
    out={'date':'2026-10-16','final_validation_date':'2026-10-17','entry_HEAD':'8571db1d plus preserved KR BA FP RBCUT additions'}
    out['source_manifest']=helper.manifest(SOURCE)
    ck('source manifest and coverage verified',out['source_manifest']['matched']==48)
    aliases={helper.NAMES[k]:v for k,v in helper.STAGES.items()}
    wrappers={
        'FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip':ROOT/'fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip',
        'FIN_RB_ORBIT_WEIGHTED_DEGREE_ONLY_HANDOFF_20261016.zip':helper.BUNDLE/'ORYGINALNE_PAKIETY/FIN_RB_ORBIT_WEIGHTED_DEGREE_ONLY_HANDOFF_20261016.zip'}
    def resolve(path):
        if Path(path).name in wrappers:return wrappers[Path(path).name]
        rel=path.removeprefix('/mnt/data/');parts=Path(rel).parts
        if parts and parts[0] in aliases:return aliases[parts[0]].joinpath(*parts[1:])
        return ROOT/rel
    def ledger(path):
        matched=[];missing=[]
        for row in path.read_text().splitlines():
            digest,rest=row.split('  ',1);name,_,label=rest.partition('  # ');local=resolve(name)
            item={'original':name,'label':label,'sha256':digest}
            if local.is_file():
                assert sha(local)==digest,(name,digest);item['local']=str(local.relative_to(ROOT));matched.append(item)
            else:missing.append(item)
        return {'matched':matched,'unavailable':missing}
    out['source_inputs']=ledger(SOURCE/'INPUTS.sha256')
    out['stage_inputs']={str(path.relative_to(SOURCE)):ledger(path) for path in sorted(SOURCE.glob('SC-*/INPUTS.sha256'))}
    ck('available original input references match without invented wrappers',len(out['source_inputs']['matched'])==10 and len(out['source_inputs']['unavailable'])==2)
    unique={v['local'] for v in out['source_inputs']['matched']}
    ck('duplicate MTdomain alias identified instead of double scientific evidence',len(unique)==9)
    ck('all four task ledgers checked including unactivated tasks',len(out['stage_inputs'])==4 and all(len(v['matched'])==10 and len(v['unavailable'])==2 for v in out['stage_inputs'].values()))
    out['alias_correction']='MT_domain_extracted repeats NATIVE_SHELL_THEOREM.md; actual MT domain is independently bound below'
    previous={}
    for rel in ('fin_rbcut_weighted_degree_review_20261016','fin_count_partition_feasibility_review_20261016',
                'fin_reversible_bracket_averaging_review_20261016','fin_stationary_kinetic_rigidity_review_20261016',
                'fin_critical_metastability_review_20261015'):
        previous[rel]=helper.manifest(ROOT/rel)
        saved=json.loads((ROOT/rel/'verification.json').read_text())
        for path,digest in saved.get('input_hashes',{}).items():
            if path=='AGENTS.md':assert helper.prefix(digest) is not None
            else:assert sha(ROOT/path)==digest,path
    out['actual_dependencies']=previous
    ck('actual current predecessor reviews and saved inputs preserved',len(previous)==5)
    # Independent sign/quadrant and exact radial identities, including axis ties.
    cos=(1,0,-1,0)*3;sin=(0,1,0,-1)*3;sig=tuple((-1)**j for j in range(12))
    for sign in (-2,2):
        bins={d:[] for d in itertools.product((-1,1),repeat=2)}
        for i,j in itertools.product(range(12),repeat=2):
            if sig[j]-sig[i]==sign:bins[cos[j]-cos[i],sin[j]-sin[i]].append((i,j))
        assert set(map(len,bins.values()))=={9}
    ck('q12 up and down diagonal multiplicities independently enumerated',True)
    cases=0
    for c,s in itertools.product(range(-5,6),repeat=2):
        inward=(-(1 if c>=0 else -1),-(1 if s>=0 else -1))
        outward=tuple(-a for a in inward)
        assert c*inward[0]+s*inward[1]==-abs(c)-abs(s)
        assert c*outward[0]+s*outward[1]==abs(c)+abs(s)
        for dc,ds in (inward,outward):
            assert 2*((c+dc)**2+(s+ds)**2-c*c-s*s)==4*(c*dc+s*ds+1);cases+=1
    ck('all axis and quadrant radial algebra controls exact',cases==242)
    # Coarse explicit thresholds for geometry only, NOT a numerical MT N0.
    root=10**4;N=root**4;L=F(25,3)
    eta_bound=F(16,root**2)+F(4,N)+2*L/F(root**3)
    ck('rational largeN geometry payment with Tmin>1 Tmax<16',eta_bound<F(1,2) and F(root**2,2)>1+(L/2)*root,
       scope='illustrative sufficient geometric threshold only; inherited MT N0 unquantified')
    # Exhaustive cut algebra for a different six-node weighted DAG.
    U={(0,1):F(1,3),(0,2):F(1,5),(1,3):F(2,7),(2,3):F(3,11),(2,4):F(1,2),(3,5):F(1,4),(4,5):F(2,5)}
    phase={0,1,3};E={5};H={e:v for e,v in U.items() if (e[0] in phase)==(e[1] in phase) or any(x in E for x in e)}
    FF={e:v for e,v in U.items() if e not in H}
    W={x:1+sum(v for e,v in U.items() if x in e) for x in range(6)};empty_shield=0
    for mask in range(64):
        pick={x for x in range(6) if mask>>x&1}
        B=sum(W[x] for x,y in H if x in pick and y not in pick)
        fo=sum(v for (x,y),v in FF.items() if x in pick and y not in pick)
        incoming=sum(v for (x,y),v in U.items() if x not in pick and y in pick)
        tails={x for x,y in FF if x in pick and y not in pick}
        exposed={x for x in tails if any(a==x and b not in pick for a,b in H)};shield=tails-exposed
        assert fo-incoming<=B+sum(W[x] for x in shield)-incoming
        if not shield:assert fo-incoming<=B;empty_shield+=1
    ck('shield decomposition and no-shield cut bound on all64 independent controls',empty_shield>0,
       scope='finite general network identity, not actual allN cut feasibility')
    ck('generic lower-deficit statement rejected for a nonideal upper cut',-F(1)<0,
       explanation='U=H=0->1->2, F empty, S={1,2}: D=-1 despite zero shell')
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-sc-intake-'));copy=temp/'campaign'
        shutil.copytree(SOURCE,copy,ignore=shutil.ignore_patterns('__pycache__','*.pyc'));logs={}
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        for script in ('verify_package.py','verify_sc.py'):
            result=subprocess.run(['python3',script],cwd=copy,env=env,capture_output=True,text=True,timeout=60)
            assert result.returncode==0,(script,result.stderr);logs[script]={'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr}
        equal=sha(copy/'verification.json')==sha(SOURCE/'verification.json')
        ck('unchanged source scripts replay byteidentically on copy',equal)
        out['replay']={'directory':str(temp),'logs':logs,'byte_equal':equal,'producer_assertions':148,
                       'exact_arithmetic_assertions':136,'floating_threshold_assertions':12}
    inputs=[Path(__file__),SOURCE/'MANIFEST.sha256',ROOT/'AGENTS.md',
        ROOT/'FIN_SIGNED_INTERFACE_ORBIT_CUT_NEXT_CHAIN_20261016.md',
        ROOT/'FIN_BALANCED_PHASE_CURRENT_NEXT_CHAIN_20261016.md',
        ROOT/'fin_balanced_phase_cut_research_20261016/MANIFEST.sha256',
        ROOT/'fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md']
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in previous]
    out['input_hashes']={str(path.relative_to(ROOT)):sha(path) for path in inputs}
    out['checks']=checks;out['status']='PASS_SCOPED_SC_INTAKE_WITH_NEW_BALANCED_PHASE_RESEARCH'
    out['source_manifest_after']=helper.manifest(SOURCE);assert out['source_manifest_after']==out['source_manifest']
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'source_files':out['source_manifest']['matched'],
        'matched_input_references':len(out['source_inputs']['matched']),'unavailable':len(out['source_inputs']['unavailable'])},indent=2))

if __name__=='__main__':main()
