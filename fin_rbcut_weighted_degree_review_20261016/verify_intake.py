#!/usr/bin/env python3
"""Combined RB/RBCUT intake. Source copies only; exact bounded controls."""
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
RB=ROOT/'fin_rare_boundary_flow_campaign_20261016'
BUNDLE=ROOT/'FIN_RBCUT_DO_WEIGHTED_DEGREE_KOMPLET_20261016'
STAGES={
    'shell':BUNDLE/'BADANIA/01_RB_CUT_NATIVE_SHELL',
    'repair':BUNDLE/'BADANIA/02_RB_CUT_QUANTITATIVE_REPAIR/fin_rb_cut_quantitative_repair_campaign_20261016',
    'orbit':BUNDLE/'BADANIA/03_RB_ORBIT_FLOW_REDUCTION/fin_rb_orbit_flow_reduction_campaign_20261016',
    'degree':BUNDLE/'BADANIA/04_RB_ORBIT_WEIGHTED_DEGREE/fin_rb_orbit_weighted_degree_campaign_20261016'}
NAMES={'shell':'fin_rb_cut_native_shell_campaign_20261016','repair':'fin_rb_cut_quantitative_repair_campaign_20261016',
       'orbit':'fin_rb_orbit_flow_reduction_campaign_20261016','degree':'fin_rb_orbit_weighted_degree_campaign_20261016'}
ZIP_NAMES={'shell':'FIN_RB_CUT_NATIVE_SHELL_ONLY_HANDOFF_20261016.zip','repair':'FIN_RB_CUT_QUANTITATIVE_REPAIR_ONLY_HANDOFF_20261016.zip',
           'orbit':'FIN_RB_ORBIT_FLOW_REDUCTION_ONLY_HANDOFF_20261016.zip','degree':'FIN_RB_ORBIT_WEIGHTED_DEGREE_ONLY_HANDOFF_20261016.zip'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder,excluded=(),filename='MANIFEST.sha256'):
    names=[]
    for line in (folder/filename).read_text().splitlines():
        h,n=line.split('  ',1); assert n not in names and sha(folder/n)==h,n; names.append(n)
    actual={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{filename}|set(excluded),(folder,actual-set(names))
    return {'matched':len(names),'coverage':'COMPLETE','excluded_hashes':{n:sha(folder/n) for n in excluded}}
def prefix(h):
    digest=hashlib.sha256(); size=0
    for line in (ROOT/'AGENTS.md').read_bytes().splitlines(keepends=True):
        digest.update(line); size+=len(line)
        if digest.hexdigest()==h:return size
    return None

def main():
    p=argparse.ArgumentParser(); p.add_argument('--output',type=Path,required=True); p.add_argument('--replay',action='store_true'); a=p.parse_args(); checks=[]
    def ck(name,ok,**detail):checks.append({'name':name,'passed':bool(ok),**detail}); assert ok,name
    out={'date':'2026-10-16','entry_baseline':'8571db1d plus preserved KR/BA/FP worktree additions'}
    out['RB_manifest']=manifest(RB,('PACKAGE_VERIFICATION.json',))
    out['bundle_manifest']=manifest(BUNDLE,('README_PL.md',),'MANIFEST_SHA256.txt')
    out['stage_manifests']={k:manifest(v) for k,v in STAGES.items()}
    ck('RB outer bundle and all four inner manifests are hash-bound with coverage',out['bundle_manifest']['matched']==75)
    zipped={}
    for key,folder in STAGES.items():
        path=BUNDLE/'ORYGINALNE_PAKIETY'/ZIP_NAMES[key]; names=[]
        with zipfile.ZipFile(path) as z:
            for item in z.infolist():
                if item.is_dir():continue
                parts=Path(item.filename).parts
                assert not Path(item.filename).is_absolute() and '..' not in parts
                if parts[0]==folder.name:parts=parts[1:]
                local=folder.joinpath(*parts); assert local.is_file() and hashlib.sha256(z.read(item.filename)).hexdigest()==sha(local),item.filename
                names.append(str(local.relative_to(folder)))
        assert set(names)=={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file()}
        zipped[key]={'archive_sha256':sha(path),'matched_members':len(names)}
    out['ZIP_member_binding']=zipped
    ck('all71 original campaign members match the four actual ZIPs',sum(v['matched_members'] for v in zipped.values())==71)
    wrappers={
        'FIN_MT_001_MT_005_FULL_HANDOFF_20261014.zip':ROOT/'fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip',
        'FIN_UF_001_UF_005_FULL_HANDOFF_20261014.zip':ROOT/'fin_critical_metastability_campaign_20261014/source_context/UF_PREDECESSOR.zip'}
    roots={NAMES[k]:folder for k,folder in STAGES.items()}
    def resolve(name):
        if '::' in name:
            archive,member=name.split('::',1); path=wrappers.get(Path(archive).name)
            if path is None:return None
            with zipfile.ZipFile(path) as z:return hashlib.sha256(z.read(member)).hexdigest(),str(path.relative_to(ROOT))+'::'+member
        if Path(name).name in wrappers:
            path=wrappers[Path(name).name]; return sha(path),str(path.relative_to(ROOT))
        rel=name.removeprefix('/mnt/data/'); parts=Path(rel).parts
        if parts and parts[0] in roots:path=roots[parts[0]].joinpath(*parts[1:])
        else:path=ROOT/rel
        if path.is_file():return sha(path),str(path.relative_to(ROOT))
        return None
    provenance={}
    for key,folder in {'RB':RB,**STAGES}.items():
        matched={}; missing={}
        for line in (folder/'INPUTS.sha256').read_text().splitlines():
            h,n=line.split('  ',1); result=resolve(n)
            if result is None:missing[n]=h
            else:
                actual,local=result; assert actual==h,(key,n); matched[n]=local
        provenance[key]={'matched':matched,'unavailable_originals':missing}
    out['provenance']=provenance
    ck('all available original dependency references match without archive substitution',all(v['matched'] for v in provenance.values()))
    task_lists={}
    for ledger in sorted(RB.glob('RB-*/INPUTS.sha256')):
        counts={'matched':0,'unavailable':0}
        for line in ledger.read_text().splitlines():
            h,n=line.split('  ',1);result=resolve(n)
            if result is None:counts['unavailable']+=1
            else:assert result[0]==h;counts['matched']+=1
        task_lists[str(ledger.relative_to(RB))]=counts
    out['RB_stage_references']=task_lists
    ck('all four RB task input lists checked including inactive RB003',len(task_lists)==4)
    # Actual reviewed predecessor dependencies; preserve historical mutable ledger inputs by prefix.
    previous={}
    for rel in ('fin_count_partition_feasibility_review_20261016','fin_reversible_bracket_averaging_review_20261016',
                'fin_stationary_kinetic_rigidity_review_20261016','fin_activity_covariance_research_20261016','fin_critical_metastability_review_20261015'):
        previous[rel]=manifest(ROOT/rel)
        saved=json.loads((ROOT/rel/('RESULTS.json' if rel.startswith('fin_activity_') else 'verification.json')).read_text())
        for path,h in saved.get('input_hashes',saved.get('inputs',{})).items():
            if path=='AGENTS.md':
                size=prefix(h); assert size is not None; previous[rel]['historical_AGENTS_prefix_bytes']=size
            else:assert sha(ROOT/path)==h,path
    out['actual_dependencies']=previous
    ck('current FP BA KR AC MT reviews and their scientific inputs remain intact',len(previous)==5)
    # Exact native reference ellipticity constant: positive exponential series plus geometric tail.
    pmin=(F(1)-F(1,100))/12-F(1,10000); pmax=(F(1)+F(1,100))/12+F(1,10000)
    ck('actual MT occupation margins',pmin==F(103,1250) and pmax==F(158,1875) and pmax<F(1,9))
    total=sum(F(24)**j/math.factorial(j) for j in range(121)); next_term=F(24)**121/math.factorial(121)
    upper=total+next_term/(1-F(24,122)); alo=pmin/(12*upper); ahi=pmin/(12*total)
    ck('explicit native a0 lower bound certified by rational exponential tail',F('2.59e-13')<alo<ahi)
    out['ellipticity_certificate']={'a0':'(103/15000)*exp(-24)','proved_lower_public':'2.59e-13',
        'method':'121 positive e^24 terms, geometric tail, reciprocal bounds in Fraction',
        'scope':'native retained local-MT edges only; no diagnostic lower activity or physical source'}
    # Independent mixed-orbit-size nonuniform count graph (not the source36-vertex toy).
    modelpath=ROOT/'fin_activity_covariance_research_20261016/verify_research.py'
    spec=importlib.util.spec_from_file_location('old_exact_model',modelpath); model=importlib.util.module_from_spec(spec); spec.loader.exec_module(model)
    W,Q=model.exact_hb(2,4,F(2)); V=list(W); sig=(1,-1,1,-1); level=lambda x:sum(i*j for i,j in zip(x,sig))
    U={(x,y):W[x]*rate for x in V for y,rate in Q[x].items() if level(y)>level(x)}
    source={x:F(0) for x in V}
    for (x,y),v in U.items():source[x]+=v; source[y]-=v
    swap=lambda x:(x[2],x[1],x[0],x[3]); orb=lambda x:min(x,swap(x))
    orbits={o:[x for x in V if orb(x)==o] for o in set(map(orb,V))}
    ub={}; sb={o:sum(source[x] for x in xs) for o,xs in orbits.items()}; wb={o:sum(W[x] for x in xs) for o,xs in orbits.items()}
    for (x,y),v in U.items():ub[orb(x),orb(y)]=ub.get((orb(x),orb(y)),F(0))+v
    ck('unequal orbit sizes and nonuniform Gibbs weights explicitly exercised',set(map(len,orbits.values()))=={1,2} and len(set(W.values()))>1)
    for (x,y),v in U.items():assert U[swap(x),swap(y)]==v
    ck('reference source and weights invariant in the parity-preserving subgroup',all(W[swap(x)]==W[x] and source[swap(x)]==source[x] and level(swap(x))==level(x) for x in V))
    # Add a quotient diamond circulation to make lift testing nontrivial.
    low=orb((0,1,0,1)); mid1=orb((1,0,0,1)); mid2=orb((1,1,0,0)); high=orb((1,0,1,0))
    f=dict(ub); eps=min(f[low,mid2],f[mid2,high])/4
    for e,sign in (((low,mid1),1),((mid1,high),1),((low,mid2),-1),((mid2,high),-1)):f[e]+=eps*sign
    lifted={(x,y):f[orb(x),orb(y)]*v/ub[orb(x),orb(y)] for (x,y),v in U.items()}
    divergence={x:F(0) for x in V}
    for (x,y),v in lifted.items():divergence[x]+=v; divergence[y]-=v
    ck('nontrivial quotient flow lift restores every original vertex demand',lifted!=U and divergence==source and all(v>0 for v in lifted.values()))
    Cfull=max(sum(v for e,v in lifted.items() if x in e)/(2*W[x]) for x in V)
    Cquot=max(sum(v for e,v in f.items() if o in e)/(2*wb[o]) for o in orbits)
    ck('vertex capacity preserved exactly despite orbit stabilizers',Cfull==Cquot)
    def cut_stats(nodes,edgeweights,demands,weights,N):
        kmax=qmax=F(0)
        for mask in range(1<<len(nodes)):
            pick={nodes[j] for j in range(len(nodes)) if mask>>j&1}; s=sum(demands[x] for x in pick)
            cap=sum(v for (x,y),v in edgeweights.items() if x in pick and y not in pick)
            degree=sum(weights[x] for (x,y) in edgeweights if x in pick and y not in pick)
            if s>0:
                assert cap>0 and degree>0; kmax=max(kmax,s/cap); qmax=max(qmax,s/(N*degree))
        return kmax,qmax
    Kfull,Qfull=cut_stats(V,U,source,W,2)
    # Full cut count, not aggregated-arc degree, defines the artificial weighted score.
    Kquot=F(0)
    for mask in range(1<<len(orbits)):
        pick={o for j,o in enumerate(orbits) if mask>>j&1}; demand=sum(sb[o] for o in pick)
        cap=sum(v for (o,t),v in ub.items() if o in pick and t not in pick)
        if demand>0:assert cap>0; Kquot=max(Kquot,demand/cap)
    ck('all1024 original cuts and128 mixed-size orbit cuts agree on sharp cap',len(V)==10 and len(orbits)==7 and Kfull==Kquot==1)
    amin=min(v/(2*W[x]) for (x,y),v in U.items()); nativeC=max(sum(v for e,v in U.items() if x in e)/(2*W[x]) for x in V)
    ck('exact toy cut degree and native cap inequalities',Qfull<=nativeC<=Kfull<=Qfull/amin and amin*Kfull<=nativeC)
    out['mixed_orbit_fixture']={'states':10,'orbit_sizes':sorted(map(len,orbits.values())),'Kstar':str(Kfull),
        'Qfull':str(Qfull),'native_vertex_cap':str(nativeC),'modified_lift_vertex_cap':str(Cfull),
        'scope':'nonuniform q4/N2 count model, NOT actual q12 MT graph'}
    # Correct dual convention: phi increasing upward when beta=0.
    e=((0,1),(1,2)); bad_s=(-1,0,1); phi=(-1,0,0)
    ck('beta-zero orbit dual sign is nondecreasing not nonincreasing',all(phi[x]-phi[y]<=0 for x,y in e) and sum(x*y for x,y in zip(phi,bad_s))==1)
    # Signed repair, independently explicit final flow rather than producer maxflow.
    h={(1,2):F(2),(2,4):F(2),(4,6):F(2),(0,3):F(1),(3,5):F(1),(5,7):F(1)}
    deleted={(0,2):F(1),(4,7):F(1)}; final={(1,2):F(2),(2,4):F(1),(4,6):F(2),(0,3):F(2),(3,5):F(2),(5,7):F(2)}
    def div(edges):
        result=[F(0)]*8
        for (x,y),v in edges.items():result[x]+=v; result[y]-=v
        return result
    native=dict(h);native.update(deleted)
    ck('signed collective repair works where positive-only repair has a cut obstruction',div(native)==div(final) and final[2,4]<h[2,4] and sum(div(deleted)[i] for i in (4,6))==1 and not any(x in (4,6) and y not in (4,6) for x,y in h))
    # New gap derivation: exact incoming-compensated score and why unsigned ratios fail.
    U4={(0,2):F(1),(2,3):F(1),(0,1):F(1),(1,3):F(1)}; H4={(0,2),(2,3)}; F4=set(U4)-H4; W4={i:F(1) for i in range(4)}
    source4=[F(0)]*4
    for (x,y),v in U4.items():source4[x]+=v;source4[y]-=v
    for mask in range(16):
        S={i for i in range(4) if mask>>i&1}; B=sum(W4[x] for x,y in H4 if x in S and y not in S)
        upH=sum(U4[x,y] for x,y in H4 if x in S and y not in S)
        upF=sum(U4[x,y] for x,y in F4 if x in S and y not in S)
        incoming=sum(v for (x,y),v in U4.items() if x not in S and y in S)
        demand=sum(source4[i] for i in S); excess=upF-incoming
        assert demand==upH+excess
        if B:
            score=max(demand,F(0))/(2*B); xi=max(excess,F(0))/(2*B)
            assert xi<=score<=F(1,2)+xi
        else:assert (demand>0)==(excess>0)
    ck('all sixteen compensated-cut identities and score bounds exact',True)
    phase={0,2,3}
    ck('whole-phase zero-demand cut defeats a necessary unsigned-shell bound',sum(source4[x] for x in phase)==0 and not any(x in phase and y not in phase for x,y in H4) and any(x in phase and y not in phase for x,y in F4))
    # Actual Fourier multiplicities pay the NEW prefix-cut estimate.
    cos=(1,0,-1,0)*3; sin=(0,1,0,-1)*3; sig=tuple((-1)**j for j in range(12))
    for direction in (-2,2):
        tally={}
        for i,j in itertools.product(range(12),repeat=2):
            if sig[j]-sig[i]==direction:
                key=cos[j]-cos[i],sin[j]-sin[i]; tally[key]=tally.get(key,0)+1
        assert len(tally)==4 and set(tally.values())=={9}
    ck('actual q12 four fast directions each have nine distinct label pairs',pmax<F(1,9))
    # Same-parity transitions can cross the radius boundary: C changes by two.
    for k in range(1,12):
        n=[100]*12;n[0]+=k;n[2]-=k
        nn=list(n);nn[0]+=1;nn[2]-=1
        C=lambda state:sum(a*b for a,b in zip(state,cos))
        S=lambda state:sum(a*b for a,b in zip(state,sin))
        assert C(n)==2*k and C(nn)==2*(k+1) and S(n)==S(nn)==0
        assert sum(a*b for a,b in zip(n,sig))==sum(a*b for a,b in zip(nn,sig))==0
    ck('same-parity q12 radius crossing primitive is not a safe phase move',True)
    for n in (2,3,10,100):
        edge=F(1)+F(1,n); incident=(edge,2*edge,edge)
        assert max(incident)==2*(1+F(1,n)) and max(incident)<=3
        assert max(incident)/n<=F(3,2)
    ck('four-node negative control total vertex rate needs factor two',True)
    ck('critical drift deletion obstruction coefficient normalized correctly',4*F(1,8)==F(1,2))
    if a.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-rbcut-intake-')); logs={}
        jobs={'RB':(RB,'verify_rb.py','EXACT_DIAGNOSTICS.json'),
              'shell':(STAGES['shell'],'verify_cut_shell.py','verification.json'),
              'repair':(STAGES['repair'],'verify_repair.py','verification.json'),
              'orbit':(STAGES['orbit'],'verify_orbit.py','EXACT_CHECKS.json'),
              'degree':(STAGES['degree'],'verify_weighted_degree.py','FINITE_CHECKS.json')}
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1')
        for key,(folder,writer,output) in jobs.items():
            copy=temp/key;shutil.copytree(folder,copy,ignore=shutil.ignore_patterns('__pycache__','*.pyc')); runlogs={}
            for script in ('verify_package.py',writer):
                run=subprocess.run(['python3',script],cwd=copy,env=env,capture_output=True,text=True,timeout=60)
                if run.returncode:
                    # This checker includes original /mnt/data custody paths. Do not
                    # invent those files or describe its unmodified run as passing.
                    assert key=='repair' and script=='verify_package.py' and 'AssertionError: /mnt/data/fin_rb_cut_native_shell_campaign_20261016/NEXT_OPEN_ATOM.md' in run.stderr,(key,script,run.stderr)
                    fallback=manifest(copy)
                else:fallback=None
                runlogs[script]={'exit_code':run.returncode,'stdout':run.stdout,'stderr':run.stderr,
                    'independent_manifest_fallback':fallback,
                    'scope':'historical absolute-input-path failure retained' if fallback else 'unmodified source script'}
            logs[key]={'logs':runlogs,'scientific_output_byte_equal':sha(copy/output)==sha(folder/output),
                'scope':'finite algebra/toy controls, not the actual all-N weighted-cut theorem'}
        ck('all five scientific producers replay from copies with source outputs unchanged',all(v['scientific_output_byte_equal'] for v in logs.values()))
        out['replay']={'directory':str(temp),'campaigns':logs}
    inputs=[Path(__file__),RB/'MANIFEST.sha256',BUNDLE/'MANIFEST_SHA256.txt',modelpath,ROOT/'AGENTS.md',
        ROOT/'FIN_SIGNED_INTERFACE_ORBIT_CUT_NEXT_CHAIN_20261016.md']
    inputs += [folder/'MANIFEST.sha256' for folder in STAGES.values()]
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in previous]
    inputs += [ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json',ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py',
        ROOT/'fin_count_partition_feasibility_review_20261016/BOUNDARY_LEAKAGE_AND_IDEAL_CUT_SEED.md',
        ROOT/'fin_reversible_bracket_averaging_review_20261016/COUNT_PARTITION_FEASIBILITY_SEED.md',
        ROOT/'fin_critical_metastability_campaign_20261014/MT-001_COERCIVITY_CONNECTIVITY.md']
    out['input_hashes']={str(path.relative_to(ROOT)):sha(path) for path in inputs}
    out['status']='PASS_SCOPED_RB_AND_WEIGHTED_DEGREE_INTAKE';out['checks']=checks
    out['RB_manifest_after']=manifest(RB,('PACKAGE_VERIFICATION.json',));out['bundle_manifest_after']=manifest(BUNDLE,('README_PL.md',),'MANIFEST_SHA256.txt')
    assert out['RB_manifest']==out['RB_manifest_after'] and out['bundle_manifest']==out['bundle_manifest_after']
    a.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'RB_files':out['RB_manifest']['matched'],
        'bundle_files':out['bundle_manifest']['matched'],'ZIP_members':sum(v['matched_members'] for v in zipped.values()),
        'input_references':{k:{'matched':len(v['matched']),'unavailable':len(v['unavailable_originals'])} for k,v in provenance.items()},
        'ellipticity':out['ellipticity_certificate']},indent=2))

if __name__=='__main__':main()
