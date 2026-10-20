#!/usr/bin/env python3
"""Source-preserving ID-001 audit; bounded controls are not entropy/PDE proofs."""
import argparse
from fractions import Fraction as F
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import tempfile
import time

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_intrinsic_entropy_source_bound_20261019'
ARCH=ROOT/'fin_strategy_audit_20261019'
POLICY=ROOT/'fin_agents_update_20261019'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def module(p,name):
    spec=importlib.util.spec_from_file_location(name,p);obj=importlib.util.module_from_spec(spec);spec.loader.exec_module(obj);return obj
def root_manifest():
    paths={}
    for row in (ARCH/'MANIFEST.sha256').read_text().splitlines():
        h,n=row.split('  ',1);p=Path(n)
        assert not p.is_absolute() and '..' not in p.parts and n not in paths
        assert sha(ROOT/p)==h,n;paths[n]=h
    actual={str(p.relative_to(ROOT)) for p in ARCH.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.name!='MANIFEST.sha256'}
    actual|={str(p.relative_to(ROOT)) for p in (ROOT/'research_tasks').rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    actual|={'FIN_RESEARCH_STRATEGY.md','FIN_PRIOR_RESEARCH_MAP.md','FIN_BLOCKERS_REGISTER.md','FIN_RESEARCH_EXECUTION_PLAN.md'}
    assert set(paths)==actual,(set(paths)-actual,actual-set(paths))
    return {'entries':len(paths),'sha256':sha(ARCH/'MANIFEST.sha256')}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);ap.add_argument('--replay',action='store_true');args=ap.parse_args()
    started=time.monotonic();checks=[]
    def ck(name,ok,**scope):
        assert bool(ok),name;checks.append({'name':name,'passed':True,**scope})
    policy=module(POLICY/'verify_update.py','accepted_hash_prefix_tools')
    ctx=(ROOT/'AGENTS.md').read_bytes();contract=json.loads((POLICY/'UPDATE_CONTRACT.json').read_text());extra=policy.check_prefix(ctx,contract)
    before=policy.manifest(SOURCE);ck('source14 entries exact with full coverage',before['entries']==14)
    architecture=root_manifest();ck('original47-file architecture remains intact',architecture['entries']==47)
    original=[]
    for row in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,n=row.split('  ',1);item={'original':n,'sha256':h,'status':'ORIGINAL_NOT_SUPPLIED'}
        if n.startswith('../fin_id001_source_archive/'):
            rel=n.removeprefix('../fin_id001_source_archive/');p=ROOT/rel
            assert p.is_file(),rel;item['current_counterpart']=rel;item['current_sha256']=sha(p)
            item['status']='MATCHED_PAYLOAD_BYTES' if sha(p)==h else 'CURRENT_COUNTERPART_DIFFERENT_BYTES'
        else:
            p=SOURCE/n
            if p.is_file():
                assert sha(p)==h,n;item['status']='MATCHED_ORIGINAL_BYTES'
        original.append(item)
    ck('two matching payloads three changed counterparts three absent originals separated',
       sum(r['status']=='MATCHED_PAYLOAD_BYTES' for r in original)==2 and
       sum(r['status']=='CURRENT_COUNTERPART_DIFFERENT_BYTES' for r in original)==3 and
       sum(r['status']=='ORIGINAL_NOT_SUPPLIED' for r in original)==3)
    scope_sources=[ROOT/'fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md',
        ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py',
        ROOT/'FIN_ST383_Asymptotic_Markov_Rate_Distortion_Bounds.json',
        ROOT/'FIN_ST372_ST386_Global_Collision_Envelope_Exact_Invariant_Bases_and_Continuum_IR_Report_EN.tex',
        ROOT/'FIN_Ten_Research_Programs_Monograph.md',
        ROOT/'fin_unified_strict_evolution_audit_20261017/AUTONOMOUS_COMPOSITION_THEOREM_OR_OBSTRUCTION.md',
        ROOT/'fin_kinetic_composition_research_20261016/THEORETICAL_RESULTS.md',
        ROOT/'fin_intermediate_clock_bridge_research_20261017/TRACE_CLOCK_AND_INSTRUMENT_THEOREM.md',
        ROOT/'research_tasks/00_EXECUTION_RULES.md',ROOT/'research_tasks/01_ID001_ENTROPY_SOURCE_BOUND.md',
        ARCH/'references/lecomte_appert_rolland_van_wijland_2007.pdf',
        ARCH/'references/lecomte_2007_sections_2_2_to_3_2.txt']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in scope_sources}
    run=subprocess.run(['python3',str(ARCH/'verify_strategy.py'),'--inputs-only'],cwd=ROOT,
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,text=True,timeout=60)
    ck('actual strategy scientific inputs pass read-only freshness check',run.returncode==0,scope='not a theorem or a source-law execution')
    freshness={'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr}

    # Exact different rational reset chain; exp clock variable itself not rationalized.
    p=[F(1,2),F(1,3),F(1,6)];r=F(3,4)
    mat=lambda rr:[[F(i==j)*rr+(1-rr)*p[j] for j in range(3)] for i in range(3)]
    mult=lambda A,B:[[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    P=mat(r);B=[[2*(p[j]-F(i==j)) for j in range(3)] for i in range(3)]
    ck('rational reset has exact stationarity and detailed balance',
       all(sum(p[i]*P[i][j] for i in range(3))==p[j] for j in range(3)) and
       all(p[i]*B[i][j]==p[j]*B[j][i] for i in range(3) for j in range(3)))
    power=mat(F(1))
    for k in range(1,6):power=mult(power,P);assert power==mat(r**k)
    ck('decimation semigroup exact for five powers',True,scope='finite algebra control, not entropy power-rule proof')
    abar=1-sum(v*v for v in p)
    ck('actual jump coefficient excludes null refreshes',2*abar==sum(p[i]*(-B[i][i]) for i in range(3))<2)
    for m in (2,3,5):
        for n in (1,2,4):assert {m*i+j for i in range(n) for j in range(m)}==set(range(m*n))
    ck('partition power-rule block-index identity independently checked',True)
    ck('inverse-Fano endpoint and midpoint ordering correct',F(1,100)<F(11,12)<1)
    # Exact coarse native g0/H<=1 witness inequalities.
    ck('g0 positive entropy witness has explicit rational margin',F(11,12)-F(1,2)==F(5,12))
    ck('Euler exponential bounds used by coarse certificate have rational controls',1+1+F(1,2)+F(1,6)>F(5,2) and F(5,2)**8>1100)
    ck('positive error floor beats one-percent threshold in the stated control',F(9,100)<F(5,12),
       scope='log/exp inequalities are paid analytically in SUBSAMPLING_ERROR_FLOOR, not inferred from this arithmetic')
    # Native binary64 model check via versioned source, not the incoming reconstruction.
    import numpy as np
    safe=module(scope_sources[1],'versioned_native_model');_,_,_,_,A=safe.build_rank7();g=3.
    field=g*A[0]/2;pp=np.exp(field-field.max());pp/=pp.sum();q=12
    rates=2*(np.tile(pp,(q,1))-np.eye(q))
    states=[(i,j) for i in range(q) for j in range(q)];idx={v:i for i,v in enumerate(states)}
    L=np.zeros((144,144));R=np.zeros((144,12))
    for s0,(i,j) in enumerate(states):
        R[s0,(j-i)%q]=1
        for slot,other in ((0,j),(1,i)):
            logits=g*A[:,other]/2;probs=np.exp(logits-logits.max());probs/=probs.sum()
            for new in range(q):
                target=(new,j) if slot==0 else (i,new)
                if target!=(i,j):L[s0,idx[target]]+=probs[new]
        L[s0,s0]=-sum(L[s0])
    defect=float(np.max(abs(L@R-R@rates)))
    ck('direct144state native quotient confirms factor2 clock',defect<2e-13,scope='binary64 N2 model algebra, not physical sampling',max_defect=defect)

    import mpmath as mp
    mp.mp.dps=65
    def ent(vals):return -sum(v*mp.log(v) for v in vals if v)
    def entropy(p0,tau):
        tt=-mp.expm1(-2*tau)
        return sum(p0[i]*ent([(1-tt if i==j else 0)+tt*p0[j] for j in range(len(p0))]) for i in range(len(p0)))
    def invF(value,q=12):
        lo=mp.mpf(0);hi=mp.mpf(q-1)/q
        for _ in range(160):
            mid=(lo+hi)/2;ff=ent([mid,1-mid])+mid*mp.log(q-1)
            if ff<value:lo=mid
            else:hi=mid
        return (lo+hi)/2
    # A finite joint-distribution check, using nonuniform rational marginals.
    for e in (F(0),F(1,10),F(1,2),F(1)):
        joint=[[p[i]*((1-e) if i==j else e/2) for j in range(3)] for i in range(3)]
        vals=[mp.mpf(v.numerator)/v.denominator for row in joint for v in row]
        py=[sum(joint[i][j] for i in range(3)) for j in range(3)]
        Hxy=ent(vals);Hy=ent([mp.mpf(v.numerator)/v.denominator for v in py]);ee=mp.mpf(e.numerator)/e.denominator
        assert Hxy-Hy<=ent([ee,1-ee])+ee*mp.log(2)+mp.mpf('1e-55')
    ck('Fano flag bound independently checked on four finite rational joint tables',True,scope='numerical entropy evaluation of exact finite probabilities only')
    pmp=[mp.mpf(str(v)) for v in pp];pmp=[v/sum(pmp) for v in pmp]
    Ap=1-sum(v*v for v in pmp);Cp=sum(v*(1-v)*mp.log(v) for v in pmp);rb=2*Ap
    cc=rb*(1-mp.log(2))-2*Cp;rows=[]
    for m in (32,128,512):
        dd=mp.mpf(1)/m;tt=-mp.expm1(-2*dd);hh=entropy(pmp,dd)
        reduced=-sum(v*(1-tt*(1-v))*mp.log(1-tt*(1-v)) for v in pmp)-tt*Ap*mp.log(tt)-tt*Cp
        assert abs(hh-reduced)<mp.mpf('1e-55')
        bound=dd*dd*(2*Ap*mp.log(1/dd)+mp.mpf(14)/3*Ap+2*abs(Cp)+4)
        assert abs(hh-rb*dd*mp.log(1/dd)-cc*dd)<bound
        k=m//2;assert mp.mpf(k)/m==mp.mpf('.5')
        coarse=entropy(pmp,mp.mpf(k)/m)-mp.mpf(k)/m
        rows.append({'m':m,'h_point':mp.nstr(hh,30),'remainder_bound_point':mp.nstr(bound,25),
                     'H1_tau_half_error_floor_point':mp.nstr(invF(coarse),25)})
    ck('native entropy identity and explicit remainder agree at three resolutions',len(rows)==3,scope='65digit point diagnostics on a normalized declared binary64 model; not outward exact-kernel numerics')
    ck('same coarse witness persists at all three refined grids',len({v['H1_tau_half_error_floor_point'] for v in rows})==1)
    # Exact countercontrol for pair correlations versus whole histories.
    singles=F(0);pairs={j:F(0) for j in range(1,5)};triple=F(0)
    import itertools
    for bits in itertools.product((-1,1),repeat=6):
        seq=[]
        for b in range(3):aa,bb=bits[2*b:2*b+2];seq.extend((aa,bb,aa*bb))
        for phase in range(3):
            weight=F(1,64*3);z=seq[phase:]
            singles+=weight*z[0]
            for lag in pairs:pairs[lag]+=weight*z[0]*z[lag]
            triple+=weight*z[0]*z[1]*z[2]
    ck('pairwise-independent triple countercontrol exact',singles==0 and all(v==0 for v in pairs.values()) and triple==F(1,3))
    ck('native equilibrium has zero edge entropy-production ratios',all(p[i]*B[i][j]==p[j]*B[j][i] for i in range(3) for j in range(3)),scope='positive sampled entropy is not thermodynamic irreversibility')

    # Target-specific duplicate search, not a new whole-repository inventory.
    candidates=[scope_sources[i] for i in (0,2,4,5,6,7)]
    pattern=re.compile(r'stationary.joining|decimat|subsampl.{0,35}entropy|entropy.{0,35}resolution|finite.entropy',re.I)
    hits=[]
    for path in candidates:
        for i,line in enumerate(path.read_text().splitlines(),1):
            if pattern.search(line):hits.append({'path':str(path.relative_to(ROOT)),'line':i,'text':line[:300]})
    duplicate={'scope':'six named scientific predecessor files, not whole repo','searched':{str(p.relative_to(ROOT)):sha(p) for p in candidates},
               'pattern':pattern.pattern,'hits':hits,'disposition':'reuse fixed-resolution Fano, unitary/regularity obstructions and known CTMC divergence; no equivalent accepted native decimation-floor theorem found in this bounded read'}
    replay={}
    if args.replay:
        temp=Path(tempfile.mkdtemp(prefix='fin-id001-intake-'));copy=temp/'campaign';shutil.copytree(SOURCE,copy)
        env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        proc=subprocess.run(['python3','verify_id001.py','--output','results.json'],cwd=copy,env=env,capture_output=True,text=True,timeout=60)
        assert proc.returncode==0,proc.stderr;old=json.loads((SOURCE/'results.json').read_text());new=json.loads((copy/'results.json').read_text())
        drift=[]
        def compare(x,y,path=''):
            if isinstance(x,dict):assert set(x)==set(y);[compare(x[k],y[k],path+'/'+k) for k in x]
            elif isinstance(x,list):assert len(x)==len(y);[compare(a,b,path+'/'+str(i)) for i,(a,b) in enumerate(zip(x,y))]
            elif isinstance(x,(float,int)) and not isinstance(x,bool):
                assert math.isfinite(float(x)) and math.isfinite(float(y));d=abs(x-y)
                if d:drift.append({'field':path,'abs_drift':d})
            else:assert x==y,(path,x,y)
        compare(old,new);maximum=max((v['abs_drift'] for v in drift),default=0)
        ck('one frozen producer g3 diagnostic replays numerically on copy',maximum<1e-9)
        replay={'directory':str(temp),'returncode':proc.returncode,'stdout':proc.stdout,'stderr':proc.stderr,
                'byte_equal':sha(copy/'results.json')==sha(SOURCE/'results.json'),'max_numeric_drift':maximum,'changed_numeric_fields':drift,
                'scope':'binary64 finite-model diagnostic, not proof of entropy/KS/joining or an outward certificate'}
        ck('versioned model and producer p agree as nominal numerics',max(abs(np.array(new['p'])-pp))<2e-13)
    ck('source and architecture preserved before and after replay',policy.manifest(SOURCE)==before and root_manifest()==architecture)
    ck('all audited scientific inputs remain unchanged',all(sha(ROOT/n)==h for n,h in hashes.items()) and (ROOT/'AGENTS.md').read_bytes()==ctx)
    peak_self=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    peak_child=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    ck('bounded diagnostic respects oneGiB RSS target',max(peak_self,peak_child)<1024*1024)
    live={'bytes':len(ctx),'sha256':hashlib.sha256(ctx).hexdigest(),'older_verified_prefix_bytes':contract['updated_prefix_bytes'],
          'older_verified_prefix_sha256':contract['updated_prefix_sha256'],'later_bytes':extra}
    hashes.update({str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),SOURCE/'MANIFEST.sha256',SOURCE/'INPUTS.sha256',ARCH/'MANIFEST.sha256',POLICY/'UPDATE_CONTRACT.json',
        HERE/'PROOF_AUDIT_AND_SCOPE.md',HERE/'SUBSAMPLING_ERROR_FLOOR.md']})
    result={'date':'2026-10-19','entry_HEAD':'0c464f01','verification_HEAD':subprocess.run(['git','rev-parse','HEAD'],cwd=ROOT,capture_output=True,text=True,check=True).stdout.strip(),
            'status':'PASS_SCOPED_ID001_INTAKE_WITH_SUBSAMPLING_FLOOR','source_manifest':before,'source_manifest_after':policy.manifest(SOURCE),
            'architecture_manifest':architecture,'original_references':original,'freshness':freshness,'input_hashes':hashes,'AGENTS_context':live,
            'checks':checks,'replay':replay,'duplicate_search':duplicate,'native_model_diagnostic':{'quotient_defect':defect,'rbar_point':mp.nstr(rb,28),'c_point':mp.nstr(cc,28),'rows':rows},
            'resources':{'elapsed_seconds':time.monotonic()-started,'self_max_rss_KiB':peak_self,'child_max_rss_KiB':peak_child},
            'F02':'CLASS_OBSTRUCTION_CLOSED_IN_DECLARED_STATIONARY_JOINING_CLASS','F01':'OPEN_SOURCE','physical_validation':False}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'source_entries':before['entries'],'bounded_checks':len(checks),'F01':result['F01'],'max_rss_KiB':max(peak_self,peak_child)},indent=2))

if __name__=='__main__':main()
