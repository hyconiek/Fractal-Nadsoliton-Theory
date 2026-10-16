#!/usr/bin/env python3
"""Independent bounded MM checks; source writes only in a temporary copy."""
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
import numpy as np
from scipy.linalg import expm, solve_continuous_lyapunov
from scipy.special import softmax
import sympy as sy
from certify_margin import build, ends, public
import mpmath as mp

HERE=Path(__file__).resolve().parent; ROOT=HERE.parent
SOURCE=ROOT/'fin_mediator_memory_observability_campaign_20261016'
SEED=ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def manifest(folder,excluded=()):
    names=[]
    for row in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=row.split('  ',1); assert n not in names and sha(folder/n)==h,n; names.append(n)
    actual={str(p.relative_to(folder)) for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==set(names)|{'MANIFEST.sha256'}|set(excluded),(folder,actual-set(names))
    return {'matched':len(names),'coverage':'COMPLETE_IN_DECLARED_SCOPE','excluded_hashes':{n:sha(folder/n) for n in excluded}}

def changes(a,b,path=''):
    if isinstance(a,dict) and isinstance(b,dict):
        assert a.keys()==b.keys(); return sum((changes(a[k],b[k],path+'/'+k) for k in a),[])
    if isinstance(a,list) and isinstance(b,list):
        assert len(a)==len(b); return sum((changes(x,y,path+'/'+str(j)) for j,(x,y) in enumerate(zip(a,b))),[])
    return [] if a==b else [{'path':path,'source':a,'replay':b}]

def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--output',type=Path,required=True); ap.add_argument('--replay',action='store_true'); args=ap.parse_args()
    checks=[]
    def ck(name,ok,**detail):
        checks.append({'name':name,'passed':bool(ok),**detail}); assert ok,name
    out={'date':'2026-10-16','baseline':'3be481a7 plus preserved NT/HM worktree changes'}
    excluded=('PACKAGE_VERIFICATION.json','verify_package.py','build_mm_documents.py')
    out['source_manifest']=manifest(SOURCE,excluded)
    records=json.loads((SOURCE/'INPUT_PROVENANCE.json').read_text())['actual_source_file_bytes']
    expected=''.join(v['sha256']+'  '+v['path']+'\n' for v in records.values())
    ck('six original input declarations and thirty stage references agree',
       (SOURCE/'INPUTS.sha256').read_text()==expected and all((SOURCE/f'MM-{i:03d}/INPUTS.sha256').read_text()==expected for i in range(1,6)))
    mapping={}
    for row in records.values():
        if row['sha256']=='5d390389192c6882b9a911c4455911d528895b47dabbf0926252e97298c5a899':
            rel='fin_intrinsic_observation_campaign_20261015/source_context/FIN_MT_FULL_HANDOFF.zip'
            assert sha(ROOT/rel)==row['sha256']; mapping[row['path']]=rel
    out['provenance']={'matched_wrappers':mapping,'unavailable_originals':[v['path'] for v in records.values() if v['path'] not in mapping],
        'stage_references':{'matched':5*len(mapping),'unavailable':5*(len(records)-len(mapping))},
        'scope':'Actual repo inputs bound separately; no equality to missing uploaded plan or remote custody inferred.'}
    ck('exactly one original wrapper currently matched',len(mapping)==1)
    prior={}
    for rel in ('fin_auxiliary_mediator_review_20261016','fin_dynamic_kernel_observability_review_20261015',
                'fin_kernel_observability_review_20261015','fin_critical_metastability_review_20261015','fin_native_transfer_review_20261016'):
        folder=ROOT/rel; r=manifest(folder); old=json.loads((folder/'verification.json').read_text())
        for path,h in old['input_hashes'].items():assert sha(ROOT/path)==h,path
        r['saved_inputs_unchanged']=len(old['input_hashes']); prior[rel]=r
    out['actual_review_dependencies']=prior
    ck('actual predecessor manifests and saved scientific inputs unchanged',len(prior)==5)
    spec=importlib.util.spec_from_file_location('versioned_seed',SEED); seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,lam,X,A=seed.build_rank7(); g=12/lam[6]; c=lam[3]/lam[6]; gamma=1-c
    features=np.sqrt(2)*np.cos(2*np.pi*np.array([1,3])[:,None]*np.arange(12)/12)
    ck('versioned null and active Fourier maps',np.linalg.norm(A@features[0])<1e-13 and np.linalg.norm(A@features[1]-lam[3]*features[1])<1e-13)
    ck('feature means variances cross terms',np.max(abs(features@features.T/12-np.eye(2)))<1e-14 and np.max(abs(features.mean(1)))<1e-14)
    # Direct labelled jump list, repeated labels, deliberately nonuniform block/field.
    labels=np.array([0,0,1,2,3,3,4,5,6,7,8,9,9,10,11,1,3,5,7,0,2,4,6,8]); N=len(labels); k=9
    a=np.array([N/k-1 if j<k else -1 for j in range(N)]); weights=np.array([np.ones(N),a])
    theta=np.array([.13,-.17,.05,.22,-.08,.16,.11]); q=softmax(X@theta)
    drift=np.zeros(4); bracket=np.zeros((4,4))
    for slot,old in enumerate(labels):
        for new in range(12):
            delta=np.array([weights[t,slot]*(features[r,new]-features[r,old])/math.sqrt(N) for r in range(2) for t in range(2)])
            drift+=q[new]*delta; bracket+=q[new]*np.outer(delta,delta)
    global_=features[:,labels].sum(1)/math.sqrt(N); private=(features[:,labels]*a).sum(1)/math.sqrt(N)
    desired=np.array([math.sqrt(N)*(features[r]@q)-global_[r] if t==0 else -private[r] for r in range(2) for t in range(2)])
    ck('actual nonuniform labelled global and private drift',np.max(abs(drift-desired))<2e-14)
    formula=np.zeros((4,4))
    for r,t,s,u in itertools.product(range(2),repeat=4):
        inner=sum(weights[t,j]*weights[u,j]*sum(q[v]*(features[r,v]-features[r,o])*(features[s,v]-features[s,o]) for v in range(12)) for j,o in enumerate(labels))/N
        formula[2*r+t,2*s+u]=inner
    ck('all sixteen joint jump brackets directly reconstructed',np.max(abs(bracket-formula))<2e-14)
    ck('mixed jump brackets not statewise zero',abs(bracket[2,3])>1e-3 and abs(bracket[0,2])>1e-3)
    ck('exact field drift and diffusion normalization',abs(g*lam[3]/12-c)<1e-15 and abs(2*4*g*lam[3]/12-8*c)<1e-14)
    # Semigroup covariance route, not source mixture evaluation.
    numeric={}
    for kap in (1,4):
        J=np.array([[-1.,1.],[kap*c,-float(kap)]]); Q=np.diag([2.,2*kap*c])
        S=solve_continuous_lyapunov(J,-Q); target=np.array([[1,c],[c,c]])/gamma
        ck(f'kappa{kap} covariance via independent Lyapunov solve',np.max(abs(S-target))<3e-13 and np.linalg.eigvalsh(S)[0]>0)
        C=np.array([(expm(J*t)@S)[0,0]+math.exp(-t) for t in range(5)]); R=C/C[0]; q0=math.exp(-1)
        D=R[2:]-q0*R[1:-1]; h=D[0]*D[2]-D[1]**2
        ck(f'kappa{kap} negative filtered matrix-semigroup score',h<0)
        numeric[f'kappa_{kap}']={'H_from_matrix_semigroup':float(h),'covariances':C.tolist()}
    q0=math.exp(-1)
    for alpha in (.5,1.):
        C=np.array([math.exp(-gamma*t)/gamma+(1-alpha)/alpha*math.exp(-t) for t in range(5)])
        R=C/C[0]; D=R[2:]-q0*R[1:-1]; h=D[0]*D[2]-D[1]**2
        ck(f'native filtered zero alpha{alpha}',abs(h)<2e-15)
        if alpha<1:ck('native ordinary Hankel false positive retained',C[1]*C[3]-C[2]**2>0)
    # Exact generic factor including the eliminated third pole.
    aa,bb,dd,x,y,z=sy.symbols('a b d x y z'); C=lambda j:aa*x**j+bb*y**j+dd*z**j
    D=lambda j:C(j+1)-z*C(j)
    ck('exact three-pole cancellation and signed factor',sy.expand(D(1)*D(3)-D(2)**2-aa*bb*x*y*(x-z)*(y-z)*(x-y)**2)==0)
    # Finite N=2 labelled reset lift: three residual atoms for exact moments, NOT a Gaussian discretization certificate.
    states=list(itertools.product(range(12),repeat=2)); index={x:i for i,x in enumerate(states)}; L=np.zeros((144,144)); lw=[]
    for row,xs in enumerate(states):
        n=np.bincount(xs,minlength=12); lw.append(g/(4)*n@A@n)
        for slot in range(2):
            p=softmax(g/2*A[:,xs[1-slot]])
            for new in range(12):
                if new!=xs[slot]:
                    target=list(xs); target[slot]=new; L[row,index[tuple(target)]]+=p[new]
    np.fill_diagonal(L,-L.sum(1)); pi=softmax(lw)
    ck('independent labelled N2 detailed balance',np.max(abs(pi[:,None]*L-(pi[:,None]*L).T))<2e-14)
    nu=np.array([.25,.5,.25]); K=np.kron(L-np.diag(np.diag(L)),np.tile(nu,(3,1)))+np.diag(np.repeat(np.diag(L),3))
    K+=np.kron(np.eye(144),1.7*(np.tile(nu,(3,1))-np.eye(3)))
    embed=np.kron(np.eye(144),np.ones((3,1)))
    ck('labelled reset lift exact finite intertwining',np.max(abs(K@embed-embed@L))<2e-14)
    piplus=np.kron(pi,nu)
    ck('labelled reset lift finite detailed balance',np.max(abs(piplus[:,None]*K-(piplus[:,None]*K).T))<2e-14)
    block=np.array([features[1,xs[0]] for xs in states]); v=np.repeat(block,3)
    ck('selected block is measurable in labelled lift',np.max(abs(K@v-np.repeat(L@block,3)))<2e-14)
    ck('block values not determined by total counts',abs(block[index[(0,1)]]-block[index[(1,0)]])>1)
    # Exact reference gate counterexample: real Gaussian mixture with arbitrarily small extra pole.
    z=F(1,3); r=F(1,2); eps=F(1,1000); R=lambda j:(1-eps)*z**j+eps*r**j
    ratio=R(2)/R(1); residual=R(3)-R(1)*ratio**2
    ck('strict reference equality is discontinuous at valid limit',residual>0 and residual==eps*(1-eps)*z*r*(z-r)**2/R(1))
    values=np.array([0.,2.,-3.,0.,4.]); original=np.sign(values[0])*np.sign(values[1:])
    # Include nonzero anchors, zeros and both gain signs.
    products=np.sign(values[1])*np.sign(values)
    invariant=all(np.array_equal(np.sign(gain*values[1])*np.sign(gain*values),products) for gain in (-3.,.2,7.))
    changing=values.copy(); changing[2]*=-1
    ck('fixed real gain invariance and changing sign negative control',invariant and not np.array_equal(np.sign(changing[1])*np.sign(changing),products))
    certificate=build(); assert certificate==json.loads((HERE/'MARGIN_CERTIFICATE.json').read_text())
    ck('source assumed gamma now paid by exact-input binding',F('0.16257283698729')<F(certificate['gamma'][0]) and F(certificate['gamma'][1])<F('0.16257283698731'))
    for key in numeric:
        low,high=map(F,certificate[key]['H']); value=F(numeric[key]['H_from_matrix_semigroup'])
        numeric[key]['binary64_point_inside_interval']=low<=value<=high
        numeric[key]['distance_to_interval']=float(max(low-value,value-high,F(0)))
    ck('certified negativity and separate approximate matrix agreement',all(F(certificate[key]['H'][1])<0 and numeric[key]['distance_to_interval']<2e-14 for key in numeric))
    wrong=mp.iv.mpf([1,2]); inward=public(wrong); inward[0]='1.1'
    ck('inward public endpoint rejected by exact rational test',not(F(inward[0])<=ends(wrong)[0]<=ends(wrong)[1]<=F(inward[1])))
    out['numeric_checks']=numeric; out['margin_certificate']=certificate
    if args.replay:
        tmp=Path(tempfile.mkdtemp(prefix='fin-mm-intake-'))/'campaign'; shutil.copytree(SOURCE,tmp,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1'); logs={}
        for script in ('verify_package.py','verify_mm.py'):
            run=subprocess.run(['python3',script],cwd=tmp,env=env,text=True,capture_output=True,timeout=60)
            assert run.returncode==0,(script,run.stdout,run.stderr); logs[script]={'exit_code':run.returncode,'stdout':run.stdout,'stderr':run.stderr}
        before=json.loads((SOURCE/'verification.json').read_text()); after=json.loads((tmp/'verification.json').read_text())
        ck('source bounded assertions replay',after['passed']==after['total'] and all(row['pass'] for row in after['checks']))
        out['replay']={'directory':str(tmp),'logs':logs,'source_flags':after['total'],
            'byte_equal':sha(SOURCE/'verification.json')==sha(tmp/'verification.json'),'changes':changes(before,after)}
    inputs=[SEED,Path(__file__),HERE/'certify_margin.py',SOURCE/'MANIFEST.sha256',ROOT/'FIN_MEDIATOR_MEMORY_OBSERVABILITY_NEXT_CHAIN_20261016.md',
        ROOT/'fin_rank7_followup/inputs/fin_handoff_audit/results.json']
    inputs += [ROOT/rel/'MANIFEST.sha256' for rel in prior]
    for rel in ('fin_auxiliary_mediator_review_20261016/PROOF_COMPLETIONS_AND_CORRECTIONS.md',
                'fin_auxiliary_mediator_review_20261016/POSTERIOR_RETENTION_AND_MM_SEED.md',
                'fin_dynamic_kernel_observability_review_20261015/PROOF_COMPLETIONS.md',
                'fin_kernel_observability_review_20261015/PROOF_COMPLETIONS_AND_ARITHMETIC.md',
                'fin_critical_metastability_review_20261015/PROOF_COMPLETIONS.md'):
        inputs.append(ROOT/rel)
    out['input_hashes']={str(p.relative_to(ROOT)):sha(p) for p in inputs}
    out['source_manifest_after']=manifest(SOURCE,excluded); assert out['source_manifest']==out['source_manifest_after']
    for rel in prior:assert manifest(ROOT/rel)=={k:v for k,v in prior[rel].items() if k!='saved_inputs_unchanged'}
    out['checks']=checks; out['status']='PASS_SCOPED_MM_INTAKE_CHECKS'
    args.output.write_text(json.dumps(out,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':out['status'],'bounded_checks':len(checks),'source_flags':out.get('replay',{}).get('source_flags'),
                      'source_files':out['source_manifest']['matched'],'certificate':certificate},indent=2))

if __name__=='__main__':main()
