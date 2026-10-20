#!/usr/bin/env python3
"""Source-preserving ES intake and independent bounded causal-filter controls."""
import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import importlib.util
import itertools
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_evolution_source_research_20261020'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def snapshot(folder):
    return {str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}


def manifest():
    names=[]
    for line in (SOURCE/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split(None,1)
        assert not Path(n).is_absolute() and '..' not in Path(n).parts and n not in names
        assert sha(SOURCE/n)==h,n
        names.append(n)
    assert set(names)|{'MANIFEST.sha256'}==set(snapshot(SOURCE))
    return len(names)


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))


def replay():
    logs={}
    with tempfile.TemporaryDirectory(prefix='fin-causal-intake-') as tmp:
        replay_root=Path(tmp)/'repo'; replay_root.mkdir()
        # The source resolves current mathematical inputs against its parent.
        # Only its own directory is copied/writable by the source scripts;
        # all other inputs are symlinked and no old writer is invoked.
        for p in ROOT.iterdir():
            if p.name!=SOURCE.name: (replay_root/p.name).symlink_to(p,target_is_directory=p.is_dir())
        cp=replay_root/SOURCE.name
        shutil.copytree(SOURCE,cp,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for cmd in [[sys.executable,'check_package.py','--check'],[sys.executable,'verify_research.py']]:
            p=subprocess.run(cmd,cwd=cp,text=True,capture_output=True,timeout=60,preexec_fn=limits,
                             env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
            assert p.returncode==0,(cmd,p.stderr)
            logs[' '.join(cmd[1:])]={'returncode':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
        logs['result_comparison']={'byte_identical':sha(cp/'RESULTS.json')==sha(SOURCE/'RESULTS.json'),
                                   'JSON_equal':json.loads((cp/'RESULTS.json').read_text())==json.loads((SOURCE/'RESULTS.json').read_text())}
    return logs


def prod(vals):
    out=F(1)
    for v in vals: out*=v
    return out


def hidden_source(eps):
    p=[F(1,2),F(1,3),F(1,6)]; h=[F(1,2),-F(1,2),-F(1,2)]; lam=F(2)
    blocks=list(itertools.product(range(3),repeat=3))
    prob={w:prod(p[i] for i in w)*(1+eps*prod(h[i] for i in w)) for w in blocks}
    states=[(w,k) for w in blocks for k in range(3)]
    R=[w[k] for w,k in states]; index={z:i for i,z in enumerate(states)}
    S=[]
    for w,k in states:
        S.append({index[(w,k+1)]:F(1)} if k<2 else {index[(v,0)]:pv for v,pv in prob.items()})
    mu=[prob[w]/3 for w,k in states]
    rates=[]
    for i,row in enumerate(S):
        a=[F(0)]*3
        for j,pv in row.items():
            if R[j]!=R[i]: a[R[j]]+=lam*pv
        rates.append(a)
    Q=[[lam*(p[j]-int(i==j)) for j in range(3)] for i in range(3)]
    return {'p':p,'h':h,'lambda':lam,'blocks':blocks,'states':states,'R':R,'S':S,'mu':mu,'a':rates,'Q':Q}


def advance(v,S,allowed=None):
    out=[F(0)]*len(v)
    for i,vi in enumerate(v):
        if not vi: continue
        for j,p in S[i].items():
            if allowed is None or j in allowed: out[j]+=vi*p
    return out


def initial(s,i):
    return [v/s['p'][i] if s['R'][z]==i else F(0) for z,v in enumerate(s['mu'])]


def killed_polynomial(v,S,allowed,x,terms):
    out=v.copy(); power=v.copy(); coef=F(1)
    for n in range(1,terms+1):
        power=advance(power,S,allowed); coef*=x/n
        out=[a+coef*b for a,b in zip(out,power)]
    return out


def normalize(v):
    total=sum(v)
    if not total: raise ValueError('zero probability observation')
    return [x/total for x in v]


def jump_update(s,eta,target):
    out=[F(0)]*len(eta)
    for i,vi in enumerate(eta):
        for j,p in s['S'][i].items():
            if s['R'][j]==target and s['R'][i]!=target: out[j]+=vi*s['lambda']*p
    return normalize(out)


def log_bounds(x,terms=24):
    x=F(x)
    assert x>0
    if x<1:
        lo,hi=log_bounds(1/x,terms); return -hi,-lo
    n=0; y=x
    while y>=2: y/=2; n+=1
    def part(v):
        z=(v-1)/(v+1)
        a=2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
        b=2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return a,a+b
    l,u=part(2); a,b=part(y)
    return n*l+a,n*u+b


def phi_bounds(a,q):
    a,q=F(a),F(q)
    if a<0 or q<=0: raise ValueError('invalid entropy support')
    if not a: return q,q
    lo,hi=log_bounds(a/q)
    return a*lo-a+q,a*hi-a+q


def controls():
    checks=[]
    def ck(name,condition):
        assert bool(condition),name
        checks.append(name)
    contexts={}
    for eps in [F(0),F(1,2)]:
        tag='iid' if eps==0 else 'memory'
        s=hidden_source(eps); contexts[tag]=s
        ck(tag+'_stationary_positive_source',sum(s['mu'])==1 and min(s['mu'])>0 and advance(s['mu'],s['S'])==s['mu'])
        ck(tag+'_generator_row_stochastic',all(sum(row.values())==1 for row in s['S']))
        for i in range(3):
            eta=initial(s,i)
            ck(tag+f'_current_conditional_mean_rates_{i}',all(sum(eta[z]*s['a'][z][j] for z in range(len(eta)))==s['Q'][i][j]
               for j in range(3) if j!=i))
        raw=sum(s['mu'][z]*sum(abs(s['a'][z][j]-s['Q'][s['R'][z]][j]) for j in range(3) if j!=s['R'][z])
                for z in range(len(s['states'])))
        ck(tag+'_raw_source_defect_exact_28_over_27',raw==F(28,27))
        # Full-history drift witness via independent sparse source powers.
        signed=[s['mu'][z]*s['h'][s['R'][z]] for z in range(len(s['states']))]
        f=[]
        for z in range(len(s['states'])):
            i=s['R'][z]
            src=sum(s['a'][z][j]*(s['h'][j]-s['h'][i]) for j in range(3) if j!=i)
            tgt=sum(s['Q'][i][j]*s['h'][j] for j in range(3))
            f.append(s['h'][i]*(src-tgt))
        v=sum(s['p'][i]*s['h'][i]**2 for i in range(3))
        for n in range(8):
            coef=sum(x*y for x,y in zip(signed,f))
            ck(tag+f'_drift_witness_coefficient_{n}',coef==(s['lambda']*eps*v**3/3 if n==1 else 0))
            signed=advance(signed,s['S'])
    s=contexts['iid']
    for i in range(3):
        eta=initial(s,i); allowed={z for z,r in enumerate(s['R']) if r==i}; power=eta
        for n in range(9):
            ck(f'iid_killed_filter_rate_coefficient_{i}_{n}',all(sum(power[z]*s['a'][z][j] for z in range(len(power)))
               ==s['Q'][i][j]*sum(power) for j in range(3) if j!=i))
            power=advance(power,s['S'],allowed)
        eta_nojump=normalize(killed_polynomial(eta,s['S'],allowed,F(1,2),10))
        for j in range(3):
            if j==i: continue
            post=jump_update(s,eta_nojump,j)
            ck(f'iid_Bayes_jump_{i}_{j}',sum(post)==1 and all(x==0 for z,x in enumerate(post) if s['R'][z]!=j))
            future=normalize(killed_polynomial(post,s['S'],{z for z,r in enumerate(s['R']) if r==j},F(1,3),6))
            ck(f'iid_actual_filter_preserves_native_rates_{i}_{j}',all(sum(future[z]*s['a'][z][k] for z in range(len(future)))
               ==s['Q'][j][k] for k in range(3) if k!=j))
    mem=contexts['memory']; i=0
    allowed={z for z,r in enumerate(mem['R']) if r==i}
    eta=initial(mem,i)
    coeff2=sum(advance(advance(eta,mem['S'],allowed),mem['S'],allowed))
    defect=coeff2-mem['p'][i]**2
    ck('waiting_law_second_refresh_coefficient',defect==F(1,2)*mem['p'][i]**2*mem['h'][i]**3/3)
    ck('waiting_law_time_quadratic_germ',defect*mem['lambda']**2/2==F(1,96))
    survival_derivative=-mem['lambda']+mem['lambda']*sum(advance(eta,mem['S'],allowed))
    ck('null_hidden_jumps_must_be_retained',survival_derivative==-mem['lambda']*(1-mem['p'][i]) and survival_derivative!=-mem['lambda'])
    bad_eta=[F(0)]*len(eta); bad_z=mem['states'].index(((0,0,0),0)); bad_eta[bad_z]=1
    try: jump_update(mem,bad_eta,1)
    except ValueError: rejected=True
    else: rejected=False
    ck('impossible_instant_jump_update_rejected',rejected)
    ck('same_source_arbitrary_start_not_native_closure',mem['a'][bad_z][1]!=mem['Q'][0][1])
    # An exact-native reset lift is a reused positive control for strong lumpability.
    p=mem['p']; q=mem['Q']; states=list(itertools.product(range(3),range(2)))
    L=[[F(0)]*6 for _ in range(6)]
    for x,(i,h) in enumerate(states):
        for y,(j,k) in enumerate(states):
            if j!=i: L[x][y]=q[i][j]/2
            elif h!=k: L[x][y]=F(3)
        L[x][x]=-sum(L[x])
    ck('reused_exact_native_lift_source_fibre_defect_zero',all(sum(L[x][y] for y,(j,k) in enumerate(states) if j==target)==q[i][target]
       for x,(i,h) in enumerate(states) for target in range(3) if target!=i))
    for a,q0 in [(F(0),F(1,3)),(F(1,6),F(1,3)),(F(2,3),F(1,3)),(F(3),F(2)),(F(1,3),F(1,3))]:
        lo,hi=phi_bounds(a,q0)
        ck('entropy_nonnegative_and_quadratic_upper_'+str(a)+'_'+str(q0),lo>=0 and hi<=(a-q0)**2/q0)
    lo,hi=phi_bounds(1,1)
    b0=phi_bounds(0,1); b2=phi_bounds(2,1)
    ck('filter_Jensen_can_strictly_improve_source_cost',hi==0 and (b0[0]+b2[0])/2>hi)
    try: phi_bounds(1,0)
    except ValueError: rejected=True
    else: rejected=False
    ck('entropy_target_zero_support_rejected',rejected)
    ck('visible_defect_rational_lower',F(1,2)*F(1,4)**3/F(6)*F(1,4)==F(1,3072))
    ck('visible_entropy_rational_lower',F(1,2)*F(1,3456)**2==F(1,23887872))
    # Versioned actual labelled N2 parent is checked numerically, not sourced anew.
    import numpy as np
    spec=importlib.util.spec_from_file_location('seed_cf',ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py')
    seed=importlib.util.module_from_spec(spec); spec.loader.exec_module(seed)
    _,_,_,_,A7=seed.build_rank7(); g=seed.G_FROZEN
    labels=list(itertools.product(range(12),repeat=2)); p12=seed.pair_difference_distribution(A7,g)
    maxerr=0.0
    for a,b in labels:
        current=(b-a)%12; rates=np.zeros(12)
        for old,other,slot in [(a,b,0),(b,a,1)]:
            f=(g/2)*A7[:,other]; w=np.exp(f-max(f)); w/=sum(w)
            for new in range(12):
                if new==old: continue
                d=(b-new)%12 if slot==0 else (new-a)%12
                if d!=current: rates[d]+=w[new]
        expected=2*p12.copy(); expected[current]=0
        maxerr=max(maxerr,float(max(abs(rates-expected))))
    ck('actual_declared_N2_parent_native_fibre_rates_nominal',maxerr<2e-14)
    return checks,{'raw_source_defect_rate':'28/27','visible_L1_lower_at_T1':'1/3072',
                   'visible_KL_lower_at_T1':'1/23887872','waiting_time_quadratic_germ':'1/96',
                   'native_N2_max_fibre_rate_error_nominal':maxerr,
                   'arithmetic':'Fraction controls; native-parent bound is numerical only'}


def main():
    p=argparse.ArgumentParser(); p.add_argument('--replay',action='store_true'); p.add_argument('--output',type=Path,default=HERE/'verification.json'); a=p.parse_args()
    before=snapshot(SOURCE); guard=sha(ROOT/'AGENTS.md'); count=manifest()
    old_inputs=[]
    for line in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,n=line.split(None,1); assert sha(ROOT/n)==h,n; old_inputs.append(n)
    logs=replay() if a.replay else {}; checks,details=controls()
    assert before==snapshot(SOURCE) and sha(ROOT/'AGENTS.md')==guard
    data={'date':'2026-10-20','status':'PASS_SCOPED_ES_INTAKE_AND_NEW_CAUSAL_TRANSFER_CONTROLS',
          'source_manifest_entries':count,'source_preserved':before,'source_file_count':len(before),
          'source_input_matches':len(old_inputs),'source_input_paths':old_inputs,'copy_only_replay':logs,
          'bounded_controls':checks,'bounded_control_count':len(checks),'details':details,'live_AGENTS_sha256':guard,
          'new_results':'conditional filter, full-path comparison and ES visible-defect witness',
          'F01':'OPEN_SOURCE','physical_law_added':False,'formal_proof':False}
    a.output.write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps({'checks':len(checks),'source_manifest_entries':count,'source_files':len(before),
                      'source_input_matches':len(old_inputs),'output':str(a.output)},indent=2))


if __name__=='__main__': main()
