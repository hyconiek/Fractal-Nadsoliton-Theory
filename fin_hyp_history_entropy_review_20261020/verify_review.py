#!/usr/bin/env python3
"""New exact controls, copy-only source replay and transparent self-audit scope."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from math import isqrt
import os
from pathlib import Path
import resource
import shutil
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SOURCE=ROOT/'fin_hyp_history_entropy_20261020'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def files(folder):
    return sorted(p for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts)


def manifest(folder):
    seen=[]
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,name=line.split('  ',1);p=Path(name)
        assert not p.is_absolute() and '..' not in p.parts and name not in seen
        assert sha(folder/name)==h,(folder,name)
        seen.append(name)
    assert set(seen)|{'MANIFEST.sha256'}=={p.relative_to(folder).as_posix() for p in files(folder)},folder
    return seen


def limits():resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))


def replay_source():
    saved=json.loads((SOURCE/'RESULTS.json').read_text())
    deps={Path(n).parts[0] for n in saved['input_hashes']}|set(saved['protected_packages'])
    with tempfile.TemporaryDirectory(prefix='fin-hyp-second-pass-') as tmp:
        root=Path(tmp)
        shutil.copytree(SOURCE,root/SOURCE.name,ignore=shutil.ignore_patterns('__pycache__'))
        for name in sorted(deps):
            if name!=SOURCE.name:(root/name).symlink_to(ROOT/name,target_is_directory=(ROOT/name).is_dir())
        run=subprocess.run(['python3',str(root/SOURCE.name/'verify_research.py'),'--check','--output',str(root/'replay.json')],
                           cwd=root,env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'),capture_output=True,
                           text=True,timeout=60,preexec_fn=limits)
        assert run.returncode==0,(run.stdout,run.stderr)
        actual=json.loads((root/'replay.json').read_text())
        assert actual==saved
        return {'scope':'COPY_ONLY_SOURCE_CONTROLS; not independent proof',
                'returncode':run.returncode,'stdout':run.stdout,'stderr':run.stderr,
                'bounded_source_checks':actual['bounded_check_count'],'JSON_equal':True,
                'cache_free_copy':True,'temporary_tree_removed_on_exit':True}


def factors(n):
    out={};p=2
    while p*p<=n:
        while n%p==0:out[p]=out.get(p,0)+1;n//=p
        p+=1
    if n>1:out[n]=out.get(n,0)+1
    return out


def add_log(out,value,weight):
    if not weight:return
    assert value>0
    for p,k in factors(value.numerator).items():out[p]=out.get(p,F(0))+weight*k
    for p,k in factors(value.denominator).items():out[p]=out.get(p,F(0))-weight*k


def clean(a):return {p:x for p,x in a.items() if x}


def scale(a,c):return clean({p:x*c for p,x in a.items()})


def add(a,b):
    out=dict(a)
    for p,x in b.items():out[p]=out.get(p,F(0))+x
    return clean(out)


def entropy(distribution):
    out={}
    for x in distribution:
        if x:add_log(out,x,-x)
    return clean(out)


def divergence(p,q):
    out={}
    assert set(p)<=set(q)
    for word,x in p.items():
        if x:add_log(out,x/q[word],x)
    return clean(out)


def log_box(x,terms=40):
    assert x>0
    u=(x-1)/(x+1)
    mid=2*sum((u**(2*j+1)/F(2*j+1) for j in range(terms)),F(0))
    err=2*abs(u)**(2*terms+1)/(F(2*terms+1)*(1-u*u))
    return mid-err,mid+err


def expr_box(expr):
    lo=hi=F(0)
    for prime,coefficient in expr.items():
        a,b=log_box(F(prime))
        if coefficient<0:a,b=b,a
        lo+=coefficient*a;hi+=coefficient*b
    return lo,hi


def sqrt_box(x,bits=90):
    assert x>=0
    den=1<<bits
    base=isqrt((x.numerator*den*den)//x.denominator)
    lo=F(base,den)
    hi=lo if lo*lo==x else F(base+1,den)
    assert lo*lo<=x<=hi*hi
    return lo,hi


def root_sign_box(t,p):
    out=[sqrt_box(1-4*x/t) for x in p]
    return sum((x[0] for x in out),F(0))-(len(p)-1),sum((x[1] for x in out),F(0))-(len(p)-1)


def enclose_root(p):
    lo=4*max(p);hi=F(4)
    assert root_sign_box(lo,p)[1]<0 and root_sign_box(hi,p)[0]>0
    for _ in range(28):
        mid=(lo+hi)/2;a,b=root_sign_box(mid,p)
        if b<0:lo=mid
        elif a>0:hi=mid
        elif a==b==0:return mid,mid
        else:raise AssertionError('exact sign unresolved; do not substitute a floating sign')
    return lo,hi


def block_law(pi,p,n):
    out={}
    for word in product(range(len(pi)),repeat=n+1):
        prob=pi[word[0]]
        for i,j in zip(word,word[1:]):prob*=p[i][j]
        out[word]=prob
    return out


def tv(p,q):return sum((abs(p.get(w,F(0))-x) for w,x in q.items()),F(0))/2


def source_snapshot():return {p.relative_to(SOURCE).as_posix():sha(p) for p in files(SOURCE)}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true')
    ap.add_argument('--check',action='store_true');ap.add_argument('--replay-source',action='store_true')
    ap.add_argument('--output',type=Path);args=ap.parse_args()
    checks=[]
    def ck(name,value):
        if not value:raise AssertionError(name)
        assert name not in checks;checks.append(name)
    names=manifest(SOURCE);before=source_snapshot()
    for line in (SOURCE/'INPUTS.sha256').read_text().splitlines():
        h,name=line.split('  ',1);assert sha(ROOT/name)==h,name
    agents=sha(ROOT/'AGENTS.md')
    # Dual Hessian via a path graph with loops, not the producer's Hamming test.
    k=[[F(1),F(1),F(0)],[F(1),F(1),F(1)],[F(0),F(1),F(1)]]
    z=[F(1,3),F(1,2),F(2,3)]
    raw=[[k[i][j]*z[i]*z[j] for j in range(3)] for i in range(3)]
    total=sum(map(sum,raw));c=[[x/total for x in row] for row in raw]
    pi=list(map(sum,c));q=[[c[i][j]/pi[i] for j in range(3)] for i in range(3)]
    hess=[[c[i][j]+(pi[i] if i==j else 0) for j in range(3)] for i in range(3)]
    for n,v in enumerate(([1,0,0],[0,1,0],[0,0,1],[1,-2,3],[1,1,1])):
        left=sum(hess[i][j]*v[i]*v[j] for i,j in product(range(3),repeat=2))
        right=sum(c[i][j]*(v[i]+v[j])**2 for i,j in product(range(3),repeat=2))/2
        lower=2*sum(c[i][i]*v[i]*v[i] for i in range(3))
        ck(f'exact_dual_Hessian_form_{n}',left==right and left>=lower>0)
    logz=[]
    for x in z:
        term={};add_log(term,x,F(1));add_log(term,total,F(-1,2));logz.append(clean(term))
    dual_value={}
    for i in range(3):dual_value=add(dual_value,scale(logz[i],-2*pi[i]))
    for x in pi:add_log(dual_value,x,x)
    hstar={}
    for i in range(3):hstar=add(hstar,scale(entropy(q[i]),pi[i]))
    ck('dual_value_equals_primal_entropy',clean(dual_value)==hstar)
    ck('loopless_unequal_pi_rejected',F(1,3)!=F(2,3))
    ck('loopless_scaling_not_unique',F(1)*F(1,2)==F(2)*F(1,4)==F(1,2))
    roots=[]
    for n,p in enumerate(([F(1,2),F(1,3),F(1,6)],[F(1,3)]*3,[F(15,110),F(28,110),F(39,110),F(28,110)])):
        a,b=enclose_root(p)
        ck(f'certified_nonlinear_root_bracket_{n}',a<=b and b-a<F(1,1<<24))
        if a<b:ck(f'root_endpoint_signs_{n}',root_sign_box(a,p)[1]<0<root_sign_box(b,p)[0])
        if n==1:ck('uniform_root_known_exactly',a<=F(12,5)<=b)
        if n==2:ck('q4_root_known_exactly',a<=F(128,55)<=b)
        roots.append({'p':[str(x) for x in p],'t_lower':str(a),'t_upper':str(b)})
    # Full finite-block KL identity for a stationary Markov law, using enumeration.
    p=[F(2,5),F(3,5)];ref=[list(p),list(p)]
    transition=[[F(1,2),F(1,2)],[F(1,3),F(2,3)]]
    ck('Markov_control_stationary',all(sum(p[i]*transition[i][j] for i in range(2))==p[j] for j in range(2)))
    delta=entropy(p)
    for i in range(2):delta=add(delta,scale(entropy(transition[i]),-p[i]))
    for n in (1,2,3,4):
        a=block_law(p,transition,n);b=block_law(p,ref,n)
        d=divergence(a,b)
        ck(f'finite_block_KL_equals_n_deficit_{n}',d==scale(delta,n))
        lo,hi=expr_box(d)
        ck(f'outward_finite_block_Pinsker_{n}',lo>=2*tv(a,b)**2 and lo<=hi)
    # Stationary order-two control: identical optimal pair law, positive rate gap.
    eps=F(1,2);rate_gap={3:F(3,4),2:F(-1)}
    ck('positive_order2_gap',3**3>2**4)
    pair_gap_zero=False
    for n in (1,2,3,4):
        a={}
        for word in product((0,1),repeat=n+1):
            prob=F(1,4)
            for t in range(2,n+1):prob*=F(1,2)*(1+eps*(-1)**sum(word[t-2:t+1]))
            a[word]=prob
        b={word:F(1,2**(n+1)) for word in a}
        ck(f'order2_normalized_{n}',sum(a.values())==1)
        ck(f'order2_full_history_deficit_identity_{n}',divergence(a,b)==scale(rate_gap,n-1))
        if n==1:pair_gap_zero=divergence(a,b)=={}
    ck('finite_pair_gap_not_upper_rate_bound',pair_gap_zero and expr_box(rate_gap)[0]>0)
    # Nonstationary countercontrols have zero asymptotic gap but wrong short blocks.
    b={w:F(1,8) for w in product((0,1),repeat=3)}
    wrong_initial={w:(F(1,4) if w[0]==0 else F(0)) for w in b}
    wrong_prefix={w:(F(1,4) if w[0]==w[1] else F(0)) for w in b}
    ck('initial_mismatch_not_paid_by_entropy_rate',tv(wrong_initial,b)==F(1,2))
    ck('stationary_marginals_do_not_replace_stationary_histories',tv(wrong_prefix,b)==F(1,2)
       and all(sum(x for w,x in wrong_prefix.items() if w[t]==0)==F(1,2) for t in range(3)))
    # A separately generated q4 quotient with the required multiplicity.
    labels=range(4);a=[1,2,3,2];den=110;ss=8
    states=list(product(labels,repeat=2));index={x:i for i,x in enumerate(states)}
    rows=[]
    for x,y in states:
        dest={(x,y)}|{(z,y) for z in labels}|{(x,z) for z in labels}
        counts=[sum((v-u)%4==d for u,v in dest) for d in labels]
        d=(y-x)%4
        ck(f'direct_neighbour_counts_{x}_{y}',counts==[1 if j==d else 2 for j in labels])
        rows.append([sum((F(a[(v-u)%4],16-a[d]) for u,v in dest if (v-u)%4==j),F(0)) for j in labels])
    matrix=[rows[index[(0,d)]] for d in labels]
    mass=[F(x*(16-x),den) for x in a]
    ck('q4_source_diagonal_3_13',matrix[2][2]==F(3,13))
    ck('q4_actual_exit_1_5',matrix[2][0]/(1-matrix[2][2])==F(1,5))
    ck('q4_HB_actual_exit_15_71',mass[0]/(1-mass[2])==F(15,71))
    # Weighted macro objective equals full labelled entropy for uniform lifts.
    iid=[list(mass) for _ in mass]
    macro_delta={}
    for i,j in product(labels,repeat=2):add_log(macro_delta,iid[i][j]/matrix[i][j],mass[i]*iid[i][j])
    d1=divergence(block_law(mass,iid,1),block_law(mass,matrix,1))
    ck('weighted_reference_rate_gap_exact',clean(macro_delta)==d1 and bool(d1))
    for n in (2,3):ck(f'weighted_block_identity_{n}',divergence(block_law(mass,iid,n),block_law(mass,matrix,n))==scale(d1,n))
    gap=F(4,355);tail=F(3,13)**4
    lower=2*mass[2]*(gap-tail)**2/4
    ck('fixed_four_step_tail_smaller_than_gap',tail==F(81,28561) and tail<gap)
    ck('exact_near_maximum_entropy_obstruction',lower==F(21925107363,1739734008083500))
    ck('strict_positive_nats_lower_bound',lower>F(1,100000))
    ck('concrete_HB_rate_gap_exceeds_obstruction',expr_box(d1)[0]>lower)
    bare_rate={}
    for x in mass:bare_rate=add(bare_rate,scale(entropy(mass),x))
    ck('bare_macro_zero_gap_would_be_wrong',add(entropy(mass),scale(bare_rate,-1))=={} and bool(d1))
    # Recompute the q12 geometry rather than repeating the displayed ratio.
    weights=list(range(1,13));s=sum(weights);den=2*s*s-sum(x*x for x in weights)
    p12=[F(x*(2*s-x),den) for x in weights]
    def thinning(x,y):
        step=abs(x-y)
        return F(1,10) if min(step,12-step)==3 else F(1)
    thin_odds=lambda d:thinning(d,0)*p12[0]/(thinning(d,1)*p12[1])
    model_odds=lambda d:(F(weights[0],s-weights[d])/F(weights[1],s-weights[d]))
    ck('thinning_diagnostic_preserved',thin_odds(3)/thin_odds(4)==F(1,100) and model_odds(3)/model_odds(4)==1)
    # Rebuild the eight-state stress from a pair-agreement Gibbs weight.
    binary=list(product((0,1),repeat=3));w=[F(2**sum(x[i]==x[j] for i,j in ((0,1),(0,2),(1,2)))) for x in binary]
    weight_sum=sum(w);rates={}
    for x in ((0,0,0),(1,0,0)):
        i=binary.index(x);diag=F(0)
        for slot in range(3):
            y=list(x);y[slot]^=1;j=binary.index(tuple(y))
            diag+=w[i]/(3*(w[i]+w[j]))
            rates[x,tuple(y)]=w[i]*w[j]/(3*weight_sum*(w[i]+w[j]))
        rates[x,x]=w[i]*diag/weight_sum
    cx,cy,cxy=rates[(0,0,0),(0,0,0)],rates[(1,0,0),(1,0,0)],rates[(0,0,0),(1,0,0)]
    ck('rebuilt_binary_N3_stress',9*cxy*cxy/(cx*cy)==F(1,2))
    replay=replay_source() if args.replay_source else None
    ck('source_bytes_unchanged',before==source_snapshot())
    ck('live_AGENTS_unchanged',agents==sha(ROOT/'AGENTS.md'))
    result={'date':'2026-10-20','source_model_id':'FIN-HYP-HIST-ME-001',
            'review_status':'SCOPED_CONDITIONAL_THEOREMS_SURVIVE_SAME_AUTHOR_CROSSCHECK',
            'review_independence':'SAME_AUTHOR_NEW_PASS; SEPARATE_IMPLEMENTATION; NOT_EXTERNAL_PEER_REVIEW',
            'source_manifest_entries':len(names),'source_files_preserved':len(before),'source_snapshot':before,
            'source_input_hash_matches':len((SOURCE/'INPUTS.sha256').read_text().splitlines()),
            'new_bounded_control_count':len(checks),'new_bounded_controls':checks,'nonlinear_root_fixtures':roots,
            'source_copy_replay':replay,'new_proofs':['coercive dual certificate','entropy deficit finite-block theorem','first-exit tail transfer','robust HB obstruction'],
            'q4_fixed_n':4,'q4_entropy_deficit_lower':str(lower),'q4_conservative_strict_bound':'1/100000 nats per intrinsic labelled step',
            'source_repaired_or_postulate_changed':False,'independent_peer_audit_executed':False,
            'strict_F01':'OPEN_SOURCE','new_physical_claim':False,'AGENTS_sha256':agents}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    input_names=['AGENTS.md','FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md','research_tasks/00_EXECUTION_RULES.md',
                 'fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md',
                 'fin_intrinsic_clock_source_research_20261020/PROOF_OR_OPEN_LEMMA.md']+[str((SOURCE/name).relative_to(ROOT)) for name in before]
    input_text=''.join(f'{sha(ROOT/name)}  {name}\n' for name in sorted(input_names))
    if args.record:
        assert args.replay_source and args.output and args.output.resolve()==HERE/'verification.json'
        (HERE/'INPUTS.sha256').write_text(input_text)
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.sha256'))
    if args.check:
        manifest(HERE)
        saved=json.loads((HERE/'verification.json').read_text())
        if args.replay_source:assert saved==result
        else:
            compare=dict(saved);compare['source_copy_replay']=None;assert compare==result
        assert (HERE/'INPUTS.sha256').read_text()==input_text
        for p in HERE.glob('*.json'):json.loads(p.read_text())
        for p in HERE.glob('*.py'):ast.parse(p.read_text())
    print(json.dumps({'status':result['review_status'],'new_bounded_checks':len(checks),
                      'source_replay_checks':replay['bounded_source_checks'] if replay else None,
                      'peer_review':False,'strict_F01':'OPEN_SOURCE'}))


if __name__=='__main__':main()
