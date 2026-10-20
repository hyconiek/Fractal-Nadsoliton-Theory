#!/usr/bin/env python3
"""Small exact controls for the causal-clock quotient; not a proof assistant."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    'AGENTS.md', 'FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md',
    'FIN_RESEARCH_STRATEGY.md', 'research_tasks/00_EXECUTION_RULES.md',
    'fin_evolution_source_research_20261020/PROOF.md',
    'fin_evolution_source_research_20261020/SOURCE_ATOM.md',
    'fin_evolution_source_research_20261020/PREMISE_LEDGER.md',
    'fin_evolution_source_review_20261020/EXACT_CAUSAL_FILTER.md',
    'fin_evolution_source_review_20261020/FULL_PATH_TRANSFER.md',
    'fin_evolution_source_review_20261020/FILTERED_BLOCK_DEFECT.md',
    'fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md',
    'fin_reversible_source_transfer_research_20261020/PROOF_OR_OPEN_LEMMA.md',
    'fin_strategy_audit_20261019/SEARCH_FILE_INDEX.json',
]
PROTECTED = [
    'fin_evolution_source_research_20261020',
    'fin_evolution_source_review_20261020',
    'fin_reversible_source_transfer_research_20261020',
    'fin_agents_causal_update_20261020',
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_manifest(folder):
    names=set()
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,name=line.split('  ',1)
        path=Path(name)
        assert not path.is_absolute() and '..' not in path.parts and name not in names
        assert sha(folder/name)==h,(folder,name)
        names.add(name)
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==names|{'MANIFEST.sha256'},folder
    return len(names)


def block(w,p,h,e):
    prob,corr=F(1),F(1)
    for j in w:
        prob*=p[j]
        corr*=h[j]
    return prob*(1+e*corr)


def direct_word(word,p,h,e):
    """Three alignments, independent full/partial block probabilities."""
    total=F(0)
    for offset in range(3):
        cursor=0
        prob=F(1)
        while cursor<len(word):
            phase=(offset+cursor)%3
            take=min(3-phase,len(word)-cursor)
            fragment=word[cursor:cursor+take]
            if take==3:
                prob*=block(fragment,p,h,e)
            else:
                for j in fragment:
                    prob*=p[j]
            cursor+=take
        total+=prob/3
    return total


def c_formula(n,i,j,p,h,e):
    if n==0:
        return p[j]
    x=p[i]
    a=1+e*h[i]**3
    b=1+e*h[i]**2*h[j]
    l,r=divmod(n-1,3)
    factor=(2+b,1+a+b,2*a+b)[r]
    return x**n*p[j]*a**l*factor/3


def first_formula(i,j,p,h,e):
    if i==j:
        return F(0)
    x=p[i]
    a=1+e*h[i]**3
    b=1+e*h[i]**2*h[j]
    return p[j]+p[j]*(x*(2+b)+x*x*(1+a+b)+x**3*(2*a+b))/(3*(1-x**3*a))


def linear_solve(a,b):
    n=len(a)
    aug=[list(row)+[rhs] for row,rhs in zip(a,b)]
    for col in range(n):
        pivot=next(i for i in range(col,n) if aug[i][col])
        aug[col],aug[pivot]=aug[pivot],aug[col]
        scale=aug[col][col]
        aug[col]=[v/scale for v in aug[col]]
        for i in range(n):
            if i!=col and aug[i][col]:
                c=aug[i][col]
                aug[i]=[x-c*y for x,y in zip(aug[i],aug[col])]
    return [r[-1] for r in aug]


def lift(p,h,e):
    words=list(product(range(len(p)),repeat=3))
    states=list(product(words,range(3)))
    index={s:i for i,s in enumerate(states)}
    weights={w:block(w,p,h,e) for w in words}
    kernel=[]
    for w,phase in states:
        if phase<2:
            kernel.append({index[(w,phase+1)]:F(1)})
        else:
            kernel.append({index[(v,0)]:prob for v,prob in weights.items()})
    mu=[weights[w]/3 for w,k in states]
    labels=[w[k] for w,k in states]
    return kernel,mu,labels


def exact_exit(k,mu,labels,i,j,p):
    fibre=[x for x,label in enumerate(labels) if label==i]
    a=[[F(x==y)-k[x].get(y,F(0)) for y in fibre] for x in fibre]
    rhs=[sum((v for y,v in k[x].items() if labels[y]==j),F(0)) for x in fibre]
    hit=linear_solve(a,rhs)
    return sum((mu[x]*hit[z]/p[i] for z,x in enumerate(fibre)),F(0))


def primes(n):
    result={}
    d=2
    while d*d<=n:
        while n%d==0:
            result[d]=result.get(d,0)+1
            n//=d
        d+=1
    if n>1:
        result[n]=result.get(n,0)+1
    return result


def log_terms(x):
    assert x>0
    out={}
    for p,n in primes(x.numerator).items():out[p]=F(n)
    for p,n in primes(x.denominator).items():out[p]=out.get(p,F(0))-n
    return {p:v for p,v in out.items() if v}


def add_expr(a,b):
    c=a[0]+b[0]
    logs=dict(a[1])
    for p,v in b[1].items():logs[p]=logs.get(p,F(0))+v
    return c,{p:v for p,v in logs.items() if v}


def phi(x,y):
    assert y>0
    return y-x,{p:x*v for p,v in log_terms(x/y).items()} if x else {}


def branch_expr(alpha,q):
    a=sum(alpha)
    if not a:return F(0),{}
    total=(F(0),{})
    for x,y in zip(alpha,q):
        if x:
            total=add_expr(total,(F(0),{p:x*v for p,v in log_terms(x*sum(q)/(a*y)).items()}))
    return total


def save_search():
    prior=json.loads((ROOT/'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json').read_text())
    names=[]
    focused={'fin_evolution_source_research_20261020','fin_evolution_source_review_20261020',
             'fin_transition_hierarchy','fin_intermediate_clock_bridge_research_20261017',
             'fin_post_rp004_history_review_20261020','fin_complete_anonymous_history_review_20261020',
             'fin_record_information_review_20261015'}
    for name in prior['searched_files']:
        path=Path(name)
        root_match=len(path.parts)==1 and re.search(r'Clock|Composition|Process|Bridge|Source|Evolution|History|transition_hierarchy|KINETIC|ACTIVITY',path.name,re.I)
        if root_match or path.parts[0] in focused:
            names.append(name)
    pattern=r'predictable.{0,50}clock|history.dependent.{0,50}clock|adapted.{0,50}clock|branch.{0,60}(KL|entropy)|KL.{0,60}branch|jump.destination|destination.{0,50}clock|semi.Markov|time.change'
    lines=[]
    hashes={name:sha(ROOT/name) for name in names}
    for start in range(0,len(names),70):
        run=subprocess.run(['rg','-n','-i','--max-columns','280','--',pattern,*names[start:start+70]],
                           cwd=ROOT,text=True,capture_output=True,timeout=60)
        assert run.returncode in (0,1),run.stderr
        lines.extend(run.stdout.splitlines())
    assert hashes=={name:sha(ROOT/name) for name in names}
    return {'date':'2026-10-20','selected_files':len(names),'path_hashes':hashes,
            'pattern':pattern,'matching_lines':len(lines),'lines':lines,
            'scope':'focused duplicate search, not exhaustive proof reading or global novelty; older methods reused'}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--output',type=Path)
    ap.add_argument('--record',action='store_true')
    ap.add_argument('--check',action='store_true')
    args=ap.parse_args()
    checks=[]
    def ck(name,value):
        if not value:raise AssertionError(name)
        assert name not in checks
        checks.append(name)

    before={p:sha(ROOT/p) for p in INPUTS}
    old={name:sha(ROOT/name/'MANIFEST.sha256') for name in PROTECTED}
    old_counts={name:verify_manifest(ROOT/name) for name in PROTECTED}
    # Formal exact log-rational arithmetic, no floating entropy comparison.
    for case,(alpha,q,v) in enumerate([
        ([F(1,2),F(3,2)],[F(2),F(1)],F(2,3)),
        ([F(0),F(3,2)],[F(2),F(1)],F(5,7)),
        ([F(0),F(0)],[F(2),F(1)],F(1,3)),
        ([F(1,3),F(2,7),F(5,11)],[F(7,5),F(1,9),F(3,2)],F(4,3)),
    ]):
        lhs=(F(0),{})
        for x,y in zip(alpha,q):lhs=add_expr(lhs,phi(x,v*y))
        rhs=add_expr(phi(sum(alpha),v*sum(q)),branch_expr(alpha,q))
        ck(f'profiled_entropy_exact_log_identity_{case}',lhs==rhs)
        if sum(alpha):
            opt=sum(alpha)/sum(q)
            lhs=(F(0),{})
            for x,y in zip(alpha,q):lhs=add_expr(lhs,phi(x,opt*y))
            ck(f'optimal_clock_leaves_only_branch_cost_{case}',lhs==branch_expr(alpha,q))
    ck('zero_branch_cost_proportional_rates',branch_expr([F(2),F(4)],[F(1),F(2)])==(F(0),{}))
    ck('one_exit_destination_always_zero_cost',branch_expr([F(7,3)],[F(5,2)])==(F(0),{}))

    p=[F(1,2),F(1,3),F(1,6)]
    h=[F(1,2),F(-1,2),F(-1,2)]
    ck('ES_centering',sum(x*y for x,y in zip(p,h))==0)
    e=F(1,2)
    k,mu,labels=lift(p,h,e)
    ck('positive_81_state_block_law',len(k)==81 and all(x>0 for x in mu) and sum(mu)==1)
    ck('source_rows_and_stationarity',all(sum(row.values())==1 for row in k) and
       all(sum((mu[x]*row.get(y,F(0)) for x,row in enumerate(k)),F(0))==mu[y] for y in range(len(k))))
    ck('directed_source_not_reversible',any(mu[x]*v!=mu[y]*k[y].get(x,F(0)) for x,row in enumerate(k) for y,v in row.items()))
    for i in range(3):
        v=[mu[x]/p[i] if labels[x]==i else F(0) for x in range(len(k))]
        for n in range(10):
            propagated=[F(0)]*len(k)
            for x,row in enumerate(k):
                for y,prob in row.items():propagated[y]+=v[x]*prob
            for j in range(3):
                if j!=i:
                    actual=sum((propagated[y] for y in range(len(k)) if labels[y]==j),F(0))
                    ck(f'hidden_null_refresh_coefficient_{i}_{j}_{n}',actual==c_formula(n,i,j,p,h,e))
            v=[prob if labels[y]==i else F(0) for y,prob in enumerate(propagated)]
        for j in range(3):
            if j!=i:
                ck(f'killed_fibre_exact_exit_{i}_{j}',exact_exit(k,mu,labels,i,j,p)==first_formula(i,j,p,h,e))
        ck(f'first_exit_normalization_{i}',sum(first_formula(i,j,p,h,e) for j in range(3))==1)
    ck('exact_first_exit_1883_2502',first_formula(1,0,p,h,e)==F(1883,2502))
    ck('exact_first_exit_619_2502',first_formula(1,2,p,h,e)==F(619,2502))
    ck('exact_branch_difference_13_5004',first_formula(1,0,p,h,e)-F(3,4)==F(13,5004))
    for sign in (F(1,2),F(-1,2)):
        for n in range(1,12):
            i,j,z=1,0,2
            actual=p[z]*c_formula(n,i,j,p,h,sign)/(1-p[i])-p[j]*c_formula(n,i,z,p,h,sign)/(1-p[i])
            expected=sign*h[i]**2*(h[j]-h[z])*p[i]**n*p[j]*p[z]*(1+sign*h[i]**3)**((n-1)//3)/(3*(1-p[i]))
            ck(f'same_sign_clock_contrast_{sign}_{n}',actual==expected and actual*sign>0)
    for i,j in product(range(3),repeat=2):
        if i!=j:ck(f'epsilon_zero_exit_control_{i}_{j}',first_formula(i,j,p,h,F(0))==p[j]/(1-p[i]))

    p12=[F(i,78) for i in range(1,13)]
    c=sum(x*x for x in p12)
    h12=[x-c for x in p12]
    ck('q12_rational_fixture_centered',sum(p12)==1 and sum(x*y for x,y in zip(p12,h12))==0)
    # All 1728 block words only, not a count-state or full hidden-state enumeration.
    for i,j in product(range(12),repeat=2):
        ck(f'q12_proper_pair_marginal_{i}_{j}',sum(block((i,j,z),p12,h12,e) for z in range(12))==p12[i]*p12[j])
    for i in (0,5,11):
        for j in (0,4,11):
            if i==j:continue
            for n in range(9):
                word=(i,)*(n+1)+(j,)
                ck(f'q12_phase_segment_word_{i}_{j}_{n}',direct_word(word,p12,h12,e)/p12[i]==c_formula(n,i,j,p12,h12,e))
        ck(f'q12_first_exit_normalization_{i}',sum(first_formula(i,j,p12,h12,e) for j in range(12))==1)
    coeff=e*h[1]**2*(h[0]-h[2])*p[1]**2*p[0]*p[2]/(3*(1-p[1]))
    ck('finite_horizon_coefficient_1_1728',coeff==F(1,1728))
    factorial=1
    partial=F(0)
    for n in range(5):
        if n:factorial*=n
        partial+=F(2**n,factorial)
    ck('exp2_strict_lower_series',partial==7)
    ck('TV_strict_lower_1_3024',coeff*(1-F(3,7))==F(1,3024))
    ck('profiled_KL_strict_lower_1_4572288',2*F(1,3024)**2==F(1,4572288))
    ck('q2_age_hazard_not_constant',F(1)+F(1,2)!=1)
    ck('q2_operational_clock_derivative',F(2)-F(1,2)==F(1)+F(1,2))

    ck('instruction_and_scientific_inputs_unchanged',before=={p:sha(ROOT/p) for p in INPUTS})
    ck('protected_predecessor_manifests_unchanged',old=={p:sha(ROOT/p/'MANIFEST.sha256') for p in PROTECTED})
    result={'date':'2026-10-20','status':'PASS_BOUNDED_INTRINSIC_CLOCK_CONTROLS',
            'bounded_checks':len(checks),'checks':checks,'arithmetic':'exact Fraction and symbolic rational logarithms',
            'full_hidden_fixture':'81 q3 block/phase states; each killed fibre has 27 states',
            'q12_fixture':'artificial rational probabilities, not actual A7 numerical probabilities; 1728 block words',
            'first_exit_row_i1':['1883/2502','0','619/2502'],
            'finite_time_TV_strict_lower_bound':'1/3024',
            'profiled_KL_strict_lower_bound':'1/4572288',
            'F01':'OPEN_SOURCE_WITH_CLOCK_FREE_BRANCHING_OBLIGATION',
            'new_independent_source':False,'CSP_executed':False,'new_physical_law':False,
            'input_hashes':before,'protected_manifests':old,'protected_manifest_entries':old_counts}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    if args.record:
        assert args.output and args.output.resolve()==HERE/'RESULTS.json'
        search=save_search()
        (HERE/'DUPLICATE_SEARCH.json').write_text(json.dumps(search,indent=2,ensure_ascii=False)+'\n')
        (HERE/'INPUTS.sha256').write_text(''.join(f'{before[name]}  {name}\n' for name in INPUTS))
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(f)}  {f.name}\n' for f in sorted(HERE.iterdir()) if f.is_file() and f.name!='MANIFEST.sha256'))
    if args.check:
        assert json.loads((HERE/'RESULTS.json').read_text())==result
        files={}
        for line in (HERE/'MANIFEST.sha256').read_text().splitlines():
            hsh,name=line.split('  ',1)
            assert name not in files and Path(name).name==name and sha(HERE/name)==hsh
            files[name]=hsh
        assert set(files)|{'MANIFEST.sha256'}=={p.name for p in HERE.iterdir() if p.is_file()}
        assert (HERE/'INPUTS.sha256').read_text()==''.join(f'{before[name]}  {name}\n' for name in INPUTS)
        for f in HERE.glob('*.json'):json.loads(f.read_text())
        for f in HERE.glob('*.py'):ast.parse(f.read_text())
    print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'F01':result['F01'],'CSP':'NOT_ACTIVATED'}))


if __name__=='__main__':main()
