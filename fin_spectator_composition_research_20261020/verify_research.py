#!/usr/bin/env python3
"""Exact finite controls for the unchanged hypothesis's composition test."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from math import comb, factorial
from pathlib import Path
import re
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
INPUTS=[
    'AGENTS.md','research_tasks/00_EXECUTION_RULES.md',
    'fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md',
    'fin_hyp_history_entropy_20261020/MODEL_CONTRACT.json',
    'fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md',
    'fin_reference_transport_research_20261020/PROOF_OR_OPEN_LEMMA.md',
    'fin_reference_transport_research_20261020/NEXT_OPEN_ATOM.md',
    'fin_reversible_source_transfer_research_20261020/PROOF_OR_OPEN_LEMMA.md',
    'fin_replication_consistency/README.md','fin_transition_hierarchy/REPORT.md',
    'FIN_ST2532_ST2621_Configuration_Generator_and_Dynamic_RG_Report_EN.md',
    'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json',
]
PROTECTED=['fin_hyp_history_entropy_20261020','fin_reference_transport_research_20261020',
           'fin_hyp_history_entropy_review_20261020','fin_reversible_source_transfer_research_20261020']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen=set()
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);p=Path(n)
        assert n not in seen and not p.is_absolute() and '..' not in p.parts
        assert sha(folder/n)==h,(folder,n);seen.add(n)
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==seen|{'MANIFEST.sha256'},folder
    return len(seen)


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0)) for j in range(len(b[0]))] for i in range(len(a))]


def check_markdown(folder):
    for path in folder.glob('*.md'):
        text=path.read_text()
        assert all(line==line.rstrip() for line in text.splitlines()),path
        prose=[];fenced=False
        for line in text.splitlines():
            if line.lstrip().startswith('```'):fenced=not fenced
            elif not fenced and not line.startswith(('    ','\t')):prose.append(line)
        prose=re.sub(r'`[^`\n]*`','', '\n'.join(prose))
        for link in re.findall(r'\[[^\]\n]+\]\(([^)\n]+)\)',prose):
            link=link.strip('<>').split('#')[0]
            if link and not link.startswith(('http:','https:')):
                assert (path.parent/link).exists(),(path,link)


def eye(n):return [[F(i==j) for j in range(n)] for i in range(n)]


def compress(p,pi,labels):
    q=max(labels)+1
    mass=[sum((pi[x] for x in range(len(pi)) if labels[x]==i),F(0)) for i in range(q)]
    rows=[[sum((p[x][y] for y in range(len(pi)) if labels[y]==j),F(0)) for j in range(q)] for x in range(len(pi))]
    out=[[sum((pi[x]*rows[x][j] for x in range(len(pi)) if labels[x]==i),F(0))/mass[i] for j in range(q)] for i in range(q)]
    return mass,rows,out


def solve(a,b):
    n=len(a);aug=[list(row)+[v] for row,v in zip(a,b)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if aug[i][j])
        aug[j],aug[pivot]=aug[pivot],aug[j]
        v=aug[j][j];aug[j]=[x/v for x in aug[j]]
        for i in range(n):
            if i!=j and aug[i][j]:
                c=aug[i][j];aug[i]=[x-c*y for x,y in zip(aug[i],aug[j])]
    return [row[-1] for row in aug]


def product_fixture(biased_coordinate):
    q=3;states=list(product(range(q),repeat=2));a=[1,2,3]
    c=[[F(a[x[biased_coordinate]]*a[y[biased_coordinate]],192)
        if sum(u!=v for u,v in zip(x,y))<=1 else F(0) for y in states] for x in states]
    pi=list(map(sum,c));p=[[c[i][j]/pi[i] for j in range(9)] for i in range(9)]
    return states,c,pi,p


def count_kernel(N,q):
    states=[x for x in product(range(N+1),repeat=q) if sum(x)==N]
    index={x:i for i,x in enumerate(states)};d=1+N*(q-1)
    p=[[F(0) for _ in states] for _ in states]
    for i,n in enumerate(states):
        p[i][i]=F(1,d)
        for a in range(q):
            if n[a]:
                for b in range(q):
                    if a!=b:
                        m=list(n);m[a]-=1;m[b]+=1
                        p[i][index[tuple(m)]]+=F(n[a],d)
    pi=[]
    for n in states:
        denominator=q**N
        for x in n:denominator*=factorial(x)
        pi.append(F(factorial(N),denominator))
    return states,p,pi


def search():
    old=json.loads((ROOT/'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json').read_text())
    paths=[]
    for name in old['searched_files']:
        p=Path(name)
        if (len(p.parts)==1 and re.search(r'(Composition|Replication|Projectiv|Process|Generator|Markov|Kinetic|Entropy|Family|Hierarchy)',p.name,re.I)) or p.parts[0] in {'fin_transition_hierarchy','fin_replication_consistency'}:
            paths.append(name)
    for d in ('fin_hyp_history_entropy_20261020','fin_hyp_history_entropy_review_20261020','fin_reference_transport_research_20261020'):
        paths.extend(str(p.relative_to(ROOT)) for p in (ROOT/d).glob('*.md'))
    paths=sorted(set(paths));hits=[]
    pattern=r'projectiv.{0,65}(kinetic|Markov|source|generator)|spectator|product.{0,35}(equilibrium|stationary|measure)|maximum.entropy.{0,80}(composition|product|autonom)|delet.{0,60}(copy|slot|Markov|intertwin)'
    for i in range(0,len(paths),100):
        r=subprocess.run(['rg','-n','-i','--max-columns','330','--',pattern,*paths[i:i+100]],cwd=ROOT,capture_output=True,text=True,timeout=60)
        assert r.returncode in (0,1);hits.extend(r.stdout.splitlines())
    return {'selected_files':len(paths),'pattern':pattern,'matching_lines':len(hits),'lines':hits,
            'scope':'bounded duplicate search, not every proof read or global novelty',
            'delta':'test the law selected by the authorized variational hypothesis; older packages assume or construct compatible families'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true')
    ap.add_argument('--check',action='store_true');ap.add_argument('--output',type=Path)
    args=ap.parse_args();checks=[]
    def ck(name,value):
        if not value:raise AssertionError(name)
        assert name not in checks;checks.append(name)
    inputs={p:sha(ROOT/p) for p in INPUTS}
    old={p:{'entries':manifest(ROOT/p),'manifest_sha256':sha(ROOT/p/'MANIFEST.sha256')} for p in PROTECTED}
    nu=[F(1,8),F(5,16),F(9,16)];uniform=[F(1,3)]*3
    st,c,pi,p=product_fixture(1);labels=[x[0] for x in st]
    ck('hidden_bias_fixture_product_static_law',pi==[uniform[x]*nu[a] for x,a in st])
    ck('hidden_bias_source_stochastic',all(sum(row)==1 for row in p))
    ck('hidden_bias_source_reversible',all(c[i][j]==c[j][i] for i,j in product(range(9),repeat=2)))
    ck('hidden_bias_optimal_scaling_certificate',all(c[i][j]**2==c[i][i]*c[j][j] for i,j in product(range(9),repeat=2) if c[i][j]))
    mass,rows,c1=compress(p,pi,labels);_,_,c2=compress(mm(p,p),pi,labels)
    square=mm(c1,c1);defect=sum(c2[i][i]-square[i][i] for i in range(3))
    ck('visible_uniform_stationary_marginal',mass==uniform)
    ck('explicit_visible_kernel',c1==[[F(9,16) if i==j else F(7,32) for j in range(3)] for i in range(3)])
    beta=[rows[st.index((0,a))][1] for a in range(3)]
    mean=sum(nu[i]*beta[i] for i in range(3))
    variance=sum(nu[i]*(beta[i]-mean)**2 for i in range(3))
    ck('hidden_activity_exact',[str(x) for x in beta]==['1/8','1/5','1/4'] and mean==F(7,32))
    ck('hidden_activity_variance',variance==F(9,5120))
    ck('memory_Gram_defect_81_2560',defect==18*variance==F(81,2560)>0)
    ck('visible_not_strongly_lumpable',rows[st.index((0,0))]!=rows[st.index((0,2))])
    for i in range(3):
        fibre=[x for x,label in enumerate(labels) if label==i]
        a=[[F(x==y)-p[x][y] for y in fibre] for x in fibre]
        for j in range(3):
            if j!=i:
                rhs=[sum(p[x][y] for y,label in enumerate(labels) if label==j) for x in fibre]
                ck(f'killed_fibre_uniform_jump_order_{i}_{j}',solve(a,rhs)==[F(1,2)]*3)
    ck('exit_tail_guard',1-2*min(beta)==F(3,4))
    st2,c2,pi2,p2=product_fixture(0)
    mass2,rows2,coarse=compress(p2,pi2,[x[0] for x in st2])
    expected=[[F(3,8),F(1,4),F(3,8)],
              [F(1,10),F(3,5),F(3,10)],
              [F(1,12),F(1,6),F(3,4)]]
    ck('visible_bias_fixture_product_static_law',pi2==[nu[x]*uniform[a] for x,a in st2])
    ck('visible_bias_projection_strongly_Markov',all(rows2[i]==expected[x] for i,(x,a) in enumerate(st2)) and coarse==expected)
    ck('visible_bias_stationary_measure',mass2==nu)
    j_joint=coarse[2][0]/(1-coarse[2][2]);j_isolated=nu[0]/(1-nu[2])
    ck('clock_free_autonomy_failure',j_joint==F(1,3) and j_isolated==F(2,7) and j_joint-j_isolated==F(1,21))
    ck('unpostselected_initial_exit_record_bound',nu[2]*(j_joint-j_isolated)==F(3,112))
    # The transported reference is K1+2I, certified independently at psi=a/8.
    a=[1,2,3];transported=[[F(3 if i==j else 1) for j in range(3)] for i in range(3)]
    flow=[[transported[i][j]*F(a[i]*a[j],64) for j in range(3)] for i in range(3)]
    ck('transported_reference_matches_actual_projection',list(map(sum,flow))==nu and [[flow[i][j]/nu[i] for j in range(3)] for i in range(3)]==coarse)
    ck('isolated_source_is_iid_certificate',sum(nu)==1 and [[nu[i]*nu[j]/nu[i] for j in range(3)] for i in range(3)]==[nu]*3)
    # q2 warns that an odd-cycle argument cannot be used on every graph.
    pbin=[F(2,7),F(5,7)];isolated=[list(pbin),list(pbin)]
    induced=[[F(1,2),F(1,2)],[F(1,5),F(4,5)]]
    ck('binary_nonuniform_countercontrol_stationary',all(sum(pbin[i]*induced[i][j] for i in range(2))==pbin[j] for j in range(2)))
    ck('binary_jump_orders_always_alternate',all(induced[i][1-i]/(1-induced[i][i])==isolated[i][1-i]/(1-isolated[i][i])==1 for i in range(2)))
    ck('binary_nonuniform_lazy_equivalence',induced==[[F(i==j)+F(7,10)*(isolated[i][j]-F(i==j)) for j in range(2)] for i in range(2)])
    # Uniform N=3 -> M=2, q=3: direct labelled-source rows, not just count means.
    q=3;N=3;M=2;dN=7;dM=5
    full=list(product(range(q),repeat=N));small=list(product(range(q),repeat=M))
    row_fixture=[]
    for x in full:
        neighbours={x}
        for slot in range(N):
            for val in range(q):
                y=list(x);y[slot]=val;neighbours.add(tuple(y))
        row_fixture.append([F(sum(y[:M]==v for y in neighbours),dN) for v in small])
    pm=[[F(sum(a!=b for a,b in zip(x,y))<=1,dM) for y in small] for x in small]
    lazy=[[F(i==j)+F(dM,dN)*(pm[i][j]-F(i==j)) for j in range(len(small))] for i in range(len(small))]
    ck('uniform_labelled_deletion_actual_kernel',all(row_fixture[i]==lazy[small.index(x[:M])] for i,x in enumerate(full)))
    ck('uniform_degree_counts',all(sum(row)==1 for row in lazy) and len(full)==27 and len(small)==9)
    ns,pn,pin=count_kernel(N,q);ms,pmc,pim=count_kernel(M,q)
    delete=[]
    for n in ns:
        row=[]
        for m in ms:
            numerator=1
            for na,ma in zip(n,m):
                numerator*=comb(na,ma) if ma<=na else 0
            row.append(F(numerator,comb(N,M)))
        delete.append(row)
    small_lazy=[[F(i==j)+F(dM,dN)*(pmc[i][j]-F(i==j)) for j in range(len(ms))] for i in range(len(ms))]
    ck('hypergeometric_delete_rows',all(sum(row)==1 for row in delete) and len(ns)==10 and len(ms)==6)
    ck('uniform_count_intertwining',mm(pn,delete)==mm(delete,small_lazy))
    ck('uniform_equilibrium_deletion_consistent',all(sum(pin[i]*delete[i][j] for i in range(len(ns)))==pim[j] for j in range(len(ms))))
    ck('source_count_kernels_stationary',all(sum(pin[i]*pn[i][j] for i in range(len(ns)))==pin[j] for j in range(len(ns)))
       and all(sum(pim[i]*pmc[i][j] for i in range(len(ms)))==pim[j] for j in range(len(ms))))
    ck('all_inputs_preserved',inputs=={p:sha(ROOT/p) for p in INPUTS})
    ck('old_packages_preserved',old=={p:{'entries':manifest(ROOT/p),'manifest_sha256':sha(ROOT/p/'MANIFEST.sha256')} for p in PROTECTED})
    result={'date':'2026-10-20','model_id':'FIN-HYP-HIST-ME-001_UNCHANGED',
            'status':'PASS_BOUNDED_SPECTATOR_COMPOSITION_CONTROLS','bounded_checks':len(checks),'checks':checks,
            'arithmetic':'exact Fraction; nine-state and 27-state labelled controls, 10-to-6 count deletion',
            'memory_defect':'81/2560','clock_free_branch_difference':'1/21',
            'product_measures':'diagnostic inputs, not actual positive-g A7 equilibrium',
            'source_inputs':inputs,'protected_packages':old,
            'new_postulate_or_clock':False,'dynamic_autonomy_adopted_as_physical_axiom':False,
            'strict_F01':'OPEN_SOURCE','primary_FIN_changed':False,'independent_peer_review':False,
            'goal_objective':'badaj teorię','goal_finished':False}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    if args.record:
        assert args.output and args.output.resolve()==HERE/'RESULTS.json'
        (HERE/'DUPLICATE_SEARCH.json').write_text(json.dumps(search(),indent=2)+'\n')
        (HERE/'INPUTS.sha256').write_text(''.join(f'{h}  {n}\n' for n,h in inputs.items()))
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.sha256'))
    if args.check:
        manifest(HERE)
        assert json.loads((HERE/'RESULTS.json').read_text())==result
        assert (HERE/'INPUTS.sha256').read_text()==''.join(f'{h}  {n}\n' for n,h in inputs.items())
        for p in HERE.glob('*.json'):json.loads(p.read_text())
        for p in HERE.glob('*.py'):ast.parse(p.read_text())
        check_markdown(HERE)
    print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'strict_F01':'OPEN_SOURCE','new_postulate':False}))


if __name__=='__main__':main()
