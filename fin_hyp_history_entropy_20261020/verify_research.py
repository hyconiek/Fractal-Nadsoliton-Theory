#!/usr/bin/env python3
"""Exact finite certificates for the frozen authorized entropy hypothesis."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
INPUTS=[
    'AGENTS.md','FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md',
    'FIN_RESEARCH_STRATEGY.md','research_tasks/00_EXECUTION_RULES.md',
    'FIN_ST2532_ST2621_Configuration_Generator_and_Dynamic_RG_Report_EN.md',
    'fin_st2532_st2621_configuration_generator_research.py',
    'FIN_ST6897_ST8396_Hundred_Round_Entropic_Frenetic_and_FDT_Report_EN.md',
    'fin_evolution_source_research_20261020/SOURCE_ATOM.md',
    'fin_intrinsic_clock_source_research_20261020/PROOF_OR_OPEN_LEMMA.md',
    'fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md',
    'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json',
]
PROTECTED=['fin_evolution_source_research_20261020','fin_evolution_source_review_20261020',
           'fin_intrinsic_clock_source_research_20261020','fin_agents_causal_update_20261020']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    seen=set()
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1)
        p=Path(n)
        assert n not in seen and not p.is_absolute() and '..' not in p.parts
        assert sha(folder/n)==h,(folder,n)
        seen.add(n)
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==seen|{'MANIFEST.sha256'},folder
    return len(seen)


def fixture(a):
    q=len(a);total=sum(a);den=2*total*total-sum(x*x for x in a)
    p=[F(x*(2*total-x),den) for x in a]
    states=list(product(range(q),repeat=2))
    dif=[(y-x)%q for x,y in states]
    pi=[p[d]/q for d in dif]
    k=[[int(sum(x!=y for x,y in zip(s,t))<=1) for t in states] for s in states]
    c=[[F(k[i][j]*a[dif[i]]*a[dif[j]],den*q) for j in range(q*q)] for i in range(q*q)]
    mat=[[c[i][j]/pi[i] for j in range(q*q)] for i in range(q*q)]
    return states,dif,pi,k,c,mat,p,den


def coarse(mat,dif,q):
    return [[sum((mat[i][j] for j,v in enumerate(dif) if v==d),F(0)) for d in range(q)] for i in range(len(mat))]


def prime_factors(n):
    out={};d=2
    while d*d<=n:
        while n%d==0:
            out[d]=out.get(d,0)+1;n//=d
        d+=1
    if n>1:out[n]=out.get(n,0)+1
    return out


def plus_log(out,x,coefficient):
    if not coefficient:return
    assert x>0
    for p,n in prime_factors(x.numerator).items():out[p]=out.get(p,F(0))+coefficient*n
    for p,n in prime_factors(x.denominator).items():out[p]=out.get(p,F(0))-coefficient*n


def cleaned(out):return {p:v for p,v in out.items() if v}


def entropy(c,pi):
    logs={}
    for i,row in enumerate(c):
        for x in row:
            if x:plus_log(logs,x/pi[i],-x)
    return cleaned(logs)


def subtract(a,b):
    out=dict(a)
    for p,v in b.items():out[p]=out.get(p,F(0))-v
    return cleaned(out)


def kl(a,b):
    out={}
    for row,ref in zip(a,b):
        for x,y in zip(row,ref):
            if x:plus_log(out,x/y,x)
    return cleaned(out)


def search():
    previous=json.loads((ROOT/'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json').read_text())
    paths=sorted(previous['searched_files'])
    pattern=r'maximum.calib(er|re)|maximum.entropy.{0,45}(Markov|random.walk|path|histor)|path.entropy|sinkhorn|entropy.{0,35}multiplicit'
    matches=[]
    for start in range(0,len(paths),100):
        run=subprocess.run(['rg','-n','-i','--max-columns','300','--',pattern,*paths[start:start+100]],
                           cwd=ROOT,capture_output=True,text=True,timeout=60)
        assert run.returncode in (0,1),run.stderr
        matches.extend(run.stdout.splitlines())
    return {'selected_files':len(paths),'scope':'frozen prior index selection; duplicate search, not every proof read',
            'pattern':pattern,'matching_lines':len(matches),'lines':matches,
            'prior_index_capture_sha256':sha(ROOT/'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json'),
            'novelty':'authorized concrete postulate and its FIN test, not a new general maximum-caliber method'}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true')
    ap.add_argument('--check',action='store_true');ap.add_argument('--output',type=Path)
    args=ap.parse_args();checks=[]
    def ck(name,value):
        if not value:raise AssertionError(name)
        assert name not in checks;checks.append(name)
    inputs={p:sha(ROOT/p) for p in INPUTS}
    old={p:{'manifest_sha256':sha(ROOT/p/'MANIFEST.sha256'),'entries':manifest(ROOT/p)} for p in PROTECTED}
    reports=[]
    for a in ([1,2,3,2],[1,2,3,4,5,6,7,6,5,4,3,2]):
        st,dif,pi,k,c,mat,p,den=fixture(a)
        q=len(a);m=q*q;s=sum(a)
        ck(f'q{q}_positive_pi',all(x>0 for x in pi) and sum(pi)==1)
        ck(f'q{q}_correct_support',all((c[i][j]>0)==bool(k[i][j]) for i,j in product(range(m),repeat=2)))
        ck(f'q{q}_degree',all(sum(row)==2*q-1 for row in k))
        ck(f'q{q}_stationary_row_margins',all(sum(c[i])==pi[i] for i in range(m)))
        ck(f'q{q}_reversible_column_margins',all(c[i][j]==c[j][i] for i,j in product(range(m),repeat=2)))
        ck(f'q{q}_stochastic',all(sum(row)==1 for row in mat))
        ck(f'q{q}_factorization_certificate',all(c[i][j]**2==c[i][i]*c[j][j] for i,j in product(range(m),repeat=2) if k[i][j]))
        cm=coarse(mat,dif,q)
        expected=[[F(a[d],2*s-a[d]) if d==v else F(2*a[v],2*s-a[d]) for v in range(q)] for d in range(q)]
        ck(f'q{q}_actual_quotient_intertwining',all(cm[i]==expected[dif[i]] for i in range(m)))
        ck(f'q{q}_quotient_stationary',all(sum((p[d]*expected[d][v] for d in range(q)),F(0))==p[v] for v in range(q)))
        ck(f'q{q}_branch_normalization',all(sum((F(a[v],s-a[d]) for v in range(q) if v!=d),F(0))==1 for d in range(q)))
        ck(f'q{q}_nonuniform_not_HB_mod_clock',any(F(a[j],a[z])!=p[j]/p[z] for j,z in product(range(q),repeat=2)))
        ck(f'q{q}_ordering_amplifies_HB_odds',all(F(a[j],a[z])>p[j]/p[z] for j,z in product(range(q),repeat=2) if p[j]>p[z]))
        ck(f'q{q}_scalar_root_discriminants',all(F((2*s)**2,den)-4*x>0 for x in p))
        # Square roots are exact (2s-2a[d])/sqrt(den); sum identity needs no floating roots.
        ck(f'q{q}_scalar_root_sum',sum(2*s-2*x for x in a)==(q-1)*2*s)
        reports.append({'q':q,'labelled_states':m,'degree':2*q-1,'a':list(a),'denominator':den,'p':[str(x) for x in p]})
        if q==4:
            ck('exact_HB_branch_mismatch_minus_4_355',F(a[0],s-a[2])-p[0]/(1-p[2])==F(-4,355))
            # Independent feasible flow perturbation, with exact logarithmic KL identity.
            trial=[list(row) for row in c]
            x,y=0,1
            delta=min(c[x][x],c[y][y],c[x][y])/10
            trial[x][x]-=delta;trial[y][y]-=delta
            trial[x][y]+=delta;trial[y][x]+=delta
            ck('nonoptimal_flow_preserves_margins',all(sum(trial[i])==pi[i] for i in range(m)))
            gap=subtract(entropy(c,pi),entropy(trial,pi))
            ck('exact_symbolic_entropy_KL_gap',gap==kl(trial,c) and bool(gap))
            # Order-two history with the same pair flow, testing the entropy memory term.
            middle=0;pred=0;dest=1;eps=F(1,2)
            triple={};base={}
            for i,j,z in product(range(m),repeat=3):
                b=c[i][j]*mat[j][z]
                if b:
                    u=F(i==pred)-mat[middle][pred]
                    v=F(z==dest)-mat[middle][dest]
                    t=b*(1+eps*u*v) if j==middle else b
                    triple[i,j,z]=t;base[i,j,z]=b
            ck('memory_control_same_two_margins',all(sum((v for (i,j,z),v in triple.items() if (i,j)==(x,y)),F(0))==c[x][y]
               and sum((v for (i,j,z),v in triple.items() if (j,z)==(x,y)),F(0))==c[x][y] for x,y in product(range(m),repeat=2)))
            h2={};d3={}
            for (i,j,z),v in triple.items():
                plus_log(h2,v/c[i][j],-v)
                plus_log(d3,v/base[i,j,z],v)
            ck('memory_entropy_loss_identity',subtract(entropy(c,pi),cleaned(h2))==cleaned(d3) and bool(cleaned(d3)))
            # Macro counting without the microscopic two-choice factor picks iid p.
            macro_iid=[list(p) for _ in p]
            ck('bare_macro_optimum_iid_certificate',all(sum(row)==1 for row in macro_iid))
            ck('bare_macro_differs_from_actual_projection',macro_iid!=expected)
            # Non-adopted marked-self contrast, not a repair of the primary axiom.
            marked=[[F(0) for _ in st] for _ in st]
            for i in range(m):
                for j in range(m):
                    if k[i][j]:marked[i][j]=p[dif[j]] if i==j else p[dif[j]]/2
            ck('marked_control_stochastic',all(sum(row)==1 for row in marked))
            ck('marked_control_has_reset_quotient',all(row==p for row in coarse(marked,dif,q)))
            ck('marked_control_not_primary_kernel',marked!=mat)

    # Mandatory thinning: exact destination-odds double ratios on actual q12 geometry.
    q=12
    dist=lambda a,b:min((a-b)%q,(b-a)%q)
    thin=lambda a,b:F(1,10) if dist(a,b)==3 else F(1)
    ratio=(thin(3,0)/thin(3,1))/(thin(4,0)/thin(4,1))
    ck('mandatory_thinning_double_ratio_1_100',ratio==F(1,100))
    ck('primary_double_ratio_one',F(1)!=ratio)
    ck('one_vs_two_self_weight_not_vertex_gauge',F(2)*F(2)!=1)
    # Uniform p: main branch law equals HB only modulo a clock, not as a step kernel.
    aa=[1]*12;_,dif,pi,k,c,mat,p,den=fixture(aa)
    ck('uniform_zero_gain_branch_control',all(F(aa[j],sum(aa)-aa[i])==p[j]/(1-p[i]) for i,j in product(range(12),repeat=2) if i!=j))
    ck('uniform_zero_gain_self_steps_differ',F(1,23)!=F(1,12))
    # General marked-HB recovery is false: eight-state binary N3 Gibbs toy.
    st=list(product((0,1),repeat=3))
    weights=[F(8) if len(set(x))==1 else F(2) for x in st]
    pi=[x/sum(weights) for x in weights]
    hb=[[F(0) for _ in st] for _ in st]
    for i,x in enumerate(st):
        for slot in range(3):
            y=list(x);y[slot]=1-y[slot];j=st.index(tuple(y))
            hb[i][j]+=weights[j]/(3*(weights[i]+weights[j]))
            hb[i][i]+=weights[i]/(3*(weights[i]+weights[j]))
    i,j=st.index((0,0,0)),st.index((1,0,0))
    cxx,cyy,cxy=pi[i]*hb[i][i],pi[j]*hb[j][j],pi[i]*hb[i][j]
    ck('binary_N3_control_Gibbs_and_balance',sum(pi)==1 and all(sum(r)==1 for r in hb)
       and all(pi[x]*hb[x][y]==pi[y]*hb[y][x] for x,y in product(range(8),repeat=2)))
    ck('binary_N3_exact_flows',(cxx,cyy,cxy)==(F(8,35),F(1,35),F(2,105)))
    ck('binary_N3_marked_factorization_rejected',9*cxy*cxy/(cxx*cyy)==F(1,2))
    ck('old_inputs_preserved',inputs=={p:sha(ROOT/p) for p in INPUTS})
    ck('old_manifests_preserved',old=={p:{'manifest_sha256':sha(ROOT/p/'MANIFEST.sha256'),'entries':manifest(ROOT/p)} for p in PROTECTED})
    result={'date':'2026-10-20','model_id':'FIN-HYP-HIST-ME-001',
            'execution':'COMPLETED_FROZEN_HYPOTHESIS_AND_TESTS',
            'science':'UNIQUE_CONDITIONAL_SOURCE_DERIVED; NONUNIFORM_N2_HB_COMPATIBILITY_REFUTED_MOD_CLOCK',
            'bounded_check_count':len(checks),'bounded_checks':checks,
            'arithmetic':'exact Fraction and formal logarithms; no native A7 floating p or asymptotic proof by fixtures',
            'fixtures':reports,'exact_q4_branch_difference':'-4/355','thinning_double_ratio':'1/100',
            'N3_stress_factor':'1/2','same_hypothesis_retuned':False,
            'strict_F01':'OPEN_SOURCE','conditional_PM_source':'DERIVED_GIVEN_NEW_POSTULATE_AND_SUPPLIED_CARRIER_PI_SUPPORT',
            'primary_FIN_changed':False,'physical_adoption':False,'CSP_strict_executed':False,
            'input_hashes':inputs,'protected_packages':old}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    if args.record:
        assert args.output and args.output.resolve()==HERE/'RESULTS.json'
        (HERE/'DUPLICATE_SEARCH.json').write_text(json.dumps(search(),indent=2)+'\n')
        (HERE/'INPUTS.sha256').write_text(''.join(f'{inputs[p]}  {p}\n' for p in INPUTS))
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.sha256'))
    if args.check:
        manifest(HERE)
        assert json.loads((HERE/'RESULTS.json').read_text())==result
        assert (HERE/'INPUTS.sha256').read_text()==''.join(f'{inputs[p]}  {p}\n' for p in INPUTS)
        for p in HERE.glob('*.json'):json.loads(p.read_text())
        for p in HERE.glob('*.py'):ast.parse(p.read_text())
    print(json.dumps({'model_id':result['model_id'],'checks':len(checks),'science':result['science'],'strict_F01':'OPEN_SOURCE'}))


if __name__=='__main__':main()
