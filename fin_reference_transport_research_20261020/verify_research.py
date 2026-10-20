#!/usr/bin/env python3
"""Exact controls for reference transport; no imported producer helpers."""
import argparse
import ast
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path
import subprocess

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
INPUTS=[
    'AGENTS.md','research_tasks/00_EXECUTION_RULES.md',
    'fin_hyp_history_entropy_20261020/AUTHORIZATION_AND_MODEL.md',
    'fin_hyp_history_entropy_20261020/MODEL_CONTRACT.json',
    'fin_hyp_history_entropy_20261020/PROOF_OR_OPEN_LEMMA.md',
    'fin_hyp_history_entropy_review_20261020/PROOF_AUDIT_AND_DUAL.md',
    'fin_hyp_history_entropy_review_20261020/NEXT_DECISION.md',
    'fin_reversible_source_transfer_research_20261020/PROOF_OR_OPEN_LEMMA.md',
    'fin_generative_operational_review_20261014/EXACT_N2_PROCESS_AND_SHARP_SEPARATION.md',
    'fin_rank7_mathphysics_next/proofs/MP7-031_finite_N_variational_model.md',
    'fin_kinetic_invariants_campaign_20261014/KI-002/REVERSIBILITY_AND_MOBILITY.md',
    'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json',
]
PROTECTED=['fin_hyp_history_entropy_20261020','fin_hyp_history_entropy_review_20261020',
           'fin_intrinsic_clock_source_research_20261020','fin_reversible_source_transfer_research_20261020']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def manifest(folder):
    names=set()
    for line in (folder/'MANIFEST.sha256').read_text().splitlines():
        h,n=line.split('  ',1);p=Path(n)
        assert n not in names and not p.is_absolute() and '..' not in p.parts
        assert sha(folder/n)==h,(folder,n)
        names.add(n)
    actual={p.relative_to(folder).as_posix() for p in folder.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
    assert actual==names|{'MANIFEST.sha256'},folder
    return len(names)


def mm(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))),F(0)) for j in range(len(b[0]))] for i in range(len(a))]


def transpose(a):return [list(x) for x in zip(*a)]


def equiv(k,r):
    n=len(k)
    if any(k[i][i]<=0 or r[i][i]<=0 for i in range(n)):return False
    for i,j in product(range(n),repeat=2):
        if (k[i][j]>0)!=(r[i][j]>0):return False
        if k[i][j]>0 and k[i][j]**2*r[i][i]*r[j][j]!=r[i][j]**2*k[i][i]*k[j][j]:return False
    return True


def aggregate(p,groups):
    return [[sum((p[x][y] for y in block),F(0)) for block in groups] for x in range(len(p))]


@dataclass(frozen=True,eq=False)
class Tau:
    """a+b*t with t^2=1-t and the positive root 3/5<t<5/8."""
    a:F=F(0)
    b:F=F(0)
    def __post_init__(self):
        object.__setattr__(self,'a',F(self.a));object.__setattr__(self,'b',F(self.b))
    @staticmethod
    def make(x):return x if isinstance(x,Tau) else Tau(x)
    def __add__(self,other):
        y=Tau.make(other);return Tau(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self):return Tau(-self.a,-self.b)
    def __sub__(self,other):return self+-Tau.make(other)
    def __rsub__(self,other):return Tau.make(other)+-self
    def __mul__(self,other):
        y=Tau.make(other)
        return Tau(self.a*y.a+self.b*y.b,self.a*y.b+self.b*y.a-self.b*y.b)
    __rmul__=__mul__
    def inverse(self):
        d=self.a*self.a-self.a*self.b-self.b*self.b
        if not d:raise ZeroDivisionError
        return Tau((self.a-self.b)/d,-self.b/d)
    def __truediv__(self,other):return self*Tau.make(other).inverse()
    def __rtruediv__(self,other):return Tau.make(other)*self.inverse()
    def __pow__(self,n):
        if n<0:return self.inverse()**(-n)
        out=Tau(1)
        for _ in range(n):out=out*self
        return out
    def __eq__(self,other):
        y=Tau.make(other);return self.a==y.a and self.b==y.b
    def box(self):
        lo,hi=F(3,5),F(5,8)
        if self.b<0:lo,hi=hi,lo
        return self.a+self.b*lo,self.a+self.b*hi


def duplicate_search():
    source=json.loads((ROOT/'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json').read_text())
    names=set(source['searched_files'])
    for d in ('fin_hyp_history_entropy_20261020','fin_hyp_history_entropy_review_20261020'):
        names.update(str(p.relative_to(ROOT)) for p in (ROOT/d).glob('*.md'))
    names=sorted(names)
    pattern=r'equitable.{0,70}(entropy|scal|reference)|entropy.{0,70}equitable|reference.{0,50}gauge|gauge.{0,50}reference|nested.{0,50}(multiplicity|counting)|quotient.{0,50}multiplicity'
    matches=[]
    for start in range(0,len(names),100):
        r=subprocess.run(['rg','-n','-i','--max-columns','320','--',pattern,*names[start:start+100]],cwd=ROOT,text=True,capture_output=True,timeout=60)
        assert r.returncode in (0,1),r.stderr;matches.extend(r.stdout.splitlines())
    return {'selected_files':len(names),'pattern':pattern,'matching_lines':len(matches),'lines':matches,
            'scope':'bounded indexed duplicate search; not all proofs read or a global novelty certificate',
            'prior_scope_sha256':sha(ROOT/'fin_reversible_source_transfer_research_20261020/SEARCH_RESULTS.json')}


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--record',action='store_true')
    ap.add_argument('--check',action='store_true');ap.add_argument('--output',type=Path)
    args=ap.parse_args();checks=[]
    def ck(name,value):
        if not value:raise AssertionError(name)
        assert name not in checks;checks.append(name)
    before={name:sha(ROOT/name) for name in INPUTS}
    protected={name:{'entries':manifest(ROOT/name),'sha256':sha(ROOT/name/'MANIFEST.sha256')} for name in PROTECTED}
    k=[[F(1),F(1),F(0)],[F(1),F(1),F(1)],[F(0),F(1),F(1)]]
    b=[2,3,5];r=[[b[i]*k[i][j]*b[j] for j in range(3)] for i in range(3)]
    ck('gauge_equivalence_positive',equiv(k,r))
    wrong=[list(row) for row in k];wrong[0][0]=2
    ck('self_multiplicity_not_gauge',not equiv(k,wrong))
    missing=[list(row) for row in k];missing[0][1]=missing[1][0]=0
    ck('support_change_rejected',not equiv(k,missing))
    # Loop-only fibres: raw gauge changes rows, canonical normalization restores them.
    h=[[F(1),F(0),F(1),F(2)],[F(0),F(1),F(2),F(1)],
       [F(1),F(2),F(1),F(0)],[F(2),F(1),F(0),F(1)]]
    bb=[1,2,3,4];raw=[[h[i][j]*bb[i]*bb[j] for j in range(4)] for i in range(4)]
    a=[1,1,2,2];flow=[[h[i][j]*a[i]*a[j]/34 for j in range(4)] for i in range(4)]
    pi=list(map(sum,flow));p=[[flow[i][j]/pi[i] for j in range(4)] for i in range(4)]
    raw_flow=[[raw[i][j]*F(a[i],bb[i])*F(a[j],bb[j])/34 for j in range(4)] for i in range(4)]
    ck('raw_reference_gauge_preserves_source_flow',flow==raw_flow and equiv(h,raw))
    groups=[[0,1],[2,3]]
    cp=aggregate(p,groups)
    ck('loop_only_strong_quotient',cp[0]==cp[1] and cp[2]==cp[3])
    raw_counts=aggregate(raw,groups)
    ck('raw_equitable_not_necessary',raw_counts[0]!=raw_counts[1])
    canonical=[[raw[i][j]/(bb[i]*bb[j]) for j in range(4)] for i in range(4)]
    ck('canonical_restores_block_condition',canonical==h and aggregate(h,groups)[0]==aggregate(h,groups)[1])
    # Canonical balance at the base law is H/16, independent of the raw gauge.
    balanced_raw=[[raw[i][j]*F(1,4*bb[i])*F(1,4*bb[j]) for j in range(4)] for i in range(4)]
    ck('balanced_reference_is_gauge_canonical',balanced_raw==[[x/16 for x in row] for row in h])
    broken=[[F(1),F(0),F(1),F(1)],[F(0),F(1),F(1),F(0)],
            [F(1),F(1),F(1),F(0)],[F(1),F(0),F(0),F(1)]]
    ck('loop_only_nonequitable_reference_detected',aggregate(broken,groups)[0]!=aggregate(broken,groups)[1])
    # Fifty-state nested fixture. All embeddings and sources are rational.
    sizes=[9,16,25];cuts=[0,9,25,50]
    block=[0]*9+[1]*16+[2]*25
    coef=[[F(25,9),F(1),F(1)],[F(1),F(2),F(1)],[F(1),F(1),F(3)]]
    fine=[[coef[block[x]][block[y]] for y in range(50)] for x in range(50)]
    v=[[F(1,[3,4,5][i]) if block[x]==i else F(0) for i in range(3)] for x in range(50)]
    u=[[F(3,5),F(0)],[F(4,5),F(0)],[F(0),F(1)]]
    direct=[[F(1,5) if ((x<25)==(j==0)) else F(0) for j in range(2)] for x in range(50)]
    identity3=[[F(i==j) for j in range(3)] for i in range(3)]
    identity2=[[F(i==j) for j in range(2)] for i in range(2)]
    k1=mm(mm(transpose(v),fine),v)
    k2=mm(mm(transpose(u),k1),u)
    ck('first_embedding_isometry',mm(transpose(v),v)==identity3)
    ck('second_embedding_isometry',mm(transpose(u),u)==identity2)
    ck('nested_embedding_equals_direct',mm(v,u)==direct)
    ck('first_exact_transported_reference',k1==[[25,12,15],[12,32,20],[15,20,75]])
    ck('first_block_invariance',mm(fine,v)==mm(v,k1))
    ck('second_block_invariance',mm(k1,u)==mm(u,k2))
    ck('direct_staged_references_equal',k2==mm(mm(transpose(direct),fine),direct)==[[41,25],[25,75]])
    rows=list(map(sum,fine));pi50=[x/4150 for x in rows]
    pfine=[[fine[x][y]/rows[x] for y in range(50)] for x in range(50)]
    coarse_pi=[sum(pi50[:25]),sum(pi50[25:])]
    ck('nested_stationary_masses',coarse_pi==[F(33,83),F(50,83)] and sum(pi50)==1)
    ck('nested_fine_source_reversible',all(pi50[x]*pfine[x][y]==pi50[y]*pfine[y][x] for x,y in product(range(50),repeat=2)))
    first=aggregate(pfine,[list(range(cuts[i],cuts[i+1])) for i in range(3)])
    middle=[first[cuts[i]] for i in range(3)]
    second=aggregate(middle,[[0,1],[2]])
    final=aggregate(pfine,[list(range(25)),list(range(25,50))])
    expected=[[F(41,66),F(25,66)],[F(1,4),F(3,4)]]
    ck('nested_actual_kernel_intertwining',all(final[x]==expected[int(x>=25)] for x in range(50))
       and second[0]==second[1]==expected[0] and second[2]==expected[1])
    correct_invariant=k2[0][1]**2/(k2[0][0]*k2[1][1])
    wrong_invariant=F(1225,2)/(F(81,2)*75)
    ck('correct_nested_edge_invariant',correct_invariant==F(25,123))
    ck('flattened_mass_changes_gauge_class',wrong_invariant==F(49,243) and wrong_invariant!=correct_invariant)
    ck('intermediate_mass_not_uniform',sum(pi50[:9])/sum(pi50[9:25])==F(9,16))
    # Actual q12 count geometry and a separate exact symmetric probability fixture.
    q=12;states=list(product(range(q),repeat=2));where={x:i for i,x in enumerate(states)}
    diff=[(y-x)%q for x,y in states];radial=[min(d,q-d) for d in diff]
    ell=[1,2,2,2,2,2,1];av=[1,2,3,4,5,6,7,6,5,4,3,2]
    den=2*sum(av)**2-sum(x*x for x in av);mass=[F(x*(96-x),den) for x in av]
    p_rad=[ell[i]*mass[i] for i in range(7)]
    phase_rows=[];rad_rows=[]
    for x,y in states:
        d=(y-x)%q;r=min(d,q-d)
        dest={(x,y)}|{(z,y) for z in range(q)}|{(x,z) for z in range(q)}
        counts_d=[sum((b-a)%q==v for a,b in dest) for v in range(q)]
        counts_r=[sum(min((b-a)%q,q-(b-a)%q)==s for a,b in dest) for s in range(7)]
        phase_rows.append(counts_d);rad_rows.append(counts_r)
    ck('q12_difference_counts_all_states',all(row==[1 if v==diff[i] else 2 for v in range(q)] for i,row in enumerate(phase_rows)))
    ck('q12_radial_counts_all_states',all(row==[2*ell[s]-(s==radial[i]) for s in range(7)] for i,row in enumerate(rad_rows)))
    ck('q12_nested_fibre_sizes',[radial.count(i) for i in range(7)]==[12*x for x in ell])
    counts=[[F(2*ell[j]-(i==j)) for j in range(7)] for i in range(7)]
    ck('q12_reference_symmetry_in_weighted_frame',all(ell[i]*counts[i][j]==ell[j]*counts[j][i] for i,j in product(range(7),repeat=2)))
    ck('q12_self_counts_retain_invisible_real_moves',[counts[i][i] for i in range(7)]==[1,3,3,3,3,3,1])
    ck('q12_radial_reference_not_bare_counting',counts[0][1]*counts[1][0]/(counts[0][0]*counts[1][1])==F(8,3))
    p2=[[counts[i][j]*av[j]/(96-av[i]) for j in range(7)] for i in range(7)]
    ck('q12_fixture_exact_denominator',den==4378)
    ck('q12_radial_stochastic',all(sum(row)==1 for row in p2))
    ck('q12_radial_stationary',all(sum(p_rad[i]*p2[i][j] for i in range(7))==p_rad[j] for j in range(7)))
    for r in (0,2,6):
        x,y=0,r;dest={(x,y)}|{(z,y) for z in range(q)}|{(x,z) for z in range(q)}
        actual=[sum((F(av[(b-a)%q],96-av[r]) for a,b in dest if min((b-a)%q,q-(b-a)%q)==j),F(0)) for j in range(7)]
        ck(f'q12_actual_radial_quotient_{r}',actual==p2[r])
    jhyp=p2[2][0]/(1-p2[2][2]);jhb=p_rad[0]/(1-p_rad[2])
    ck('q12_fixture_radial_branch_1_42',jhyp==F(1,42))
    ck('q12_fixture_radial_HB_branch_19_764',jhb==F(19,764))
    ck('q12_fixture_orientation_free_mismatch',jhyp-jhb==F(-17,16044))
    # A second actual route: labelled states -> unordered counts -> radial classes.
    count_key=lambda state:tuple(state.count(j) for j in range(q))
    types=sorted({count_key(x) for x in states});type_index={x:i for i,x in enumerate(types)}
    orbit=[[i for i,x in enumerate(states) if count_key(x)==n] for n in types]
    mult=list(map(len,orbit));type_rad=[radial[g[0]] for g in orbit]
    count_rows=[];count_dynamics=[]
    for n,g in zip(types,orbit):
        x,y=states[g[0]];d=(y-x)%q
        dest={(x,y)}|{(z,y) for z in range(q)}|{(x,z) for z in range(q)}
        by_type=[sum(count_key(t)==v for t in dest) for v in types]
        dynamics=[sum((F(av[(b-a)%q],96-av[d]) for a,b in dest if count_key((a,b))==v),F(0)) for v in types]
        count_rows.append(by_type);count_dynamics.append(dynamics)
    ck('q12_count_orbits_and_total_mass',len(types)==78 and sum(mult)==144 and set(mult)=={1,2})
    valid=True
    for i,n in enumerate(types):
        for j,np in enumerate(types):
            if i==j:
                valid=valid and count_rows[i][j]==1
            elif count_rows[i][j]:
                old=[a for a in range(q) if np[a]-n[a]==-1]
                new=[a for a in range(q) if np[a]-n[a]==1]
                valid=valid and len(old)==len(new)==1
                if old and new:
                    valid=valid and count_rows[i][j]==n[old[0]]
                    valid=valid and count_rows[i][j]*count_rows[j][i]==n[old[0]]*(n[new[0]]+1)
                valid=valid and mult[i]*count_rows[i][j]==mult[j]*count_rows[j][i]
    ck('count_reference_sqrt_occupation_identity',valid)
    ck('count_to_radial_reference_matches_difference_route',all(
        [sum(count_rows[i][j] for j,rj in enumerate(type_rad) if rj==r) for r in range(7)]
        ==[2*ell[r]-(r==type_rad[i]) for r in range(7)] for i in range(len(types))))
    ck('count_to_radial_source_matches_difference_route',all(
        [sum((count_dynamics[i][j] for j,rj in enumerate(type_rad) if rj==r),F(0)) for r in range(7)]
        ==p2[type_rad[i]] for i in range(len(types))))
    ck('count_mass_transport_is_direct_radial_mass',[sum(mult[i] for i,r in enumerate(type_rad) if r==j) for j in range(7)]==[12*x for x in ell])
    ck('count_auxiliary_base_balance_and_degree',all(sum(row)==23 for row in count_rows)
       and all(F(mult[i],144)*F(count_rows[i][j],23)==F(mult[j],144)*F(count_rows[j][i],23)
               for i,j in product(range(len(types)),repeat=2)))
    # A genuine unweighted non-Markov projection, with exact Q(sqrt(5)) arithmetic.
    t=Tau(0,1);zero=Tau(0)
    poly=lambda x:x*x+x-1
    ck('positive_tau_root_enclosure',poly(F(3,5))<0<poly(F(5,8)))
    toy=[[t,1-t,zero],[1-t,2*t-1,1-t],[zero,1-t,t]]
    ck('golden_path_stochastic_positive',all(sum(row)==1 for row in toy)
       and all(x==0 or x.box()[0]>0 for row in toy for x in row))
    ck('golden_path_source_scaling_certificate',toy[0][1]**2==toy[0][0]*toy[1][1] and toy[1][2]**2==toy[1][1]*toy[2][2])
    def compressed(matrix):
        return [[sum((matrix[x][y] for x in g for y in h),Tau(0))/len(g) for h in ([0,1],[2])] for g in ([0,1],[2])]
    c1=compressed(toy);c2=compressed(mm(toy,toy));square=mm(c1,c1)
    defect=sum((c2[i][i]-square[i][i] for i in range(2)),Tau(0))
    expected_defect=F(3,4)*(2-3*t)
    ck('golden_path_compressed_memory_exact',defect==expected_defect and defect.box()[0]>0)
    ck('golden_balanced_reference_fails_block_invariance',toy[0][2]/3!=toy[1][2]/3)
    # Internal edges defeat a universal necessity claim for canonical invariance.
    inside=[[F(3,5),F(1,10),F(1,5),F(1,10)],
            [F(1,10),F(1,2),F(3,10),F(1,10)],
            [F(1,5),F(3,10),F(2,5),F(1,10)],
            [F(1,10),F(1,10),F(1,10),F(7,10)]]
    ips=aggregate(inside,[[0,1,2],[3]])
    ck('internal_edge_source_is_reversible',all(sum(row)==1 for row in inside)
       and all(inside[i][j]==inside[j][i] for i,j in product(range(4),repeat=2)))
    ck('internal_edge_output_is_Markov',ips[0]==ips[1]==ips[2] and inside[0][1]>0)
    canon_squares=[inside[i][3]**2/(inside[i][i]*inside[3][3]) for i in range(3)]
    ck('internal_edge_necessity_countercontrol',canon_squares==[F(1,42),F(1,35),F(1,28)] and len(set(canon_squares))==3)
    # The complete balanced criterion, unlike the diagonal shortcut, passes.
    balanced_inside=[[x/4 for x in row] for row in inside]
    iblocks=aggregate(balanced_inside,[[0,1,2],[3]])
    ck('internal_edge_balanced_reference_passes',iblocks[0]==iblocks[1]==iblocks[2])
    tilt=[F(1),F(1),F(1),F(2)]
    tilted_flow=[[inside[i][j]*tilt[i]*tilt[j] for j in range(4)] for i in range(4)]
    tilted_mass=sum(map(sum,tilted_flow))
    tilted_pi=[sum(row)/tilted_mass for row in tilted_flow]
    tilted_kernel=[[x/sum(tilted_flow[i]) for x in row] for i,row in enumerate(tilted_flow)]
    tilted_blocks=aggregate(tilted_kernel,[[0,1,2],[3]])
    ck('one_marginal_closure_extends_to_new_marginal_fixture',tilted_pi==[F(11,67)]*3+[F(34,67)]
       and tilted_blocks[0]==tilted_blocks[1]==tilted_blocks[2])
    ck('all_inputs_unchanged',before=={name:sha(ROOT/name) for name in INPUTS})
    ck('prior_packages_unchanged',protected=={name:{'entries':manifest(ROOT/name),'sha256':sha(ROOT/name/'MANIFEST.sha256')} for name in PROTECTED})
    result={'date':'2026-10-20','model_id':'FIN-HYP-HIST-ME-001_UNCHANGED',
            'status':'PASS_BOUNDED_REFERENCE_TRANSPORT_CONTROLS',
            'new_bounded_checks':len(checks),'checks':checks,
            'arithmetic':'exact Fraction and Q(sqrt(5)); no floating theorem certificate',
            'q12_probability_fixture':'constructed exact even p, not actual A7 numerical probabilities; actual A7 application is qualitative',
            'new_theoretical_objects':['complete reference gauge test','weighted reference transport','canonical balanced all-marginal closure criterion','loop-only diagonal shortcut','nested mass-aware composition','finite-N count reference','q12 two-route inversion-orbit application'],
            'new_axiom':False,'primary_model_retuned':False,'physical_source_claim':False,
            'q12_fixture_radial_difference':str(jhyp-jhb),'nested_correct_invariant':str(correct_invariant),'nested_flattened_invariant':str(wrong_invariant),
            'golden_memory_defect_coefficients':[str(defect.a),str(defect.b)],
            'source_inputs':before,'protected_packages':protected,'strict_F01':'OPEN_SOURCE',
            'independent_peer_review':False,'goal_scope':'progress toward ongoing theory research, not completion'}
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    if args.record:
        assert args.output and args.output.resolve()==HERE/'RESULTS.json'
        (HERE/'DUPLICATE_SEARCH.json').write_text(json.dumps(duplicate_search(),indent=2)+'\n')
        (HERE/'INPUTS.sha256').write_text(''.join(f'{h}  {name}\n' for name,h in before.items()))
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in sorted(HERE.iterdir()) if p.is_file() and p.name!='MANIFEST.sha256'))
    if args.check:
        manifest(HERE)
        assert json.loads((HERE/'RESULTS.json').read_text())==result
        assert (HERE/'INPUTS.sha256').read_text()==''.join(f'{h}  {name}\n' for name,h in before.items())
        for p in HERE.glob('*.json'):json.loads(p.read_text())
        for p in HERE.glob('*.py'):ast.parse(p.read_text())
    print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'new_axiom':False,'strict_F01':'OPEN_SOURCE'}))


if __name__=='__main__':main()
