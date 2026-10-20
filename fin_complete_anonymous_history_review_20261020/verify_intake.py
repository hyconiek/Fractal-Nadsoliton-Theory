#!/usr/bin/env python3
"""Bounded source intake and independent rational AH controls, not a theorem prover."""
import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT/'fin_complete_anonymous_history_campaign_20261020'
HIST = ROOT/'FIN_POST_RP004_RESEARCH_ONLY_20261019'


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def snapshot(folder):
    return {str(p.relative_to(folder)):sha(p) for p in sorted(folder.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}


def manifest():
    rows = []
    for line in (SOURCE/'MANIFEST.sha256').read_text().splitlines():
        digest, name = line.split(None, 1)
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        assert name not in rows and sha(SOURCE/name)==digest, name
        rows.append(name)
    assert set(rows)|{'MANIFEST.sha256'} == set(snapshot(SOURCE))
    return {'entries':len(rows),'complete_coverage':True}


def provenance():
    mapping = {
        'archived_AGENTS': ROOT/'AGENTS.md',
        'archived_GEN001_956': ROOT/'GEN001-gen956.md',
        'predecessor_E3_E4': HIST/'01_RESEARCH_POST_RP004/fin_anonymous_four_event_hierarchy_20261019/PAPER.md',
        'predecessor_E2': HIST/'01_RESEARCH_POST_RP004/fin_relational_event_history_20261019/PAPER.md',
        'predecessor_H1': HIST/'01_RESEARCH_POST_RP004/fin_six_mode_history_reconstruction_20261019/PAPER.md',
        'master_handoff': HIST/'00_MASTER_HANDOFF/HANDOFF_MASTER.md',
        'current_instructions': ROOT/'Wklejony tekst(2).txt',
    }
    rows = []
    for r in json.loads((SOURCE/'INPUTS.json').read_text())['inputs']:
        p = mapping[r['alias']]
        rows.append({**r, 'current_counterpart':str(p.relative_to(ROOT)),
                     'actual_sha256':sha(p) if p.is_file() else None,
                     'actual_bytes':p.stat().st_size if p.is_file() else None,
                     'state':'MISSING_ORIGINAL' if not p.is_file() else
                     'MATCH' if sha(p)==r['sha256'] and p.stat().st_size==r['size_bytes'] else 'DIFFERENT_COUNTERPART'})
    return rows


def limits():
    resource.setrlimit(resource.RLIMIT_CPU,(60,60))
    resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))


def replay():
    logs = {}
    with tempfile.TemporaryDirectory(prefix='fin-ah-intake-') as temp:
        cp = Path(temp)/'source'
        shutil.copytree(SOURCE,cp,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        for script in ['verify_package.py','verify.py']:
            r = subprocess.run([sys.executable,script],cwd=cp,text=True,capture_output=True,
                               timeout=60,preexec_fn=limits,
                               env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1'))
            assert r.returncode==0,(script,r.stderr)
            logs[script]={'returncode':r.returncode,'stdout':r.stdout,'stderr':r.stderr}
        old = json.loads((SOURCE/'RESULTS.json').read_text())
        new = json.loads((cp/'RESULTS.json').read_text())
        logs['result_comparison']={'byte_identical':sha(SOURCE/'RESULTS.json')==sha(cp/'RESULTS.json'),
                                   'JSON_equal':old==new}
    return logs


def ceilq(x): return -((-x.numerator)//x.denominator)


def log_bounds(x, terms):
    """Positive atanh series plus exact geometric remainder, with binary reduction."""
    x = F(x); assert x>=1 and terms>0
    exponent = 0
    y = x
    while y>=2: y/=2; exponent+=1
    def small(v):
        z = (v-1)/(v+1)
        lower = 2*sum((z**(2*k+1)/F(2*k+1) for k in range(terms)),F(0))
        remainder = 2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
        return lower,lower+remainder
    l2,u2 = small(F(2)); ly,uy = small(y)
    return exponent*l2+ly,exponent*u2+uy


def make_certificate(a=F(1,20),delta=F(1,10),eta=F(1,20),terms=32):
    assert 0<a<=F(1,11) and 0<delta<=1 and 0<eta<1
    b = 10*a*a
    logs = [log_bounds(24/eta,terms),log_bounds(528/eta,terms)]
    thresholds = [(2/b,2/b),tuple(8*x/(b*b) for x in logs[0]),
                  tuple(384*x/(delta*delta*b*b) for x in logs[1])]
    lower = max(t[0] for t in thresholds); upper = max(t[1] for t in thresholds)
    assert ceilq(lower)==ceilq(upper),'Refine the analytical interval before claiming exact ceil'
    m = ceilq(upper)
    return {'scope':'FORMULA_ARITHMETIC_ONLY_NOT_A_PROBABILITY_SIMULATION',
            'a':str(a),'delta':str(delta),'eta':str(eta),'b':str(b),'series_terms':terms,
            'log_intervals':[[str(x) for x in pair] for pair in logs],
            'threshold_intervals':[[str(x) for x in pair] for pair in thresholds],
            'max_threshold_interval':[str(lower),str(upper)],'m':m,'L':2*m+2,
            'not_licensed_for_weights_1_to_6_with_this_a':True}


def check_certificate(c):
    expected = make_certificate(F(c['a']),F(c['delta']),F(c['eta']),int(c['series_terms']))
    assert c==expected
    assert F(c['m']-1)<F(c['max_threshold_interval'][0])<=F(c['max_threshold_interval'][1])<=c['m']


def radial(weights):
    s = 2*sum(weights[:5])+weights[5]
    return [[F(0) if i==j else F(weights[min((i-j)%12,(j-i)%12)-1],s)
             for j in range(12)] for i in range(12)]


def multiply(a,b):
    n = len(a)
    return [[sum(a[i][k]*b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def encode(path):
    labels = {}; out = []
    for x in path:
        if x not in labels: labels[x]=len(labels)
        out.append(labels[x])
    return out


def estimate(code,nstates=12):
    rows = [0]*nstates; edges = [[0]*nstates for _ in range(nstates)]
    for a,b in zip(code,code[1:]):
        rows[a]+=1; edges[a][b]+=1
    mat = [[F(edges[i][j],rows[i]) if rows[i] else F(0) for j in range(nstates)] for i in range(nstates)]
    return rows,edges,mat


def scientific_controls():
    checks = []
    def ck(name,condition):
        assert condition,name
        checks.append(name)
    fixtures = [[1,2,3,4,5,6],[1,1,1,1,1,1],[1,1,1,1,1,10],[2,3,2,4,2,3]]
    for ix,w in enumerate(fixtures):
        j = radial(w)
        a = min(j[i][k] for i in range(12) for k in range(12) if i!=k)
        b = 10*a*a; beta = 12*b; j2 = multiply(j,j)
        ck(f'fixture_{ix}_stochastic_symmetric_no_self',all(sum(row)==1 for row in j) and
           all(j[i][k]==j[k][i] for i in range(12) for k in range(12)) and all(j[i][i]==0 for i in range(12)))
        ck(f'fixture_{ix}_two_step_lower_all_entries',min(v for row in j2 for v in row)>=b)
        ck(f'fixture_{ix}_proper_Doeblin_coefficient',0<beta<1)
        rem = [[(j2[i][k]-b)/(1-beta) for k in range(12)] for i in range(12)]
        ck(f'fixture_{ix}_Doeblin_residual_stochastic',min(x for row in rem for x in row)>=0 and all(sum(row)==1 for row in rem))
    path = [7,2,7,5,2]
    code = encode(path)
    ck('first_appearance_generic_fixture',code==[0,1,0,2,1])
    perm = [(5*x+3)%12 for x in range(12)]
    ck('code_invariance_under_nontrivial_permutation',encode([perm[x] for x in path])==code)
    rows,edges,hat = estimate(code)
    ck('terminal_state_not_an_extra_departure',rows[:3]==[2,1,1] and sum(rows)==len(path)-1)
    ck('edge_and_departure_totals_agree',all(sum(edges[i])==rows[i] for i in range(12)))
    ck('unvisited_rows_are_zero_padded',all(x==0 for row in hat[3:] for x in row))
    covered_only_at_end = estimate(encode(list(range(12))))[0]
    ck('all_names_seen_does_not_imply_all_departure_denominators',sum(v>0 for v in covered_only_at_end)==11)
    # Deterministic Euler tour: tests exact reconstruction, not sampling probability.
    available = {i:[j for j in range(12) if j!=i] for i in range(12)}
    stack = [0]; tour = []
    while stack:
        i = stack[-1]
        if available[i]: stack.append(available[i].pop())
        else: tour.append(stack.pop())
    tour.reverse()
    rt,et,ht = estimate(encode(tour))
    ck('Euler_fixture_uses_each_legal_edge_once',len(tour)==133 and all(rt[i]==11 for i in range(12)) and
       all(et[i][j]==int(i!=j) for i in range(12) for j in range(12)))
    flat = radial([1]*6)
    ck('exact_anonymous_estimator_Euler_flat_fixture',ht==flat)
    ck('Euler_anonymous_invariance_after_renaming',encode([perm[x] for x in tour])==encode(tour))
    for m in [1,2,3,17]:
        n = 2*m+1; L = 2*m+2
        ck(f'index_payment_m_{m}',n==L-1 and list(range(2,2*m+1,2))[-1]==n-1)
        # u=delta*b*m/8; u^2/(2n)>=delta^2*b^2*m/384.
        ck(f'exponent_payment_m_{m}',F(m*m,128*n)>=F(m,384))
    ck('ordered_edge_union_count',len([(i,j) for i in range(12) for j in range(12) if i!=j])==132)
    ck('Azuma_two_sided_union_coefficient',2*132==264)
    ck('ratio_error_payment',F(1,8)/F(1,2)==F(1,4))
    ck('nearest_orbit_decision_has_strict_margin',1-F(1,4)==F(3,4)>F(1,4))
    j0 = radial([1,2,3,4,5,6]); j1 = radial([1,1,1,1,1,10])
    def accepts_common_a(a, matrices):
        return 0<a<=F(1,11) and all(q[i][j]>=a for q in matrices for i in range(12) for j in range(12) if i!=j)
    ck('example_a_1_20_is_not_valid_for_1_to_6_fixture',not accepts_common_a(F(1,20),[j0,flat]))
    ck('common_lower_bound_must_hold_for_both_hypotheses',not accepts_common_a(F(1,15),[flat,j0]))
    ck('illustrative_budget_parameters_have_admissible_pair',accepts_common_a(F(1,20),[flat,j1]))
    # Flat J is invariant under every permutation, so this exact max is d_iso.
    delta_iso = max(abs(j1[i][k]-flat[i][k]) for i in range(12) for k in range(12))
    ck('admissible_pair_orbit_separation_no_factorial_search',delta_iso==F(9,22)>F(1,10))
    cert = make_certificate(); check_certificate(cert)
    ck('sample_budget_exact_ceil',cert['m']==569231067 and cert['L']==1138462136)
    bad = deepcopy(cert); bad['m']-=1; bad['L']=2*bad['m']+2
    try: check_certificate(bad)
    except AssertionError: rejected=True
    else: rejected=False
    ck('undersized_integer_budget_rejected',rejected)
    for x in [F(1),F(2),F(3,2),F(480),F(10560)]:
        lo,hi=log_bounds(x,32)
        ck('exact_log_interval_order_'+str(x),lo<=hi and (lo==hi==0 if x==1 else lo>0))
    ck('coverage_integer_threshold_payment',F(cert['m'])*F(cert['b'])/2>=1)
    return checks,cert


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--replay',action='store_true')
    parser.add_argument('--output',type=Path,default=HERE/'verification.json'); args = parser.parse_args()
    before=snapshot(SOURCE); agents=sha(ROOT/'AGENTS.md')
    manifest_report=manifest(); prov=provenance(); logs=replay() if args.replay else {}
    checks,cert=scientific_controls()
    assert before==snapshot(SOURCE) and sha(ROOT/'AGENTS.md')==agents
    result={'date':'2026-10-20','status':'PASS_SCOPED_AH_INTAKE_WITH_EXACT_ARITHMETIC_SUPPORT',
            'source_preserved':before,'source_manifest':manifest_report,'source_file_count':len(before),
            'original_inputs':prov,'copy_only_replay':logs,'bounded_checks':checks,'bounded_check_count':len(checks),
            'example_budget_certificate':cert,'live_AGENTS_sha256':agents,
            'independent_external_peer_review':False,'formal_proof_assistant':False,
            'RP004_audit_executed':False,'F01':'OPEN_SOURCE','F03':'OPEN_SOURCE','F06':'OPEN_SOURCE'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'bounded_checks':len(checks),'source_files':len(before),'manifest':manifest_report,
                      'input_rows':len(prov),'input_matches':sum(r['state']=='MATCH' for r in prov),
                      'missing_originals':[r['filename'] for r in prov if r['state']=='MISSING_ORIGINAL'],
                      'budget_m':cert['m'],'budget_L':cert['L'],'output':str(args.output)},indent=2))


if __name__=='__main__': main()
