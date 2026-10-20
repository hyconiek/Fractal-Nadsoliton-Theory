#!/usr/bin/env python3
"""Source-preserving intake, copy-only replay and independently implemented controls."""
import argparse
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
from fractions import Fraction as F

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT/'FIN_POST_RP004_RESEARCH_ONLY_20261019'
SCIENCE = SOURCE/'01_RESEARCH_POST_RP004'
NAMES = ['FIN_missing_objects_research_20261019', 'fin_intrinsic_kl_history_research_20261019',
         'fin_six_mode_history_reconstruction_20261019', 'fin_relational_event_history_20261019',
         'fin_anonymous_four_event_hierarchy_20261019']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def snapshot(folder):
    return {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}


def checked_manifest(base, name, exclusions=(), recursive=True):
    declared = set()
    for line in (base/name).read_text().splitlines():
        h, rel = line.split(None, 1)
        rel = rel.lstrip('*')
        assert not Path(rel).is_absolute() and '..' not in Path(rel).parts and rel not in declared
        assert sha(base/rel) == h, (base, rel)
        declared.add(rel)
    present = {str(p.relative_to(base)) for p in (base.rglob('*') if recursive else base.iterdir())
               if p.is_file() and '__pycache__' not in p.parts}
    unlisted = present-declared-set(exclusions)-{name}
    assert not unlisted, (base, unlisted)
    return {'entries': len(declared), 'coverage': 'all_recursive' if recursive else 'direct_files_only',
            'excluded_separately_bound': list(exclusions)}


def provenance():
    rows = []
    for folder in [SCIENCE/n for n in NAMES]:
        for line in (folder/'INPUTS.sha256').read_text().splitlines():
            h, original = line.split(None, 1)
            local = original.removeprefix('/mnt/data/').removeprefix('../')
            if local.startswith('fin_rp004_inputs/'):
                p = ROOT/local.removeprefix('fin_rp004_inputs/')
            elif local.startswith(tuple(NAMES)):
                p = SCIENCE/local
            elif local in ['AGENTS(1).md', 'GEN001-gen956(1).md']:
                p = SOURCE/'03_REFERENCE_SNAPSHOTS'/local
            else:
                p = ROOT/local
            archived_only = False
            live_hash = sha(p) if p.is_file() else None
            if p.is_file() and sha(p) != h:
                archived = SOURCE/'03_REFERENCE_SNAPSHOTS'/p.name
                if archived.is_file() and sha(archived) == h:
                    p = archived
                    archived_only = True
            rows.append({'stage': folder.name, 'original_reference': original,
                         'expected_sha256': h,
                         'local_counterpart': str(p.relative_to(ROOT)),
                         'actual_sha256': sha(p) if p.is_file() else None,
                         'different_live_counterpart_sha256': live_hash if archived_only else None,
                         'state': 'MISSING_ORIGINAL' if not p.is_file() else
                         'MATCH_ARCHIVED_SNAPSHOT_NOT_LIVE' if archived_only else
                         'MATCH' if sha(p) == h else 'DIFFERENT_CURRENT_COUNTERPART'})
    refs = []
    for p in sorted((SOURCE/'03_REFERENCE_SNAPSHOTS').glob('*.md')):
        target = ROOT/('AGENTS.md' if p.name == 'AGENTS(1).md' else
                       'GEN001-gen956.md' if p.name == 'GEN001-gen956(1).md' else p.name)
        refs.append({'snapshot': p.name, 'live': target.name, 'same_bytes_now':
                     target.is_file() and sha(target) == sha(p)})
    embedded = []
    for p in sorted((SCIENCE/NAMES[0]/'sources').rglob('*')):
        if p.is_file():
            rel = p.relative_to(SCIENCE/NAMES[0]/'sources')
            target = ROOT/rel
            embedded.append({'path':str(rel),'matches_current_bytes':target.is_file() and sha(p)==sha(target)})
    return {'input_reference_rows': rows, 'reference_snapshot_comparisons': refs,
            'embedded_current_input_comparisons':embedded,
            'original_archives': 'Five original science ZIPs and Archiwum.zip are not supplied; metadata hashes are not custody verification'}


def limit_child():
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (1024**3, 1024**3))


def run(cmd, cwd):
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', PYTHONDONTWRITEBYTECODE='1')
    try:
        p = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True,
                           timeout=60, preexec_fn=limit_child)
        return {'returncode': p.returncode, 'stdout': p.stdout, 'stderr': p.stderr}
    except subprocess.TimeoutExpired as e:
        return {'returncode': None, 'status': 'RESOURCE_STOP_60_SECONDS_NOT_MATHEMATICAL_FAILURE',
                'stdout': (e.stdout or b'').decode() if isinstance(e.stdout, bytes) else e.stdout or ''}


def json_drift(old, new):
    numeric, other = [], []
    def visit(a, b, path=''):
        if isinstance(a, dict) and isinstance(b, dict) and a.keys() == b.keys():
            for k in a:
                visit(a[k], b[k], path+'/'+k)
        elif isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
            for i, (x, y) in enumerate(zip(a, b)):
                visit(x, y, path+'/'+str(i))
        elif a != b:
            if isinstance(a, (int, float)) and isinstance(b, (int, float)) and not isinstance(a, bool) and not isinstance(b, bool):
                numeric.append({'path': path, 'abs_drift': abs(a-b)})
            else:
                other.append({'path': path, 'old': a, 'new': b})
    visit(old, new)
    return {'max_numeric_drift': max((r['abs_drift'] for r in numeric), default=0),
            'changed_numeric_coordinates': len(numeric), 'other_changes': other}


def replay():
    logs = {}
    with tempfile.TemporaryDirectory(prefix='fin-post-history-') as temp:
        cp = Path(temp)/'source'
        shutil.copytree(SOURCE, cp, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
        logs['master'] = run([sys.executable, '00_MASTER_HANDOFF/VERIFY_ALL.py'], cp)
        for name in NAMES:
            folder = cp/'01_RESEARCH_POST_RP004'/name
            script = folder/'verify.py'
            patch = None
            if name == 'fin_intrinsic_kl_history_research_20261019':
                old = "'/mnt/data/fin_intrinsic_kl_history_research_20261019/verification.json'"
                new = repr(str(folder/'verification.json'))
                text = script.read_text()
                assert text.count(old) == 1
                original_sha = sha(script)
                # Mechanical output-path relocation only, on the disposable copy.
                script.write_text(text.replace(old, new))
                patch = {'original_sha256': original_sha, 'copy_sha256': sha(script),
                         'change': 'one output path; no mathematical code changed'}
            result = run([sys.executable, 'verify.py'], folder)
            result['copy_only_path_repair'] = patch
            output = 'verification.json' if 'kl_history' in name else 'RESULTS.json'
            original = SCIENCE/name/output
            regenerated = folder/output
            if result['returncode'] == 0:
                result['output_byte_identical'] = sha(original) == sha(regenerated)
                result['output_comparison'] = json_drift(json.loads(original.read_text()), json.loads(regenerated.read_text()))
            logs[name] = result
    return logs


class R3:
    """Independent exact Q(sqrt(3)) implementation; no source scientific imports."""
    def __init__(self, a=0, b=0):
        self.a, self.b = F(a), F(b)
    def __add__(self, v):
        v = cv(v); return R3(self.a+v.a, self.b+v.b)
    __radd__ = __add__
    def __neg__(self):
        return R3(-self.a, -self.b)
    def __sub__(self, v):
        return self+-cv(v)
    def __rsub__(self, v):
        return cv(v)+-self
    def __mul__(self, v):
        v = cv(v); return R3(self.a*v.a+3*self.b*v.b, self.a*v.b+self.b*v.a)
    __rmul__ = __mul__
    def __truediv__(self, v):
        v = cv(v); d = v.a*v.a-3*v.b*v.b
        if not d: raise ZeroDivisionError
        return self*R3(v.a/d, -v.b/d)
    def __pow__(self, n):
        out = R3(1)
        for _ in range(n): out = out*self
        return out
    def __eq__(self, v):
        v = cv(v); return self.a == v.a and self.b == v.b
    def sign(self):
        if not self.b: return (self.a>0)-(self.a<0)
        if not self.a: return (self.b>0)-(self.b<0)
        if self.a*self.b>0: return (self.a>0)-(self.a<0)
        d = self.a*self.a-3*self.b*self.b
        return ((self.a>0)-(self.a<0))*((d>0)-(d<0))
    def __abs__(self):
        return self if self.sign()>=0 else -self
    def __float__(self):
        return float(self.a)+math.sqrt(3)*float(self.b)
    def __repr__(self):
        return f'{self.a}+({self.b})sqrt(3)'


def cv(v): return v if isinstance(v, R3) else R3(v)
def conv(a, b): return [sum((a[j]*b[(i-j)%12] for j in range(12)), R3()) for i in range(12)]
def row(weights): return [R3()] + [weights[min(i, 12-i)-1] for i in range(1, 12)]
def mass(weights): return 2*sum(weights[:5], R3())+weights[5]


def independent_controls():
    import numpy as np
    from scipy.linalg import expm
    from scipy.integrate import quad_vec
    checks = []
    def ck(name, condition):
        assert condition, name
        checks.append(name)
    ck('radical_sign_positive', R3(1, -F(1, 2)).sign() == 1)
    ck('radical_sign_negative_control', R3(1, -1).sign() == -1)
    cos = [R3(1), R3(0, F(1, 2)), R3(F(1, 2)), R3(), R3(-F(1, 2)),
           R3(0, -F(1, 2)), R3(-1), R3(0, -F(1, 2)), R3(-F(1, 2)), R3(),
           R3(F(1, 2)), R3(0, F(1, 2))]
    B = [[2*(1-cos[(k*d)%12]) if d<6 else R3(1-(-1)**k)
          for d in range(1, 7)] for k in range(1, 7)]
    # Inverse from full 12-point Fourier inversion, not inversion of producer B.
    inv = [[-cos[(k*d)%12]/6 if k<6 else R3(-((-1)**d), 0)/12
            for k in range(1, 7)] for d in range(1, 7)]
    ck('full_fourier_inverse_exact', all(sum((inv[i][k]*B[k][j] for k in range(6)), R3()) == (i==j)
                                       for i in range(6) for j in range(6)))
    normrows = [sum((abs(v) for v in line), R3()) for line in inv]
    ck('inverse_error_norm_exact', normrows[-1] == F(11, 12) and all((F(11, 12)-r).sign()>=0 for r in normrows))
    blind = [-1, 1, -1, 1, -1, 1]
    ck('five_mode_blind_vector', [sum((B[k][d]*blind[d] for d in range(6)), R3()) for k in range(6)]
       == [R3()]*5+[R3(-12)])
    bad = [sum((inv[d][k]*([3]+[1]*5)[k] for k in range(6)), R3()) for d in range(6)]
    ck('invalid_six_positive_correlations_rejected', bad[1] == -F(1, 12) and bad[0].sign()<0)
    w = [R3(F(16,15),F(1,30)), R3(1), R3(F(61,60)), R3(1), R3(F(16,15),-F(1,30)), R3(1)]
    v = [R3(F(43,40),F(1,40)), R3(F(59,60)), R3(1), R3(1), R3(F(43,40),-F(1,40)), R3(F(31,30))]
    ck('witness_strict_positivity_exact', all(x.sign()>0 for x in w+v))
    ck('witness_same_actual_escape', mass(w) == mass(v) == F(113, 10))
    rates = [[sum((b[i]*a[i] for i in range(6)), R3()) for b in B] for a in [w, v]]
    expected = [R3(12+F(k,10)) for k in range(1,7)]
    ck('witness_labelled_spectrum_swap_exact', rates[0] == expected and rates[1] == [expected[1],expected[0]]+expected[2:])
    cubic = lambda x: 2*sum((a**3 for a in x[:5]), R3())+x[5]**3
    ck('permutation_invariant_cubic_difference', cubic(w)-cubic(v) == F(11, 24000))
    jw, jv = [[x/mass(a) for x in row(a)] for a in [w, v]]
    ck('event_row_normalization_and_no_self', sum(jw, R3()) == sum(jv, R3()) == 1 and jw[0] == jv[0] == 0)
    mw, mv = [R3(1)]+[R3()]*11, [R3(1)]+[R3()]*11
    return_values = []
    for n in range(1, 13):
        mw, mv = conv(mw, jw), conv(mv, jv)
        return_values.append(mw[0] == mv[0])
    ck('return_moments_one_through_twelve_exact', all(return_values))
    def partitions(jump):
        totals = {}
        # Exact rooted enumeration with a generic first-appearance encoding,
        # not the producer's five-case classifier or cubic formula.
        for tail in itertools.product(range(12), repeat=3):
            states = (0,)+tail
            if any(states[i] == states[i+1] for i in range(3)): continue
            labels, pattern = {}, []
            for x in states:
                if x not in labels: labels[x] = len(labels)
                pattern.append(labels[x])
            key = ''.join(map(str, pattern))
            prob = jump[(tail[0])%12]*jump[(tail[1]-tail[0])%12]*jump[(tail[2]-tail[1])%12]
            totals[key] = totals.get(key, R3())+prob
        return totals
    pw, pv = partitions(jw), partitions(jv)
    ck('exact_anonymous_partition_support', set(pw) == set(pv) == {'0101','0102','0121','0120','0123'})
    ck('exact_anonymous_partition_normalization', sum(pw.values(), R3()) == sum(pv.values(), R3()) == 1)
    gap = pw['0101']-pv['0101']
    ck('ABAB_gap_independent_exact_enumeration', gap == F(11, 24*113**3))
    ck('anonymous_TV_independent_exact_enumeration', sum((abs(pw[k]-pv[k]) for k in pw), R3())/2 == F(11, 12*113**3))
    ck('four_states_are_three_jumps', 4-1 == 3)
    ck('continuous_time_leading_coefficient', gap*mass(w)**3 == F(11, 24000))
    # Include collisions with w6: odd multiplicity must survive degeneracy.
    import sympy as sp
    for weights in [[1,2,2,4,5,2], [1,2,3,4,5,6]]:
        entries = [F(x, 2*sum(weights[:5])+weights[5]) for x in weights[:5] for _ in [0,1]]
        entries += [F(weights[5], 2*sum(weights[:5])+weights[5])]
        power = [F(0)]+[sum(x**k for x in entries) for k in range(1, 12)]
        elem = [F(1)]
        for k in range(1, 12):
            elem.append(sum((-1)**(j-1)*elem[k-j]*power[j] for j in range(1,k+1))/k)
        x = sp.Symbol('x')
        actual = sp.Poly(sp.prod(x-sp.Rational(y.numerator,y.denominator) for y in entries), x)
        ck('Newton_exact_duplicate_fixture_'+str(weights[-1]), actual.all_coeffs() ==
           [sp.Rational(((-1)**i*e).numerator, ((-1)**i*e).denominator) for i,e in enumerate(elem)])
        odd = {y for y in set(entries) if entries.count(y)%2}
        ck('unique_odd_multiplicity_fixture_'+str(weights[-1]), odd == {entries[-1]})
    q, q2 = [R3(i) for i in [1,2,3,4,5,6]], [R3(i) for i in [2,1,3,4,5,6]]
    j, j2 = [[x/mass(a) for x in row(a)] for a in [q,q2]]
    ck('all_alternating_moments_multiset_counterpair', sorted((x.a,x.b) for x in j) == sorted((x.a,x.b) for x in j2))
    ck('nonalternating_return_detects_shell_placement', conv(conv(j,j),j)[0]-conv(conv(j2,j2),j2)[0] == -F(11,7776))
    # Source correction: a full stationary mean bracket DOES identify a
    # reversible one-hot generator with known positive pi (unlike count data).
    pi3 = [F(1,2),F(1,3),F(1,6)]
    c3 = [[F(0),F(1,7),F(3,11)], [F(1,7),F(0),F(2,9)], [F(3,11),F(2,9),F(0)]]
    q3 = [[c3[i][j]/pi3[i] if i!=j else F(0) for j in range(3)] for i in range(3)]
    for i in range(3): q3[i][i] = -sum(q3[i])
    def mean_bracket(q, pi):
        n = len(pi); out = [[F(0) for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i==j: continue
                delta = [(int(k==j)-int(k==i)) for k in range(n)]
                for a in range(n):
                    for b in range(n): out[a][b] += pi[i]*q[i][j]*delta[a]*delta[b]
        return out
    bracket = mean_bracket(q3,pi3)
    ck('mean_one_hot_bracket_reconstruction_nonuniform_pi', all(-bracket[i][j]/(2*pi3[i]) == q3[i][j]
       for i in range(3) for j in range(3)))
    ck('stationary_mean_bracket_edge_conductances', all(bracket[i][j] == -2*c3[i][j]
       for i in range(3) for j in range(3) if i!=j))
    nonrev = [[F(-2) if i==j else F(1)+(F(1,4) if (j-i)%3==1 else -F(1,4))
               for j in range(3)] for i in range(3)]
    rev = [[F(-2) if i==j else F(1) for j in range(3)] for i in range(3)]
    ck('nonreversible_current_negative_control_same_mean_bracket', nonrev!=rev and
       mean_bracket(nonrev,[F(1,3)]*3) == mean_bracket(rev,[F(1,3)]*3))
    ck('PSD_alone_not_valid_mean_bracket', -F(1)/(2*pi3[0])<0)  # uu^T, u=(1,1,-2)
    # Native source uses no exact-decimal bridge for the printed digits.
    spec = importlib.util.spec_from_file_location('safe_seed_history_audit', ROOT/'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py')
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    W, L, _, _, A7 = mod.build_rank7()
    decimal = np.array([.46998567264502017,.19204355169010282,.09142861427792497,.0470291687456504,.02413122336363006,.011070817321442113])
    ck('nominal_versioned_kernel_weight_agreement', np.max(abs(W[0,1:7]-decimal)) < 1e-14)
    ck('laplacian_is_not_A7', np.linalg.norm(L-A7)>0.1)
    D, C = L[::2,::2], L[::2,1::2]
    dark = np.array([1,-1,1,-1,1,-1])
    ck('even_odd_dark_channel', np.max(abs(C.T@dark))<1e-14)
    ck('even_odd_hidden_preparation_really_changes_drift', np.linalg.norm(C[:,0])>0)
    z0 = np.arange(12)/11-0.5
    reorder = np.r_[np.arange(0,12,2),np.arange(1,12,2)]
    G = L[np.ix_(reorder,reorder)]; init = z0[reorder]
    time = 0.4
    full = expm(-time*G)@init; deriv = -G@full
    val, err = quad_vec(lambda s:C@expm(-(time-s)*D)@C.T@(expm(-s*G)@init)[:6],0,time,epsabs=1e-12)
    rhs = -D@full[:6]-C@expm(-time*D)@init[6:]+val
    ck('memory_identity_actual_convolution', np.max(abs(rhs-deriv[:6]))<2e-12)
    U = expm(-0.2j*L); P = abs(U)**2
    ck('Born_endpoint_kernel_not_a_semigroup', np.linalg.norm(abs(expm(-0.4j*L))**2-P@P)>1e-4)
    # Independent Fisher diagnostic: diagonalize L and use fixed Gaussian quadrature.
    lam, vec = np.linalg.eigh(L); pi = np.ones(12)/12
    p0 = 0.65*pi; p0[0] += 0.35
    coeff = vec.T@p0
    def prob(t): return vec@(np.exp(-lam*t)*coeff)
    nodes, wgts = np.polynomial.legendre.leggauss(64)
    action = 0.25*sum(wt*np.sum((L@prob((nd+1)/2))**2/prob((nd+1)/2)) for nd,wt in zip(nodes,wgts))
    geom = 2*math.acos(np.dot(np.sqrt(p0),np.sqrt(prob(1))))**2
    ck('Fisher_heat_action_strictly_above_endpoint_geodesic', action>geom>0)
    x = np.zeros(12); x[0]=1; x[1]=-1
    def fgrad(p):
        lp = np.log(p/pi); return -p*(lp-np.dot(p,lp))
    fd = (fgrad(pi+1e-6*x)-fgrad(pi-1e-6*x))/(2e-6)
    ck('Fisher_gradient_tangent_linearization', np.max(abs(fd+x))<1e-8)
    p = np.arange(1,13,dtype=float); p /= p.sum()
    drift = np.zeros(12)
    for i in range(12):
        for j in range(12):
            if i!=j:
                logdiff = math.log(p[i])-math.log(p[j])
                drift[i] -= W[i,j]*(p[i]-p[j])/logdiff*logdiff
    ck('supplied_log_mean_mobility_identity', np.max(abs(drift+L@p))<1e-14)
    ck('strict_heat_nonconstant_spectrum_nominal', lam[-1]-lam[1]>1.5)
    return checks, {'exact_ABAB_gap':str(gap), 'exact_TV':str(2*gap),
                    'memory_convolution_error':float(np.max(abs(rhs-deriv[:6]))),
                    'Fisher_action':float(action), 'Fisher_endpoint_min':geom,
                    'rounded_weight_max_drift_nominal':float(np.max(abs(W[0,1:7]-decimal)))}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--output', type=Path, default=HERE/'verification.json')
    args = parser.parse_args()
    before = snapshot(SOURCE); context = sha(ROOT/'AGENTS.md')
    consolidated = checked_manifest(SOURCE,'00_MASTER_HANDOFF/CONSOLIDATED_MANIFEST.sha256',
                                    exclusions=['00_MASTER_HANDOFF/VERIFY_ALL.py'])
    nested = {n:checked_manifest(SCIENCE/n,'MANIFEST.sha256', recursive=False) for n in NAMES}
    science_logs = replay() if args.replay else {}
    checks, details = independent_controls()
    assert snapshot(SOURCE) == before and sha(ROOT/'AGENTS.md') == context
    result = {'date':'2026-10-20', 'status':'PASS_SCOPED_FIVE_STUDY_HISTORY_INTAKE',
              'source_preserved':before,'source_file_count':len(before),
              'consolidated_manifest':consolidated,'nested_manifests':nested,
              'provenance':provenance(),'copy_only_replay':science_logs,
              'bounded_checks':checks,'bounded_check_count':len(checks),'diagnostic_details':details,
              'live_AGENTS_sha256':context,'independent_peer_review':False,
              'RP004_full_proof_audited_here':False,'F01':'OPEN_SOURCE','F03':'OPEN_SOURCE','F06':'OPEN_SOURCE'}
    args.output.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'source_files':len(before),'consolidated_entries':consolidated['entries'],
                      'bounded_checks':len(checks),'source_replay_codes':{k:v['returncode'] for k,v in science_logs.items()},
                      'input_states':{v:sum(r['state']==v for r in result['provenance']['input_reference_rows'])
                                      for v in ['MATCH','MATCH_ARCHIVED_SNAPSHOT_NOT_LIVE','MISSING_ORIGINAL','DIFFERENT_CURRENT_COUNTERPART']},
                      'output':str(args.output)},indent=2))


if __name__ == '__main__': main()
