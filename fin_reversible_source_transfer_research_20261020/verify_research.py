#!/usr/bin/env python3
"""Exact small controls for RS; never a proof of the general theorem."""
import argparse
import ast
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INPUTS = [
    'fin_strategy_audit_20261019/SEARCH_FILE_INDEX.json',
    'FIN_RESEARCH_STRATEGY.md', 'research_tasks/00_EXECUTION_RULES.md',
    'FIN_CAUSAL_SOURCE_TRANSFER_NEXT_CHAIN_20261020.md',
    'FIN_Current_State_After_Programs_267_280_EN.md',
    'fin_transition_hierarchy/REPORT.md',
    'fin_auxiliary_mediator_review_20261016/PROOF_COMPLETIONS_AND_CORRECTIONS.md',
    'fin_intermediate_clock_bridge_research_20261017/ACTIVATION_BARRIER_SOURCE_GATE.md',
    'fin_evolution_source_research_20261020/PROOF.md',
    'fin_evolution_source_review_20261020/INTAKE_20261020.md',
    'fin_evolution_source_review_20261020/EXACT_CAUSAL_FILTER.md',
    'fin_evolution_source_review_20261020/FULL_PATH_TRANSFER.md',
    'fin_evolution_source_review_20261020/FILTERED_BLOCK_DEFECT.md',
]
BASE_BYTES = 1252727
BASE_SHA = '300f3864fde7843bd46f53840baebfc8190e3501c5d19de394a7266d3c1367e4'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def mul(a, b):
    return [[sum((a[i][k]*b[k][j] for k in range(len(b))), F(0))
             for j in range(len(b[0]))] for i in range(len(a))]


def sub(a, b):
    return [[x-y for x,y in zip(r,s)] for r,s in zip(a,b)]


def compress(s, mu, labels):
    q = max(labels)+1
    pi = [sum((mu[x] for x in range(len(mu)) if labels[x] == i), F(0)) for i in range(q)]
    a = [[sum((s[x][y] for y in range(len(mu)) if labels[y] == j), F(0))
          for j in range(q)] for x in range(len(mu))]
    c = [[sum((mu[x]*a[x][j] for x in range(len(mu)) if labels[x] == i), F(0))/pi[i]
          for j in range(q)] for i in range(q)]
    return pi, a, c


def trace(a):
    return sum((a[i][i] for i in range(len(a))), F(0))


def prefix_laws(s, mu, labels, m):
    """Unnormalized hidden distributions for every complete observed word."""
    q = max(labels)+1
    cur = {(i,):[mu[x] if labels[x] == i else F(0) for x in range(len(mu))]
           for i in range(q)}
    for _ in range(m):
        nxt = {}
        for prefix, v in cur.items():
            propagated = [sum((v[x]*s[x][y] for x in range(len(mu))), F(0)) for y in range(len(mu))]
            for j in range(q):
                row = [propagated[y] if labels[y] == j else F(0) for y in range(len(mu))]
                if sum(row):
                    nxt[prefix+(j,)] = row
        cur = nxt
    return {p:sum(v) for p,v in cur.items()}


def markov_word(c, pi, word):
    p = pi[word[0]]
    for a,b in zip(word,word[1:]):
        p *= c[a][b]
    return p


def record_tv(s, mu, labels, m):
    pi, _, c = compress(s, mu, labels)
    law = prefix_laws(s, mu, labels, m)
    return sum((abs(law.get(w,F(0))-markov_word(c,pi,w))
                for w in product(range(len(pi)), repeat=m+1)), F(0))/2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path)
    ap.add_argument('--record', action='store_true')
    ap.add_argument('--check', action='store_true')
    args = ap.parse_args()
    checks = []

    def ck(name, value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    agents = (ROOT/'AGENTS.md').read_bytes()
    ck('historical_instruction_prefix', len(agents) >= BASE_BYTES and
       hashlib.sha256(agents[:BASE_BYTES]).hexdigest() == BASE_SHA)
    mu = [F(i,21) for i in range(1,7)]
    s = eye(6)
    for i in range(6):
        for j in range(i+1,6):
            conductance = F(1+(i+j)%3,10000)
            s[i][j] = conductance/mu[i]
            s[j][i] = conductance/mu[j]
            s[i][i] -= s[i][j]
            s[j][j] -= s[j][i]
    ck('rational_kernel_stochastic', all(sum(r) == 1 and min(r) >= 0 for r in s))
    ck('nonuniform_detailed_balance', all(mu[i]*s[i][j] == mu[j]*s[j][i] for i,j in product(range(6), repeat=2)))
    ck('lazy_diagonal', all(s[i][i] >= F(1,2) for i in range(6)))
    cases = []
    partitions = [(0,0,1,1,2,2),(0,1,0,1,2,2),(0,0,0,0,1,1),(0,1,2,3,4,5),(0,0,0,0,0,0)]
    for case, labels in enumerate(partitions):
        pi,a,c = compress(s,mu,labels)
        _,_,c2 = compress(mul(s,s),mu,labels)
        d = sub(c2,mul(c,c))
        residual = [[a[x][j]-c[labels[x]][j] for j in range(len(pi))] for x in range(6)]
        gram = [[sum((mu[x]*residual[x][i]*residual[x][j] for x in range(6)),F(0))
                 for j in range(len(pi))] for i in range(len(pi))]
        ck(f'weighted_Gram_identity_{case}', all(pi[i]*d[i][j] == gram[i][j]
           for i,j in product(range(len(pi)),repeat=2)))
        delta = trace(d)
        ck(f'exact_trace_norm_{case}', delta == sum((gram[j][j]/pi[j] for j in range(len(pi))),F(0)) and delta >= 0)
        ck(f'weighted_not_unweighted_symmetry_{case}', all(pi[i]*d[i][j] == pi[j]*d[j][i]
           for i,j in product(range(len(pi)),repeat=2)))
        for m in (1,2,3) if len(pi) <= 3 else (1,2):
            tv = record_tv(s,mu,labels,m)
            ck(f'whole_prefix_TV_{case}_{m}', 4*tv*tv <= m*m*delta)
        if case in (0,1,2):
            ck(f'genuine_memory_{case}', delta > 0 and record_tv(s,mu,labels,2) > 0)
        else:
            ck(f'exact_lumpability_{case}', delta == 0 and record_tv(s,mu,labels,2) == 0)
        kappa = max(pi[j]/c[i][j] for i,j in product(range(len(pi)),repeat=2))
        # Exact conditional chi-square dominates KL, with the target row positive.
        for n in (0,1,2):
            prefix = {(i,):[mu[x] if labels[x] == i else F(0) for x in range(6)] for i in range(len(pi))}
            for _ in range(n):
                nxt = {}
                for word,v in prefix.items():
                    propagated = [sum((v[x]*s[x][y] for x in range(6)),F(0)) for y in range(6)]
                    for j in range(len(pi)):
                        nxt[word+(j,)] = [propagated[y] if labels[y] == j else F(0) for y in range(6)]
                prefix = nxt
            cost = F(0)
            for word,v in prefix.items():
                mass = sum(v)
                if not mass:
                    continue
                alpha = [sum((v[x]*a[x][j] for x in range(6)),F(0))/mass for j in range(len(pi))]
                cost += mass*sum(((alpha[j]-c[word[-1]][j])**2/c[word[-1]][j] for j in range(len(pi))),F(0))
            ck(f'conditional_chi_square_{case}_{n}', cost <= kappa*delta)
        cases.append({'partition':list(labels),'delta':str(delta),'three_sample_TV':str(record_tv(s,mu,labels,2))})

    # An independently implemented ES discrete analysis lift: 8 words x 3 phases.
    p = [F(1,3),F(2,3)]
    h = [F(2,3),F(-1,3)]
    words = list(product(range(2),repeat=3))
    weight = {}
    for w in words:
        prodp,prodh = F(1),F(1)
        for j in w:
            prodp *= p[j]
            prodh *= h[j]
        weight[w] = prodp*(1+F(1,2)*prodh)
    states = list(product(words,range(3)))
    where = {x:i for i,x in enumerate(states)}
    lifted = [[F(0) for _ in states] for _ in states]
    for x,(w,phase) in enumerate(states):
        if phase < 2:
            lifted[x][where[(w,phase+1)]] = 1
        else:
            for v in words:
                lifted[x][where[(v,0)]] = weight[v]
    lifted_mu = [weight[w]/3 for w,phase in states]
    labels = [w[phase] for w,phase in states]
    ck('ES_stationarity', all(sum((lifted_mu[x]*lifted[x][y] for x in range(len(states))),F(0)) == lifted_mu[y] for y in range(len(states))))
    ck('ES_source_NOT_reversible', any(lifted_mu[x]*lifted[x][y] != lifted_mu[y]*lifted[y][x] for x,y in product(range(len(states)),repeat=2)))
    pi,_,c = compress(lifted,lifted_mu,labels)
    _,_,c2 = compress(mul(lifted,lifted),lifted_mu,labels)
    ck('ES_exact_pair_semigroup', c == [p,p] and c2 == mul(c,c))
    es_tv = record_tv(lifted,lifted_mu,labels,2)
    ck('ES_reversibility_hypothesis_countercontrol', trace(sub(c2,mul(c,c))) == 0 and es_tv == F(16,2187))

    # Independently derive the actual continuous-time three-state spectrum.
    l = [[F(x) for x in r] for r in [[-1,1,0],[1,-2,1],[0,1,-1]]]
    vectors = [[1,1,1],[1,0,-1],[1,-2,1]]
    projectors = [[[F(v[i]*v[j],sum(x*x for x in v)) for j in range(3)] for i in range(3)] for v in vectors]
    ck('CTMC_projector_resolution', [[sum((z[i][j] for z in projectors),F(0)) for j in range(3)] for i in range(3)] == eye(3))
    for k,eigen in enumerate((0,-1,-3)):
        ck(f'CTMC_eigenspace_{k}', mul(l,projectors[k]) == [[eigen*x for x in row] for row in projectors[k]])
    semigroup = [[projectors[0][i][j]+projectors[1][i][j]/2+projectors[2][i][j]/8 for j in range(3)] for i in range(3)]
    mu3 = [F(1,3)]*3
    lab3 = [0,0,1]
    pi3,_,c = compress(semigroup,mu3,lab3)
    _,_,c2 = compress(mul(semigroup,semigroup),mu3,lab3)
    ck('CTMC_exact_log2_defect', trace(sub(c2,mul(c,c))) == F(27,1024))
    _,_,qbar = compress(l,mu3,lab3)
    ck('CTMC_mean_generator', qbar == [[F(-1,2),F(1,2)],[F(1),F(-1)]] and pi3 == [F(2,3),F(1,3)])
    _,a,_ = compress(l,mu3,lab3)
    v = sum((mu3[x]*sum(((a[x][j]-qbar[lab3[x]][j])**2/pi3[j] for j in range(2)),F(0)) for x in range(3)),F(0))
    ck('CTMC_infinitesimal_defect', v == F(3,4))
    ck('sqrt2_outward_rational', F(707,500)**2 < 2 < F(283,200)**2)
    ck('continuous_path_counterexample_lower_bound', F(4,9)*(F(13,32)-F(283,800)) == F(7,300))
    # Strict sqrt(2) inequality makes the actual lower bound strictly greater.
    ck('fast_family_defect_scales_quadratically', all(k*k*v == F(3*k*k,4) for k in (1,2,10)))

    for pth in HERE.glob('*.py'):
        ast.parse(pth.read_text())
    result = {'date':'2026-10-20','status':'PASS_BOUNDED_RS_CONTROLS',
              'bounded_checks':len(checks),'checks':checks,'exact_reversible_cases':cases,
              'ES_triple_TV':str(es_tv),'CTMC_log2_delta':str(F(27,1024)),
              'continuous_path_TV_strict_lower_bound':str(F(7,300)),
              'arithmetic':'stdlib Fraction; no numerical asymptotic proof',
              'general_theorem':'analytical self-audit, not mechanized or independently peer reviewed',
              'CSP_executed':False,'F01':'OPEN_SOURCE'}
    if args.output:
        args.output.write_text(json.dumps(result,indent=2)+'\n')
    if args.record:
        assert args.output and args.output.resolve() == HERE/'RESULTS.json'
        (HERE/'INPUTS.sha256').write_text(''.join(f'{sha(ROOT/p)}  {p}\n' for p in INPUTS))
        (HERE/'INSTRUCTION_CONTEXT.json').write_text(json.dumps({'historical_prefix_bytes':BASE_BYTES,'historical_prefix_sha256':BASE_SHA,
            'scope':'actual instruction context before the user-authorized append; follow live AGENTS thereafter'},indent=2)+'\n')
        claims = {'date':'2026-10-20','RS001':'DERIVED_CONDITIONAL_FINITE_LAG_GRAM_IDENTITY',
                  'RS002':'DERIVED_CONDITIONAL_WHOLE_SAMPLED_RECORD_TV_AND_KL',
                  'RS003':'DERIVED_EXACT_FINITE_REVERSIBLE_CTMC_CLOSURE_AND_ES_REPRESENTATION_EXCLUSION',
                  'RS004':'DERIVED_SCOPE_COUNTERCONTROLS_NO_COARSE_LAG_FULL_PATH_BOUND',
                  'reused':['HM generator compression','CF causal path inequalities','ES block counterfamily','P276 supplied-matrix consistency'],
                  'independent_provider':None,'CSP001_003':'NOT_ACTIVATED','F01_F03_F06':'OPEN_SOURCE',
                  'formalization':False,'independent_peer_review':False,'physical_promotion':False}
        (HERE/'CLAIM_REGISTER.json').write_text(json.dumps(claims,indent=2)+'\n')
        files = sorted(p for p in HERE.iterdir() if p.is_file() and p.name != 'MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in files))
    if args.check:
        for line in (HERE/'INPUTS.sha256').read_text().splitlines():
            h,name = line.split('  ',1)
            assert sha(ROOT/name) == h, name
        listed = set()
        for line in (HERE/'MANIFEST.sha256').read_text().splitlines():
            h,name = line.split('  ',1)
            assert name not in listed and Path(name).name == name and sha(HERE/name) == h, name
            listed.add(name)
        assert listed|{'MANIFEST.sha256'} == {p.name for p in HERE.iterdir() if p.is_file()}
        assert json.loads((HERE/'RESULTS.json').read_text()) == result
        search = json.loads((HERE/'SEARCH_RESULTS.json').read_text())
        assert all(sha(ROOT/p) == h for p,h in search['searched_files'].items())
    print(json.dumps({'status':result['status'],'bounded_checks':len(checks),'F01':result['F01'],'CSP':'NOT_ACTIVATED'}))


if __name__ == '__main__':
    main()
