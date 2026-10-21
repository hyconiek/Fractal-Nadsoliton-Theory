#!/usr/bin/env python3
"""Bounded exact algebra diagnostics, NOT a proof of PDE long-time dynamics.

Polynomials are sparse rational Laurent polynomials in z,p,c,r,w. Negative
powers occur only in the positive parameter w, never in a state coordinate.
The physical factor d^2=27*pi^2/8192 is handled as a formal coefficient:
its independent linear component is tested, not rounded to a rational.
--record authors only this new packet's metadata; --check is read-only.
"""
from __future__ import annotations
import argparse
import ast
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from math import factorial
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
NVAR = 5
ZERO = (0,) * NVAR


def poly(x=0):
    return {ZERO: F(x)} if x else {}


def mono(powers, coefficient=1):
    return {tuple(powers): F(coefficient)} if coefficient else {}


def add(*ps):
    ans = {}
    for p in ps:
        for k, v in p.items():
            ans[k] = ans.get(k, F(0)) + v
    return {k: v for k, v in ans.items() if v}


def scale(p, a):
    a = F(a)
    return {k: a*v for k, v in p.items() if a*v}


def mul(p, q):
    ans = {}
    for k, a in p.items():
        for l, b in q.items():
            m = tuple(x+y for x, y in zip(k, l))
            ans[m] = ans.get(m, F(0)) + a*b
    return {k: v for k, v in ans.items() if v}


def power(p, n):
    assert type(n) is int and n >= 0
    ans = poly(1)
    for _ in range(n):
        ans = mul(ans, p)
    return ans


def deriv(p, i):
    ans = {}
    for k, a in p.items():
        if k[i]:
            m = list(k)
            m[i] -= 1
            ans[tuple(m)] = a*k[i]
    return ans


def bracket(p, q):
    return add(mul(deriv(p, 0), deriv(q, 1)),
               scale(mul(deriv(p, 1), deriv(q, 0)), -1),
               mul(deriv(p, 2), deriv(q, 3)),
               scale(mul(deriv(p, 3), deriv(q, 2)), -1))


z, p, c, r, w = [mono(tuple(int(i == j) for i in range(5))) for j in range(5)]
Jz = add(power(p, 2), mul(w, power(z, 2)))
Jc = add(power(r, 2), scale(power(c, 2), 3))
H2 = scale(add(Jz, Jc), F(1, 2))
delta = add(poly(3), scale(w, -4))
Q = scale(delta, 3)
h = add(scale(mul(power(z, 2), c), 6), scale(power(c, 3), 4))
# B=Q*G3/d, valid for every w, not only sampled parameter fixtures.
B = add(scale(mul(add(scale(w, 2), poly(-3)), mul(power(z, 2), r)), 6),
        scale(mul(mul(z, p), c), 36),
        scale(mul(power(p, 2), r), 12),
        scale(mul(delta, mul(power(c, 2), r)), -4),
        scale(mul(delta, power(r, 3)), F(-8, 9)))


def average(a):
    ans = {}
    for ex, val in a.items():
        if any(ex[j] % 2 for j in range(4)):
            continue
        i, j, k, l = (ex[s]//2 for s in range(4))
        weight = F(factorial(2*i)*factorial(2*j),
                   4**(i+j)*factorial(i)*factorial(j)*factorial(i+j))
        weight *= F(factorial(2*k)*factorial(2*l),
                    4**(k+l)*factorial(k)*factorial(l)*factorial(k+l))
        weight *= F(1, 3**k)
        term = mul(power(Jz, i+j), power(Jc, k+l))
        term = mul(mono((0, 0, 0, 0, ex[4]-i)), term)
        ans = add(ans, scale(term, val*weight))
    return ans


def at_w(a, value):
    value = F(value)
    assert value > 0
    ans = {}
    for ex, val in a.items():
        ans = add(ans, mono(ex[:4]+(0,), val*value**ex[4]))
    return ans


V4 = [mono(ex+(0,)) for ex in product(range(5), repeat=4)
      if sum(ex) == 4 and (ex[0]+ex[1]) % 2 == 0]
KEYS = [next(iter(a)) for a in V4]


def matrix(columns):
    assert all(set(a) <= set(KEYS) for a in columns)
    return [[a.get(key, F(0)) for a in columns] for key in KEYS]


def vec(a):
    assert set(a) <= set(KEYS)
    return [a.get(k, F(0)) for k in KEYS]


def from_vec(a):
    return {k: v for k, v in zip(KEYS, a) if v}


def determinant(a):
    a = [row[:] for row in a]
    out = F(1)
    for j in range(len(a)):
        k = next((k for k in range(j, len(a)) if a[k][j]), None)
        if k is None:
            return F(0)
        if k != j:
            a[j], a[k] = a[k], a[j]
            out = -out
        pivot = a[j][j]
        out *= pivot
        for i in range(j+1, len(a)):
            fac = a[i][j]/pivot
            for k in range(j+1, len(a)):
                a[i][k] -= fac*a[j][k]
            a[i][j] = F(0)
    return out


def solve(a, b):
    a = [row[:]+[val] for row, val in zip(a, b)]
    n = len(a)
    for j in range(n):
        k = next((k for k in range(j, n) if a[k][j]), None)
        if k is None:
            raise ValueError('singular homological matrix')
        a[j], a[k] = a[k], a[j]
        pivot = a[j][j]
        a[j] = [v/pivot for v in a[j]]
        for i in range(n):
            if i != j:
                fac = a[i][j]
                a[i] = [x-fac*y for x, y in zip(a[i], a[j])]
    return [row[-1] for row in a]


def quartic_generator(value, target):
    """Actual zero-average matrix construction (formal rational component).

    The input is a homogeneous quartic in the specified 19D class. A user
    must provide its actual coefficients; this is not a parameter fitter.
    """
    if type(value) not in (int, F) or not F(5, 4) <= value <= F(7, 4):
        raise ValueError('w must be exact in the frozen interval')
    target = at_w(target, value)
    quadratic = at_w(H2, value)
    D = lambda a: bracket(a, quadratic)
    Pi = lambda a: at_w(average(a), value)
    columns = [add(D(a), Pi(a)) for a in V4]
    rhs = add(target, scale(Pi(target), -1))
    generator = from_vec(solve(matrix(columns), vec(rhs)))
    assert D(generator) == rhs and not Pi(generator)
    return generator


def science():
    checks = {}
    def ck(name, condition):
        assert name not in checks
        checks[name] = bool(condition)
        assert condition, name

    ck('symbolic_cubic_homological_equation_all_w', bracket(B, H2) == mul(Q, h))
    ck('opposite_flow_is_not_cancellation', bracket(B, H2) != scale(mul(Q, h), -1))
    ck('quartic_even_pair_dimension_19', len(V4) == 19)
    ck('symbolic_projection_idempotent_all_w', all(average(average(a))==average(a) for a in V4))
    ck('symbolic_projection_D_commutation_all_w',
       all(not average(bracket(a,H2)) and not bracket(average(a),H2) for a in V4))
    ck('clear_denominator_uniform_gap', all(abs(3*(3-4*x)) >= 6
                                          for x in (F(5,4),F(3,2),F(7,4))))
    quartic_num = scale(bracket(h, B), F(1, 2))
    alpha_num = mul(mono((0,0,0,0,-2),F(9,4)), add(scale(w,8),poly(-9)))
    beta_num = poly(24)
    gamma_num = scale(delta,F(-5,3))
    expected = add(mul(alpha_num,power(Jz,2)), mul(beta_num,mul(Jz,Jc)),
                   mul(gamma_num,power(Jc,2)))
    ck('symbolic_quartic_average_d_squared_component', average(quartic_num)==expected)
    fixed = add(scale(mul(power(z,2),power(c,2)),F(18,35)),
                scale(power(c,4),F(9,70)))
    expected_fixed = add(mul(mono((0,0,0,0,-1),F(3,70)),mul(Jz,Jc)),
                         scale(power(Jc,2),F(3,560)))
    g4_part = scale(power(z,4),F(9,70))
    ck('symbolic_quartic_average_source_fixed', average(fixed)==expected_fixed)
    ck('symbolic_quartic_average_g4', average(g4_part)==
       mul(mono((0,0,0,0,-2),F(27,560)),power(Jz,2)))
    ck('cubic_momentum_odd', all((ex[1]+ex[3])%2==1 for ex in B))
    ck('cubic_psi_sign_even', all((ex[0]+ex[1])%2==0 for ex in B))
    # Exact moment integrals with t=tanh(x), avoiding transcendental quadrature.
    I=lambda degree: F(2,degree+1) if degree%2==0 else F(0)
    ck('source_f_fourth_moment', F(9,16)*(I(0)-3*I(2)+3*I(4)-I(6))==F(18,35))
    ck('source_shape_fourth_moment', F(9,4)*(I(4)-I(6))==F(9,35))
    ck('source_mixed_fourth_moment', F(9,8)*(I(2)-2*I(4)+I(6))==F(6,35))
    # beta(3/2,5/2)=beta(5/2,3/2)=pi/16 follows by symmetry;
    # this flag checks the prefactors only, not a numerical value for pi.
    ck('cubic_shape_prefactor_relation', F(2)*F(3,2)*F(1,16)==4*F(3,4)*F(1,16))
    values=(F(5,4),F(3,2),F(7,4))
    for value in values:
        tag=str(value).replace('/','_')
        quadratic=at_w(H2,value)
        D=lambda a: bracket(a,quadratic)
        Pi=lambda a: at_w(average(a),value)
        ck('projection_idempotent_'+tag, all(Pi(Pi(a))==Pi(a) for a in V4))
        ck('projection_annihilates_D_'+tag, all(not Pi(D(a)) and not D(Pi(a)) for a in V4))
        mat=matrix([add(D(a),Pi(a)) for a in V4])
        ck('determinant_formula_'+tag,
           determinant(mat)==2**20*value**3*3**3*(value-3)**2)
        components=(fixed,g4_part,scale(at_w(quartic_num,value),1/(3*(3-4*value))))
        generators=[quartic_generator(value,a) for a in components]
        ck('all_three_quartic_components_solved_'+tag,
           all(D(g)==add(at_w(a,value),scale(Pi(at_w(a,value)),-1))
               and not Pi(g) for a,g in zip(components,generators)))
        ck('quartic_generator_reversal_and_shape_parity_'+tag,
           all((ex[1]+ex[3])%2==1 and (ex[2]+ex[3])%2==0
               for g in generators for ex in g))
        A=6*(2*value-3)/(3*(3-4*value))
        C=12/(3*(3-4*value))
        ck('initial_slip_reconstruction_'+tag,
           A+C*value==F(-6,3) and A-C*value==-6/(3-4*value))
    ck('centre_prepared_point_cubic_fixed_exactly',
       all(not {ex:co for ex,co in at_w(deriv(B,i),F(3,2)).items()
                if ex[1]==ex[2]==ex[3]==0} for i in range(4)))
    ck('centre_shape_forcing_not_zero', h.get((2,0,1,0,0))==6)
    try:
        quartic_generator(F(3),fixed)
    except ValueError:
        ck('reject_outside_interval',True)
    else:
        ck('reject_outside_interval',False)
    try:
        quartic_generator(1.5,fixed)
    except ValueError:
        ck('reject_inexact_parameter',True)
    else:
        ck('reject_inexact_parameter',False)
    # At w=ell, a genuine mixed quartic resonance occurs outside the contract.
    q3=at_w(H2,3)
    pi3=lambda a: at_w(average(a),3)
    ck('outside_equal_frequency_matrix_singular',
       determinant(matrix([add(bracket(a,q3),pi3(a)) for a in V4]))==0)
    return {'status':'PASS_BOUNDED_EXACT_ALGEBRA_NOT_LONG_TIME_PDE_PROOF',
            'checks':checks,'check_count':len(checks),
            'parameter_fixtures':[str(x) for x in values],
            'arithmetic':'Fraction; symbolic w; d^2 treated as formal coefficient',
            'scope':'No full-PDE epsilon^-2 estimate or physical-source admission'}


INPUTS = [
    'AGENTS.md',
    'research_tasks/RP004_SHARED_CONTRACT.json',
    'fin_slow_memory_modulation_20261021/PROOF_OR_OPEN_LEMMA.md',
    'fin_slow_memory_modulation_20261021/GATE_DELTA.md',
    'fin_slow_memory_modulation_20261021/MANIFEST.sha256',
    'fin_nonlinear_mode_memory_20261021/NONLINEAR_REDUCTION.md',
    'fin_nonlinear_mode_memory_20261021/SPECTRAL_MEMORY_KERNEL.md',
    'fin_nonlinear_mode_memory_20261021/NEXT_OPEN_LEMMA.md',
    'fin_nonlinear_mode_memory_20261021/MANIFEST.sha256',
    'fin_rp004_flux_review_20261019/OUTWARD_FLUX_THEOREM.md',
    'fin_rp004_flux_review_20261019/BOOTSTRAP_COMPLETION.md',
    'fin_rp004_flux_review_20261019/MANIFEST.sha256',
    'fin_rp004_full_input_self_audit_20261020/CAUSAL_AND_SPECTRAL_AUDIT.md',
    'fin_rp004_full_input_self_audit_20261020/MANIFEST.sha256',
]
PROTECTED = [
    'fin_slow_memory_modulation_20261021/MANIFEST.sha256',
    'fin_nonlinear_mode_memory_20261021/MANIFEST.sha256',
    'fin_rp004_flux_review_20261019/MANIFEST.sha256',
    'fin_rp004_full_input_self_audit_20261020/MANIFEST.sha256',
]


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def manifest_check(path, root):
    count=0
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        h, name = line.split(None,1)
        name=name.lstrip('*')
        candidates=[root/name, ROOT/name]
        target=next((p for p in candidates if p.is_file() and digest(p)==h),None)
        assert target is not None, (str(path),name)
        count+=1
    return count


def record_or_check(record):
    result=science()
    protected={rel:manifest_check(ROOT/rel,(ROOT/rel).parent) for rel in PROTECTED}
    result['protected_manifest_entry_counts']=protected
    result['protected_entry_total']=sum(protected.values())
    result['input_count']=len(INPUTS)
    for f in HERE.glob('*.py'):
        ast.parse(f.read_text())
    for f in HERE.glob('*.json'):
        json.loads(f.read_text())
    for f in HERE.glob('*.md'):
        txt=f.read_text()
        assert all(line==line.rstrip() for line in txt.splitlines()),str(f)
        prose=re.sub(r'```.*?```','',txt,flags=re.S)
        prose=re.sub(r'`[^`]*`','',prose)
        for link in re.findall(r'\[[^\]]*\]\(([^)]+)\)',prose):
            if not re.match(r'^[a-z]+://',link) and not link.startswith('#'):
                assert (f.parent/link.split('#')[0]).is_file(),(str(f),link)
    inputs=''.join(f'{digest(ROOT/rel)}  {rel}\n' for rel in INPUTS)
    payload=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if record:
        (HERE/'RESULTS.json').write_text(payload)
        (HERE/'INPUTS.sha256').write_text(inputs)
        files=sorted(p for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256')
        (HERE/'MANIFEST.sha256').write_text(''.join(f'{digest(p)}  {p.name}\n' for p in files))
    else:
        assert (HERE/'RESULTS.json').read_text()==payload,'bounded replay result changed'
        assert (HERE/'INPUTS.sha256').read_text()==inputs,'current inputs changed'
        manifest_check(HERE/'MANIFEST.sha256',HERE)
        names={line.split(None,1)[1].lstrip('*') for line in (HERE/'MANIFEST.sha256').read_text().splitlines()}
        assert names=={p.name for p in HERE.iterdir() if p.is_file() and p.name!='MANIFEST.sha256'}
    print(json.dumps({'status':result['status'],'bounded_checks':result['check_count'],
                      'inputs':len(INPUTS),'protected_entries':sum(protected.values()),
                      'mode':'record_own_packet_only' if record else 'read_only_check'},sort_keys=True))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    choice=parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--record',action='store_true')
    choice.add_argument('--check',action='store_true')
    args=parser.parse_args()
    record_or_check(args.record)
