#!/usr/bin/env python3
"""Bounded proof controls and source-preserving replay; not a PDE proof."""
import argparse
import ast
from fractions import Fraction as Q
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCE = ROOT / 'fin_rp004_flux_continuation_20261019'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(folder):
    return {str(p.relative_to(folder)): sha(p) for p in sorted(folder.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts}


def manifest(folder):
    rows = []
    declared = set()
    for line in (folder / 'MANIFEST.sha256').read_text().splitlines():
        digest, name = line.split(None, 1)
        name = name.lstrip('*')
        assert name not in declared
        assert not Path(name).is_absolute() and '..' not in Path(name).parts
        declared.add(name)
        p = folder / name
        state = 'MISSING' if not p.is_file() else (
            'MATCH' if sha(p) == digest else 'MISMATCH')
        rows.append({'path': name, 'expected_sha256': digest, 'state': state})
    # These three producer manifests are explicitly direct-file manifests.
    present = {p.name for p in folder.iterdir() if p.is_file()
               and p.name != 'MANIFEST.sha256'}
    return {'entries': len(rows), 'matching': sum(r['state'] == 'MATCH' for r in rows),
            'missing': [r['path'] for r in rows if r['state'] == 'MISSING'],
            'mismatches': [r['path'] for r in rows if r['state'] == 'MISMATCH'],
            'unlisted_direct_files': sorted(present - declared), 'rows': rows}


def original_inputs():
    ledgers = {}
    for folder in [SOURCE, SOURCE / 'RP004-A', SOURCE / 'RP004-B']:
        rows = []
        for line in (folder / 'INPUTS.sha256').read_text().splitlines():
            digest, name = line.split(None, 1)
            name = name.removeprefix('/mnt/data/fin_rp004_inputs/').removeprefix('/mnt/data/')
            p = ROOT / name
            state = 'NO_LOCAL_COUNTERPART' if not p.is_file() else (
                'MATCH' if sha(p) == digest else 'DIFFERENT_LOCAL_COUNTERPART')
            rows.append({'local_path': name, 'original_sha256': digest,
                         'current_sha256': sha(p) if p.is_file() else None,
                         'state': state})
        ledgers[str(folder.relative_to(SOURCE))] = rows
    return ledgers


def child_limits():
    resource.setrlimit(resource.RLIMIT_CPU, (60, 60))
    resource.setrlimit(resource.RLIMIT_AS, (1024 ** 3, 1024 ** 3))


def run(cmd, cwd):
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='1', MKL_NUM_THREADS='1')
    p = subprocess.run(cmd, cwd=cwd, env=env, text=True, capture_output=True,
                       timeout=60, preexec_fn=child_limits)
    assert p.returncode == 0, (cmd, p.returncode, p.stderr)
    return {'command': cmd, 'returncode': p.returncode,
            'stdout': p.stdout, 'stderr': p.stderr}


def scientific_controls():
    import sympy as s
    checks = []

    def ck(name, value):
        assert bool(value), name
        checks.append(name)

    t = s.symbols('t', real=True)
    nu, Om, q = s.symbols('nu Om q', positive=True)
    z = -(s.exp(-s.I * Om * t) - s.cos(nu * t)
          + s.I * Om / nu * s.sin(nu * t)) / (nu ** 2 - Om ** 2)
    ck('complex_retarded_oscillator_equation',
       s.simplify(s.diff(z, t, 2) + nu ** 2 * z + s.exp(-s.I * Om * t)) == 0)
    ck('complex_zero_position', z.subs(t, 0) == 0)
    ck('complex_zero_velocity', s.simplify(s.diff(z, t).subs(t, 0)) == 0)
    split = (-s.exp(-s.I * Om * t) / (nu ** 2 - Om ** 2)
             + s.exp(s.I * nu * t) / (2 * nu * (nu + Om))
             + s.exp(-s.I * nu * t) * (nu + Om) / (2 * nu * (nu ** 2 - Om ** 2)))
    ck('causal_split_before_i0', s.simplify(s.expand_complex((z - split).expand())) == 0)
    ck('real_resonance_sign', s.simplify(s.limit(s.re(z), nu, Om)
                                        + t * s.sin(Om * t) / (2 * Om)) == 0)
    # Coarse separation is exact, not a floating comparison near a front.
    ck('cone_behind_group_ray', Q(1, 3) ** 2 < Q(1, 5))
    ck('stationary_points_separated_from_poles', 4 * Q(1, 3) ** 2 /
       (1 - Q(1, 3) ** 2) == Q(1, 2) < 1)
    k = s.symbols('k', real=True)
    ck('ordinary_partial_fractions', s.simplify(1 / (k ** 2 - q ** 2) -
       (1 / (k - q) - 1 / (k + q)) / (2 * q)) == 0)
    ck('outgoing_green_derivative_jump', s.simplify(2 * s.I * q * s.I / (2 * q)) == -1)
    ck('incoming_sign_negative_control', s.simplify(2 * s.I * q * (-s.I) / (2 * q)) != -1)
    F4 = s.symbols('F4', positive=True)
    M = s.Rational(9, 4) * s.I * q * (q ** 2 + 4) * F4 / (8 * (1+s.I*q) * (2+s.I*q))
    M2 = s.simplify(s.expand_complex(M * s.conjugate(M)))
    gamma = 81*q*(q*q+4)*F4*F4 / (2048*(1+q*q))
    ck('both_roots_gamma_normalization', s.simplify(M2/(2*q)-gamma) == 0)
    ck('outgoing_amplitude_squared', s.simplify(M2/(4*q*q)-gamma/(2*q)) == 0)
    A2 = s.symbols('A2', positive=True)
    ck('plane_energy_average', s.simplify(A2*(Om**2+q**2+4)/4
       - A2*Om**2/2).subs(Om**2, q**2+4) == 0)
    ck('two_band_energy_factor', s.simplify(Om**2*(gamma/(2*q))/6
                                         - Om**2*gamma/(12*q)) == 0)
    c0 = 5 * Q(125, 73728) / 21
    ck('uniform_reference_constant', c0 == Q(625, 1548288) > 0)
    ck('nonlinear_coefficient', c0/8 == Q(625, 12386304))
    # Independent jet calculation of the full original energy balance.
    ph, ps, a, g4, g6 = s.symbols('ph ps a g4 g6', real=True)
    pt, px, st, sx, pxx, sxx, pxt, sxt = s.symbols('pt px st sx pxx sxx pxt sxt')
    V = (ph*ph-1)**2/2+(a+6*ph*ph)*ps*ps/2+g4*ps**4/4+g6*ps**6/6
    ptt, stt = pxx-s.diff(V, ph), sxx-s.diff(V, ps)
    et = pt*ptt+px*pxt+st*stt+sx*sxt+s.diff(V, ph)*pt+s.diff(V, ps)*st
    jx = -pxt*px-pt*pxx-sxt*sx-st*sxx
    ck('original_local_energy_identity', s.expand(et+jx) == 0)
    ep, Kx, ut, ux, wt, wx, vt, vx, rt, rx = s.symbols('ep Kx ut ux wt wx vt vx rt rx')
    J = -(ep**2*ut+wt)*(Kx+ep**2*ux+wx)-(ep*vt+rt)*(ep*vx+rx)
    inventory = (-ep**2*Kx*ut-ep**2*vt*vx-ep**4*ut*ux
                 -wt*(Kx+ep**2*ux)-ep**2*ut*wx-wt*wx
                 -ep*vt*rx-ep*rt*vx-rt*rx)
    ck('all_original_current_terms', s.expand(J-inventory) == 0)
    ck('omitted_kink_negative_control', s.expand(J-(inventory+ep**2*Kx*ut)) != 0)
    # epsilon=exp(-L): each error/reference ratio is polynomial times exp(-cL).
    L = s.symbols('L', positive=True)
    for power, degree in [(2, 5), (4, 7), (6, 1)]:
        ratio = s.exp(-power*L)*(1+12*L)**degree/(12*L)
        ck(f'log_window_ratio_{power}_{degree}', s.limit(ratio, L, s.oo) == 0)
    ck('bootstrap_log_smallness', s.limit(s.exp(-2*L)*(1+12*L)**3, L, s.oo) == 0)
    ck('B_leading_coefficient_not_small', Q(100)*Q(1, 1000)**4 /
       Q(1, 1000)**4 == 100)
    x = s.symbols('x', real=True)
    K = s.tanh(x)
    P = 3*K*K-3*s.I*k*K-(k*k+1)
    direct = s.exp(s.I*k*x)*P/((1-s.I*k)*(2-s.I*k))
    corr = s.exp(s.I*k*x)*(3*(K*K-1)-3*s.I*k*(K-1))/((1-s.I*k)*(2-s.I*k))
    ck('exact_jost_correction', s.simplify(direct-s.exp(s.I*k*x)-corr) == 0)
    ck('shape_projection_coefficient', s.Rational(9, 4)*(s.Rational(3, 8)-s.Rational(5, 16))
       == s.Rational(9, 64))
    return checks


def cone_diagnostics():
    """Truncated finite-time Jost integral, NOT an outward error certificate."""
    import numpy as np
    from scipy.integrate import quad_vec

    def sinc(y):
        return 1.0 if abs(y) < 1e-14 else math.sin(y)/y

    def coeff(nu, om, t):
        z = (nu-om)*t/2
        av = (nu+om)*t/2
        bh = -t*math.sin(av)*sinc(z)/(nu+om)
        bht = -(om*t*math.cos(av)*sinc(z)+math.sin(nu*t))/(nu+om)
        return -(1-math.cos(nu*t))/(nu*nu)+bh, -math.sin(nu*t)/nu+bht

    def m(k):
        f4 = 4/3 if abs(k) < 1e-12 else math.pi*k*(k*k+4)/(6*math.sinh(math.pi*k/2))
        return (9/4)*1j*k*(k*k+4)*f4/(8*(1+1j*k)*(2+1j*k))

    cases = [(-0.5, t, v) for t in [40, 80, 160] for v in [1/6, 1/3]]
    cases += [(a, 80, 1/3) for a in [-0.75, -0.25]]
    rows = []
    for a, t, v in cases:
        om = 2*math.sqrt(a+2)
        q = math.sqrt(om*om-4)
        x = t*v
        th = math.tanh(x)
        se = 1/math.cosh(x)

        def integrand(k):
            nu = math.sqrt(k*k+4)
            b, bt = coeff(nu, om, t)
            den = (1-1j*k)*(2-1j*k)
            poly = 3*th*th-3j*k*th-(k*k+1)
            mode = np.exp(1j*k*x)*poly/den
            modex = np.exp(1j*k*x)*(1j*k*poly+(6*th-3j*k)*(1-th*th))/den
            return np.array([(mode*m(k)).real*b, (mode*m(k)).real*bt,
                             (modex*m(k)).real*b])/math.pi

        val, err = quad_vec(integrand, 0, 18, epsabs=2e-8, epsrel=2e-8,
                            points=[q], limit=1000)
        b, bt = coeff(math.sqrt(3), om, t)
        proj = 9*math.pi*math.sqrt(3/2)/64
        shape = math.sqrt(3/2)*se*th
        shapex = math.sqrt(3/2)*se*(1-2*th*th)
        val += proj*np.array([shape*b, shape*bt, shapex*b])
        amp = -1j*m(q)/(2*q)*np.exp(1j*(q*x-om*t))
        outgoing = np.array([amp.real, (-1j*om*amp).real, (1j*q*amp).real])
        dif = np.abs(val-outgoing)
        assert all(math.isfinite(float(z)) for z in dif)
        rows.append({'a': a, 't': t, 'x_over_t': v,
                     'absolute_errors_u_ut_ux': dif.tolist(),
                     'sqrt_t_times_max_error': math.sqrt(t)*float(max(dif)),
                     'reported_quadrature_error_not_tail_bound': float(err)})
    return {'status': 'FINITE_POINT_FALSIFICATION_DIAGNOSTIC_ONLY',
            'momentum_cutoff': 18, 'quadrature': 'ordinary scipy quad_vec',
            'not_certified': ['tail error', 'uniform cone theorem', 'nonlinear PDE'],
            'rows': rows}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--replay', action='store_true')
    parser.add_argument('--output', type=Path, default=HERE/'verification.json')
    args = parser.parse_args()
    before = snapshot(SOURCE)
    strategy_before = snapshot(ROOT/'fin_strategy_audit_20261019')
    tasks_before = snapshot(ROOT/'research_tasks')
    guards = sha(ROOT/'AGENTS.md')
    manifests = {name: manifest(SOURCE/name) for name in ['.', 'RP004-A', 'RP004-B']}
    assert all(not v['mismatches'] and not v['unlisted_direct_files'] for v in manifests.values())
    assert manifests['.']['missing'] == ['finalize.py']
    assert not manifests['RP004-A']['missing'] and not manifests['RP004-B']['missing']
    # Later AGENTS tails are not silently equated with a historical complete input.
    freshness = run([sys.executable, 'fin_strategy_audit_20261019/verify_strategy.py',
                     '--inputs-only'], ROOT)
    logs = []
    if args.replay:
        with tempfile.TemporaryDirectory(prefix='fin-rp004-review-') as temp:
            cp = Path(temp)/'source'
            shutil.copytree(SOURCE, cp, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
            for script in ['verify_rp004.py', 'check_spectral_energy.py']:
                logs.append(run([sys.executable, script], cp))
    checks = scientific_controls()
    diagnostic = cone_diagnostics()
    assert snapshot(SOURCE) == before, 'Incoming source was modified'
    assert snapshot(ROOT/'fin_strategy_audit_20261019') == strategy_before
    assert snapshot(ROOT/'research_tasks') == tasks_before
    assert sha(ROOT/'AGENTS.md') == guards
    tree = ast.parse((SOURCE/'verify_rp004.py').read_text())
    result = {'date': '2026-10-19', 'status': 'PASS_BOUNDED_CONTROLS_WITH_DISCLOSED_SOURCE_OMISSION',
              'manifest_checks': manifests, 'original_inputs': original_inputs(),
              'source_preserved': before, 'strategy_preserved': strategy_before,
              'tasks_preserved': tasks_before, 'live_AGENTS_sha256': guards,
              'strategy_freshness': freshness, 'copy_only_source_replay': logs,
              'producer_symbolic_assert_statements': sum(isinstance(n, ast.Assert) for n in ast.walk(tree)),
              'intake_bounded_check_count': len(checks), 'intake_bounded_checks': checks,
              'new_finite_time_spectral_diagnostic': diagnostic,
              'analytic_proofs_not_established_by_flags': True,
              'independent_peer_review': False, 'physical_validation': False}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'checks': len(checks), 'source_files_preserved': len(before),
                      'manifests': {k: {n: v[n] for n in ['entries', 'matching', 'missing']}
                                    for k, v in manifests.items()},
                      'output': str(args.output)}, indent=2))


if __name__ == '__main__':
    main()
