#!/usr/bin/env python3
"""Bounded local checks for post-326 preparation/physical-bridge methodology.

Reuses the versioned small-N reconstruction; does not claim an independent
implementation of its Gibbs generator or a physical hardware replay.
"""
from pathlib import Path
import hashlib
import importlib.util
import itertools
import json
import types
import numpy as np

ROOT = Path(__file__).resolve().parents[1]

def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    obj = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(obj)
    return obj

checks = {}
for name in ('fin_physical_bridge_campaign_01', 'fin_physical_bridge_campaign_02', 'fin_physical_bridge_campaign_03'):
    root = ROOT / name
    entries = [line.split(None, 1) for line in (root / 'MANIFEST.sha256').read_text().splitlines() if line.strip()]
    bad = [p for h, p in entries if not (root / p).is_file() or hashlib.sha256((root / p).read_bytes()).hexdigest() != h]
    assert not bad, (name, bad)
    checks[name] = {'manifest_entries': len(entries), 'bad': bad}

seed_path = ROOT / 'fin_physical_bridge_campaign_01/PHYS-002/safe_seed.py'
seed = module('post326_seed', seed_path)
tests = module('post326_existing_tests', seed_path.with_name('test_safe_seed.py'))
passed = []
for name in sorted(x for x in vars(tests) if x.startswith('test_')):
    getattr(tests, name)()
    passed.append(name)
checks['smallN_existing_assertion_tests'] = passed
A = seed.build_rank7()[-1]
lam = np.fft.fft(A[0]).real
errors = {}
for N in (1, 2, 3):
    states = seed.compositions(N)
    eps = 1e-5
    plus = seed.count_pi(states, eps, A)
    minus = seed.count_pi(states, -eps, A)
    err = []
    for k in range(1, 7):
        slope = (seed.static_S(states, plus, N, k) - seed.static_S(states, minus, N, k)) / (2 * eps)
        err.append(abs(slope - (N - 1) * lam[k] / (12 * N)))
    errors[str(N)] = max(err)
    assert max(err) < 1e-8
checks['static_slope_fd_max_error'] = errors

g = seed.G_FROZEN
contrast = 0.0
for a, b, c in itertools.product((0, 1), (0, 2), (0, 4)):
    n = np.bincount([a, b, c], minlength=12)
    field = (g / 4) * (A @ n)
    log_weight = (g / 8) * (n @ A @ n) + field.max() + np.log(np.exp(field - field.max()).sum())
    contrast += (-1) ** ((a == 1) + (b == 2) + (c == 4)) * log_weight
assert abs(contrast) > 0.2
checks['marginal_four_to_three'] = {'contrast': contrast, 'sign_convention': 'plus for (0,0,0); a pairwise log weight has zero alternating third difference'}

path = ROOT / 'fin_physical_bridge_campaign_03/PHYS-013/ctmc_replay.py'
source = path.read_text()
start = source.index('# Requires')
end = source.index('\ng=3.0', start)
source = source[:start] + source[end:]
# Explicit binding replaces a missing historical fin_core import; prints muted.
env = {'fc': types.SimpleNamespace(A7_default=lambda: (A, None)), '__file__': str(path), 'print': lambda *args, **kwargs: None}
exec(compile(source, str(path), 'exec'), env)
checks['PHYS013_nominal_CTMC'] = {'tv': float(env['tv']), 'source_sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bound_A7_source_sha256': hashlib.sha256(seed_path.read_bytes()).hexdigest(), 'physical_calibration_recomputed': False}
assert abs(env['tv'] - 0.0004127833002122317) < 1e-12
print(json.dumps(checks, indent=2))
