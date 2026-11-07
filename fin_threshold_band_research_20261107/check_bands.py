"""Finite controls for the threshold budget and Darboux endpoint identity."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
checks = []


def check(name, value):
    checks.append({'name': name, 'pass': bool(value)})
    if not value:
        raise AssertionError(name)


def spectral_integral(lower, upper, exponent):
    return (upper ** (exponent + 1) - lower ** (exponent + 1)) / (exponent + 1)


budgets = []
for exponent in [0.0, 0.5, 1.0]:
    for epsilon in [1 / 16, 1 / 64, 1 / 256, 1 / 1024]:
        core = epsilon ** (2 * exponent) / (2 * exponent + 1)
        if exponent == 0:
            intermediate = np.log(1 / epsilon)
        else:
            intermediate = (1 - epsilon ** (2 * exponent)) / (2 * exponent)
        dyadic = core
        left = epsilon
        while left < 1:
            right = min(2 * left, 1)
            dyadic += spectral_integral(left, right, 2 * exponent) / left
            left = right
        continuous = core + intermediate
        check(f'dyadic equivalence alpha={exponent} epsilon={epsilon}', continuous <= dyadic <= 2 * continuous)
        # T=epsilon^-2 and the pointwise forcing envelope is epsilon^3 |k|^alpha.
        squared_budget = 2 * epsilon ** 2 * continuous
        budgets.append({'alpha': exponent, 'epsilon': epsilon, 'squared_budget': float(squared_budget), 'budget_divided_by_epsilon_squared': float(squared_budget / epsilon ** 2), 'dyadic_to_integral_ratio': float(dyadic / continuous)})

snapshot_results = []
for count in [2048, 4096]:
    half_domain = 20.0
    dx = 2 * half_domain / count
    x = np.arange(count) * dx - half_domain
    k = 2 * np.pi * np.fft.fftfreq(count, dx)
    kink = np.tanh(x)
    sech_squared = 1 / np.cosh(x) ** 2
    bound_even = np.sqrt(3 / 4) * sech_squared
    bound_odd = np.sqrt(3 / 2) / np.cosh(x) * kink
    u = np.exp(-x * x / 2)
    v = x * np.exp(-x * x / 2)
    u -= bound_even * dx * np.dot(bound_even, u)
    v -= bound_odd * dx * np.dot(bound_odd, v)
    source = -6 * v * v * u - u ** 3
    spectrum = np.fft.fft(source)
    source_x = np.fft.ifft(1j * k * spectrum).real
    source_xx = np.fft.ifft(-k * k * spectrum).real
    denominator = np.sqrt((1 + k * k) * (4 + k * k))
    direct = np.fft.ifft(np.fft.fft(source_xx + 3 * kink * source_x + 2 * source) / denominator).real
    commuted_spectrum = (-k * k * spectrum + 3j * k * np.fft.fft(kink * source) + np.fft.fft((2 - 3 * sech_squared) * source)) / denominator
    commuted = np.fft.ifft(commuted_spectrum).real
    relative_identity_error = float(np.linalg.norm(direct - commuted) / np.linalg.norm(direct))
    threshold_by_transform = float(dx * np.sum(direct) / np.sqrt(2 * np.pi))
    r0 = (3 * kink * kink - 1) / 2
    threshold_by_pairing = float(dx * np.dot(r0, source) / np.sqrt(2 * np.pi))
    check(f'projected even input {count}', abs(dx * np.dot(bound_even, u)) < 1e-12)
    check(f'projected odd input {count}', abs(dx * np.dot(bound_odd, v)) < 1e-12)
    check(f'commuted Darboux identity {count}', relative_identity_error < 1e-10)
    check(f'threshold adjoint identity {count}', abs(threshold_by_transform - threshold_by_pairing) < 1e-11)
    check(f'parity permits nonzero threshold {count}', abs(threshold_by_pairing) > 1e-6)
    snapshot_results.append({'grid_points': count, 'relative_identity_error': relative_identity_error, 'threshold_by_transform': threshold_by_transform, 'threshold_by_pairing': threshold_by_pairing})

for source_file in json.loads((HERE / 'INPUTS.json').read_text())['files']:
    raw = (HERE.parent / source_file['path']).read_bytes()
    check('source hash ' + source_file['path'], hashlib.sha256(raw).hexdigest() == source_file['sha256'])

results = {'scope': 'Budget identities and a projected comparison snapshot. Analytical uniform bounds are in the proof.', 'numpy_version': np.__version__, 'budgets': budgets, 'snapshot_results': snapshot_results, 'passed': len(checks), 'failed': 0, 'checks': checks}
(HERE / 'verification.json').write_text(json.dumps(results, indent=2) + '\n')
print(json.dumps({'passed': len(checks), 'failed': 0, 'snapshot_results': snapshot_results, 'smallest_epsilon_budgets': [row for row in budgets if row['epsilon'] == 1 / 1024]}, indent=2))
