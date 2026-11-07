"""Resonant reference history and spectral-modulus comparison diagnostics."""
from pathlib import Path
import hashlib
import json
import numpy as np

HERE = Path(__file__).resolve().parent
checks = []


def check(name, condition):
    checks.append({'name': name, 'pass': bool(condition)})
    if not condition:
        raise AssertionError(name)


def gauss_interval(lower, upper, count):
    nodes, weights = np.polynomial.legendre.leggauss(count)
    return (lower + upper) / 2 + (upper - lower) * nodes / 2, weights * (upper - lower) / 2


def resonant_history(epsilon, time_step, frequency_nodes, separation=0.0):
    # The diagnostic uses a fixed Gaussian momentum profile and carrier m.
    # The spectral oscillator solution is exact, avoiding a PDE time solver.
    mass = np.sqrt(16 / 3)
    final_time = epsilon ** -2
    steps = round(final_time / time_step)
    time_step = final_time / steps
    wave_numbers, frequency_weights = gauss_interval(0, 8, frequency_nodes)
    observation_x, observation_weights = gauss_interval(-1, 1, 64)
    observation_spectrum = np.cos(np.outer(wave_numbers, observation_x)) @ (observation_weights / np.sqrt(1 + observation_x ** 2)) / np.sqrt(2 * np.pi)
    frequency = np.sqrt(mass ** 2 + wave_numbers ** 2)
    delta = wave_numbers ** 2 / (frequency + mass)
    spatial_weights = 2 * frequency_weights * observation_spectrum * np.exp(-wave_numbers ** 2 / 2) * np.cos(wave_numbers * separation)
    squared_trace = 0.0
    for first in range(0, steps + 1, 128):
        indices = np.arange(first, min(first + 128, steps + 1))
        time = indices * time_step
        # (cos(m t)-cos(nu t))/(nu^2-m^2), evaluated stably near k=0.
        oscillator = time[:, None] * np.sin(np.outer(time, (frequency + mass) / 2))
        oscillator *= np.sinc(np.outer(time, delta) / (2 * np.pi)) / (frequency + mass)
        trace = epsilon ** 3 * (oscillator @ spatial_weights)
        weights = np.ones(len(indices))
        weights[indices == 0] = 0.5
        weights[indices == steps] = 0.5
        squared_trace += time_step * np.dot(weights, trace ** 2)
    return {'separation': separation, 'epsilon': epsilon, 'T': final_time, 'dt': time_step, 'frequency_nodes': frequency_nodes, 'squared_trace': float(squared_trace), 'trace_divided_by_epsilon_squared': float(squared_trace / epsilon ** 2)}


histories = []
for config in [(1 / 8, 0.1, 384), (1 / 16, 0.1, 768), (1 / 32, 0.1, 1536), (1 / 32, 0.05, 1536), (1 / 32, 0.1, 2048)]:
    result = resonant_history(*config)
    histories.append(result)
    print(result, flush=True)
check('resonant reference finite scale comparison', all(0 < row['trace_divided_by_epsilon_squared'] < 2 for row in histories))
check('time quadrature refinement', abs(histories[3]['squared_trace'] / histories[2]['squared_trace'] - 1) < 0.001)
check('frequency quadrature refinement', abs(histories[4]['squared_trace'] / histories[2]['squared_trace'] - 1) < 0.001)


def modulus_integral(separation, step):
    # y=kL, k in [0,1/2]. The value at zero is its continuous limit zero.
    count = int(np.ceil(separation / (2 * step)))
    y = np.linspace(0, separation / 2, count + 1)
    difference = np.cos(y) * np.exp(-0.5 * (y / separation) ** 2) - 1
    integrand = np.zeros_like(y)
    integrand[1:] = difference[1:] ** 2 / y[1:]
    return float(np.trapezoid(integrand, y))


translated_histories = [resonant_history(1 / 32, 0.1, 1536, separation) for separation in [16.0, 64.0, 256.0]]
for row in translated_histories:
    check(f'translated resonant history {row["separation"]}', 0 <= row['trace_divided_by_epsilon_squared'] < 2)

moduli = []
for separation in [16, 32, 64, 128, 256]:
    coarse = modulus_integral(separation, 0.02)
    fine = modulus_integral(separation, 0.01)
    check(f'modulus quadrature L={separation}', abs(fine - coarse) < 1e-4)
    moduli.append({'L': separation, 'weighted_difference_integral': fine, 'value_minus_three_halves_log_L': float(fine - 1.5 * np.log(separation))})
last_slope = (moduli[-1]['weighted_difference_integral'] - moduli[-2]['weighted_difference_integral']) / np.log(2)
check('logarithmic comparison slope', 1.3 < last_slope < 1.7)
for source in json.loads((HERE / 'INPUTS.json').read_text())['files']:
    raw = (HERE.parent / source['path']).read_bytes()
    check('source hash ' + source['path'], hashlib.sha256(raw).hexdigest() == source['sha256'])
output = {'scope': 'Finite resonant history and modulus diagnostics. Uniform bounds are proved analytically.', 'numpy_version': np.__version__, 'resonant_histories': histories, 'translated_histories': translated_histories, 'translated_profiles': moduli, 'last_logarithmic_slope': float(last_slope), 'passed': len(checks), 'failed': 0, 'checks': checks}
(HERE / 'verification.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps({'passed': len(checks), 'failed': 0, 'translated_histories': translated_histories, 'translated_profiles': moduli, 'last_logarithmic_slope': float(last_slope)}, indent=2))
