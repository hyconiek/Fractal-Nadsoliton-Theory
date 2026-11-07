"""Finite quadrature diagnostics for the causal and spectral interfaces."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
checks = []


def check(name, condition):
    checks.append({'name': name, 'pass': bool(condition)})
    if not condition:
        raise AssertionError(name)


def quadrature(lower, upper, count):
    nodes, weights = np.polynomial.legendre.leggauss(count)
    return (lower + upper) / 2 + (upper - lower) * nodes / 2, weights * (upper - lower) / 2


def source_spectrum(k):
    # Unitary transform of W-adjoint applied to f cubed.
    transform = np.pi * k * (k * k + 4) * (k * k + 16) / (120 * np.sinh(np.pi * k / 2))
    transform *= (3 / 4) ** 1.5 / np.sqrt(2 * np.pi)
    return -4 / 7 * np.sqrt((1 + k * k) / (4 + k * k)) * transform


spectrum_errors = []
for count in [512, 1024]:
    x, weights = quadrature(-20, 20, count)
    sech_squared = 1 / np.cosh(x) ** 2
    phi = (3 / 4) ** 1.5 * sech_squared ** 3
    differential = 20 * phi - 24 * sech_squared * phi
    wave_numbers = np.array([0.5, 0.75, 1.0, 2.0])
    numerical = (np.cos(np.outer(wave_numbers, x)) @ (weights * differential))
    numerical /= np.sqrt(2 * np.pi * (1 + wave_numbers ** 2) * (4 + wave_numbers ** 2))
    error = float(np.max(np.abs(numerical - source_spectrum(wave_numbers))))
    spectrum_errors.append(error)
    check(f'Darboux spectral identity with {count} nodes', error < 1e-10)

annular_norms = []
for count in [64, 128]:
    k, weights = quadrature(0.5, 1, count)
    annular_norms.append(float(np.sqrt(2 * np.dot(weights, source_spectrum(k) ** 2))))
check('annulus quadrature refinement', abs(annular_norms[1] - annular_norms[0]) < 1e-12)
check('nonzero source annulus', annular_norms[1] > 0.05)

mass_squared = 16 / 3
annular_mass_bounds = []
for width in [4, 8, 16, 32]:
    atom_bound = np.sqrt((mass_squared + 1) * 1.5 * width / np.sqrt(np.pi)) * np.exp(-width ** 2 / 8)
    necessary_mass = annular_norms[1] / (2 * atom_bound)
    annular_mass_bounds.append({'R': width, 'atom_annular_upper_bound': float(atom_bound), 'necessary_mass_for_half_annular_error': float(necessary_mass)})
check('atomic mass obstruction increases', all(a['necessary_mass_for_half_annular_error'] < b['necessary_mass_for_half_annular_error'] for a, b in zip(annular_mass_bounds, annular_mass_bounds[1:])))

# An exact forced free spectral oscillator checks the source-time integral.
# Momentum forcing is cos(m t) Omega h_R and the birth atom is sin(Omega age) h_R.
k = np.linspace(-2, 2, 401)
frequency = np.sqrt(k * k + mass_squared)
carrier = np.sqrt(mass_squared)
final_time = 8.0
width = 4
packet_spectrum = np.sqrt(width) * np.pi ** (-0.25) * np.exp(-(width * k) ** 2 / 2)
exact = np.empty_like(k)
nonzero = np.abs(k) > 1e-14
exact[nonzero] = frequency[nonzero] * packet_spectrum[nonzero] * (np.cos(carrier * final_time) - np.cos(frequency[nonzero] * final_time)) / k[nonzero] ** 2
exact[~nonzero] = packet_spectrum[~nonzero] * final_time * np.sin(carrier * final_time) / 2
causal_runs = []
for steps in [256, 512, 1024]:
    source_time = np.linspace(0, final_time, steps + 1)
    time_weights = np.ones(steps + 1) * final_time / steps
    time_weights[[0, -1]] *= 0.5
    delayed_atoms = np.sin(np.outer(frequency, final_time - source_time)) * packet_spectrum[:, None]
    causal = delayed_atoms @ (time_weights * np.cos(carrier * source_time))
    relative_error = float(np.linalg.norm(causal - exact) / np.linalg.norm(exact))
    coefficient_mass = float(np.dot(time_weights, np.abs(np.cos(carrier * source_time))))
    causal_runs.append({'steps': steps, 'relative_spectral_error': relative_error, 'source_time_coefficient_mass': coefficient_mass})
    check(f'causal quadrature {steps}', relative_error < 0.001)
for coarse, fine in zip(causal_runs, causal_runs[1:]):
    check(f'causal second-order refinement {fine["steps"]}', 3.8 < coarse['relative_spectral_error'] / fine['relative_spectral_error'] < 4.2)

for source in json.loads((HERE / 'INPUTS.json').read_text())['files']:
    raw = (ROOT / source['path']).read_bytes()
    check('source hash ' + source['path'], hashlib.sha256(raw).hexdigest() == source['sha256'])

output = {
    'scope': 'Finite identity and quadrature diagnostics, with source provenance',
    'numpy_version': np.__version__,
    'spectrum_absolute_errors': spectrum_errors,
    'annular_source_norms': annular_norms,
    'atomic_mass_lower_bounds': annular_mass_bounds,
    'causal_quadrature': causal_runs,
    'passed': len(checks),
    'failed': 0,
    'checks': checks,
}
(HERE / 'verification.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps({key: output[key] for key in ['annular_source_norms', 'atomic_mass_lower_bounds', 'causal_quadrature', 'passed', 'failed']}, indent=2))
