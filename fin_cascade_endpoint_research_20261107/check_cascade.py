"""Full signed cascade on projected comparison snapshots, with grid refinement."""
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


def snapshot(count, half_domain=20.0):
    dx = 2 * half_domain / count
    x = np.arange(count) * dx - half_domain
    k = 2 * np.pi * np.fft.fftfreq(count, dx)
    ordered_k = np.fft.fftshift(k)
    mass_squared = 16 / 3
    omega_phi = np.sqrt(k * k + 4)
    omega_psi = np.sqrt(k * k + mass_squared)
    kink = np.tanh(x)
    r0 = (3 * kink * kink - 1) / 2
    f = np.sqrt(3 / 4) / np.cosh(x) ** 2
    g = np.sqrt(3 / 2) * kink / np.cosh(x)

    def project(field, bound):
        return field - bound * dx * np.dot(bound, field)

    gaussian = np.exp(-x * x / 2)
    u = project(gaussian, f)
    v = project(x * gaussian, g)
    pu = project(0.3 * x * x * gaussian, f)
    pv = project(0.2 * x * gaussian, g)
    force_phi = -6 * kink * (v * v + u * u)
    force_psi = -12 * kink * v * u

    def inverse_omega(field, omega):
        return np.fft.ifft(np.fft.fft(field) / omega)

    h_phi = {sign: v + sign * 1j * inverse_omega(pv, omega_phi) for sign in [-1, 1]}
    h_psi = {sign: u + sign * 1j * inverse_omega(pu, omega_psi) for sign in [-1, 1]}
    fg_phi = {sign: sign * 1j * inverse_omega(force_phi, omega_phi) for sign in [-1, 1]}
    fg_psi = {sign: sign * 1j * inverse_omega(force_psi, omega_psi) for sign in [-1, 1]}
    frequencies_p = ordered_k[:, None]
    frequencies_q = ordered_k[None, :]
    output_frequency = frequencies_p + frequencies_q
    output_omega = np.sqrt(output_frequency ** 2 + mass_squared)
    input_phi_omega = np.sqrt(frequencies_p ** 2 + 4)
    input_psi_omega = np.sqrt(frequencies_q ** 2 + mass_squared)
    sum_indices = (np.arange(count)[:, None] + np.arange(count)[None, :]).ravel()

    def bilinear(first, second, symbol):
        first_hat = np.fft.fftshift(np.fft.fft(first))
        second_hat = np.fft.fftshift(np.fft.fft(second))
        # Exclude the unpaired Nyquist mode to preserve conjugate symmetry.
        first_hat[0] = 0
        second_hat[0] = 0
        products = (symbol * first_hat[:, None] * second_hat[None, :]).ravel() / count
        # Linear convolution followed by spectral truncation avoids wraparound aliasing.
        full = np.bincount(sum_indices, weights=products.real, minlength=2 * count - 1)
        full = full + 1j * np.bincount(sum_indices, weights=products.imag, minlength=2 * count - 1)
        retained = full[count // 2:3 * count // 2].copy()
        retained[0] = 0
        return np.fft.ifft(np.fft.ifftshift(retained))

    cascade = {}
    for sigma in [-1, 1]:
        total = np.zeros(count, dtype=complex)
        for tau in [-1, 1]:
            for rho in [-1, 1]:
                phase = sigma * output_omega - tau * input_phi_omega - rho * input_psi_omega
                symbol = 1 / (output_omega * phase)
                total += bilinear(fg_phi[tau], h_psi[rho], symbol)
                total += bilinear(h_phi[tau], fg_psi[rho], symbol)
        cascade[sigma] = total

    coefficient = -3
    position = -coefficient / 2 * kink * (cascade[1] - cascade[-1])
    momentum_preimage = -coefficient / (2j) * kink * (cascade[1] + cascade[-1])
    momentum = np.fft.ifft(omega_psi * np.fft.fft(momentum_preimage))
    omega_r0 = np.fft.ifft(omega_psi * np.fft.fft(r0)).real
    gamma1 = dx * np.dot(r0, position) / np.sqrt(2 * np.pi)
    gamma2_cascade = dx * np.dot(omega_r0, momentum_preimage) / np.sqrt(2 * np.pi)
    gamma2_direct_pair = dx * np.dot(r0, momentum) / np.sqrt(2 * np.pi)
    slope = -dx * np.dot(r0, u ** 3) / np.sqrt(2 * np.pi)
    gamma2_values = []
    for g4 in [0.0, 1.0, 2.0]:
        direct = -6 * v * v * u - g4 * u ** 3
        gamma2 = dx * np.dot(r0, direct) / np.sqrt(2 * np.pi) + gamma2_cascade
        gamma2_values.append(float(gamma2.real))
    reality_error = float(max(np.max(np.abs(position.imag)), np.max(np.abs(momentum.imag))))
    check(f'physical reality {count}', reality_error < 1e-10)
    check(f'momentum duality {count}', abs(gamma2_direct_pair - gamma2_cascade) < 1e-11)
    check(f'fixed-snapshot g4 derivative {count}', abs(gamma2_values[2] - gamma2_values[1] - slope) < 1e-12)
    check(f'nonzero g4 slope {count}', abs(slope) > 1e-6)
    check(f'continuous projections {count}', max(abs(dx * np.dot(f, u)), abs(dx * np.dot(f, pu)), abs(dx * np.dot(g, v)), abs(dx * np.dot(g, pv))) < 1e-12)
    return {'grid_points': count, 'half_domain': half_domain, 'gamma1': float(gamma1.real), 'gamma2_cascade': float(gamma2_cascade.real), 'gamma2_at_g4_zero_one_two': gamma2_values, 'fixed_snapshot_g4_slope': float(slope), 'reality_error': reality_error}


results = [snapshot(256), snapshot(512), snapshot(768, 30.0)]
for index in [1, 2]:
    reference = results[1] if index == 2 else results[0]
    difference = max(abs(results[index]['gamma1'] - reference['gamma1']), abs(results[index]['gamma2_cascade'] - reference['gamma2_cascade']))
    check(f'grid or domain refinement {index}', difference < 1e-7)
for source in json.loads((HERE / 'INPUTS.json').read_text())['files']:
    raw = (HERE.parent / source['path']).read_bytes()
    check('source hash ' + source['path'], hashlib.sha256(raw).hexdigest() == source['sha256'])
output = {'scope': 'Full signed cascade on fixed comparison snapshots, with grid/domain refinement and source provenance', 'numpy_version': np.__version__, 'results': results, 'passed': len(checks), 'failed': 0, 'checks': checks}
(HERE / 'verification.json').write_text(json.dumps(output, indent=2) + '\n')
print(json.dumps({'results': results, 'passed': len(checks), 'failed': 0}, indent=2))
