#!/usr/bin/env python3
"""Independent finite checks for the FIN missing-principle monograph.

No repository interpretation is imported.  The script uses only the frozen
C12 profiles and elementary finite-dimensional linear algebra.
"""

from __future__ import annotations

import json
import math
from typing import Callable

import numpy as np


N = 12


def cyclic_distance(i: int, j: int) -> int:
    delta = abs(i - j)
    return min(delta, N - delta)


def strict_profile(d: int) -> float:
    if d == 0:
        return 0.0
    return math.cos(0.18575 * d + 0.1625) / (1.0 + d**1.8)


def legacy_profile(d: int) -> float:
    if d == 0:
        return 0.0
    return 4.0 * math.log(2.0) * math.cos(math.pi * d / 4.0 + math.pi / 6.0) / (1.0 + 0.01 * d)


def radial_matrix(profile: Callable[[int], float]) -> np.ndarray:
    return np.array([[profile(cyclic_distance(i, j)) for j in range(N)] for i in range(N)])


def cycle_laplacian() -> np.ndarray:
    shift = np.roll(np.eye(N), 1, axis=1)
    return 2.0 * np.eye(N) - shift - shift.T


def grouped_symbol(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the seven distinct C12 Laplacian nodes and radial-circulant symbol."""
    first_row = matrix[0]
    symbol = np.fft.fft(first_row).real
    lap_nodes = 2.0 - 2.0 * np.cos(2.0 * np.pi * np.arange(7) / N)
    return lap_nodes, symbol[:7]


def golden_minimize(function: Callable[[float], float], lo: float, hi: float, iterations: int = 180) -> float:
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    c = hi - (hi - lo) / phi
    d = lo + (hi - lo) / phi
    fc, fd = function(c), function(d)
    for _ in range(iterations):
        if fc < fd:
            hi, d, fd = d, c, fc
            c = hi - (hi - lo) / phi
            fc = function(c)
        else:
            lo, c, fc = c, d, fd
            d = lo + (hi - lo) / phi
            fd = function(d)
    return (lo + hi) / 2.0


def best_resolvent_fit(mu: np.ndarray, target: np.ndarray) -> dict[str, float]:
    """Fit target_j = A/(mu_j+m2)+C, including C as contact shift."""

    def solve(log_m2: float) -> tuple[float, float, float, np.ndarray]:
        m2 = math.exp(log_m2)
        design = np.column_stack((1.0 / (mu + m2), np.ones_like(mu)))
        (amplitude, contact), *_ = np.linalg.lstsq(design, target, rcond=None)
        fitted = design @ np.array([amplitude, contact])
        error = float(np.linalg.norm(target - fitted) / np.linalg.norm(target))
        return error, m2, float(amplitude), float(contact), fitted

    optimum = golden_minimize(lambda x: solve(x)[0], math.log(1e-5), math.log(1e3))
    error, m2, amplitude, contact, _ = solve(optimum)
    return {"nrmse_7mode": error, "m2": m2, "A": amplitude, "contact": contact}


def heat_dimension(mu_full: np.ndarray, times: np.ndarray) -> np.ndarray:
    weights = np.exp(-np.outer(times, mu_full))
    return 2.0 * times * ((weights @ mu_full) / weights.sum(axis=1))


def main() -> None:
    laplacian = cycle_laplacian()
    strict = radial_matrix(strict_profile)
    legacy = radial_matrix(legacy_profile)
    mu, strict_symbol = grouped_symbol(strict)
    _, legacy_symbol = grouped_symbol(legacy)

    # Exact graph-filter interpolation on the seven distinct Laplacian nodes.
    strict_coeff = np.polynomial.polynomial.polyfit(mu, strict_symbol, 6)
    legacy_coeff = np.polynomial.polynomial.polyfit(mu, legacy_symbol, 6)
    strict_poly = sum(strict_coeff[k] * np.linalg.matrix_power(laplacian, k) for k in range(7))
    legacy_poly = sum(legacy_coeff[k] * np.linalg.matrix_power(laplacian, k) for k in range(7))
    strict_poly_residual = np.linalg.norm(strict_poly - strict) / np.linalg.norm(strict)
    legacy_poly_residual = np.linalg.norm(legacy_poly - legacy) / np.linalg.norm(legacy)

    resolvent_fit = best_resolvent_fit(mu, strict_symbol)

    # A finite boundary density has no atom at the legacy shell.
    legacy_shell_energy = 2.0 - 2.0 * math.cos(math.pi / 4.0)
    shell_gap = float(np.min(np.abs(mu - legacy_shell_energy)))

    # Explicit first-order-condition witness for a=e_0, b=e_0, i=0, j=1.
    first_order_witness = float(strict[0, 1])

    # Heat spectral dimension of one finite C12 graph: zero at both ends.
    full_mu = np.linalg.eigvalsh(laplacian)
    times = np.logspace(-5, 5, 4001)
    dimensions = heat_dimension(full_mu, times)
    max_index = int(np.argmax(dimensions))

    # Modular flow on a finite commutative diagonal algebra is trivial.
    probabilities = np.arange(1.0, N + 1.0)
    probabilities /= probabilities.sum()
    rho = np.diag(probabilities)
    observable = np.diag(np.linspace(-1.0, 1.0, N))
    modular_commutator = float(np.linalg.norm(rho @ observable - observable @ rho))

    # A spectrally constructed state cannot break reflection symmetry.
    reflection = np.zeros((N, N))
    for i in range(N):
        reflection[(-i) % N, i] = 1.0
    spectral_state = np.linalg.matrix_power(np.eye(N), 1)
    evals, evecs = np.linalg.eigh(strict)
    boltzmann = np.exp(-evals - np.max(-evals))
    spectral_state = (evecs * (boltzmann / boltzmann.sum())) @ evecs.T
    reflection_residual = float(np.linalg.norm(reflection @ spectral_state @ reflection.T - spectral_state))

    # Same Shannon information, arbitrarily different physical mean energy.
    bit_probabilities = np.array([0.5, 0.5])
    shannon_nats = float(-np.sum(bit_probabilities * np.log(bit_probabilities)))
    mean_energy_1 = float(bit_probabilities @ np.array([0.0, 1.0]))
    mean_energy_100 = float(bit_probabilities @ np.array([0.0, 100.0]))

    # Symmetric bifurcation: an exact symmetric seed stays unselected; signs of
    # infinitesimal perturbations select opposite attractors.
    def double_well_flow(x0: float, dt: float = 1e-3, steps: int = 20000) -> float:
        x = x0
        for _ in range(steps):
            x += dt * (-4.0 * x * (x * x - 1.0))
        return x

    bifurcation = {
        "seed_0": double_well_flow(0.0),
        "seed_plus": double_well_flow(1e-9),
        "seed_minus": double_well_flow(-1e-9),
    }

    results = {
        "exact_filter_algebra": {
            "strict_relative_residual": float(strict_poly_residual),
            "legacy_relative_residual": float(legacy_poly_residual),
            "dimension": 7,
        },
        "strict_screened_resolvent_fit": resolvent_fit,
        "legacy_shell": {
            "energy_2_minus_sqrt2": legacy_shell_energy,
            "distance_to_finite_spectrum": shell_gap,
            "boundary_density_off_spectrum": 0.0,
        },
        "real_spectral_triple_first_order_witness": first_order_witness,
        "finite_heat_spectral_dimension": {
            "d_at_t_1e_minus_5": float(dimensions[0]),
            "maximum": float(dimensions[max_index]),
            "time_of_maximum": float(times[max_index]),
            "d_at_t_1e5": float(dimensions[-1]),
        },
        "commutative_modular_flow_commutator_norm": modular_commutator,
        "spectral_state_reflection_residual": reflection_residual,
        "same_information_different_energy": {
            "entropy_nats": shannon_nats,
            "mean_energy_scale_1": mean_energy_1,
            "mean_energy_scale_100": mean_energy_100,
        },
        "symmetric_double_well": bifurcation,
    }
    print(json.dumps(results, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
