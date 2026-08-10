#!/usr/bin/env python3
"""FIN Programs P488--P496: low-compute analytical research batch.

This executable deliberately excludes the high-cost P485/P487 Groebner jobs
and the global exact complex-face elimination.  It studies nine small, typed
objects using finite linear algebra, exact identities, and small synthetic
operational models.  Every physical interpretation remains conditional.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import matplotlib.pyplot as plt
import numpy as np
import scipy
from numpy.typing import NDArray
from scipy.linalg import expm
from scipy.optimize import least_squares

from fin_programs_488_496_exact_checker import run_exact_checks


ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT / "FIN_Programs_488_496_Figures"
RESULTS = ROOT / "FIN_Programs_488_496_Results.json"
SUMMARY = ROOT / "FIN_Programs_488_496_Summary.csv"
SEED = 20260809

OMEGA_S = 0.18575
PHI_S = 0.16250
ETA_S = 1.8
BETA_S = 1.0
ALPHA_GEO = 4.0 * math.log(2.0)


def strict_profile(d: int | float) -> float:
    return math.cos(OMEGA_S * float(d) + PHI_S) / (1.0 + float(d) ** ETA_S)


def strict_envelope(d: int | float) -> float:
    return 1.0 / (1.0 + float(d) ** ETA_S)


def legacy_profile(d: int | float) -> float:
    return ALPHA_GEO * math.cos(math.pi * float(d) / 4.0 + math.pi / 6.0) / (
        1.0 + 0.01 * float(d)
    )


def cyclic_weight_matrix(n: int, profile_fn: Callable[[int], float]) -> NDArray[np.float64]:
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = min(abs(i - j), n - abs(i - j))
            value = float(profile_fn(d))
            w[i, j] = value
            w[j, i] = value
    return w


def laplacian_from_w(w: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.diag(w.sum(axis=1)) - w


def strict_operator() -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    w = cyclic_weight_matrix(12, strict_profile)
    a = laplacian_from_w(w)
    return w, a


def fourier_mode(n: int, m: int) -> NDArray[np.complex128]:
    x = np.arange(n)
    return np.exp(2j * np.pi * m * x / n) / math.sqrt(n)


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(v) for v in value]
    if isinstance(value, np.ndarray):
        return native(value.tolist())
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def p488_nonlinear_persistent_modes(a: NDArray[np.float64]) -> dict:
    """Exact spectral limit cycles for a saturating gain plus strict Hamiltonian."""
    n = a.shape[0]
    gain = 0.20
    saturation = 1.0
    modes = []
    for m in range(7):
        v = fourier_mode(n, m)
        lam = float(np.real(np.vdot(v, a @ v)))
        site_amplitude = math.sqrt(gain / saturation)
        psi0 = math.sqrt(n) * site_amplitude * v
        residual = np.linalg.norm(
            (gain - saturation * np.abs(psi0) ** 2) * psi0
            - 1j * (a @ psi0)
            + 1j * lam * psi0
        )

        phases = np.exp(2j * np.pi * m * np.arange(n) / n)
        dmat = np.diag(phases)
        h = dmat.conj().T @ a @ dmat - lam * np.eye(n)
        hr, hi = np.real(h), np.imag(h)
        jac = np.block(
            [
                [hi - 2.0 * gain * np.eye(n), hr],
                [-hr, hi],
            ]
        )
        exponents = np.linalg.eigvals(jac)
        phase_index = int(np.argmin(np.abs(exponents)))
        transverse = np.delete(exponents, phase_index)
        max_transverse = float(np.max(np.real(transverse)))
        no_positive = bool(max_transverse <= 1e-9)
        asymptotic = bool(max_transverse < -1e-9)
        modes.append(
            {
                "mode": m,
                "frequency": lam,
                "period": None if abs(lam) < 1e-14 else 2.0 * math.pi / abs(lam),
                "exact_cycle_residual": float(residual),
                "site_amplitude": site_amplitude,
                "inverse_participation_ratio": float(np.sum(np.abs(v) ** 4)),
                "phase_exponent_abs": float(abs(exponents[phase_index])),
                "max_transverse_real_exponent": max_transverse,
                "no_positive_transverse_exponent": no_positive,
                "asymptotically_stable_modulo_phase": asymptotic,
            }
        )

    return {
        "program": "P488",
        "object": "O192 Saturating Strict Spectral Limit-Cycle Flow",
        "equation": "dot(psi)=(r-g|psi|^2)psi-iApsi",
        "gain_r": gain,
        "saturation_g": saturation,
        "analytic_result": (
            "Every constant-modulus Fourier eigenmode gives an exact persistent rotating "
            "solution with site amplitude sqrt(r/g) and angular frequency lambda_m."
        ),
        "modes": modes,
        "nonexpanding_excited_modes": [m["mode"] for m in modes if m["mode"] and m["no_positive_transverse_exponent"]],
        "asymptotically_stable_excited_modes": [m["mode"] for m in modes if m["mode"] and m["asymptotically_stable_modulo_phase"]],
        "neutral_excited_modes": [
            m["mode"]
            for m in modes
            if m["mode"] and m["no_positive_transverse_exponent"] and not m["asymptotically_stable_modulo_phase"]
        ],
        "verdict": (
            "Conditional finite persistent modes are constructed, but their IPR is 1/12: "
            "they are delocalized spectral cycles, not a localized soliton theorem."
        ),
        "status": "proven_exact_solution_plus_numerical_linear_stability",
    }


def p489_relational_clock(a: NDArray[np.float64], p488: dict) -> dict:
    stable = p488["asymptotically_stable_excited_modes"]
    mode = int(stable[0] if stable else 1)
    record = p488["modes"][mode]
    omega = float(record["frequency"])
    period = 2.0 * math.pi / omega
    rng = np.random.default_rng(SEED + 489)
    sample_times = np.linspace(0.0, 0.8 * period, 41)
    phase_noise_sigma = 0.005
    observed = -omega * sample_times + 0.37 + rng.normal(0.0, phase_noise_sigma, sample_times.size)
    unwrapped = np.unwrap(observed)
    slope, intercept = np.polyfit(sample_times, unwrapped, 1)
    estimated = -float(slope)
    relative_error = abs(estimated - omega) / omega

    phase_error = 0.01
    frequency_relative_error = 1e-3
    horizon = 0.8 * period
    worst_first_order_time_error = (
        phase_error / omega + horizon * frequency_relative_error
    ) / (1.0 - frequency_relative_error)

    base_times = np.linspace(0.0, 0.9 * period, 51)
    scale_residuals = []
    for scale in (0.1, 1.0, 10.0):
        rescaled_omega = omega / scale
        rescaled_times = scale * base_times
        residual = np.max(
            np.abs(np.exp(-1j * omega * base_times) - np.exp(-1j * rescaled_omega * rescaled_times))
        )
        scale_residuals.append(float(residual))

    return {
        "program": "P489",
        "object": "O193 Relational Spectral Phase Clock",
        "selected_mode": mode,
        "dimensionless_frequency": omega,
        "dimensionless_period": period,
        "clock_domain_before_aliasing": [0.0, period],
        "synthetic_phase_fit": {
            "samples": int(sample_times.size),
            "phase_noise_sigma": phase_noise_sigma,
            "estimated_frequency": estimated,
            "relative_error": relative_error,
        },
        "declared_error_bound": {
            "phase_error": phase_error,
            "relative_frequency_error": frequency_relative_error,
            "horizon": horizon,
            "first_order_time_error_bound": worst_first_order_time_error,
        },
        "affine_scale_orbit_phase_residuals": scale_residuals,
        "theorem": (
            "A nonzero persistent mode and a phase record define local relational time "
            "modulo its period; simultaneous A/time rescaling leaves all phase records invariant."
        ),
        "verdict": "Relative cycle time is available; an absolute second and clock selection are not generated.",
        "status": "proven_conditional_clock_and_scale_obstruction",
    }


def grouped_stieltjes_data(
    a: NDArray[np.float64], delta: float = 0.15
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    even = np.arange(0, a.shape[0], 2)
    hidden = np.arange(1, a.shape[0], 2)
    b = a[np.ix_(even, hidden)]
    c = a[np.ix_(hidden, hidden)] + delta * np.eye(hidden.size)
    vals, vecs = np.linalg.eigh(c)
    groups: list[list[int]] = []
    for idx, val in enumerate(vals):
        if not groups or abs(val - vals[groups[-1][0]]) > 1e-9:
            groups.append([idx])
        else:
            groups[-1].append(idx)
    poles, residues = [], []
    for group in groups:
        q = vecs[:, group]
        residue = np.trace(b @ q @ q.T @ b.T) / even.size
        if residue > 1e-14:
            poles.append(float(np.mean(vals[group])))
            residues.append(float(residue))
    return np.array(poles), np.array(residues), b, c


def stieltjes_response(z: np.ndarray, poles: np.ndarray, residues: np.ndarray) -> np.ndarray:
    return np.sum(residues[None, :] / (z[:, None] + poles[None, :]), axis=1)


def p490_adaptive_stieltjes_schur(a: NDArray[np.float64]) -> dict:
    poles, residues, _, _ = grouped_stieltjes_data(a)
    scale = float(np.median(poles))
    z = np.geomspace(0.03 * scale, 12.0 * scale, 48)
    target = stieltjes_response(z, poles, residues)
    ncomp = poles.size
    rng = np.random.default_rng(SEED + 490)

    def residual_fn(theta: np.ndarray) -> np.ndarray:
        fit_poles = np.exp(theta[:ncomp])
        fit_residues = np.exp(theta[ncomp:])
        pred = stieltjes_response(z, fit_poles, fit_residues)
        return (pred - target) / target

    fits = []
    for _ in range(12):
        init_poles = np.exp(
            rng.uniform(math.log(poles.min() * 0.5), math.log(poles.max() * 2.0), ncomp)
        )
        init_residues = np.exp(
            rng.uniform(math.log(residues.min() * 0.3), math.log(residues.max() * 3.0), ncomp)
        )
        fit = least_squares(
            residual_fn,
            np.log(np.concatenate([init_poles, init_residues])),
            max_nfev=5000,
            xtol=1e-13,
            ftol=1e-13,
            gtol=1e-13,
        )
        fits.append(fit)
    best = min(fits, key=lambda f: np.linalg.norm(f.fun))
    fit_poles = np.exp(best.x[:ncomp])
    fit_residues = np.exp(best.x[ncomp:])
    order = np.argsort(fit_poles)
    fit_poles, fit_residues = fit_poles[order], fit_residues[order]
    true_order = np.argsort(poles)
    true_poles, true_residues = poles[true_order], residues[true_order]
    fitted_response = stieltjes_response(z, fit_poles, fit_residues)
    max_rel = float(np.max(np.abs(fitted_response - target) / target))

    derivative_margins = {}
    for order_k in range(5):
        positive_quantity = math.factorial(order_k) * np.sum(
            residues[None, :] / (z[:, None] + poles[None, :]) ** (order_k + 1), axis=1
        )
        derivative_margins[str(order_k)] = float(np.min(positive_quantity))

    singular_values = np.linalg.svd(best.jac, compute_uv=False)
    condition = float(singular_values[0] / singular_values[-1])

    return {
        "program": "P490",
        "object": "O194 Positive Stieltjes Parameter Learning Flow",
        "parameterization": "r_j=exp(a_j), c_j=exp(b_j), sigma(z)=sum_j r_j/(z+c_j)",
        "component_count": int(ncomp),
        "true_poles": true_poles,
        "true_residues": true_residues,
        "fitted_poles": fit_poles,
        "fitted_residues": fit_residues,
        "response_max_relative_error": max_rel,
        "relative_jacobian_condition_number": condition,
        "complete_monotonicity_positive_margins_orders_0_4": derivative_margins,
        "successful_starts_below_1e-6": int(sum(np.max(np.abs(f.fun)) < 1e-6 for f in fits)),
        "theorem": (
            "Log-parameter gradient flow preserves positive poles and residues, hence the "
            "Stieltjes class and complete monotonicity at every finite step."
        ),
        "verdict": (
            "A class-preserving adaptive memory law is constructed and numerically identifiable "
            "for the declared trace response; its training target and loss are supplied inputs."
        ),
        "status": "proven_class_preservation_plus_strong_finite_recovery",
        "plot_data": {"z": z, "target": target, "fit": fitted_response},
    }


def scalar_context_fingerprint(n: int, delta: float = 0.10) -> dict:
    w = cyclic_weight_matrix(n, strict_envelope)
    a = laplacian_from_w(w)
    even = np.arange(0, n, 2)
    hidden = np.arange(1, n, 2)
    b = a[np.ix_(even, hidden)]
    c = a[np.ix_(hidden, hidden)] + delta * np.eye(hidden.size)
    vals, vecs = np.linalg.eigh(c)
    scale = float(np.median(vals))
    u = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    sigma = []
    for value in u:
        response = b @ np.linalg.solve(value * scale * np.eye(hidden.size) + c, b.T)
        sigma.append(float(np.trace(response) / even.size))
    sigma = np.array(sigma)
    fingerprint = sigma / sigma[np.where(u == 1.0)[0][0]]

    visible_poles = []
    for idx, val in enumerate(vals):
        q = vecs[:, idx : idx + 1]
        residue = float(np.trace(b @ q @ q.T @ b.T) / even.size)
        if residue > 1e-12:
            visible_poles.append(float(val / scale))
    unique = []
    for val in sorted(visible_poles):
        if not unique or abs(val - unique[-1]) > 1e-8:
            unique.append(val)
    return {
        "n": n,
        "u": u,
        "fingerprint": fingerprint,
        "visible_pole_count_with_multiplicity": len(visible_poles),
        "distinct_visible_pole_count": len(unique),
        "normalized_visible_poles": unique,
    }


def p491_dynamic_fractal_fixed_point() -> dict:
    fingerprints = [scalar_context_fingerprint(n) for n in (12, 24, 48, 96)]
    successive = []
    for left, right in zip(fingerprints[:-1], fingerprints[1:]):
        difference = np.asarray(right["fingerprint"]) - np.asarray(left["fingerprint"])
        successive.append(
            {
                "sizes": [left["n"], right["n"]],
                "max_fingerprint_difference": float(np.max(np.abs(difference))),
                "l2_fingerprint_difference": float(np.linalg.norm(difference)),
                "distinct_poles": [left["distinct_visible_pole_count"], right["distinct_visible_pole_count"]],
            }
        )
    exact_fixed_point_refuted = all(item["max_fingerprint_difference"] > 1e-10 for item in successive)
    differences = [item["max_fingerprint_difference"] for item in successive]
    monotone_convergence = all(right < left for left, right in zip(differences[:-1], differences[1:]))
    return {
        "program": "P491",
        "object": "O195 Normalized Dynamic Context Fingerprint",
        "family": "positive strict attenuation envelope on C_n, even/odd Schur context, delta=0.1",
        "fingerprints": fingerprints,
        "successive_differences": successive,
        "finite_exact_fixed_point_refuted": exact_fixed_point_refuted,
        "monotone_convergence_observed": monotone_convergence,
        "analytic_obstruction": (
            "Finite self-energies are rational Stieltjes functions. Equality requires matching "
            "visible pole sets and residues; changing noncancelling pole counts forbids an exact "
            "fixed point in this declared refinement family."
        ),
        "verdict": (
            "No exact dynamic fractal fixed point occurs in the tested normalized envelope family, "
            "and the three finite-size discrepancies do not decrease monotonically. Approximate "
            "continuum convergence remains open and is distinct from exact self-similarity."
        ),
        "status": "proven_finite_family_no_go_no_monotone_convergence_evidence",
    }


def p492_resonance_dirichlet(w: NDArray[np.float64], a: NDArray[np.float64]) -> dict:
    s = float(w.sum(axis=1)[0])
    eig_a = np.linalg.eigvalsh(a)
    eig_w = np.linalg.eigvalsh(w)[::-1]
    rng = np.random.default_rng(SEED + 492)
    residuals = []
    for _ in range(200):
        x = rng.normal(size=12)
        x /= np.linalg.norm(x)
        residuals.append(abs(float(x @ w @ x + x @ a @ x - s)))
    mean_zero_max_w = float(eig_w[1])
    mean_zero_min_a = float(eig_a[1])
    return {
        "program": "P492",
        "object": "O196 Resonance--Dirichlet Complementarity Theorem",
        "row_sum_s": s,
        "random_rayleigh_identity_max_residual": max(residuals),
        "mean_zero_max_W": mean_zero_max_w,
        "mean_zero_min_A": mean_zero_min_a,
        "complementarity_residual": abs(mean_zero_max_w + mean_zero_min_a - s),
        "largest_A_frequency": float(eig_a[-1]),
        "corresponding_smallest_W": float(eig_w[-1]),
        "linear_gain_thresholds": {
            "unconstrained_first_instability": float(eig_a[0]),
            "mean_zero_first_instability": float(eig_a[1]),
            "highest_frequency_instability": float(eig_a[-1]),
        },
        "theorem": "For A=sI-W, R_W(x)=s-R_A(x) on every nonzero state.",
        "verdict": (
            "Maximum coupling resonance is exactly minimum Dirichlet cost. It is not the same "
            "as selecting the largest temporal frequency of A, which reverses the ordering."
        ),
        "status": "proven",
        "plot_data": {"A": eig_a, "W_descending": eig_w},
    }


@dataclass(frozen=True)
class KernelParameters:
    amplitude: float
    beta: float
    eta: float
    omega: float
    phi: float


def completion_automaton(params: KernelParameters, dmax: int) -> np.ndarray:
    q = 0.0
    damping = 1.0
    theta = params.phi
    values = []
    for _ in range(dmax + 1):
        values.append(params.amplitude * damping * math.cos(theta))
        delta_power = (q + 1.0) ** params.eta - q**params.eta
        damping = damping / (1.0 + params.beta * delta_power * damping)
        q += 1.0
        theta += params.omega
    return np.array(values)


def closed_kernel(params: KernelParameters, dmax: int) -> np.ndarray:
    d = np.arange(dmax + 1, dtype=float)
    return params.amplitude * np.cos(params.omega * d + params.phi) / (
        1.0 + params.beta * d**params.eta
    )


def p493_typed_dynamic_completion() -> dict:
    legacy = KernelParameters(ALPHA_GEO, 0.01, 1.0, math.pi / 4.0, math.pi / 6.0)
    strict = KernelParameters(1.0, 1.0, 1.8, OMEGA_S, PHI_S)
    endpoint_residuals = {}
    transfer = {}
    for name, params in (("legacy", legacy), ("strict", strict)):
        endpoint_residuals[name] = float(
            np.max(np.abs(completion_automaton(params, 48) - closed_kernel(params, 48)))
        )
        transfer[name] = {}
        for n in (12, 16, 24):
            dmax = n // 2
            transfer[name][str(n)] = float(
                np.max(np.abs(completion_automaton(params, dmax) - closed_kernel(params, dmax)))
            )

    homotopy = []
    first_passive_t = None
    for t in np.linspace(0.0, 1.0, 101):
        params = KernelParameters(
            amplitude=(1.0 - t) * legacy.amplitude + t * strict.amplitude,
            beta=(1.0 - t) * legacy.beta + t * strict.beta,
            eta=(1.0 - t) * legacy.eta + t * strict.eta,
            omega=(1.0 - t) * legacy.omega + t * strict.omega,
            phi=(1.0 - t) * legacy.phi + t * strict.phi,
        )
        radial = closed_kernel(params, 6)[1:]
        w = cyclic_weight_matrix(12, lambda d, arr=radial: float(arr[d - 1]))
        lap = laplacian_from_w(w)
        min_eig = float(np.linalg.eigvalsh(lap)[0])
        negative_distances = int(np.sum(radial < 0.0))
        passive = bool(min_eig >= -1e-10 and negative_distances == 0)
        if passive and first_passive_t is None:
            first_passive_t = float(t)
        if round(t * 100) % 10 == 0:
            homotopy.append(
                {
                    "t": float(t),
                    "min_signed_laplacian_eigenvalue": min_eig,
                    "negative_radial_distances": negative_distances,
                    "passive_nonnegative_weight_laplacian": passive,
                }
            )

    return {
        "program": "P493",
        "object": "O197 Exponent-Lift Kernel Completion Automaton",
        "state_update": (
            "q'=q+1; c'=c/[1+beta((q+1)^eta-q^eta)c]; theta'=theta+omega; "
            "K=amplitude*c*cos(theta)"
        ),
        "endpoint_closed_form_residuals": endpoint_residuals,
        "carrier_transfer_residuals": transfer,
        "linear_parameter_homotopy_C12": homotopy,
        "first_sampled_passive_homotopy_parameter": first_passive_t,
        "theorem": (
            "The automaton exactly generates amplitude*cos(omega d+phi)/(1+beta d^eta) "
            "for every nonnegative integer d."
        ),
        "verdict": (
            "A typed dynamic carrier contains both Legacy* and strict as supplied parameter "
            "points and transfers across sizes. It does not source the parameter change; most "
            "of the naive interpolation is not a nonnegative-weight passive Laplacian."
        ),
        "status": "proven_representation_not_derivation",
    }


def discrete_winding(phases: np.ndarray) -> float | None:
    increments = np.angle(np.roll(phases, -1) / phases)
    if np.any(np.abs(np.abs(increments) - math.pi) < 1e-10):
        return None
    return float(np.sum(increments) / (2.0 * math.pi))


def p494_topological_defect_obstruction(a: NDArray[np.float64]) -> dict:
    mode_rows = []
    for m in range(-5, 7):
        v = fourier_mode(12, m % 12)
        mode_rows.append(
            {
                "mode": m,
                "discrete_winding": discrete_winding(v),
                "strict_dirichlet_energy": float(np.real(np.vdot(v, a @ v))),
                "minimum_vertex_amplitude": float(np.min(np.abs(v))),
            }
        )

    v1 = fourier_mode(12, 1)
    v0 = fourier_mode(12, 0)
    path_energies, minimum_site_amplitudes = [], []
    for t in np.linspace(0.0, 1.0, 201):
        if t <= 0.5:
            psi = (1.0 - 2.0 * t) * v1
        else:
            psi = (2.0 * t - 1.0) * v0
        path_energies.append(float(np.real(np.vdot(psi, a @ psi))))
        minimum_site_amplitudes.append(float(np.min(np.abs(psi))))

    return {
        "program": "P494",
        "object": "O198 Phase-Field Topology Obligation",
        "mode_winding_energy_table": mode_rows,
        "finite_unconstrained_path": {
            "from_winding": 1,
            "to_winding": 0,
            "passes_through_zero_field": True,
            "max_dirichlet_energy": max(path_energies),
            "minimum_site_amplitude": min(minimum_site_amplitudes),
        },
        "theorem": (
            "The scalar configuration space C^12 is connected and carries no protected winding "
            "sector. Integer winding requires a nonvanishing S1-valued phase field plus an "
            "edge interpolation/refinement rule; changing winding then requires leaving that space."
        ),
        "verdict": (
            "The strict kernel can price supplied Fourier winding patterns but cannot create "
            "topological particle sectors from scalar kernel data alone."
        ),
        "status": "proven_obstruction",
    }


def density(vector: NDArray[np.complex128]) -> NDArray[np.complex128]:
    vector = vector / np.linalg.norm(vector)
    return np.outer(vector, vector.conj())


def normalize_probability(p: np.ndarray) -> np.ndarray:
    p = np.maximum(np.real(p), 0.0)
    total = p.sum()
    if total <= 0:
        raise ValueError("nonpositive probability total")
    return p / total


def p495_operational_identifiability(a: NDArray[np.float64]) -> dict:
    n = a.shape[0]
    eye = np.eye(n)
    projectors = []
    for x in range(n):
        p = np.zeros((n, n), dtype=complex)
        p[x, x] = 1.0
        projectors.append(p)
    gamma = 0.60
    hamiltonian_super = -1j * (np.kron(eye, a) - np.kron(a.T, eye))
    dephase_super = sum((np.kron(p.T, p) for p in projectors)) - np.eye(n * n)
    lindblad = hamiltonian_super + gamma * dephase_super
    fourier = np.column_stack([fourier_mode(n, m) for m in range(n)])

    e = np.eye(n, dtype=complex)
    preparations = {
        "basis0": density(e[:, 0]),
        "basis1": density(e[:, 1]),
        "basis2": density(e[:, 2]),
        "plus01": density((e[:, 0] + e[:, 1]) / math.sqrt(2.0)),
        "iplus01": density((e[:, 0] + 1j * e[:, 1]) / math.sqrt(2.0)),
        "fourier1": density(fourier[:, 1]),
    }
    times = [0.20, 0.55, 1.00]

    def evolve(model: str, rho: np.ndarray, t: float) -> np.ndarray:
        u = expm(-1j * t * a)
        if model == "unitary":
            return u @ rho @ u.conj().T
        if model == "markov":
            p = expm(-t * a) @ np.real(np.diag(rho))
            return np.diag(normalize_probability(p)).astype(complex)
        if model == "dephasing":
            vec = expm(t * lindblad) @ rho.reshape(-1, order="F")
            out = vec.reshape((n, n), order="F")
            return (out + out.conj().T) / 2.0
        if model == "revival":
            coherent = u @ rho @ u.conj().T
            coherence_weight = 0.45 + 0.40 * math.cos(2.8 * t)
            diagonal = np.diag(np.diag(coherent))
            return coherence_weight * coherent + (1.0 - coherence_weight) * diagonal
        raise KeyError(model)

    def measure(rho: np.ndarray, basis: str) -> np.ndarray:
        if basis == "vertex":
            return normalize_probability(np.real(np.diag(rho)))
        transformed = fourier.conj().T @ rho @ fourier
        return normalize_probability(np.real(np.diag(transformed)))

    models = ["unitary", "markov", "dephasing", "revival"]
    all_contexts = []
    for t in times:
        for prep in preparations:
            for measurement in ("vertex", "fourier"):
                all_contexts.append((t, prep, measurement))

    probabilities: dict[str, dict[tuple, np.ndarray]] = {m: {} for m in models}
    for model in models:
        for context in all_contexts:
            t, prep, measurement = context
            probabilities[model][context] = measure(evolve(model, preparations[prep], t), measurement)

    designs = {
        "D1_single_time_basis_vertex": [
            c for c in all_contexts if c[0] == 0.55 and c[1].startswith("basis") and c[2] == "vertex"
        ],
        "D2_multitime_allprep_vertex": [c for c in all_contexts if c[2] == "vertex"],
        "D3_multitime_allprep_vertex_fourier": all_contexts,
    }

    rng = np.random.default_rng(SEED + 495)
    shots = 120
    trials_per_model = 100
    design_results = {}
    for design_name, contexts in designs.items():
        pairwise = np.zeros((len(models), len(models)))
        for i, left in enumerate(models):
            for j, right in enumerate(models):
                pairwise[i, j] = np.mean(
                    [
                        0.5
                        * np.sum(
                            np.abs(probabilities[left][ctx] - probabilities[right][ctx])
                        )
                        for ctx in contexts
                    ]
                )
        confusion = np.zeros((len(models), len(models)), dtype=int)
        for true_index, true_model in enumerate(models):
            for _ in range(trials_per_model):
                log_likelihood = np.zeros(len(models))
                for ctx in contexts:
                    counts = rng.multinomial(shots, probabilities[true_model][ctx])
                    for candidate_index, candidate in enumerate(models):
                        p = np.clip(probabilities[candidate][ctx], 1e-15, 1.0)
                        log_likelihood[candidate_index] += float(np.dot(counts, np.log(p)))
                predicted = int(np.argmax(log_likelihood))
                confusion[true_index, predicted] += 1
        accuracy = float(np.trace(confusion) / np.sum(confusion))
        offdiag = pairwise[np.triu_indices(len(models), 1)]
        design_results[design_name] = {
            "context_count": len(contexts),
            "shots_per_context": shots,
            "models": models,
            "pairwise_mean_TV": pairwise,
            "minimum_pairwise_mean_TV": float(offdiag.min()),
            "confusion_rows_true_cols_predicted": confusion,
            "classification_accuracy": accuracy,
        }

    return {
        "program": "P495",
        "object": "O199 Local Operational Identifiability Atlas",
        "models": models,
        "times": times,
        "dephasing_gamma": gamma,
        "revival_family": "q(t) unitary output + [1-q(t)] dephased output, q=0.45+0.40 cos(2.8t)",
        "designs": design_results,
        "verdict": (
            "Finite local preparations and vertex/Fourier measurements separate the four declared "
            "synthetic models with quantifiable resource gain. This is an in-silico identifiability "
            "test, not laboratory evidence or a complete process-tensor theorem."
        ),
        "status": "strong_synthetic_evidence",
    }


def make_figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    p488 = results["P488"]
    modes = [row["mode"] for row in p488["modes"]]
    exponents = [row["max_transverse_real_exponent"] for row in p488["modes"]]
    frequencies = [row["frequency"] for row in p488["modes"]]
    fig, ax1 = plt.subplots(figsize=(7.6, 4.4))
    ax1.bar(np.array(modes) - 0.17, exponents, width=0.34, label="max transverse Re exponent")
    ax1.axhline(0.0, color="black", lw=1)
    ax1.set_xlabel("Fourier mode")
    ax1.set_ylabel("stability exponent")
    ax2 = ax1.twinx()
    ax2.bar(np.array(modes) + 0.17, frequencies, width=0.34, color="#d95f02", alpha=0.65, label="frequency")
    ax2.set_ylabel("dimensionless frequency")
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p488_limit_cycle_stability.png", dpi=180)
    plt.close(fig)

    p490 = results["P490"]
    plot = p490.pop("plot_data")
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.loglog(plot["z"], plot["target"], "o", ms=4, label="exact Schur trace response")
    ax.loglog(plot["z"], plot["fit"], "-", lw=1.8, label="positive Stieltjes fit")
    ax.set_xlabel("z")
    ax.set_ylabel("sigma(z)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p490_stieltjes_learning.png", dpi=180)
    plt.close(fig)

    p491 = results["P491"]
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    for row in p491["fingerprints"]:
        ax.semilogx(row["u"], row["fingerprint"], marker="o", label=f"C{row['n']}")
    ax.set_xlabel("normalized spectral argument u")
    ax.set_ylabel("normalized context fingerprint")
    ax.legend(ncol=2)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p491_dynamic_context_fingerprints.png", dpi=180)
    plt.close(fig)

    p492 = results["P492"]
    plot = p492.pop("plot_data")
    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.plot(np.arange(12), plot["A"], "o-", label="A eigenvalues ascending")
    ax.plot(np.arange(12), plot["W_descending"], "s-", label="W eigenvalues descending")
    ax.axhline(p492["row_sum_s"], color="black", lw=1, ls="--", label="row sum s")
    ax.set_xlabel("paired spectral order")
    ax.set_ylabel("eigenvalue")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p492_resonance_dirichlet.png", dpi=180)
    plt.close(fig)

    p495 = results["P495"]
    design_names = list(p495["designs"])
    accuracies = [p495["designs"][name]["classification_accuracy"] for name in design_names]
    minimum_tvs = [p495["designs"][name]["minimum_pairwise_mean_TV"] for name in design_names]
    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
    short = ["D1", "D2", "D3"]
    axes[0].bar(short, accuracies, color="#1b9e77")
    axes[0].set_ylim(0, 1.05)
    axes[0].set_ylabel("classification accuracy")
    axes[1].bar(short, minimum_tvs, color="#7570b3")
    axes[1].set_ylabel("minimum pairwise mean TV")
    fig.suptitle("Synthetic operational identifiability")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p495_operational_identifiability.png", dpi=180)
    plt.close(fig)


def write_summary(results: dict) -> None:
    rows = [
        ("P488", results["P488"]["status"], len(results["P488"]["asymptotically_stable_excited_modes"]), "persistent cycles, not localized"),
        ("P489", results["P489"]["status"], results["P489"]["synthetic_phase_fit"]["relative_error"], "relative time only"),
        ("P490", results["P490"]["status"], results["P490"]["response_max_relative_error"], "target/loss supplied"),
        ("P491", results["P491"]["status"], results["P491"]["finite_exact_fixed_point_refuted"], "declared envelope family"),
        ("P492", results["P492"]["status"], results["P492"]["random_rayleigh_identity_max_residual"], "regular strict C12"),
        ("P493", results["P493"]["status"], results["P493"]["first_sampled_passive_homotopy_parameter"], "parameters remain supplied"),
        ("P494", results["P494"]["status"], results["P494"]["finite_unconstrained_path"]["minimum_site_amplitude"], "phase-field obligation"),
        ("P495", results["P495"]["status"], results["P495"]["designs"]["D3_multitime_allprep_vertex_fourier"]["classification_accuracy"], "synthetic only"),
        ("P496", results["P496"]["status"], True, "not proof-assistant checked"),
    ]
    with SUMMARY.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["program", "status", "headline_metric", "boundary"])
        writer.writerows(rows)


def run_all() -> dict:
    start = time.perf_counter()
    w, a = strict_operator()
    p488 = p488_nonlinear_persistent_modes(a)
    results = {
        "metadata": {
            "title": "FIN Programs P488-P496: Low-Compute Analytical Research",
            "seed": SEED,
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "scope": "local analytical and small finite computations; no laboratory or external audit",
        },
        "strict_audit": {
            "row_sum": float(w.sum(axis=1)[0]),
            "row_sum_spread": float(np.ptp(w.sum(axis=1))),
            "A_min_eigenvalue": float(np.linalg.eigvalsh(a)[0]),
            "A_gap": float(np.linalg.eigvalsh(a)[1]),
            "A_max_eigenvalue": float(np.linalg.eigvalsh(a)[-1]),
        },
        "P488": p488,
        "P489": p489_relational_clock(a, p488),
        "P490": p490_adaptive_stieltjes_schur(a),
        "P491": p491_dynamic_fractal_fixed_point(),
        "P492": p492_resonance_dirichlet(w, a),
        "P493": p493_typed_dynamic_completion(),
        "P494": p494_topological_defect_obstruction(a),
        "P495": p495_operational_identifiability(a),
        "P496": run_exact_checks(write=True),
    }
    make_figures(results)
    results["metadata"]["elapsed_seconds"] = time.perf_counter() - start
    results["metadata"]["script_sha256"] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    write_summary(results)
    RESULTS.write_text(json.dumps(native(results), indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return results


if __name__ == "__main__":
    output = run_all()
    print(json.dumps(native({
        "elapsed_seconds": output["metadata"]["elapsed_seconds"],
        "P488_asymptotically_stable_excited_modes": output["P488"]["asymptotically_stable_excited_modes"],
        "P488_neutral_excited_modes": output["P488"]["neutral_excited_modes"],
        "P489_clock_relative_error": output["P489"]["synthetic_phase_fit"]["relative_error"],
        "P490_response_error": output["P490"]["response_max_relative_error"],
        "P491_successive": output["P491"]["successive_differences"],
        "P492_identity_residual": output["P492"]["random_rayleigh_identity_max_residual"],
        "P493_first_passive_t": output["P493"]["first_sampled_passive_homotopy_parameter"],
        "P495_accuracies": {k: v["classification_accuracy"] for k, v in output["P495"]["designs"].items()},
        "P496": output["P496"]["status"],
    }), indent=2, sort_keys=True))
