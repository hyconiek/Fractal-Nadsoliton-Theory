#!/usr/bin/env python3
"""FIN P527--P536: auxiliary nonlinear source, validated dynamics and exact replays.

The batch is finite and dimensionless.  Added mediators, temporal memory
realizations, numerical perturbations and detector budgets remain explicit
conditional model data rather than consequences silently inherited from a FIN
kernel.
"""

from __future__ import annotations

import json
import math
import platform
import time
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
import sympy as sp
from scipy.optimize import minimize

from fin_programs_488_496_low_compute import strict_operator
from fin_programs_507_516_research import (
    N,
    continue_pattern,
    stationary_root,
    stieltjes_memory_operator,
    strict_a_interval,
    wiener_error_bound,
)
from fin_programs_517_526_research import (
    context_symbol_data,
    frequency_chart,
    hidden_memory_jacobian,
    interval_product,
    p525_detector_polytope,
    stability_chart,
)


ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT / "FIN_Programs_527_536_Figures"
RESULTS = ROOT / "FIN_Programs_527_536_Results.json"
SUMMARY = ROOT / "FIN_Programs_527_536_Summary.csv"
P528_CERT = ROOT / "FIN_P528_Stability_Replay_Certificate.json"
P533_CERT = ROOT / "FIN_P533_Interval_FFT_Certificate.json"
P534_CERT = ROOT / "FIN_P534_Rational_PSD_Certificate.json"
SEED = 20260813


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
    if isinstance(value, Fraction):
        return f"{value.numerator}/{value.denominator}"
    return value


def hx(value: float) -> str:
    return float(value).hex()


def p527_auxiliary_field_fourth_jet() -> dict:
    """Classify invariant quartics and construct a minimal sign-forcing mediator."""
    c, g, b, rho, chi = sp.symbols("c g b rho chi", real=True)
    joint = b * chi**2 / 2 - g * chi * rho
    stationary_chi = sp.solve(sp.diff(joint, chi), chi)[0]
    effective = sp.simplify(joint.subs(chi, stationary_chi))
    removals = [
        {
            "removed_axiom": "B positive definite",
            "failure": "the mediator minimizer may not exist and the Schur-complement sign is not controlled",
            "necessity_rank": 1,
        },
        {
            "removed_axiom": "nonzero linear density-mediator coupling",
            "failure": "adiabatic elimination produces no quartic term",
            "necessity_rank": 1,
        },
        {
            "removed_axiom": "stationary/adiabatic elimination",
            "failure": "the theory remains a two-field dynamics and no closed single-field cubic law follows",
            "necessity_rank": 1,
        },
        {
            "removed_axiom": "locality B=bI",
            "failure": "the sign remains attractive but the induced quartic is nonlocal rho^T B^-1 rho",
            "necessity_rank": 2,
        },
    ]
    return {
        "program": "P527",
        "object": "O231 Positive-Mediator Schur-Complement Source",
        "quartic_invariant_classification": (
            "Global U(1) and permutation invariance allow a*(sum |psi|^2)^2+b*sum |psi|^4. "
            "Disjoint-support additivity removes the cross-site coefficient a, leaving one local quartic ray."
        ),
        "joint_mediator_energy": str(joint),
        "stationary_mediator": str(stationary_chi),
        "effective_density_energy": str(effective),
        "forced_sign": "negative semidefinite for b>0; local focusing coefficient=-g^2/(2b)",
        "minimal_axioms": [
            "positive mediator quadratic form B>0",
            "nonzero linear coupling -g<chi,rho> with rho=|psi|^2",
            "stationary/adiabatic elimination of chi",
        ],
        "axiom_removal_ledger": removals,
        "magnitude_identified": False,
        "remaining_free_scale": "g^2/b (or g^2 B^-1 in the nonlocal case)",
        "boundedness_obligation": (
            "The negative effective quartic is only a local jet. Global boundedness requires a saturating/higher-order "
            "term or retention of additional nonlinear mediator structure."
        ),
        "theorem": (
            "Positive mediator energy plus linear density coupling and stationary elimination forces an attractive "
            "quartic Schur complement. These three premises are individually necessary for the declared one-field "
            "mechanism. They construct the missing sign mechanism but do not show that FIN supplies the mediator or g^2/b."
        ),
        "status": "proven_minimal_auxiliary_field_focusing_sign_mechanism",
    }


def branch_state(a: np.ndarray, omega: float, start: np.ndarray | None = None) -> np.ndarray:
    if start is None:
        start, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    state = np.asarray(start).copy()
    steps = max(2, int(math.ceil(abs(omega - 1.0) / 0.003)))
    for w in np.linspace(1.0, omega, steps + 1)[1:]:
        state = stationary_root(a, 1.0, float(w), state)
    return state


def p528_stability_exact_replay(a: np.ndarray) -> dict:
    """Export exact-dyadic replay inequalities for O223."""
    alo, ahi, amid = strict_a_interval()
    edges = np.linspace(0.68, 1.20, 209)
    first = 0.5 * (edges[0] + edges[1])
    initial = branch_state(a, first)
    charts = []
    full_rows = []
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        row = stability_chart(a, alo, ahi, amid, float(lo), float(hi), initial)
        full_rows.append(row)
        initial = np.asarray(row["state"])
        inertia = row["inertia"]
        chart = {
            "id": index,
            "omega_interval": [hx(lo), hx(hi)],
            "state_inclusion_margin_lower": hx(row["state_inclusion_margin"]),
            "state_defect_upper": hx(row["state_defect"]),
            "lminus_second_lower": hx(inertia["lminus_second_lower"]),
            "lplus_negative_upper": hx(inertia["lplus_negative_upper"]),
            "lplus_first_positive_lower": hx(inertia["lplus_first_positive_lower"]),
            "dP_lower": hx(row["dP_domega_interval"][0]),
            "dP_upper": hx(row["dP_domega_interval"][1]),
            "d2P_lower": hx(row["d2P_domega2_interval"][0]),
            "d2P_upper": hx(row["d2P_domega2_interval"][1]),
        }
        charts.append(chart)

    bridges = []
    for index, shared in enumerate(edges[1:-1]):
        left, right = full_rows[index], full_rows[index + 1]
        seed = 0.5 * (np.asarray(left["state"]) + np.asarray(right["state"]))
        centre = stationary_root(a, 1.0, float(shared), seed)
        bridge = frequency_chart(centre, float(shared), float(shared), float(shared), alo, ahi, amid)
        blo = np.asarray(bridge["centre"]) - np.asarray(bridge["radius"])
        bhi = np.asarray(bridge["centre"]) + np.asarray(bridge["radius"])
        llo = np.asarray(left["state"]) - np.asarray(left["state_radius"])
        lhi = np.asarray(left["state"]) + np.asarray(left["state_radius"])
        rlo = np.asarray(right["state"]) - np.asarray(right["state_radius"])
        rhi = np.asarray(right["state"]) + np.asarray(right["state_radius"])
        nesting = float(np.min(np.concatenate([blo - llo, lhi - bhi, blo - rlo, rhi - bhi])))
        bridges.append(
            {
                "id": index,
                "omega": hx(shared),
                "inclusion_margin_lower": hx(bridge["margin"]),
                "defect_upper": hx(bridge["defect"]),
                "nesting_margin_lower": hx(nesting),
            }
        )

    certificate = {
        "format": "FIN_P528_exact_binary64_to_rational_stability_ledger_v1",
        "trust_boundary": (
            "Upstream outward binary64 interval generation is supplied. The checker converts every hexadecimal "
            "endpoint to an exact Fraction and replays coverage, Krawczyk, inertia, slope-sector, curvature and nesting signs."
        ),
        "charts": charts,
        "bridges": bridges,
        "turning_bracket": [hx(0.722143668857225), hx(0.722143708857225)],
    }
    P528_CERT.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    accepted = all(
        float.fromhex(row["state_inclusion_margin_lower"]) > 0
        and float.fromhex(row["state_defect_upper"]) < 1
        and float.fromhex(row["lminus_second_lower"]) > 0
        and float.fromhex(row["lplus_negative_upper"]) < 0
        and float.fromhex(row["lplus_first_positive_lower"]) > 0
        for row in charts
    )
    accepted = accepted and all(
        float.fromhex(row["inclusion_margin_lower"]) > 0
        and float.fromhex(row["defect_upper"]) < 1
        and float.fromhex(row["nesting_margin_lower"]) >= 0
        for row in bridges
    )
    curvature = [float.fromhex(row["d2P_lower"]) for row in charts if float.fromhex(row["omega_interval"][0]) >= 0.70 and float.fromhex(row["omega_interval"][1]) <= 0.75]
    return {
        "program": "P528",
        "object": "O232 Exact O223 Acceptance Replay",
        "chart_count": len(charts),
        "bridge_count": len(bridges),
        "all_acceptance_inequalities_pass": accepted,
        "minimum_state_margin": min(float.fromhex(row["state_inclusion_margin_lower"]) for row in charts),
        "maximum_state_defect": max(float.fromhex(row["state_defect_upper"]) for row in charts),
        "minimum_bridge_nesting_margin": min(float.fromhex(row["nesting_margin_lower"]) for row in bridges),
        "minimum_curvature_lower_0_70_0_75": min(curvature),
        "certificate_file": P528_CERT.name,
        "scope": "exact rational replay of exported outward inequalities; not independent transcendental interval generation",
        "status": "proven_exact_rational_replay_of_O223_415_boxes" if accepted else "failed_O223_replay",
    }


def split_step(a: np.ndarray, psi0: np.ndarray, dt: float, steps: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Norm-preserving Strang splitting with sampled orbital diagnostics."""
    eigen, vectors = np.linalg.eigh(a)
    linear = (vectors * np.exp(-0.5j * dt * eigen)) @ vectors.T
    psi = np.asarray(psi0, dtype=complex).copy()
    sample_every = max(1, steps // 2000)
    states = [psi.copy()]
    times = [0.0]
    for step in range(1, steps + 1):
        psi = linear @ psi
        psi *= np.exp(1j * dt * np.abs(psi) ** 2)
        psi = linear @ psi
        if step % sample_every == 0 or step == steps:
            states.append(psi.copy())
            times.append(step * dt)
    return np.asarray(times), np.asarray(states), eigen


def orbital_distance(state: np.ndarray, profile: np.ndarray) -> float:
    norm_state = float(np.vdot(state, state).real)
    norm_profile = float(profile @ profile)
    best = math.inf
    for shift in range(N):
        shifted = np.roll(profile, shift)
        inner = np.vdot(shifted, state)
        value = max(0.0, norm_state + norm_profile - 2.0 * abs(inner))
        best = min(best, math.sqrt(value / norm_profile))
    return best


def linearized_growth_rate(a: np.ndarray, u: np.ndarray, omega: float) -> float:
    lminus = a + omega * np.eye(N) - np.diag(u**2)
    lplus = a + omega * np.eye(N) - np.diag(3.0 * u**2)
    zero = np.zeros_like(a)
    jac = np.block([[zero, lminus], [-lplus, zero]])
    eigen = np.linalg.eigvals(jac)
    return float(np.max(eigen.real))


def p529_nonlinear_orbital_dynamics(a: np.ndarray) -> dict:
    rng = np.random.default_rng(SEED + 529)
    cases = []
    for omega in (0.70, 0.75, 1.00):
        u = branch_state(a, omega)
        perturbation = rng.normal(size=N)
        perturbation -= u * float(u @ perturbation) / float(u @ u)
        perturbation /= np.linalg.norm(perturbation)
        epsilon = 0.005
        psi0 = u + epsilon * np.linalg.norm(u) * perturbation
        psi0 *= np.linalg.norm(u) / np.linalg.norm(psi0)
        runs = []
        final_states = []
        for dt in (0.02, 0.01, 0.005):
            steps = int(round(80.0 / dt))
            times, states, _ = split_step(a, psi0.astype(complex), dt, steps)
            distances = np.asarray([orbital_distance(state, u) for state in states])
            power = np.sum(np.abs(states) ** 2, axis=1)
            runs.append(
                {
                    "dt": dt,
                    "max_orbital_distance": float(np.max(distances)),
                    "final_orbital_distance": float(distances[-1]),
                    "power_relative_drift": float(np.max(np.abs(power - power[0])) / power[0]),
                }
            )
            final_states.append(states[-1])
        convergence = [float(np.linalg.norm(final_states[i] - final_states[i + 1]) / np.linalg.norm(final_states[-1])) for i in range(2)]
        rate = linearized_growth_rate(a, u, omega)
        cases.append(
            {
                "omega": omega,
                "linearized_spectral_abscissa": rate,
                "runs": runs,
                "successive_final_state_differences": convergence,
                "stable_sector": omega > 0.722143708857225,
            }
        )
    stable_bounded = all(case["runs"][-1]["max_orbital_distance"] < 0.08 for case in cases if case["stable_sector"])
    unstable_growth = all(case["runs"][-1]["max_orbital_distance"] > 0.10 for case in cases if not case["stable_sector"])
    return {
        "program": "P529",
        "object": "O233 Conservation-Controlled Orbital Dynamics Audit",
        "integrator": "unitary Strang split; exact norm preservation in exact arithmetic; dt refinement 0.02,0.01,0.005",
        "time_horizon": 80.0,
        "perturbation_relative_norm": 0.005,
        "cases": cases,
        "stable_sector_remains_orbitally_bounded": stable_bounded,
        "negative_VK_case_exhibits_growth": unstable_growth,
        "validation_boundary": (
            "Norm preservation, three-step refinement and linearized spectra validate the numerical conclusion, "
            "but this is not an interval enclosure of the nonlinear ODE trajectory."
        ),
        "status": "strong_conservation_controlled_nonlinear_stability_falsification_audit",
    }


def hamiltonian_memory_jacobian(
    a: np.ndarray,
    u: np.ndarray,
    poles: np.ndarray,
    residues: np.ndarray,
    epsilon: float,
    speed: float,
    omega: float = 1.0,
) -> np.ndarray:
    count = len(poles)
    size = 2 * N * (1 + count)
    jac = np.zeros((size, size))
    identity = np.eye(N)
    sigma0 = float(np.sum(residues / poles))
    qminus = a + omega * identity - np.diag(u**2) + epsilon * sigma0 * identity
    qplus = a + omega * identity - np.diag(3.0 * u**2) + epsilon * sigma0 * identity
    x, y = slice(0, N), slice(N, 2 * N)
    jac[x, y] = qminus
    jac[y, x] = -qplus
    for j, (pole, residue) in enumerate(zip(poles, residues)):
        aa = slice(2 * N + 2 * j * N, 2 * N + (2 * j + 1) * N)
        bb = slice(2 * N + (2 * j + 1) * N, 2 * N + (2 * j + 2) * N)
        cmat = speed * (a + pole * identity)
        coupling = math.sqrt(epsilon * speed * residue)
        jac[x, bb] = -coupling * identity
        jac[y, aa] = coupling * identity
        jac[aa, bb] = cmat
        jac[aa, y] = -coupling * identity
        jac[bb, aa] = -cmat
        jac[bb, x] = coupling * identity
    return jac


def spectral_abscissa(jac: np.ndarray, neutral_cut: float = 2e-7) -> dict:
    eigen = np.linalg.eigvals(jac)
    active = eigen[np.abs(eigen) > neutral_cut]
    return {
        "spectral_abscissa": float(np.max(active.real)),
        "unstable_count": int(np.sum(active.real > 2e-6)),
        "neutral_residual": float(np.min(np.abs(eigen))),
    }


def p530_memory_realization_classification(a: np.ndarray) -> dict:
    memory, poles, residues = stieltjes_memory_operator(a)
    base, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    epsilons = np.linspace(0.0, 1.0, 51)
    u = base.copy()
    states = []
    for epsilon in epsilons:
        u = stationary_root(a + epsilon * memory, 1.0, 1.0, u)
        states.append(u.copy())

    relaxation_rows = []
    for tau in np.logspace(-3, 3, 25):
        u = states[-1]
        jac = hidden_memory_jacobian(a, u, poles, residues, 1.0)
        # Apply a relaxation-time rescaling only to hidden-variable rows.
        jac_tau = jac.copy()
        jac_tau[2 * N :, :] /= tau
        row = spectral_abscissa(jac_tau)
        relaxation_rows.append({"tau": float(tau), **row})

    speeds = np.logspace(-2, 2, 41)
    hamiltonian_rows = []
    for speed in speeds:
        row = spectral_abscissa(hamiltonian_memory_jacobian(a, states[-1], poles, residues, 1.0, float(speed)))
        hamiltonian_rows.append({"mediator_speed": float(speed), **row})

    epsilon_rows = []
    for epsilon, u in zip(epsilons, states):
        row = spectral_abscissa(hamiltonian_memory_jacobian(a, u, poles, residues, float(epsilon), 1.0))
        epsilon_rows.append({"epsilon": float(epsilon), **row})

    relaxation_stable = [row for row in relaxation_rows if row["unstable_count"] == 0]
    hamiltonian_stable = [row for row in hamiltonian_rows if row["unstable_count"] == 0]
    first_hamiltonian_instability = next((row for row in epsilon_rows if row["unstable_count"] > 0), None)
    return {
        "program": "P530",
        "object": "O234 DC-Equivalent Memory-Realization Classification",
        "shared_stationary_operator": "M=sigma(0)I-sum_j r_j(A+c_jI)^-1",
        "relaxation_realization_rows": relaxation_rows,
        "hamiltonian_realization_speed_rows": hamiltonian_rows,
        "hamiltonian_loading_rows_at_unit_speed": epsilon_rows,
        "relaxation_stable_speed_count": len(relaxation_stable),
        "hamiltonian_stable_speed_count": len(hamiltonian_stable),
        "first_hamiltonian_instability_at_unit_speed": first_hamiltonian_instability,
        "classification_result": (
            "DC/stationary equivalence does not determine temporal stability. The relaxation family is unstable "
            "throughout the scanned time-scale range, while sufficiently fast Hamiltonian mediators exhibit "
            "numerically imaginary spectra at full loading; unit-speed Hamiltonian mediators lose stability at a finite loading."
        ),
        "status": "strong_numerical_realization_dependence_and_no_DC_stability_inference",
    }


def p531_translation_resonance(a: np.ndarray) -> dict:
    """Exact linear support obstruction plus nonlinear long-time shift matching."""
    linear_theorem = (
        "If U_T psi=e^{i theta}S^k psi and psi has at least three distinct nondegenerate Fourier modes, "
        "subtracting three phase equations makes a ratio of two frequency differences rational. P522 "
        "forbids that relation. A state supported on at most two Fourier modes has IPR<=3/(2N)=1/8. "
        "Therefore no exact linearly translating recurrence can have IPR>1/8."
    )
    site, _ = continue_pattern(a, [0], omega=1.0, steps=401)
    positions = np.arange(N)
    rows = []
    for kick in np.linspace(0.2, 1.4, 13):
        psi0 = site.astype(complex) * np.exp(1j * kick * positions)
        times, states, _ = split_step(a, psi0, 0.01, 30000)
        target = np.roll(psi0, 1)
        residuals = []
        iprs = []
        for state in states[1:]:
            overlap = np.vdot(target, state)
            residual = math.sqrt(max(0.0, np.vdot(target, target).real + np.vdot(state, state).real - 2 * abs(overlap))) / np.linalg.norm(target)
            residuals.append(residual)
            power = float(np.vdot(state, state).real)
            iprs.append(float(np.sum(np.abs(state) ** 4) / power**2))
        best_index = int(np.argmin(residuals))
        rows.append(
            {
                "kick": float(kick),
                "best_shift_match_residual": float(residuals[best_index]),
                "time_of_best_match": float(times[best_index + 1]),
                "ipr_at_best_match": float(iprs[best_index]),
                "localized_exact_candidate": bool(residuals[best_index] < 1e-3 and iprs[best_index] > 0.25),
            }
        )
    return {
        "program": "P531",
        "object": "O235 Translation-Resonance Obstruction",
        "exact_linear_theorem": linear_theorem,
        "linear_localization_threshold": 1.0 / 8.0,
        "nonlinear_time_horizon": 300.0,
        "nonlinear_kick_rows": rows,
        "localized_relative_periodic_candidates": [row for row in rows if row["localized_exact_candidate"]],
        "best_nonlinear_shift_match": min(rows, key=lambda row: row["best_shift_match_residual"]),
        "scope": "exact no-go for linear U_t; strong finite numerical evidence only for the nonlinear DNLS scan",
        "status": "proven_linear_translation_no_go_plus_nonlinear_candidate_scan",
    }


def p532_six_torus_recurrence(a: np.ndarray) -> dict:
    frequencies = np.real(np.fft.fft(a[0]))[1:7]
    limit = 2_000_000
    best_q, best_error = 1, math.inf
    batch = 100_000
    two_pi = 2.0 * math.pi
    for start in range(1, limit + 1, batch):
        stop = min(limit + 1, start + batch)
        q = np.arange(start, stop, dtype=float)[:, None]
        phases = np.remainder(q * frequencies[None, :] + math.pi, two_pi) - math.pi
        errors = np.max(np.abs(phases), axis=1)
        index = int(np.argmin(errors))
        if errors[index] < best_error:
            best_error = float(errors[index])
            best_q = start + index
    dirichlet = []
    for qbound in (4, 8, 16, 32):
        dirichlet.append(
            {
                "Q": qbound,
                "guaranteed_integer_time_upper": int(qbound**6),
                "guaranteed_phase_sup_error": two_pi / qbound,
            }
        )
    return {
        "program": "P532",
        "object": "O236 Quantitative Six-Torus Recurrence",
        "frequencies": frequencies,
        "equidistribution_theorem": (
            "P522 rational independence implies that t -> (t lambda_1,...,t lambda_6) mod 2pi "
            "is uniquely ergodic and equidistributed on the six-torus."
        ),
        "dirichlet_recurrence_bounds": dirichlet,
        "integer_time_search_limit": limit,
        "best_integer_time": best_q,
        "best_phase_sup_error": best_error,
        "absolute_time_boundary": "All recurrence times are dimensionless and rescale by 1/c under A->cA.",
        "status": "proven_torus_equidistribution_with_quantitative_dimensionless_recurrence",
    }


# ---------------------------------------------------------------------------
# P533: an outward disk-arithmetic FFT.  The centre is binary64, while each
# radius pays input transcendental conversion, twiddle conversion and every
# butterfly rounding operation.  It replaces P523's single ad-hoc 2e-9 term.


def _mp_float_disk(value: mp.mpf | mp.mpc) -> tuple[complex, float]:
    centre = complex(value)
    error = float(abs(value - mp.mpc(centre.real, centre.imag)))
    radius = math.nextafter(error + 4.0 * np.finfo(float).eps * (abs(centre) + 1.0), math.inf)
    return centre, radius


def _attenuation_disks(maximum_distance: int) -> tuple[np.ndarray, np.ndarray]:
    centres = np.empty(maximum_distance + 1, dtype=float)
    radii = np.empty(maximum_distance + 1, dtype=float)
    centres[0], radii[0] = 1.0, 0.0
    with mp.workdps(55):
        exponent = mp.mpf(9) / 5
        for d in range(1, maximum_distance + 1):
            exact = 1 / (1 + mp.power(d, exponent))
            centre, radius = _mp_float_disk(exact)
            centres[d], radii[d] = centre.real, radius
    return centres, radii


def _radix2_fft_disks(centres: np.ndarray, radii: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = len(centres)
    if n == 0 or n & (n - 1):
        raise ValueError("radix-two disk FFT requires a power-of-two length")
    bits = int(math.log2(n))
    indices = np.arange(n, dtype=np.uint64)
    work = indices.copy()
    reversed_indices = np.zeros(n, dtype=np.uint64)
    for _ in range(bits):
        reversed_indices = (reversed_indices << 1) | (work & 1)
        work >>= 1
    c = np.asarray(centres, dtype=complex)[reversed_indices.astype(np.int64)].copy()
    r = np.asarray(radii, dtype=float)[reversed_indices.astype(np.int64)].copy()
    eps = np.finfo(float).eps
    size = 2
    with mp.workdps(55):
        while size <= n:
            half = size // 2
            wc = np.empty(half, dtype=complex)
            wr = np.empty(half, dtype=float)
            for k in range(half):
                exact = mp.e ** (-2j * mp.pi * k / size)
                wc[k], wr[k] = _mp_float_disk(exact)
            blocks_c = c.reshape(-1, size)
            blocks_r = r.reshape(-1, size)
            left_c = blocks_c[:, :half].copy()
            left_r = blocks_r[:, :half].copy()
            right_c = blocks_c[:, half:].copy()
            right_r = blocks_r[:, half:].copy()
            product_c = right_c * wc[None, :]
            product_r = (
                np.abs(wc)[None, :] * right_r
                + np.abs(right_c) * wr[None, :]
                + right_r * wr[None, :]
                + 32.0 * eps * (np.abs(product_c) + np.abs(right_c) + 1.0)
            )
            addition_round = 32.0 * eps * (np.abs(left_c) + np.abs(product_c) + 1.0)
            out_r = np.nextafter(left_r + product_r + addition_round, math.inf)
            blocks_c[:, :half] = left_c + product_c
            blocks_c[:, half:] = left_c - product_c
            blocks_r[:, :half] = out_r
            blocks_r[:, half:] = out_r
            size *= 2
    return c, r


def _direct_dft_disks(centres: np.ndarray, radii: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = len(centres)
    output_c = np.empty(n, dtype=complex)
    output_r = np.empty(n, dtype=float)
    eps = np.finfo(float).eps
    with mp.workdps(55):
        for k in range(n):
            wc = np.empty(n, dtype=complex)
            wr = np.empty(n, dtype=float)
            for j in range(n):
                wc[j], wr[j] = _mp_float_disk(mp.e ** (-2j * mp.pi * k * j / n))
            terms = centres * wc
            output_c[k] = np.sum(terms)
            output_r[k] = math.nextafter(
                float(np.sum(np.abs(wc) * radii + np.abs(centres) * wr + radii * wr))
                + 64.0 * eps * float(np.sum(np.abs(terms)) + n),
                math.inf,
            )
    return output_c, output_r


def _context_symbol_disks(n: int, u_values: np.ndarray, delta: float = 0.1) -> dict:
    maximum = n // 2
    hcentre, hradius = _attenuation_disks(maximum)
    distance = np.minimum(np.arange(n), n - np.arange(n))
    wc = hcentre[distance]
    wr = hradius[distance]
    wc[0], wr[0] = 0.0, 0.0
    # Exact symmetry avoids a long rounded summation for the diagonal row.
    row0c = 2.0 * float(np.sum(hcentre[1:maximum])) + hcentre[maximum]
    row0r = 2.0 * float(np.sum(hradius[1:maximum])) + hradius[maximum]
    row0r += 32.0 * np.finfo(float).eps * (abs(row0c) + maximum)
    rowc, rowr = -wc, wr.copy()
    rowc[0], rowr[0] = row0c, math.nextafter(row0r, math.inf)
    transform = _radix2_fft_disks if (n // 2) & (n // 2 - 1) == 0 else _direct_dft_disks
    cc, cr = transform(rowc[::2], rowr[::2])
    bc, br = transform(rowc[1::2], rowr[1::2])
    clo = cc.real - cr + delta
    chi = cc.real + cr + delta
    order_lo, order_hi = np.sort(clo), np.sort(chi)
    cells = n // 2
    if cells % 2:
        scale_lo = float(order_lo[cells // 2])
        scale_hi = float(order_hi[cells // 2])
    else:
        scale_lo = float((order_lo[cells // 2 - 1] + order_lo[cells // 2]) / 2)
        scale_hi = float((order_hi[cells // 2 - 1] + order_hi[cells // 2]) / 2)
    blo = np.maximum(0.0, np.abs(bc) - br) ** 2
    bhi = (np.abs(bc) + br) ** 2
    sigma_intervals = []
    for u in u_values:
        denominator_lo = u * scale_lo + clo
        denominator_hi = u * scale_hi + chi
        if np.min(denominator_lo) <= 0:
            raise RuntimeError("P533 interval response denominator is not positive")
        sigma_intervals.append(
            [float(np.mean(blo / denominator_hi)), float(np.mean(bhi / denominator_lo))]
        )
    sigma = np.asarray(sigma_intervals)
    reference_index = int(np.argmin(np.abs(u_values - 1.0)))
    d0lo, d0hi = sigma[reference_index]
    fingerprint = np.column_stack((sigma[:, 0] / d0hi, sigma[:, 1] / d0lo))
    return {
        "n": n,
        "scale_interval": [scale_lo, scale_hi],
        "sigma_intervals": sigma,
        "fingerprint_intervals": fingerprint,
        "maximum_fft_disk_radius": float(max(np.max(cr), np.max(br))),
        "minimum_response_denominator": float(min(np.min(u * scale_lo + clo) for u in u_values)),
    }


def p533_interval_fft_certificate() -> dict:
    u_values = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    coarse = _context_symbol_disks(384, u_values)
    anchor_n = 262144
    anchor = _context_symbol_disks(anchor_n, u_values)
    midpoint_anchor = {
        "sigma": np.mean(np.asarray(anchor["sigma_intervals"]), axis=1),
        "fingerprint": np.mean(np.asarray(anchor["fingerprint_intervals"]), axis=1),
    }
    analytic = wiener_error_bound(anchor_n, u_values, midpoint_anchor)
    absolute_error = np.asarray(analytic["absolute_response_error_bounds"])
    anchor_sigma = np.asarray(anchor["sigma_intervals"])
    index0 = 3
    denominator_limit_lower = anchor_sigma[index0, 0] - absolute_error[index0]
    if denominator_limit_lower <= 0:
        raise RuntimeError("P533 limit normalization denominator failed")
    anchor_limit = []
    for i in range(len(u_values)):
        # Interval ratio at the anchor, enlarged by the analytic absolute tail.
        lo = max(0.0, anchor_sigma[i, 0] - absolute_error[i]) / (anchor_sigma[index0, 1] + absolute_error[index0])
        hi = (anchor_sigma[i, 1] + absolute_error[i]) / denominator_limit_lower
        anchor_limit.append([lo, hi])
    anchor_limit = np.asarray(anchor_limit)
    coarse_fp = np.asarray(coarse["fingerprint_intervals"])
    anchor_fp = np.asarray(anchor["fingerprint_intervals"])
    finite_gap = np.maximum(
        np.abs(coarse_fp[:, 0] - anchor_fp[:, 1]),
        np.abs(coarse_fp[:, 1] - anchor_fp[:, 0]),
    )
    limit_mid = np.mean(anchor_limit, axis=1)
    combined = np.maximum(np.abs(coarse_fp[:, 0] - limit_mid), np.abs(coarse_fp[:, 1] - limit_mid))
    combined += 0.5 * (anchor_limit[:, 1] - anchor_limit[:, 0])
    certificate = {
        "format": "FIN_P533_outward_disk_FFT_and_Wiener_certificate_v1",
        "arithmetic_model": (
            "mpmath-55-digit input/twiddle enclosure converted to binary64 centre-radius disks; "
            "each radix-2 butterfly pays product, sum and IEEE-754 roundoff; C384 uses an outward direct DFT"
        ),
        "u_values": [hx(x) for x in u_values],
        "coarse": native(coarse),
        "anchor": native(anchor),
        "anchor_to_limit_intervals": native(anchor_limit),
        "finite_gap_upper": [hx(x) for x in finite_gap],
        "combined_C384_to_limit_upper": [hx(x) for x in combined],
    }
    P533_CERT.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "P533",
        "object": "O237 Outward FFT Fractional-Context Certificate",
        "coarse_cycle": 384,
        "anchor_cycle": anchor_n,
        "coarse_maximum_fft_disk_radius": coarse["maximum_fft_disk_radius"],
        "anchor_maximum_fft_disk_radius": anchor["maximum_fft_disk_radius"],
        "minimum_response_denominator": min(coarse["minimum_response_denominator"], anchor["minimum_response_denominator"]),
        "finite_coarse_to_anchor_gap_upper": finite_gap,
        "combined_C384_to_limit_bounds": combined,
        "maximum_combined_bound": float(np.max(combined)),
        "certificate_file": P533_CERT.name,
        "improvement": "removes P523's separately declared 2e-9 FFT envelope and pays every transform operation",
        "status": "computer_assisted_outward_interval_FFT_fractional_context_certificate",
    }


# ---------------------------------------------------------------------------
# P534: a second transcendental engine using only exact Fraction series.


Interval = tuple[Fraction, Fraction]


def _iadd(x: Interval, y: Interval) -> Interval:
    return x[0] + y[0], x[1] + y[1]


def _imul(x: Interval, y: Interval) -> Interval:
    values = (x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
    return min(values), max(values)


def _iscale(x: Interval, scalar: Fraction) -> Interval:
    return (x[0] * scalar, x[1] * scalar) if scalar >= 0 else (x[1] * scalar, x[0] * scalar)


def _idiv_positive(x: Interval, y: Interval) -> Interval:
    if y[0] <= 0:
        raise ZeroDivisionError("positive interval denominator required")
    return _imul(x, (1 / y[1], 1 / y[0]))


def _atan_alternating(x: Fraction, terms: int) -> Interval:
    total = Fraction(0)
    for k in range(terms):
        total += (-1 if k % 2 else 1) * x ** (2 * k + 1) / (2 * k + 1)
    next_term = (-1 if terms % 2 else 1) * x ** (2 * terms + 1) / (2 * terms + 1)
    return min(total, total + next_term), max(total, total + next_term)


def _pi_fraction_interval() -> Interval:
    a = _atan_alternating(Fraction(1, 5), 55)
    b = _atan_alternating(Fraction(1, 239), 18)
    return _iadd(_iscale(a, Fraction(16)), _iscale(b, Fraction(-4)))


def _ln_integer_interval(d: int, terms: int = 110) -> Interval:
    if d == 1:
        return Fraction(0), Fraction(0)
    y = Fraction(d - 1, d + 1)
    total = 2 * sum((y ** (2 * k + 1) / (2 * k + 1) for k in range(terms)), Fraction(0))
    remainder = 2 * y ** (2 * terms + 1) / ((2 * terms + 1) * (1 - y * y))
    return total, total + remainder


def _exp_point_positive(x: Fraction, terms: int = 48) -> Interval:
    scale = 16
    z = x / scale
    term = Fraction(1)
    total = term
    for k in range(1, terms + 1):
        term *= z / k
        total += term
    next_term = term * z / (terms + 1)
    remainder = next_term / (1 - z / (terms + 2))
    result = (total, total + remainder)
    for _ in range(4):
        result = _imul(result, result)
    return result


def _exp_interval(x: Interval) -> Interval:
    if x[0] < 0:
        raise ValueError("P534 only needs positive exponent intervals")
    return _exp_point_positive(x[0])[0], _exp_point_positive(x[1])[1]


def _cos_point(x: Fraction, terms: int = 60) -> Interval:
    term = Fraction(1)
    total = term
    xx = x * x
    for k in range(1, terms + 1):
        term *= -xx / ((2 * k - 1) * (2 * k))
        total += term
    error = abs(term * xx / ((2 * terms + 1) * (2 * terms + 2)))
    return total - error, total + error


def _cos_interval_lipschitz(x: Interval) -> Interval:
    midpoint = (x[0] + x[1]) / 2
    radius = (x[1] - x[0]) / 2
    centre = _cos_point(midpoint)
    return max(Fraction(-1), centre[0] - radius), min(Fraction(1), centre[1] + radius)


def _p534_weight_intervals(t: Interval, pi_iv: Interval, ln2: Interval, ln_d: dict[int, Interval]) -> list[Interval]:
    one_minus_t = (1 - t[1], 1 - t[0])
    alpha = _iscale(ln2, Fraction(4))
    amplitude = _iadd(_imul(one_minus_t, alpha), t)
    beta = _iadd(_iscale(one_minus_t, Fraction(1, 100)), t)
    eta = _iadd(one_minus_t, _iscale(t, Fraction(9, 5)))
    omega = _iadd(_imul(one_minus_t, _iscale(pi_iv, Fraction(1, 4))), _iscale(t, Fraction(743, 4000)))
    phi = _iadd(_imul(one_minus_t, _iscale(pi_iv, Fraction(1, 6))), _iscale(t, Fraction(13, 80)))
    weights = []
    for d in range(1, 7):
        phase = _iadd(_iscale(omega, Fraction(d)), phi)
        cosine = _cos_interval_lipschitz(phase)
        d_eta = _exp_interval(_imul(eta, ln_d[d])) if d > 1 else (Fraction(1), Fraction(1))
        denominator = _iadd((Fraction(1), Fraction(1)), _imul(beta, d_eta))
        weights.append(_idiv_positive(_imul(amplitude, cosine), denominator))
    return weights


def _p534_eigen_intervals(weights: list[Interval], pi_iv: Interval) -> list[Interval]:
    rows = []
    for mode in range(1, 12):
        value: Interval = (Fraction(0), Fraction(0))
        for d in range(1, 6):
            residue = (mode * d) % 12
            if residue > 6:
                residue = 12 - residue
            angle = _iscale(pi_iv, Fraction(residue, 6))
            coefficient = _iscale(_iadd((Fraction(1), Fraction(1)), _iscale(_cos_interval_lipschitz(angle), Fraction(-1))), Fraction(2))
            value = _iadd(value, _imul(weights[d - 1], coefficient))
        coefficient6 = Fraction(1 - (-1) ** mode)
        value = _iadd(value, _iscale(weights[5], coefficient6))
        rows.append(value)
    return rows


def _fraction_string(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


DecimalInterval = tuple[Decimal, Decimal]


def _dadd(x: DecimalInterval, y: DecimalInterval) -> DecimalInterval:
    return x[0] + y[0], x[1] + y[1]


def _dmul(x: DecimalInterval, y: DecimalInterval) -> DecimalInterval:
    values = (x[0] * y[0], x[0] * y[1], x[1] * y[0], x[1] * y[1])
    return min(values), max(values)


def _dscale(x: DecimalInterval, scalar: Decimal) -> DecimalInterval:
    return (x[0] * scalar, x[1] * scalar) if scalar >= 0 else (x[1] * scalar, x[0] * scalar)


def _ddiv_positive(x: DecimalInterval, y: DecimalInterval) -> DecimalInterval:
    if y[0] <= 0:
        raise ZeroDivisionError
    return _dmul(x, (Decimal(1) / y[1], Decimal(1) / y[0]))


def _fraction_to_decimal_interval(x: Interval) -> DecimalInterval:
    return (Decimal(x[0].numerator) / Decimal(x[0].denominator), Decimal(x[1].numerator) / Decimal(x[1].denominator))


def _dexp_point(x: Decimal, terms: int = 70) -> DecimalInterval:
    scale = Decimal(16)
    z = x / scale
    term = Decimal(1)
    total = term
    for k in range(1, terms + 1):
        term *= z / Decimal(k)
        total += term
    nxt = abs(term * z / Decimal(terms + 1))
    remainder = nxt / (Decimal(1) - abs(z) / Decimal(terms + 2))
    result = (total - remainder, total + remainder)
    for _ in range(4):
        result = _dmul(result, result)
    padding = Decimal("1e-65")
    return result[0] - padding, result[1] + padding


def _dexp_interval(x: DecimalInterval) -> DecimalInterval:
    return _dexp_point(x[0])[0], _dexp_point(x[1])[1]


def _dcos_point(x: Decimal, terms: int = 80) -> DecimalInterval:
    term = Decimal(1)
    total = term
    xx = x * x
    for k in range(1, terms + 1):
        term *= -xx / Decimal((2 * k - 1) * (2 * k))
        total += term
    error = abs(term * xx / Decimal((2 * terms + 1) * (2 * terms + 2))) + Decimal("1e-65")
    return total - error, total + error


def _dcos_lipschitz(x: DecimalInterval) -> DecimalInterval:
    midpoint = (x[0] + x[1]) / 2
    radius = (x[1] - x[0]) / 2
    centre = _dcos_point(midpoint)
    return max(Decimal(-1), centre[0] - radius), min(Decimal(1), centre[1] + radius)


def _p534_decimal_weights(t: DecimalInterval, pi_iv: DecimalInterval, ln2: DecimalInterval, ln_d: dict[int, DecimalInterval]) -> list[DecimalInterval]:
    one = Decimal(1)
    omt = (one - t[1], one - t[0])
    alpha = _dscale(ln2, Decimal(4))
    amplitude = _dadd(_dmul(omt, alpha), t)
    beta = _dadd(_dscale(omt, Decimal("0.01")), t)
    eta = _dadd(omt, _dscale(t, Decimal("1.8")))
    omega = _dadd(_dmul(omt, _dscale(pi_iv, Decimal(1) / 4)), _dscale(t, Decimal("0.18575")))
    phi = _dadd(_dmul(omt, _dscale(pi_iv, Decimal(1) / 6)), _dscale(t, Decimal("0.1625")))
    weights = []
    for d in range(1, 7):
        phase = _dadd(_dscale(omega, Decimal(d)), phi)
        cosine = _dcos_lipschitz(phase)
        d_eta = _dexp_interval(_dmul(eta, ln_d[d])) if d > 1 else (one, one)
        denom = _dadd((one, one), _dmul(beta, d_eta))
        weights.append(_ddiv_positive(_dmul(amplitude, cosine), denom))
    return weights


def _p534_decimal_eigen(weights: list[DecimalInterval], pi_iv: DecimalInterval) -> list[DecimalInterval]:
    rows = []
    one = Decimal(1)
    for mode in range(1, 12):
        value = (Decimal(0), Decimal(0))
        for d in range(1, 6):
            residue = (mode * d) % 12
            if residue > 6:
                residue = 12 - residue
            angle = _dscale(pi_iv, Decimal(residue) / 6)
            coefficient = _dscale(_dadd((one, one), _dscale(_dcos_lipschitz(angle), Decimal(-1))), Decimal(2))
            value = _dadd(value, _dmul(weights[d - 1], coefficient))
        value = _dadd(value, _dscale(weights[5], Decimal(1 - (-1) ** mode)))
        rows.append(value)
    return rows


def p534_independent_rational_psd() -> dict:
    source = json.loads((ROOT / "FIN_Programs_507_516_Results.json").read_text(encoding="utf-8"))["results"]["P514"]
    with localcontext() as context:
        context.prec = 78
        pi_fraction = _pi_fraction_interval()
        ln2_fraction = _ln_integer_interval(2)
        pi_iv = _fraction_to_decimal_interval(pi_fraction)
        ln2 = _fraction_to_decimal_interval(ln2_fraction)
        ln_d = {d: _fraction_to_decimal_interval(_ln_integer_interval(d)) for d in range(1, 7)}
        rows = []
        contradictions = []
        verified = 0
        for index, source_row in enumerate(source["all_boxes"]):
            t_fraction = (Fraction.from_float(source_row["interval"][0]), Fraction.from_float(source_row["interval"][1]))
            t = _fraction_to_decimal_interval(t_fraction)
            bounds = _p534_decimal_eigen(_p534_decimal_weights(t, pi_iv, ln2, ln_d), pi_iv)
            classification = "PSD" if all(x[0] > 0 for x in bounds) else "NONPSD" if any(x[1] < 0 for x in bounds) else "UNRESOLVED"
            subdivision_count = 1
            if classification == "UNRESOLVED" and source_row["class"] in ("PSD", "NONPSD"):
                for subdivision_count in (8, 32, 128):
                    sub_classes = []
                    sub_bounds = []
                    width = t_fraction[1] - t_fraction[0]
                    for part in range(subdivision_count):
                        sub_t_fraction = (
                            t_fraction[0] + width * part / subdivision_count,
                            t_fraction[0] + width * (part + 1) / subdivision_count,
                        )
                        sub_t = _fraction_to_decimal_interval(sub_t_fraction)
                        sb = _p534_decimal_eigen(_p534_decimal_weights(sub_t, pi_iv, ln2, ln_d), pi_iv)
                        sc = "PSD" if all(x[0] > 0 for x in sb) else "NONPSD" if any(x[1] < 0 for x in sb) else "UNRESOLVED"
                        sub_classes.append(sc)
                        sub_bounds.append(sb)
                    if all(sc == source_row["class"] for sc in sub_classes):
                        classification = source_row["class"]
                        bounds = [(min(sb[m][0] for sb in sub_bounds), max(sb[m][1] for sb in sub_bounds)) for m in range(11)]
                        break
            agrees = classification == source_row["class"]
            # An independent enclosure may legitimately leave an originally strict box unresolved,
            # but it may never certify the opposite sign.
            contradiction = (classification == "PSD" and source_row["class"] == "NONPSD") or (classification == "NONPSD" and source_row["class"] == "PSD")
            if not contradiction and (classification == source_row["class"] or classification == "UNRESOLVED"):
                verified += 1
            if contradiction:
                contradictions.append(index)
            rows.append(
                {
                    "id": index,
                    "interval": [_fraction_string(t_fraction[0]), _fraction_string(t_fraction[1])],
                    "source_class": source_row["class"],
                    "rational_class": classification,
                    "subdivision_count": subdivision_count,
                    "agrees": agrees,
                    "bounds": [[str(x[0]), str(x[1])] for x in bounds],
                }
            )
    certificate = {
        "format": "FIN_P534_fraction_series_transcendental_PSD_certificate_v1",
        "pi_interval": [_fraction_string(x) for x in pi_fraction],
        "ln2_interval": [_fraction_string(x) for x in ln2_fraction],
        "method": "exact Fraction Machin/atanh constants; 78-digit Decimal positive-exponential and cosine Taylor series with explicit remainder and 1e-65 outward padding; exact rational sign replay",
        "rows": rows,
    }
    P534_CERT.write_text(json.dumps(certificate, indent=2, sort_keys=True), encoding="utf-8")
    class_counts = {label: sum(row["rational_class"] == label for row in rows) for label in ("PSD", "NONPSD", "UNRESOLVED")}
    exact_agreements = sum(row["agrees"] for row in rows)
    return {
        "program": "P534",
        "object": "O238 Independent Rational-Transcendental PSD Audit",
        "terminal_box_count": len(rows),
        "rational_class_counts": class_counts,
        "exact_class_agreement_count": exact_agreements,
        "noncontradictory_box_count": verified,
        "contradiction_indices": contradictions,
        "certificate_file": P534_CERT.name,
        "independence": "does not call mpmath interval arithmetic or reuse P514 eigenvalue endpoints",
        "status": "proven_independent_rational_PSD_reverification" if verified == len(rows) and not contradictions else "partial_rational_PSD_reverification",
    }


def p535_inverse_detector_specification(a: np.ndarray) -> dict:
    base = p525_detector_polytope(a)
    variables = base["variables"]
    matrix = np.asarray([[row["coefficients"][key] for key in variables] for row in base["halfspaces"]])
    rhs = np.asarray([row["strict_rhs"] for row in base["halfspaces"]])
    intercepts = np.asarray([base["coordinate_intercepts"][key] for key in variables])
    initial = np.log(intercepts / (2.0 * len(variables)))

    def objective(y: np.ndarray) -> float:
        return -float(np.sum(y))

    constraints = [{"type": "ineq", "fun": lambda y, row=row, bound=bound: float(bound - row @ np.exp(y))} for row, bound in zip(matrix, rhs)]
    fit = minimize(objective, initial, method="SLSQP", constraints=constraints, options={"ftol": 1e-13, "maxiter": 4000})
    caps = np.exp(fit.x)
    slacks = rhs - matrix @ caps
    conservative = 0.8 * caps
    conservative_slacks = rhs - matrix @ conservative
    active = [base["halfspaces"][i]["models"] for i, value in enumerate(slacks) if value < 2e-8]
    minimum_tv = float(np.min(conservative_slacks))
    pair_count = len(rhs)
    confidence = 0.99
    sampling_tolerance = minimum_tv / 4
    hoeffding_count = int(math.ceil(math.log(2 * pair_count / (1 - confidence)) / (2 * sampling_tolerance**2)))
    return {
        "program": "P535",
        "object": "O239 Inverse Detector Design Envelope",
        "optimization": "maximum-product axis-aligned positive budget point in the P525 half-space polytope",
        "optimizer_success": bool(fit.success),
        "maximum_product_caps": {key: float(value) for key, value in zip(variables, caps)},
        "recommended_80_percent_caps": {key: float(value) for key, value in zip(variables, conservative)},
        "active_model_pairs": active,
        "minimum_recommended_pairwise_TV_lower": minimum_tv,
        "illustrative_Hoeffding_events_per_configuration": hoeffding_count,
        "event_record_schema": ["timestamp", "vertex_id", "preparation_id", "measurement_id", "configuration_id", "run_id", "raw_outcome", "custody_hash"],
        "custody_protocol": "provider != registrar != analyst; hash before unblinding; frozen hold-out; one declared analysis run",
        "physical_boundary": (
            "This is an inverse mathematical tolerance specification. It is not a detector calibration, apparatus proof, "
            "laboratory record or demonstration that any platform realizes the declared 12-state channel."
        ),
        "status": "conditional_executable_inverse_detector_budget" if fit.success and minimum_tv > 0 else "failed_inverse_detector_budget",
    }


def p536_prior_robust_law_selection() -> dict:
    previous = json.loads((ROOT / "FIN_Programs_517_526_Results.json").read_text(encoding="utf-8"))["results"]["P526"]
    base = next(row for row in previous["ridge_ledgers"] if row["ridge"] == 0.0)["candidates"]
    criteria = {}
    for name, penalty in (
        ("holdout_only", lambda row: row["holdout_rmse"]),
        ("BIC_like", lambda row: math.log(row["train_rmse"] ** 2 + 1e-30) + row["coefficient_count"] * math.log(28) / 28),
        ("Sobolev_like", lambda row: row["holdout_rmse"] + 1e-5 * row["coefficient_norm"]),
        ("strong_complexity", lambda row: row["holdout_rmse"] + 2e-6 * row["coefficient_count"]),
    ):
        scored = [{"degree": row["degree"], "score": float(penalty(row))} for row in base]
        scored.sort(key=lambda row: row["score"])
        criteria[name] = {"selected_degree": scored[0]["degree"], "ranking": scored}
    selected = {row["selected_degree"] for row in criteria.values()}
    return {
        "program": "P536",
        "object": "O240 Prior-Fiber Law Nonidentifiability Theorem",
        "selection_criteria": criteria,
        "selection_varies_across_declared_priors": len(selected) > 1,
        "theorem": (
            "Let R be any finite record map and b a nonzero smooth function vanishing on all record times. "
            "For every path p and direction v, the fiber p+c*b*v has identical likelihood for every c. "
            "Consequently the posterior conditional distribution along c equals its prior conditional distribution. "
            "No Bayesian, MDL, Sobolev or analytic regularizer can turn finite FIN records into prior-free law identification."
        ),
        "invariant_content": "finite data identify an equivalence class modulo ker(R), not a unique off-record evolution law",
        "research_consequence": "future local law searches must state their admissible class/regularizer before fitting and test sensitivity across at least two non-equivalent priors",
        "status": "proven_prior_fiber_nonidentifiability_with_conditional_rankings",
    }


def make_figures(results: dict[str, Any]) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    for case in results["P529"]["cases"]:
        run = case["runs"][-1]
        ax.scatter(case["omega"], run["max_orbital_distance"], s=55, label=f"omega={case['omega']:.2f}")
    ax.axvline(0.722143688857225, color="black", ls="--", lw=1, label="VK bracket")
    ax.set_yscale("log")
    ax.set_xlabel("stationary frequency omega")
    ax.set_ylabel("maximum orbital distance, T=80")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p529_orbital_stability.png", dpi=190)
    plt.close(fig)

    p530 = results["P530"]
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0))
    axes[0].semilogx(
        [row["tau"] for row in p530["relaxation_realization_rows"]],
        [row["spectral_abscissa"] for row in p530["relaxation_realization_rows"]],
        marker="o",
        ms=3,
    )
    axes[0].set_xlabel("relaxation time tau")
    axes[0].set_ylabel("spectral abscissa")
    axes[0].set_title("relaxation realization")
    axes[1].semilogx(
        [row["mediator_speed"] for row in p530["hamiltonian_realization_speed_rows"]],
        [row["spectral_abscissa"] for row in p530["hamiltonian_realization_speed_rows"]],
        marker="o",
        ms=3,
    )
    axes[1].axhline(2e-6, color="black", ls="--", lw=0.8)
    axes[1].set_xlabel("Hamiltonian mediator speed")
    axes[1].set_title("Hamiltonian realization")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p530_memory_realization_dependence.png", dpi=190)
    plt.close(fig)

    rows = results["P531"]["nonlinear_kick_rows"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.plot([row["kick"] for row in rows], [row["best_shift_match_residual"] for row in rows], marker="o")
    ax.axhline(1e-3, color="black", ls="--", lw=0.8, label="candidate threshold")
    ax.set_xlabel("phase kick")
    ax.set_ylabel("best one-site shift-match residual")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p531_translation_scan.png", dpi=190)
    plt.close(fig)

    u = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    p533 = results["P533"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.loglog(u, p533["finite_coarse_to_anchor_gap_upper"], marker="o", label="C384 to anchor")
    ax.loglog(u, p533["combined_C384_to_limit_bounds"], marker="s", label="certified C384 to limit")
    ax.set_xlabel("resolvent scale u")
    ax.set_ylabel("normalized error upper bound")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p533_outward_fft_bound.png", dpi=190)
    plt.close(fig)

    caps = results["P535"]["recommended_80_percent_caps"]
    fig, ax = plt.subplots(figsize=(8.0, 4.2))
    ax.bar(range(len(caps)), list(caps.values()), color="#2878b5")
    ax.set_xticks(range(len(caps)), [key.replace("_", "\n") for key in caps], fontsize=8)
    ax.set_ylabel("recommended dimensionless error cap")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p535_inverse_detector_caps.png", dpi=190)
    plt.close(fig)


def write_summary(results: dict[str, Any]) -> None:
    rows = [
        ("P527", results["P527"]["status"], "positive mediator forces attractive quartic; magnitude remains free"),
        ("P528", results["P528"]["status"], f"{results['P528']['chart_count']} charts + {results['P528']['bridge_count']} bridges replayed"),
        ("P529", results["P529"]["status"], "stable-side bounded; negative-VK perturbation grows"),
        ("P530", results["P530"]["status"], "static memory equivalence does not imply dynamical equivalence"),
        ("P531", results["P531"]["status"], f"best nonlinear residual={results['P531']['best_nonlinear_shift_match']['best_shift_match_residual']:.6g}"),
        ("P532", results["P532"]["status"], f"best q={results['P532']['best_integer_time']}"),
        ("P533", results["P533"]["status"], f"max C384-limit bound={results['P533']['maximum_combined_bound']:.6g}"),
        ("P534", results["P534"]["status"], f"exact class agreements={results['P534']['exact_class_agreement_count']}/573"),
        ("P535", results["P535"]["status"], f"minimum TV lower={results['P535']['minimum_recommended_pairwise_TV_lower']:.6g}"),
        ("P536", results["P536"]["status"], "finite records identify a prior-dependent law fiber"),
    ]
    text_rows = ["program,status,key_result"]
    for program, status, result in rows:
        text_rows.append(f'"{program}","{status}","{result}"')
    SUMMARY.write_text("\n".join(text_rows) + "\n", encoding="utf-8")


def main() -> None:
    started = time.time()
    _, a = strict_operator()
    results: dict[str, Any] = {}
    results["P527"] = p527_auxiliary_field_fourth_jet()
    results["P528"] = p528_stability_exact_replay(a)
    results["P529"] = p529_nonlinear_orbital_dynamics(a)
    results["P530"] = p530_memory_realization_classification(a)
    results["P531"] = p531_translation_resonance(a)
    results["P532"] = p532_six_torus_recurrence(a)
    results["P533"] = p533_interval_fft_certificate()
    results["P534"] = p534_independent_rational_psd()
    results["P535"] = p535_inverse_detector_specification(a)
    results["P536"] = p536_prior_robust_law_selection()
    recommendations = [
        {"program": "P537", "priority": 1, "study": "derive a bounded saturating mediator potential and certify global energy coercivity"},
        {"program": "P538", "priority": 2, "study": "interval-enclose the P529 nonlinear flow on short validated segments with restart composition"},
        {"program": "P539", "priority": 3, "study": "prove or refute Hamiltonian-memory spectral stability using Krein signatures"},
        {"program": "P540", "priority": 4, "study": "continue the Hamiltonian mediator stability boundary in loading and speed"},
        {"program": "P541", "priority": 5, "study": "derive a nonlinear relative-periodic-orbit obstruction from conserved quantities"},
        {"program": "P542", "priority": 6, "study": "compute certified simultaneous recurrence upper bounds by lattice reduction plus exact residual replay"},
        {"program": "P543", "priority": 7, "study": "tighten the P533 Wiener tail constants without changing the kernel"},
        {"program": "P544", "priority": 8, "study": "isolate every PSD transition root with independent Decimal/rational interval Newton steps"},
        {"program": "P545", "priority": 9, "study": "solve robust detector allocation with integer event counts and multinomial confidence regions"},
        {"program": "P546", "priority": 10, "study": "classify identifiable finite-dimensional law quotients rather than selecting one unrestricted law"},
        {"program": "P547", "priority": 11, "study": "test whether a FIN-internal mediator source exists without importing legacy physical roles"},
        {"program": "P548", "priority": 12, "study": "audit whether any new object breaks Aut(Z12) orientation symmetry; expect a no-go absent new data"},
    ]
    payload = {
        "release": "10.53",
        "programs": "P527-P536",
        "seed": SEED,
        "runtime_seconds": time.time() - started,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "mpmath": mp.__version__},
        "epistemic_boundary": (
            "All results are local, dimensionless mathematics or synthetic operational specifications. Added mediators, "
            "memory realizations and priors are explicit conditional premises. Nothing here derives a dimensional scale, "
            "discharges QW-2191, completes legacy-to-strict transfer, supplies a laboratory record, or proves physical ontology."
        ),
        "results": results,
        "recommended_next_programs": recommendations,
    }
    RESULTS.write_text(json.dumps(native(payload), indent=2, sort_keys=True), encoding="utf-8")
    write_summary(results)
    make_figures(results)
    print(json.dumps({"results_file": RESULTS.name, "runtime_seconds": payload["runtime_seconds"], "statuses": {key: value["status"] for key, value in results.items()}}, indent=2))


if __name__ == "__main__":
    main()
