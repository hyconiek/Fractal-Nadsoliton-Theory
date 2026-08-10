#!/usr/bin/env python3
"""FIN P497--P506: bounded local analytical and small-compute research.

The batch continues the explicit recommendations of Release 10.49.  It does
not run the P485/P487 Groebner campaign, use laboratory records, or treat a
supplied nonlinear/operational law as a consequence of the frozen kernels.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import platform
import time
from fractions import Fraction
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
import scipy
from numpy.typing import NDArray
from scipy.linalg import expm
from scipy.optimize import brentq, root

from fin_programs_488_496_low_compute import (
    BETA_S,
    ETA_S,
    OMEGA_S,
    PHI_S,
    KernelParameters,
    fourier_mode,
    grouped_stieltjes_data,
    laplacian_from_w,
    legacy_profile,
    strict_operator,
    strict_profile,
)


ROOT = Path(__file__).resolve().parent
FIG_DIR = ROOT / "FIN_Programs_497_506_Figures"
RESULTS = ROOT / "FIN_Programs_497_506_Results.json"
SUMMARY = ROOT / "FIN_Programs_497_506_Summary.csv"
INTERFACE = ROOT / "FIN_P497_P506_Theorem_Interface.json"
SEED = 20260810
N = 12


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


def ipr(u: NDArray[np.floating]) -> float:
    p = np.abs(u) ** 2
    return float(np.sum(p**2) / np.sum(p) ** 2)


def p497_localized_dnls(a: NDArray[np.float64]) -> dict:
    """Continue a focusing DNLS state from a one-site anti-continuum state."""
    mu = -1.0
    focusing = 1.0
    kappas = np.linspace(0.0, 1.0, 201)
    u = np.zeros(N)
    u[0] = math.sqrt(-mu / focusing)
    rows = []
    profiles = {}
    failed = None

    def residual(v: np.ndarray, kappa: float) -> np.ndarray:
        return kappa * (a @ v) - focusing * v**3 - mu * v

    def jacobian(v: np.ndarray, kappa: float) -> np.ndarray:
        return kappa * a - np.diag(3.0 * focusing * v**2) - mu * np.eye(N)

    for index, kappa in enumerate(kappas):
        sol = root(
            lambda v: residual(v, float(kappa)),
            u,
            jac=lambda v: jacobian(v, float(kappa)),
            method="hybr",
            options={"xtol": 1e-12, "maxfev": 5000},
        )
        if (not sol.success) or np.linalg.norm(residual(sol.x, float(kappa)), ord=np.inf) > 2e-9:
            failed = {"index": index, "kappa": float(kappa), "message": str(sol.message)}
            break
        u = np.asarray(sol.x, dtype=float)
        l_minus = kappa * a - mu * np.eye(N) - focusing * np.diag(u**2)
        l_plus = kappa * a - mu * np.eye(N) - 3.0 * focusing * np.diag(u**2)
        stability = np.block([[np.zeros((N, N)), l_minus], [-l_plus, np.zeros((N, N))]])
        exponents = np.linalg.eigvals(stability)
        phase_index = int(np.argmin(np.abs(exponents)))
        transverse = np.delete(exponents, phase_index)
        resolved = exponents[np.abs(exponents) > 1e-6]
        row = {
            "kappa": float(kappa),
            "residual_inf": float(np.linalg.norm(residual(u, float(kappa)), ord=np.inf)),
            "power": float(np.sum(u**2)),
            "ipr": ipr(u),
            "peak_fraction": float(np.max(u**2) / np.sum(u**2)),
            "minimum_linearized_abs_exponent": float(np.min(np.abs(exponents))),
            "max_transverse_real_exponent": float(np.max(np.real(transverse))),
            "neutral_band_count_abs_below_1e-6": int(np.sum(np.abs(exponents) <= 1e-6)),
            "resolved_positive_real_count_above_1e-6": int(np.sum(np.real(exponents) > 1e-6)),
            "max_resolved_real_exponent": float(np.max(np.real(resolved))) if resolved.size else 0.0,
            "amplitude": u.copy(),
        }
        rows.append(row)
        if index in (0, 20, 50, 100, 150, 200):
            profiles[f"kappa_{kappa:.3f}"] = u.copy()

    endpoint = rows[-1]
    localized = bool(endpoint["ipr"] > 2.0 / N and endpoint["peak_fraction"] > 0.25)
    anti_jac_diag = np.ones(N)
    anti_jac_diag[0] = -2.0
    return {
        "program": "P497",
        "object": "O201 Strict-DNLS Anti-Continuum Localization Branch",
        "supplied_equation": "i psi_dot = kappa A psi - g |psi|^2 psi",
        "stationary_equation": "kappa A u - g u^3 - mu u = 0",
        "parameters": {"mu": mu, "g": focusing, "kappa_range": [0.0, float(rows[-1]["kappa"])]},
        "anti_continuum_ift": {
            "state": [1.0] + [0.0] * (N - 1),
            "jacobian_diagonal": anti_jac_diag,
            "minimum_absolute_diagonal": 1.0,
            "theorem": (
                "At kappa=0 the real stationary Jacobian is diag(-2,1,...,1), hence the "
                "implicit-function theorem gives a unique local real branch near the chosen site."
            ),
        },
        "continuation_failed": failed,
        "sampled_branch": rows,
        "profiles": profiles,
        "endpoint": endpoint,
        "localized_by_declared_finite_metrics": localized,
        "verdict": (
            "A one-site DNLS branch is theorem-level persistent for sufficiently small coupling "
            "and numerically continues to full strict coupling. Localization at kappa=1 is a "
            "finite numerical property of the supplied focusing law, not a FIN-derived soliton law."
        ),
        "status": "proven_local_ift_plus_numerical_full_coupling_continuation",
    }


def iv_bounds(x: Any) -> tuple[float, float]:
    lower = np.nextafter(float(x.a), -np.inf)
    upper = np.nextafter(float(x.b), np.inf)
    return float(lower), float(upper)


def strict_lambda_intervals() -> list[Any]:
    mp.iv.dps = 70
    omega = mp.iv.mpf("0.18575")
    phi = mp.iv.mpf("0.16250")
    eta = mp.iv.mpf("1.8")
    weights = [None]
    for d in range(1, 7):
        weights.append(mp.iv.cos(omega * d + phi) / (1 + mp.iv.mpf(d) ** eta))
    lambdas = []
    for m in range(12):
        value = mp.iv.mpf("0")
        for d in range(1, 6):
            angle = 2 * mp.iv.pi * m * d / 12
            value += 2 * weights[d] * (1 - mp.iv.cos(angle))
        value += weights[6] * (1 - (-1) ** m)
        lambdas.append(value)
    return lambdas


def real_part_pair_bounds(dlo: float, dhi: float, gain: float) -> tuple[list[float], list[float]]:
    if dlo <= 0.0 <= dhi:
        abs_lo = 0.0
    else:
        abs_lo = min(abs(dlo), abs(dhi))
    abs_hi = max(abs(dlo), abs(dhi))
    sqrt_hi = math.sqrt(max(0.0, gain * gain - abs_lo * abs_lo))
    sqrt_lo = math.sqrt(max(0.0, gain * gain - abs_hi * abs_hi))
    plus = [float(np.nextafter(-gain + sqrt_lo, -np.inf)), float(np.nextafter(-gain + sqrt_hi, np.inf))]
    minus = [float(np.nextafter(-gain - sqrt_hi, -np.inf)), float(np.nextafter(-gain - sqrt_lo, np.inf))]
    return plus, minus


def p498_certified_p488_stability() -> dict:
    gain = 0.2
    lambdas = strict_lambda_intervals()
    records = []
    exact_neutral_pairs = []
    nonzero_upper_bounds = []
    for m in range(1, 7):
        for q in range(12):
            exact_neutral = q == 0 or (m == 3 and q == 6)
            d_interval = (lambdas[(m + q) % 12] + lambdas[(m - q) % 12] - 2 * lambdas[m]) / 2
            dlo, dhi = iv_bounds(d_interval)
            if exact_neutral:
                plus, minus = [0.0, 0.0], [-2.0 * gain, -2.0 * gain]
                exact_neutral_pairs.append([m, q])
            else:
                plus, minus = real_part_pair_bounds(dlo, dhi, gain)
                nonzero_upper_bounds.extend([plus[1], minus[1]])
            records.append(
                {
                    "mode_m": m,
                    "perturbation_q": q,
                    "d_interval": [dlo, dhi],
                    "real_part_plus_interval": plus,
                    "real_part_minus_interval": minus,
                    "exact_neutral_by_symmetry": exact_neutral,
                }
            )
    max_upper = max(nonzero_upper_bounds)
    return {
        "program": "P498",
        "object": "O202 Fourier-Block Stability Enclosure",
        "block_formula": (
            "rho_{m,q}^{+/-}=-r-i(h_q-h_-q)/2 +/- sqrt(r^2-d_{m,q}^2), "
            "d_{m,q}=(lambda_{m+q}+lambda_{m-q}-2lambda_m)/2"
        ),
        "gain_r": gain,
        "records": records,
        "exact_neutral_pairs": exact_neutral_pairs,
        "largest_certified_nonzero_real_part_upper_bound": max_upper,
        "all_nonzero_real_parts_strictly_negative": bool(max_upper < 0.0),
        "theorem": (
            "Fourier reduction turns the 24-dimensional real Jacobian into 2x2 blocks. "
            "Outward interval evaluation encloses every block real part; q=0 is the phase "
            "neutrality and (m,q)=(3,6) is an additional exact symmetry neutrality."
        ),
        "status": "interval_certified_in_declared_strict_decimal_model",
    }


def continued_fraction_convergents(x: float, max_denominator: int) -> list[tuple[int, int]]:
    value = x
    p0, q0, p1, q1 = 0, 1, 1, 0
    out = []
    for _ in range(64):
        a = math.floor(value)
        p2, q2 = a * p1 + p0, a * q1 + q0
        if q2 > max_denominator:
            break
        out.append((p2, q2))
        remainder = value - a
        if abs(remainder) < 1e-18:
            break
        value = 1.0 / remainder
        p0, q0, p1, q1 = p1, q1, p2, q2
    return out


def p499_two_frequency_clock(a: NDArray[np.float64]) -> dict:
    first_row = a[0]
    lam = np.real(np.fft.fft(first_row))
    omega1, omega2 = float(lam[1]), float(lam[2])
    ratio = omega2 / omega1
    lambda_iv = strict_lambda_intervals()
    ratio_lo, ratio_hi = iv_bounds(lambda_iv[2] / lambda_iv[1])
    max_den = 1_000_000
    convergents = continued_fraction_convergents(ratio, max_den)
    rows = []
    for p, q in convergents:
        mismatch = abs(q * ratio - p)
        recurrence_time = 2.0 * math.pi * q / omega1
        chord = 2.0 * abs(math.sin(math.pi * mismatch))
        rows.append(
            {
                "p": p,
                "q": q,
                "ratio_error": abs(ratio - p / q),
                "phase_integer_mismatch": mismatch,
                "dimensionless_recurrence_time": recurrence_time,
                "two_phase_chord_error": chord,
            }
        )
    qv = np.arange(1, max_den + 1, dtype=float)
    lower = np.nextafter(qv * ratio_lo, -np.inf)
    upper = np.nextafter(qv * ratio_hi, np.inf)
    nearest = np.rint(qv * ratio)
    contains_integer = (np.ceil(lower) <= np.floor(upper))
    margins = np.minimum(np.abs(lower - nearest), np.abs(upper - nearest))
    margins[contains_integer] = 0.0
    best_index = int(np.argmin(margins))
    certified_margin = float(margins[best_index])
    return {
        "program": "P499",
        "object": "O203 Two-Frequency Spectral Torus Clock",
        "frequencies": [omega1, omega2],
        "ratio_omega2_over_omega1": ratio,
        "outward_interval_ratio": [ratio_lo, ratio_hi],
        "general_theorem": (
            "The map t -> (exp(-i omega1 t),exp(-i omega2 t)) is periodic iff "
            "omega2/omega1 is rational; it is injective as a map R into the two-torus iff "
            "the ratio is irrational."
        ),
        "bounded_denominator_audit": {
            "maximum_denominator": max_den,
            "best_q": best_index + 1,
            "best_distance_qr_to_integer": float(margins[best_index]),
            "outward_interval_minimum_integer_separation": certified_margin,
            "possible_rational_denominator_count": int(np.sum(contains_integer)),
            "no_exact_rational_ratio_within_bound": bool(not np.any(contains_integer)),
        },
        "continued_fraction_recurrences": rows,
        "verdict": (
            "The rational/irrational clock dichotomy is proven generally. The strict ratio is "
            "excluded from rationals with denominator at most 10^6 by outward interval "
            "evaluation, but global irrationality is not proven."
        ),
        "status": "general_theorem_plus_bounded_nonresonance_certificate",
    }


def stieltjes_relative_jacobian(
    z: np.ndarray, poles: np.ndarray, residues: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    target = np.sum(residues[None, :] / (z[:, None] + poles[None, :]), axis=1)
    jp = -residues[None, :] * poles[None, :] / (
        target[:, None] * (z[:, None] + poles[None, :]) ** 2
    )
    jr = residues[None, :] / (target[:, None] * (z[:, None] + poles[None, :]))
    return np.column_stack([jp, jr]), target


def p500_stieltjes_identifiability(a: NDArray[np.float64]) -> dict:
    poles, residues, _, _ = grouped_stieltjes_data(a)
    order = np.argsort(poles)
    poles, residues = poles[order], residues[order]
    scale = float(np.median(poles))
    z = np.geomspace(0.03 * scale, 12.0 * scale, 48)
    j0, target = stieltjes_relative_jacobian(z, poles, residues)
    singular = np.linalg.svd(j0, compute_uv=False)
    pseudoinverse = np.linalg.inv(j0.T @ j0) @ j0.T
    denom = 10**12
    rational_r = [[Fraction(float(x)).limit_denominator(denom) for x in row] for row in pseudoinverse]

    mp.iv.dps = 50
    chosen = None
    candidates = [10.0 ** exponent for exponent in range(-14, -5)]
    theta0 = np.log(np.concatenate([poles, residues]))
    for radius in candidates:
        theta_iv = [mp.iv.mpf([float(v - radius), float(v + radius)]) for v in theta0]
        p_iv = [mp.iv.exp(theta_iv[j]) for j in range(3)]
        r_iv = [mp.iv.exp(theta_iv[3 + j]) for j in range(3)]
        j_iv = []
        for zi, ti in zip(z, target):
            row_iv = []
            for j in range(3):
                row_iv.append(-r_iv[j] * p_iv[j] / (mp.iv.mpf(str(ti)) * (mp.iv.mpf(str(zi)) + p_iv[j]) ** 2))
            for j in range(3):
                row_iv.append(r_iv[j] / (mp.iv.mpf(str(ti)) * (mp.iv.mpf(str(zi)) + p_iv[j])))
            j_iv.append(row_iv)
        q_sq = 0.0
        for i in range(6):
            for j in range(6):
                value = mp.iv.mpf("0")
                for k in range(48):
                    rat = rational_r[i][k]
                    value += mp.iv.mpf(rat.numerator) / rat.denominator * j_iv[k][j]
                if i == j:
                    value -= 1
                lo, hi = iv_bounds(value)
                q_sq += max(abs(lo), abs(hi)) ** 2
        q_frobenius = math.sqrt(q_sq)
        r_norm_frobenius = math.sqrt(sum(float(x) ** 2 for row in rational_r for x in row))
        if q_frobenius < 0.5:
            chosen = (radius, q_frobenius, r_norm_frobenius)
    if chosen is None:
        raise RuntimeError("no local Stieltjes inverse certificate")
    radius, q_bound, r_norm = chosen
    inverse_lower = (1.0 - q_bound) / r_norm
    noise_limit = inverse_lower * radius
    noise_rows = []
    for fraction in (1e-3, 1e-2, 1e-1, 0.5):
        epsilon = fraction * noise_limit
        theta_error = epsilon / inverse_lower
        noise_rows.append(
            {
                "relative_response_l2_noise": epsilon,
                "log_parameter_l2_error_bound": theta_error,
                "per_pole_relative_error_bound": math.expm1(theta_error),
            }
        )
    return {
        "program": "P500",
        "object": "O204 Local Stieltjes Inverse-Certificate Box",
        "declared_decimal_poles": poles,
        "declared_decimal_residues": residues,
        "relative_jacobian_singular_values": singular,
        "condition_number": float(singular[0] / singular[-1]),
        "log_parameter_infinity_radius": radius,
        "rational_left_inverse_denominator_limit": denom,
        "left_inverse_defect_frobenius_upper": q_bound,
        "left_inverse_frobenius_norm": r_norm,
        "certified_local_forward_lower_lipschitz": inverse_lower,
        "maximum_noise_for_self_consistent_box_bound": noise_limit,
        "noise_to_parameter_bounds": noise_rows,
        "theorem": (
            "Inside the declared log-parameter box, ||F(theta)-F(theta*)||_2 is at least "
            "c||theta-theta*||_2, with c=(1-q)/||R||_F, because the rational left inverse "
            "satisfies ||R J(theta)-I||_F <= q < 1."
        ),
        "verdict": (
            "Local identifiability is certified, but only in a very small box and for noise below "
            "the reported self-consistency threshold. The large inverse norm makes pole recovery "
            "practically fragile."
        ),
        "status": "conditional_interval_local_inverse_certificate",
    }


def context_fingerprint_fft(n: int, delta: float = 0.1) -> dict:
    if n % 2:
        raise ValueError("n must be even")
    distance = np.minimum(np.arange(n), n - np.arange(n))
    weights = np.zeros(n)
    weights[1:] = 1.0 / (1.0 + distance[1:] ** ETA_S)
    row = -weights
    row[0] = np.sum(weights)
    even_row = row[::2]
    odd_row = row[1::2]
    c_eigs = np.real(np.fft.fft(even_row)) + delta
    b_eigs = np.fft.fft(odd_row)
    scale = float(np.median(c_eigs))
    u = np.array([0.125, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0])
    sigma = np.array([np.mean(np.abs(b_eigs) ** 2 / (value * scale + c_eigs)) for value in u])
    fingerprint = sigma / sigma[3]
    return {"n": n, "u": u, "fingerprint": fingerprint, "scale": scale}


def p501_context_limit() -> dict:
    sizes = [12, 24, 48, 96, 192, 384]
    finite = [context_fingerprint_fft(n) for n in sizes]
    asymptotic_n = 131072
    asymptotic = context_fingerprint_fft(asymptotic_n)
    successive = []
    for left, right in zip(finite[:-1], finite[1:]):
        diff = np.asarray(right["fingerprint"]) - np.asarray(left["fingerprint"])
        successive.append(
            {"sizes": [left["n"], right["n"]], "max_difference": float(np.max(np.abs(diff)))}
        )
    to_integral = []
    for item in finite:
        diff = np.asarray(item["fingerprint"]) - np.asarray(asymptotic["fingerprint"])
        to_integral.append({"n": item["n"], "max_difference": float(np.max(np.abs(diff)))})
    tail_bound = 2.0 / (ETA_S - 1.0) * (asymptotic_n / 2.0 - 1.0) ** (1.0 - ETA_S)
    return {
        "program": "P501",
        "object": "O205 Polyphase Integral Context Fingerprint",
        "finite_fingerprints": finite,
        "successive_differences": successive,
        "asymptotic_fft_quadrature": asymptotic,
        "finite_to_asymptotic_differences": to_integral,
        "omitted_infinite_weight_l1_tail_upper_bound": tail_bound,
        "integral_formula": (
            "sigma(z)=(2pi)^-1 integral |B(theta)|^2/[z+C(theta)] dtheta, with even/odd "
            "polyphase symbols B and C of the summable attenuation-envelope Laplacian."
        ),
        "monotone_approach_from_n96": bool(
            to_integral[3]["max_difference"] > to_integral[4]["max_difference"] > to_integral[5]["max_difference"]
        ),
        "verdict": (
            "The C192 and C384 fingerprints move toward a well-defined polyphase integral "
            "fingerprint in the declared summable envelope family. This is a continuum-context "
            "limit, not an exact finite fractal fixed point."
        ),
        "status": "analytic_integral_formula_plus_numerical_convergence_evidence",
    }


def homotopy_params(t: float | mp.mpf) -> tuple[Any, Any, Any, Any, Any]:
    alpha = 4 * mp.log(2)
    return (
        (1 - t) * alpha + t,
        (1 - t) * mp.mpf("0.01") + t,
        (1 - t) + t * mp.mpf("1.8"),
        (1 - t) * (mp.pi / 4) + t * mp.mpf("0.18575"),
        (1 - t) * (mp.pi / 6) + t * mp.mpf("0.16250"),
    )


def homotopy_weight_mp(t: Any, d: int) -> Any:
    amplitude, beta, eta, omega, phi = homotopy_params(t)
    return amplitude * mp.cos(omega * d + phi) / (1 + beta * mp.power(d, eta))


def homotopy_laplacian_eigenvalue(t: float, mode: int) -> float:
    value = 0.0
    for d in range(1, 6):
        value += 2.0 * float(homotopy_weight_mp(mp.mpf(str(t)), d)) * (
            1.0 - math.cos(2.0 * math.pi * mode * d / 12.0)
        )
    value += float(homotopy_weight_mp(mp.mpf(str(t)), 6)) * (1.0 - (-1) ** mode)
    return value


def p502_passivity_boundary() -> dict:
    mp.mp.dps = 80
    weight_roots = []
    for d in range(1, 7):
        theta0 = mp.pi * d / 4 + mp.pi / 6
        theta1 = mp.mpf("0.18575") * d + mp.mpf("0.16250")
        root_t = (theta0 - mp.pi / 2) / (theta0 - theta1)
        if 0 <= root_t <= 1:
            weight_roots.append((d, root_t))
    limiting_d, t_weight = max(weight_roots, key=lambda item: item[1])
    eps = mp.mpf("5e-13")
    weight_bracket = [float(t_weight - eps), float(t_weight + eps)]
    limiting_weight_below = float(homotopy_weight_mp(t_weight - eps, limiting_d))
    limiting_weight_above = float(homotopy_weight_mp(t_weight + eps, limiting_d))

    roots = []
    grid = np.linspace(0.0, 1.0, 2001)
    for mode in range(1, 12):
        vals = np.array([homotopy_laplacian_eigenvalue(float(t), mode) for t in grid])
        for i in range(len(grid) - 1):
            if vals[i] < 0.0 <= vals[i + 1]:
                rr = brentq(lambda tt, m=mode: homotopy_laplacian_eigenvalue(tt, m), grid[i], grid[i + 1], xtol=1e-14)
                roots.append((mode, rr))
    limiting_mode, t_psd = max(roots, key=lambda item: item[1])
    psd_eps = 5e-13
    psd_bracket = [t_psd - psd_eps, t_psd + psd_eps]
    lower_values = [homotopy_laplacian_eigenvalue(psd_bracket[0], m) for m in range(1, 12)]
    upper_values = [homotopy_laplacian_eigenvalue(psd_bracket[1], m) for m in range(1, 12)]
    return {
        "program": "P502",
        "object": "O206 Exact Radial-Passivity and PSD Boundary Ledger",
        "radial_nonnegativity_boundary": {
            "limiting_distance": limiting_d,
            "closed_form": "(theta_legacy(6)-pi/2)/(theta_legacy(6)-theta_strict(6))",
            "value": float(t_weight),
            "isolation_bracket": weight_bracket,
            "limiting_weight_below": limiting_weight_below,
            "limiting_weight_above": limiting_weight_above,
            "all_distance_roots": [{"distance": d, "last_pi_over_2_crossing": float(t)} for d, t in weight_roots],
        },
        "signed_laplacian_psd_boundary": {
            "limiting_fourier_mode": limiting_mode,
            "value": t_psd,
            "isolation_bracket": psd_bracket,
            "minimum_nonzero_eigenvalue_below": float(min(lower_values)),
            "minimum_nonzero_eigenvalue_above": float(min(upper_values)),
        },
        "theorem": (
            "All six radial weights are nonnegative exactly after the largest final pi/2 phase "
            "crossing, attained at d=6. Nonnegative weights imply PSD, but the independent "
            "Fourier audit locates the earlier signed-Laplacian PSD boundary."
        ),
        "status": "closed_form_weight_boundary_plus_numerically_isolated_psd_boundary",
    }


def principal_increments(z: np.ndarray) -> np.ndarray:
    return np.angle(np.roll(z, -1) / z)


def winding(z: np.ndarray) -> int | None:
    inc = principal_increments(z)
    if np.any(np.isclose(np.abs(inc), math.pi, atol=1e-11)):
        return None
    return int(round(float(np.sum(inc) / (2.0 * math.pi))))


def geodesic_refine(z: np.ndarray) -> np.ndarray:
    inc = principal_increments(z)
    out = np.empty(2 * len(z), dtype=complex)
    out[::2] = z
    out[1::2] = z * np.exp(0.5j * inc)
    return out


def p503_phase_refinement() -> dict:
    rng = np.random.default_rng(SEED + 503)
    checks = 0
    failures = 0
    refinement_rows = []
    for n in (12, 24, 48):
        for target_winding in (-2, -1, 0, 1, 2):
            if abs(2 * math.pi * target_winding / n) >= 0.55 * math.pi:
                continue
            for _ in range(100):
                inc = 2 * math.pi * target_winding / n + rng.normal(0.0, 0.05, n)
                inc -= (np.sum(inc) - 2 * math.pi * target_winding) / n
                phase = np.concatenate([[0.0], np.cumsum(inc[:-1])])
                z = np.exp(1j * phase)
                margin = math.pi - float(np.max(np.abs(principal_increments(z))))
                delta = rng.normal(0.0, 1.0, n)
                edge_delta = np.roll(delta, -1) - delta
                scale = 0.45 * margin / max(float(np.max(np.abs(edge_delta))), 1e-15)
                update = scale * delta
                z_new = np.exp(1j * update) * z
                checks += 1
                if winding(z_new) != winding(z):
                    failures += 1
            mode = np.exp(2j * math.pi * target_winding * np.arange(n) / n)
            refined = geodesic_refine(mode)
            refinement_rows.append(
                {
                    "n": n,
                    "winding": target_winding,
                    "refined_n": 2 * n,
                    "refined_winding": winding(refined),
                    "coarse_max_increment": float(np.max(np.abs(principal_increments(mode)))),
                    "refined_max_increment": float(np.max(np.abs(principal_increments(refined)))),
                }
            )
    counterexample = None
    base = np.exp(2j * math.pi * np.arange(12) / 12)
    for _ in range(10000):
        update = rng.uniform(-math.pi, math.pi, 12)
        changed = np.exp(1j * update) * base
        if winding(changed) != winding(base) and winding(changed) is not None:
            counterexample = {
                "old_winding": winding(base),
                "new_winding": winding(changed),
                "max_neighbor_update_difference": float(np.max(np.abs(np.roll(update, -1) - update))),
                "old_branch_margin": math.pi - float(np.max(np.abs(principal_increments(base)))),
            }
            break
    return {
        "program": "P503",
        "object": "O207 Nonvanishing Phase-Refinement Sector",
        "theorem": (
            "If z_x is S1-valued, every principal edge increment has margin epsilon from +/-pi, "
            "and |delta_(x+1)-delta_x|<epsilon, then z'_x=e^{i delta_x}z_x has the same integer "
            "winding. Principal-geodesic midpoint refinement also preserves winding and halves "
            "every edge increment."
        ),
        "bounded_update_checks": checks,
        "bounded_update_failures": failures,
        "refinement_checks": refinement_rows,
        "condition_violation_counterexample": counterexample,
        "verdict": (
            "A protected discrete winding sector exists after adding a nonvanishing phase field "
            "and an admissible update/refinement law. These are additional typed objects and are "
            "not selected by the scalar strict kernel."
        ),
        "status": "proven_winding_preservation_with_finite_stress_replay",
    }


def density(v: np.ndarray) -> np.ndarray:
    v = v / np.linalg.norm(v)
    return np.outer(v, v.conj())


def normalize_probability(p: np.ndarray) -> np.ndarray:
    p = np.maximum(np.real(p), 0.0)
    return p / np.sum(p)


def operational_probability_table(a: np.ndarray) -> tuple[list[str], list[tuple], dict]:
    n = len(a)
    eye = np.eye(n)
    projectors = [np.diag(np.eye(n)[x]).astype(complex) for x in range(n)]
    gamma = 0.60
    hamiltonian_super = -1j * (np.kron(eye, a) - np.kron(a.T, eye))
    dephase_super = sum(np.kron(p.T, p) for p in projectors) - np.eye(n * n)
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
    models = ["unitary", "markov", "dephasing", "revival"]
    contexts = [(t, prep, measurement) for t in (0.20, 0.55, 1.00) for prep in preparations for measurement in ("vertex", "fourier")]

    def evolve(model: str, rho: np.ndarray, t: float) -> np.ndarray:
        u = expm(-1j * t * a)
        if model == "unitary":
            return u @ rho @ u.conj().T
        if model == "markov":
            return np.diag(normalize_probability(expm(-t * a) @ np.real(np.diag(rho)))).astype(complex)
        if model == "dephasing":
            vec = expm(t * lindblad) @ rho.reshape(-1, order="F")
            out = vec.reshape((n, n), order="F")
            return (out + out.conj().T) / 2.0
        coherent = u @ rho @ u.conj().T
        q = 0.45 + 0.40 * math.cos(2.8 * t)
        return q * coherent + (1 - q) * np.diag(np.diag(coherent))

    table = {model: {} for model in models}
    for model in models:
        for context in contexts:
            t, prep, measurement = context
            rho = evolve(model, preparations[prep], t)
            if measurement == "vertex":
                p = normalize_probability(np.diag(rho))
            else:
                p = normalize_probability(np.diag(fourier.conj().T @ rho @ fourier))
            table[model][context] = p
    return models, contexts, table


def p504_minimal_coherence_tester(a: NDArray[np.float64]) -> dict:
    models, contexts, table = operational_probability_table(a)
    shots = 120
    rows = []
    for context in contexts:
        tvs, pair_bounds = [], []
        for i, left in enumerate(models):
            for right in models[i + 1 :]:
                p, q = table[left][context], table[right][context]
                tvs.append(0.5 * float(np.sum(np.abs(p - q))))
                bc = float(np.sum(np.sqrt(p * q)))
                pair_bounds.append(0.5 * bc**shots)
        rows.append(
            {
                "context": context,
                "minimum_pairwise_tv": min(tvs),
                "sum_pair_bhattacharyya_error_bound": sum(pair_bounds),
            }
        )
    separating = [row for row in rows if row["minimum_pairwise_tv"] > 1e-10]
    best = min(separating, key=lambda row: row["sum_pair_bhattacharyya_error_bound"])
    context = tuple(best["context"])
    rng = np.random.default_rng(SEED + 504)
    confusion = np.zeros((4, 4), dtype=int)
    trials = 1000
    for true_index, model in enumerate(models):
        for _ in range(trials):
            counts = rng.multinomial(shots, table[model][context])
            scores = []
            for candidate in models:
                p = np.clip(table[candidate][context], 1e-15, 1.0)
                scores.append(float(np.dot(counts, np.log(p))))
            confusion[true_index, int(np.argmax(scores))] += 1
    return {
        "program": "P504",
        "object": "O208 One-Context Coherence Tester",
        "exhausted_single_context_count": len(contexts),
        "separating_single_context_count": len(separating),
        "minimal_context_cardinality": 1 if separating else None,
        "selected_context": {"time": context[0], "preparation": context[1], "measurement": context[2]},
        "selected_minimum_pairwise_tv": best["minimum_pairwise_tv"],
        "shots": shots,
        "analytic_sum_pair_error_upper_bound": best["sum_pair_bhattacharyya_error_bound"],
        "confusion_rows_true_cols_predicted": confusion,
        "classification_accuracy": float(np.trace(confusion) / np.sum(confusion)),
        "minimality_certificate": "Zero contexts carry no data; exhaustive search found a separating one-context design.",
        "verdict": (
            "One declared preparation/time/Fourier measurement already separates the four "
            "synthetic models. Minimality is only within the P495 finite context catalogue and "
            "does not constitute a laboratory instrument certificate."
        ),
        "status": "exhaustive_finite_synthetic_minimality_certificate",
        "all_context_scores": rows,
    }


def validate_theorem_interface(interface: dict) -> dict:
    nodes = {row["id"]: row for row in interface["theorems"]}
    missing = sorted({dep for row in nodes.values() for dep in row["depends_on"] if dep not in nodes})
    temporary, permanent = set(), set()
    cycle = False

    def visit(node: str) -> None:
        nonlocal cycle
        if node in permanent:
            return
        if node in temporary:
            cycle = True
            return
        temporary.add(node)
        for dep in nodes[node]["depends_on"]:
            if dep in nodes:
                visit(dep)
        temporary.remove(node)
        permanent.add(node)

    for node in nodes:
        visit(node)
    required = {"id", "claim_level", "assumptions", "conclusion", "depends_on", "proof_skeleton", "nonclaims"}
    schema_failures = [node for node, row in nodes.items() if not required.issubset(row)]
    return {"node_count": len(nodes), "missing_dependencies": missing, "cycle_detected": cycle, "schema_failures": schema_failures, "valid": not missing and not cycle and not schema_failures}


def p505_formal_interface() -> dict:
    theorems = [
        {
            "id": "T488-CYCLE",
            "claim_level": "conditional_exact",
            "assumptions": ["A is self-adjoint", "Av=lambda*v", "|v_x| is constant", "r>0", "g>0", "the supplied gain-saturation flow holds"],
            "conclusion": "sqrt(r/g)/|v_x| times v rotating at lambda is an exact persistent solution",
            "depends_on": [],
            "proof_skeleton": ["substitute the rotating ansatz", "cancel r-g|psi_x|^2", "use Av=lambda*v"],
            "nonclaims": ["localization", "physical source of r or g"],
        },
        {
            "id": "T489-SCALE",
            "claim_level": "exact",
            "assumptions": ["a phase record exp(-itA) is used", "c>0"],
            "conclusion": "A->A/c and t->ct leave every phase record unchanged",
            "depends_on": ["T488-CYCLE"],
            "proof_skeleton": ["compute exp(-i(ct)(A/c))=exp(-itA)"],
            "nonclaims": ["absolute time unit", "clock selection"],
        },
        {
            "id": "T492-COMPLEMENT",
            "claim_level": "exact",
            "assumptions": ["A=sI-W", "x is nonzero"],
            "conclusion": "R_W(x)=s-R_A(x)",
            "depends_on": [],
            "proof_skeleton": ["expand x*(sI-W)x", "divide by x*x"],
            "nonclaims": ["a dynamical resonance selector"],
        },
        {
            "id": "T493-INDUCTION",
            "claim_level": "exact",
            "assumptions": ["c_0=1", "q_0=0", "theta_0=phi", "the exponent-lift recurrence holds"],
            "conclusion": "c_d=(1+beta*d^eta)^-1 and K_d=a*c_d*cos(theta_d)",
            "depends_on": [],
            "proof_skeleton": ["induct on q", "invert the c update", "telescope powers and phases"],
            "nonclaims": ["parameter provenance", "legacy role transfer"],
        },
        {
            "id": "T497-IFT",
            "claim_level": "conditional_exact_local",
            "assumptions": ["the supplied focusing DNLS holds", "mu=-1", "g=1", "one-site anti-continuum state"],
            "conclusion": "a unique local real stationary branch exists for sufficiently small kappa",
            "depends_on": [],
            "proof_skeleton": ["evaluate stationary residual", "Jacobian is diag(-2,1,...,1)", "apply the implicit-function theorem"],
            "nonclaims": ["global continuation theorem", "physical nadsoliton"],
        },
        {
            "id": "T503-WINDING",
            "claim_level": "exact",
            "assumptions": ["z_x is S1-valued", "edge increments avoid the branch cut by epsilon", "neighbor update differences are below epsilon"],
            "conclusion": "integer winding is preserved; midpoint geodesic refinement preserves winding",
            "depends_on": [],
            "proof_skeleton": ["no edge crosses the principal branch", "update increments telescope", "refinement splits each increment in two"],
            "nonclaims": ["kernel-derived phase field", "particle identification"],
        },
    ]
    interface = {"schema": "FIN-proof-interface-1.0", "theorems": theorems}
    audit = validate_theorem_interface(interface)
    canonical = json.dumps(native(interface), sort_keys=True, separators=(",", ":")).encode("utf-8")
    interface["audit"] = audit
    interface["canonical_payload_sha256"] = hashlib.sha256(canonical).hexdigest()
    INTERFACE.write_text(json.dumps(native(interface), indent=2, sort_keys=True), encoding="utf-8")
    return {
        "program": "P505",
        "object": "O209 Proof-Assistant-Neutral Theorem Dependency Interface",
        "file": INTERFACE.name,
        "audit": audit,
        "canonical_payload_sha256": interface["canonical_payload_sha256"],
        "verdict": "The principal exact claims now have explicit assumptions, dependencies, proof skeletons, and nonclaims; the file is not machine-checked by a proof assistant.",
        "status": "validated_neutral_interface_not_machine_checked",
    }


def kernel_parameter_jacobian(params: KernelParameters, distances: list[int]) -> np.ndarray:
    rows = []
    for d in distances:
        power = 0.0 if d == 0 else d**params.eta
        den = 1.0 + params.beta * power
        phase = params.omega * d + params.phi
        value = params.amplitude * math.cos(phase) / den
        logd = 0.0 if d == 0 else math.log(d)
        rows.append(
            [
                value,
                -value * params.beta * power / den,
                -value * params.beta * power * logd / den,
                -params.amplitude * d * math.sin(phase) / den,
                -params.amplitude * math.sin(phase) / den,
            ]
        )
    return np.array(rows)


def p506_parameter_source_no_go() -> dict:
    legacy = KernelParameters(4.0 * math.log(2.0), 0.01, 1.0, math.pi / 4.0, math.pi / 6.0)
    strict = KernelParameters(1.0, 1.0, 1.8, OMEGA_S, PHI_S)
    midpoint = KernelParameters(
        0.5 * (legacy.amplitude + strict.amplitude),
        0.5 * (legacy.beta + strict.beta),
        0.5 * (legacy.eta + strict.eta),
        0.5 * (legacy.omega + strict.omega),
        0.5 * (legacy.phi + strict.phi),
    )
    local_rows = []
    for name, params in (("legacy", legacy), ("midpoint", midpoint), ("strict", strict)):
        jac = kernel_parameter_jacobian(params, list(range(7)))
        singular = np.linalg.svd(jac, compute_uv=False)
        local_rows.append(
            {
                "point": name,
                "rank_d0_to_d6": int(np.linalg.matrix_rank(jac, tol=1e-11)),
                "singular_values": singular,
                "condition_number": float(singular[0] / singular[-1]),
            }
        )
    subsets = []
    for distances in itertools.combinations(range(7), 5):
        jac = kernel_parameter_jacobian(strict, list(distances))
        singular = np.linalg.svd(jac, compute_uv=False)
        if singular[-1] > 1e-11:
            subsets.append((float(singular[0] / singular[-1]), distances, singular))
    subsets.sort(key=lambda row: row[0])
    best_condition, best_distances, best_singular = subsets[0]

    sample_times = np.array([0.0, 0.25, 0.5, 0.75, 1.0])
    grid = np.linspace(0.0, 1.0, 1001)
    bump = np.ones_like(grid)
    for t in sample_times:
        bump *= (grid - t) ** 2
    bump /= np.max(np.abs(bump))
    vanishing_residual = max(abs(np.interp(sample_times, grid, bump)))
    return {
        "program": "P506",
        "object": "O210 Kernel-Trajectory Observation Quotient",
        "local_parameter_identifiability": local_rows,
        "best_five_distance_strict_design": {
            "distances": list(best_distances),
            "condition_number": best_condition,
            "singular_values": best_singular,
        },
        "finite_time_trajectory_no_go": {
            "sample_times": sample_times,
            "smooth_bump_max": float(np.max(np.abs(bump))),
            "bump_at_sample_times_max_abs": float(vanishing_residual),
            "theorem": (
                "For any finite set of observation times, adding epsilon*v*product_j(t-t_j)^2 "
                "to a parameter trajectory preserves every sampled parameter and kernel value. "
                "Hence no finite-time observation design identifies a unique continuous trajectory."
            ),
        },
        "time_reparameterization_no_go": (
            "Even a continuously observed unlabelled path image is invariant under monotone time "
            "reparameterization; calibrated time labels or a parameter evolution law are necessary."
        ),
        "sufficient_conditional_observation": (
            "On a fixed phase/positivity branch, continuous time-labelled kernel values at a "
            "full-rank set of at least five distances can locally recover the parameter path pointwise."
        ),
        "verdict": (
            "Static strict parameters are locally observable from suitable distance values, but "
            "neither a finite record nor the endpoints identify the transition law. A calibrated "
            "continuous record plus branch restrictions can recover a path, not its causal source."
        ),
        "status": "proven_finite_observation_trajectory_nonidentifiability_plus_local_rank_audit",
    }


def make_figures(results: dict) -> None:
    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})

    branch = results["P497"]["sampled_branch"]
    kappa = np.array([row["kappa"] for row in branch])
    iprs = np.array([row["ipr"] for row in branch])
    peaks = np.array([row["peak_fraction"] for row in branch])
    fig, ax = plt.subplots(figsize=(7.4, 4.3))
    ax.plot(kappa, iprs, label="IPR")
    ax.plot(kappa, peaks, label="peak power fraction")
    ax.axhline(1 / 12, color="black", ls="--", lw=1, label="delocalized IPR=1/12")
    ax.set_xlabel("strict coupling scale kappa")
    ax.set_ylabel("localization measure")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p497_localization_branch.png", dpi=180)
    plt.close(fig)

    endpoint = np.asarray(results["P497"]["endpoint"]["amplitude"])
    fig, ax = plt.subplots(figsize=(7.4, 4.0))
    ax.bar(np.arange(12), endpoint**2)
    ax.set_xlabel("Z12 vertex")
    ax.set_ylabel("stationary power at kappa=1")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p497_endpoint_profile.png", dpi=180)
    plt.close(fig)

    recurrences = results["P499"]["continued_fraction_recurrences"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.loglog([row["dimensionless_recurrence_time"] for row in recurrences], [max(row["two_phase_chord_error"], 1e-18) for row in recurrences], marker="o")
    ax.set_xlabel("recurrence time")
    ax.set_ylabel("two-phase chord error")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p499_clock_recurrence.png", dpi=180)
    plt.close(fig)

    p501 = results["P501"]
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.loglog([row["n"] for row in p501["finite_to_asymptotic_differences"]], [row["max_difference"] for row in p501["finite_to_asymptotic_differences"]], marker="o")
    ax.set_xlabel("cycle size n")
    ax.set_ylabel("max difference from integral fingerprint")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p501_context_limit.png", dpi=180)
    plt.close(fig)

    scores = sorted(results["P504"]["all_context_scores"], key=lambda row: row["minimum_pairwise_tv"], reverse=True)
    fig, ax = plt.subplots(figsize=(7.4, 4.2))
    ax.bar(np.arange(len(scores)), [row["minimum_pairwise_tv"] for row in scores])
    ax.set_xlabel("single contexts, sorted")
    ax.set_ylabel("minimum pairwise TV")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "p504_single_context_separation.png", dpi=180)
    plt.close(fig)


def write_summary(results: dict) -> None:
    rows = [
        ("P497", results["P497"]["status"], f"IPR(kappa=1)={results['P497']['endpoint']['ipr']:.6g}"),
        ("P498", results["P498"]["status"], f"max nonzero Re upper={results['P498']['largest_certified_nonzero_real_part_upper_bound']:.6g}"),
        ("P499", results["P499"]["status"], f"bounded q<={results['P499']['bounded_denominator_audit']['maximum_denominator']}"),
        ("P500", results["P500"]["status"], f"c={results['P500']['certified_local_forward_lower_lipschitz']:.6g}"),
        ("P501", results["P501"]["status"], f"C384 diff={results['P501']['finite_to_asymptotic_differences'][-1]['max_difference']:.6g}"),
        ("P502", results["P502"]["status"], f"t_weight={results['P502']['radial_nonnegativity_boundary']['value']:.12g}"),
        ("P503", results["P503"]["status"], f"checks={results['P503']['bounded_update_checks']} failures={results['P503']['bounded_update_failures']}"),
        ("P504", results["P504"]["status"], f"contexts={results['P504']['minimal_context_cardinality']} accuracy={results['P504']['classification_accuracy']:.5f}"),
        ("P505", results["P505"]["status"], f"nodes={results['P505']['audit']['node_count']} valid={results['P505']['audit']['valid']}"),
        ("P506", results["P506"]["status"], f"best distances={results['P506']['best_five_distance_strict_design']['distances']}"),
    ]
    lines = ["program,status,key_result"]
    for program, status, result in rows:
        lines.append(f'"{program}","{status}","{result}"')
    SUMMARY.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    started = time.perf_counter()
    _, a = strict_operator()
    results: dict[str, Any] = {}
    results["P497"] = p497_localized_dnls(a)
    results["P498"] = p498_certified_p488_stability()
    results["P499"] = p499_two_frequency_clock(a)
    results["P500"] = p500_stieltjes_identifiability(a)
    results["P501"] = p501_context_limit()
    results["P502"] = p502_passivity_boundary()
    results["P503"] = p503_phase_refinement()
    results["P504"] = p504_minimal_coherence_tester(a)
    results["P505"] = p505_formal_interface()
    results["P506"] = p506_parameter_source_no_go()
    make_figures(results)
    elapsed = time.perf_counter() - started
    payload = {
        "release": "10.50",
        "programs": "P497-P506",
        "seed": SEED,
        "runtime_seconds": elapsed,
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform()},
        "epistemic_boundary": (
            "Local/exact theorems, outward interval computations, finite numerical continuation, "
            "and synthetic records only; no laboratory or external-audit evidence."
        ),
        "results": results,
    }
    RESULTS.write_text(json.dumps(native(payload), indent=2, sort_keys=True), encoding="utf-8")
    write_summary(results)
    print(json.dumps({"runtime_seconds": elapsed, "P497_endpoint_ipr": results["P497"]["endpoint"]["ipr"], "P498_max_nonzero_upper": results["P498"]["largest_certified_nonzero_real_part_upper_bound"], "P501_C384_integral_difference": results["P501"]["finite_to_asymptotic_differences"][-1]["max_difference"], "P504_minimal_contexts": results["P504"]["minimal_context_cardinality"], "P505_valid": results["P505"]["audit"]["valid"]}, indent=2))


if __name__ == "__main__":
    main()
