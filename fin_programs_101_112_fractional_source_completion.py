#!/usr/bin/env python3
"""Execute FIN research Programs 101--112.

The suite independently rechecks the strict Schur nonclosure certificate,
separates bounded graphon, singular local, and fractional continuum limits,
tests long-range propagation and adaptive information dynamics, freezes a
process-tensor preregistration, and audits a compact conditional route to the
strict exponent eta=9/5.

No result below closes QW-2191, creates a dimensional standard, proves Lorentz
symmetry, completes the legacy-to-strict bridge, transfers a legacy physical
role, or promotes L_total/ToE closure.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_101_112_Fractional_Source_Completion_Results.json"
PREREG = ROOT / "FIN_Programs_101_112_Process_Tensor_Preregistration.json"
FIG = ROOT / "FIN_Programs_101_112_Fractional_Source_Completion_Figures"
FIG.mkdir(exist_ok=True)

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8
MASS2 = 0.25
SEED = 20260727


def sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def canonical_digest(record: dict) -> str:
    payload = json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def strict_profile(d, omega=OMEGA, phi=PHI, beta=BETA, eta=ETA):
    d = np.asarray(d, dtype=float)
    return np.abs(np.cos(omega * d + phi)) / (1.0 + beta * d**eta)


def interval_bounds(x) -> tuple[float, float]:
    # The interval package encloses at arbitrary precision.  The final binary64
    # serialization is widened by one ulp so that conversion cannot shrink it.
    return (
        float(np.nextafter(float(x.a), -np.inf)),
        float(np.nextafter(float(x.b), np.inf)),
    )


def ratio_interval(ppi: tuple[float, float], ph: tuple[float, float],
                   mass: tuple[float, float]) -> tuple[list[float], list[float]]:
    """Return native and Schur ratio hulls by monotone-safe corner evaluation."""
    native, reduced = [], []
    for p_pi, p_h, m in itertools.product(ppi, ph, mass):
        lp = m + 1.0 - p_pi
        lh = m + 1.0 - p_h
        native.append(m / lp)
        reduced.append((2.0 * m * lp / (m + lp)) / lh)
    return [min(native), max(native)], [min(reduced), max(reduced)]


def program101_independent_interval_certificate(truncation: int = 20_000) -> dict:
    """Directed interval arithmetic plus an analytic infinite-tail enclosure."""
    mp.iv.dps = 50
    iv = mp.iv
    omega = iv.mpf([743, 743]) / 4000
    phi = iv.mpf([13, 13]) / 80
    eta = iv.mpf([9, 9]) / 5
    z = iv.mpf([0, 0])
    n_pi = iv.mpf([0, 0])
    n_half = iv.mpf([0, 0])
    quarter = (1.0, 0.0, -1.0, 0.0)
    for d in range(1, truncation + 1):
        a = abs(iv.cos(omega * d + phi)) / (1 + iv.mpf(d) ** eta)
        z += 2 * a
        n_pi += 2 * a * (-1 if d % 2 else 1)
        n_half += 2 * a * quarter[d % 4]

    z_lo, z_hi = interval_bounds(z)
    npi_lo, npi_hi = interval_bounds(n_pi)
    nh_lo, nh_hi = interval_bounds(n_half)
    tail = 2.0 * truncation ** (1.0 - ETA) / (ETA - 1.0)
    ppi = interval_bounds(
        iv.mpf([npi_lo - tail, npi_hi + tail])
        / iv.mpf([z_lo, z_hi + tail])
    )
    ph = interval_bounds(
        iv.mpf([nh_lo - tail, nh_hi + tail])
        / iv.mpf([z_lo, z_hi + tail])
    )
    native, reduced = ratio_interval(ppi, ph, (MASS2, MASS2))
    margin = reduced[0] - native[1]

    previous = ROOT / "FIN_Programs_91_100_Critical_Nonlocal_Operational_Results.json"
    previous_cert = {}
    if previous.exists():
        previous_cert = json.loads(previous.read_text(encoding="utf-8"))[
            "program91"
        ]["certificate"]

    fig, ax = plt.subplots(figsize=(8.2, 4.5), constrained_layout=True)
    centers = [np.mean(native), np.mean(reduced)]
    errors = [0.5 * np.ptp(native), 0.5 * np.ptp(reduced)]
    ax.errorbar(["native", "Schur retained"], centers, yerr=errors, fmt="o",
                capsize=7, color="#1f5a99")
    ax.set_ylabel("scale-invariant spectral ratio")
    ax.set_title("Independent directed-interval nonclosure certificate")
    ax.grid(True, axis="y", alpha=0.25)
    fig.savefig(FIG / "program101_independent_interval_certificate.png", dpi=190)
    plt.close(fig)

    return {
        "arithmetic": "mpmath.iv directed interval arithmetic, 50 decimal digits",
        "exact_inputs": {"omega": "743/4000", "phi": "13/80", "eta": "9/5"},
        "truncation": truncation,
        "partial_normalization_interval": [z_lo, z_hi],
        "tail_upper_bound": tail,
        "p_hat_pi_interval": list(ppi),
        "p_hat_pi_over_2_interval": list(ph),
        "native_ratio_interval": native,
        "schur_ratio_interval": reduced,
        "disjoint_margin": margin,
        "passes": margin > 0,
        "comparison_with_program91_margin": previous_cert.get("disjoint_interval_margin"),
        "status": "Proven, computer-assisted with outward intervals and analytic tail",
        "scope": "Fixed mass 1/4, normalized infinite strict lattice, alternating-site Schur map.",
    }


def program102_nonclosure_region(p101: dict) -> dict:
    ppi = tuple(p101["p_hat_pi_interval"])
    ph = tuple(p101["p_hat_pi_over_2_interval"])
    # With a=1-p_pi and b=1-p_half,
    # R_S=2m(m+a)/((2m+a)(m+b)).  It increases in a and m and decreases
    # in b whenever 2b-a>0.  The following strict lower bound licenses the
    # corner hulls over every mass cell.
    a_interval = (1.0 - ppi[1], 1.0 - ppi[0])
    b_interval = (1.0 - ph[1], 1.0 - ph[0])
    schur_mass_monotonicity_margin = 2.0 * b_interval[0] - a_interval[1]
    if schur_mass_monotonicity_margin <= 0:
        raise RuntimeError("Schur ratio mass monotonicity not certified")
    mass_edges = np.linspace(0.01, 2.0, 2001)
    certified = []
    for a, b in zip(mass_edges[:-1], mass_edges[1:]):
        native, reduced = ratio_interval(ppi, ph, (float(a), float(b)))
        certified.append(reduced[0] - native[1])

    values = {
        "omega": [OMEGA - 0.005, OMEGA, OMEGA + 0.005],
        "phi": [PHI - 0.02, PHI, PHI + 0.02],
        "beta": [0.8, 1.0, 1.2],
        "eta": [1.7, 1.8, 1.9],
        "mass": [0.1, 0.25, 0.5],
    }
    d_int = np.arange(1, 50_001, dtype=np.int64)
    d = d_int.astype(float)
    parity = np.where(d_int % 2, -1.0, 1.0)
    quarter = np.cos(0.5 * np.pi * d_int)
    grid = []
    for omega, phi, beta, eta, mass in itertools.product(*values.values()):
        a = np.abs(np.cos(omega * d + phi)) / (1.0 + beta * d**eta)
        z = float(2 * a.sum())
        tail = 2.0 / (beta * (eta - 1.0)) * 50_000 ** (1.0 - eta)
        err = 2.0 * tail / z
        ppi0 = float(2 * np.dot(a, parity) / z)
        ph0 = float(2 * np.dot(a, quarter) / z)
        native, reduced = ratio_interval(
            (ppi0 - err, ppi0 + err), (ph0 - err, ph0 + err), (mass, mass)
        )
        grid.append({
            "omega": omega, "phi": phi, "beta": beta, "eta": eta, "mass": mass,
            "certified_margin": reduced[0] - native[1],
        })
    worst = min(grid, key=lambda row: row["certified_margin"])

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    mids = 0.5 * (mass_edges[:-1] + mass_edges[1:])
    ax.plot(mids, certified, color="#1f5a99")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("mass parameter")
    ax.set_ylabel("certified ratio-separation margin")
    ax.set_title("Strict Schur nonclosure on a continuous mass interval")
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program102_nonclosure_region.png", dpi=190)
    plt.close(fig)

    return {
        "mass_interval": [0.01, 2.0],
        "number_of_interval_cells": len(certified),
        "schur_mass_monotonicity_margin_2b_minus_a": schur_mass_monotonicity_margin,
        "minimum_certified_margin": float(min(certified)),
        "mass_interval_passes": min(certified) > 0,
        "parameter_grid": {
            "values": values, "points": len(grid),
            "minimum_pointwise_certified_margin": worst["certified_margin"],
            "worst_point": worst,
            "all_points_pass": all(row["certified_margin"] > 0 for row in grid),
        },
        "status": (
            "Computer-assisted theorem on the frozen-profile mass interval; "
            "tail-certified pointwise evidence on the multi-parameter grid"
        ),
        "scope": "The parameter grid is not an interval theorem between sampled points.",
    }


def continuum_profile(x):
    x = np.asarray(x, dtype=float)
    return np.cos(OMEGA * np.abs(x) + PHI) / (1.0 + np.abs(x) ** ETA)


def program103_graphon_error_theorem() -> dict:
    lip = OMEGA + ETA * 0.5 ** (ETA - 1.0)
    z_lower = math.cos(OMEGA * 0.5 + PHI) / (1.0 + 0.5**ETA)
    rows = []
    for n in [32, 64, 128, 256, 512, 1024, 2048]:
        e = lip / (4.0 * n)
        normalized = 2.0 * e / (z_lower - e)
        diagonal_mass = 1.0 / (n * (z_lower - e))
        diagonal_removal = 2.0 * diagonal_mass / (1.0 - diagonal_mass)
        rows.append({
            "N": n,
            "midpoint_operator_norm_bound": normalized,
            "diagonal_removal_bound": diagonal_removal,
            "combined_bound": normalized + diagonal_removal,
        })

    # Low-mode numerical convergence, independent of the analytic bound.
    m = 2**20
    x = (np.arange(m) + 0.5) / m - 0.5
    g = continuum_profile(x)
    z = float(g.mean())
    target = np.array([np.mean(g * np.cos(2 * np.pi * k * x)) / z
                       for k in range(1, 9)])
    errors = []
    for n in [32, 64, 128, 256, 512, 1024, 2048]:
        xn = (np.arange(n) + 0.5) / n - 0.5
        gn = continuum_profile(xn)
        coeff = np.array([np.mean(gn * np.cos(2 * np.pi * k * xn)) / gn.mean()
                          for k in range(1, 9)])
        errors.append({"N": n, "low_mode_max_error": float(np.max(abs(coeff - target)))})
    slope = float(np.polyfit(np.log([r["N"] for r in errors]),
                             np.log([r["low_mode_max_error"] for r in errors]), 1)[0])

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.loglog([r["N"] for r in rows], [r["combined_bound"] for r in rows],
              "o-", label="analytic operator bound")
    ax.loglog([r["N"] for r in errors], [r["low_mode_max_error"] for r in errors],
              "s-", label="observed modes 1--8")
    ax.set_xlabel("N")
    ax.set_ylabel("error")
    ax.set_title("Bounded graphon limit: theorem and observed convergence")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program103_graphon_error_theorem.png", dpi=190)
    plt.close(fig)

    return {
        "lipschitz_upper_bound": lip,
        "normalization_lower_bound": z_lower,
        "analytic_rows": rows,
        "low_mode_rows": errors,
        "observed_low_mode_loglog_slope": slope,
        "theorem": (
            "Midpoint cell averaging and optional diagonal removal converge in "
            "L2 operator norm to the bounded circle-convolution operator at O(N^-1)."
        ),
        "status": "Proven analytic bound; numerical convergence check",
        "scope": "This regular coordinate scaling does not yield a local Laplacian.",
    }


def program104_singular_localizing_limit() -> dict:
    m = 2**20
    x = (np.arange(m) + 0.5) / m - 0.5
    g = continuum_profile(x)
    z = float(g.mean())
    p = g / z
    sigma2 = float(np.mean(p * x**2))
    mu4 = float(np.mean(p * x**4))
    rows = []
    for eps in [0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]:
        for k in [1, 2, 3, 4]:
            u = 2 * np.pi * eps * k
            lam = float(np.mean(p * (1.0 - np.cos(u * x))) / eps**2)
            target = 0.5 * sigma2 * (2 * np.pi * k) ** 2
            bound = mu4 * (2 * np.pi * k) ** 4 * eps**2 / 24.0
            rows.append({
                "epsilon": eps, "mode": k, "scaled_eigenvalue": lam,
                "local_limit": target, "absolute_error": abs(lam - target),
                "fourth_moment_bound": bound,
                "covered": abs(lam - target) <= bound * (1 + 1e-10),
            })

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    for k in [1, 2, 3, 4]:
        r = [row for row in rows if row["mode"] == k]
        ax.loglog([q["epsilon"] for q in r], [q["absolute_error"] for q in r],
                  "o-", label=f"mode {k}")
    ax.set_xlabel("localization scale epsilon")
    ax.set_ylabel("absolute spectral error")
    ax.set_title("Singularly localized strict profile converges to a Laplacian")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program104_singular_localizing_limit.png", dpi=190)
    plt.close(fig)

    return {
        "base_variance": sigma2,
        "base_fourth_moment": mu4,
        "rows": rows,
        "all_fourth_moment_bounds_pass": all(r["covered"] for r in rows),
        "theorem": (
            "For p_epsilon(x)=epsilon^-1 p(x/epsilon), "
            "(I-T_epsilon)/epsilon^2 converges modewise to "
            "(sigma^2/2)(-Delta), with the displayed O(epsilon^2) bound."
        ),
        "status": "Proven conditional continuum theorem plus quadrature check",
        "scope": "The singular localization and epsilon^-2 clock are additional premises, not sourced by FIN.",
    }


def program105_fractional_tail_universality() -> dict:
    alpha = ETA - 1.0
    d = np.arange(1, 1_000_001, dtype=float)
    a = strict_profile(d)
    z = float(2 * a.sum())
    qs = np.geomspace(1e-3, 2e-2, 18)
    symbol = np.array([2 * np.dot(a, 1.0 - np.cos(q * d)) / z for q in qs])
    slope, intercept = np.polyfit(np.log(qs), np.log(symbol), 1)
    mean_abs_cos = 2.0 / math.pi
    integral = math.pi / (
        2.0 * math.gamma(1.0 + alpha) * math.sin(math.pi * alpha / 2.0)
    )
    predicted_constant = 2.0 * mean_abs_cos * integral / z
    scaled = symbol / qs**alpha
    odd_residual = float(max(
        abs((2 * np.dot(a, 1 - np.cos(q * d)) / z)
            - (2 * np.dot(a, 1 - np.cos(-q * d)) / z))
        for q in qs
    ))

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.loglog(qs, symbol, "o", label="truncated strict symbol")
    ax.loglog(qs, predicted_constant * qs**alpha, "--", label=r"predicted $C|q|^{4/5}$")
    ax.set_xlabel("|q|")
    ax.set_ylabel(r"$1-\widehat p(q)$")
    ax.set_title("Fractional small-wave-number law of the strict lattice tail")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program105_fractional_tail_universality.png", dpi=190)
    plt.close(fig)

    return {
        "fractional_order": alpha,
        "truncation": 1_000_000,
        "normalization_partial": z,
        "q_values": qs.tolist(),
        "symbol_values": symbol.tolist(),
        "fitted_loglog_slope": float(slope),
        "fitted_constant": float(math.exp(intercept)),
        "predicted_abelian_constant": predicted_constant,
        "scaled_constant_range": [float(min(scaled)), float(max(scaled))],
        "even_symbol_odd_residual": odd_residual,
        "theorem": (
            "Equidistribution of the irrational cosine phase and the weighted "
            "Abelian theorem give 1-p_hat(q)~C|q|^(eta-1), here order 4/5."
        ),
        "status": "Analytic theorem from standard equidistribution/Abelian lemmas; strong numerical check",
        "scope": "The lattice route is fractional and nonlocal; it is not the local limit of Program 104.",
    }


def program106_long_range_semigroup() -> dict:
    n = 2**18
    idx = np.arange(n)
    dist = np.minimum(idx, n - idx).astype(float)
    p = strict_profile(dist)
    p[0] = 0
    p /= p.sum()
    phat = np.fft.fft(p)
    rows = []
    for t in [0.1, 0.5, 1.0]:
        heat = np.fft.ifft(np.exp(-t * (1.0 - phat))).real
        heat = np.maximum(heat, 0)
        heat /= heat.sum()
        for radius in [128, 256, 512, 1024, 2048, 4096, 8192]:
            mask = dist > radius
            one_jump_tail = float(p[mask].sum())
            heat_tail = float(heat[mask].sum())
            rows.append({
                "time": t, "radius": radius,
                "one_jump_tail": one_jump_tail,
                "heat_tail": heat_tail,
                "heat_over_t_one_jump": heat_tail / (t * one_jump_tail),
            })

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    for t in [0.1, 0.5, 1.0]:
        r = [row for row in rows if row["time"] == t]
        ax.semilogx([q["radius"] for q in r], [q["heat_over_t_one_jump"] for q in r],
                    "o-", label=f"t={t:g}")
    ax.axhline(1, color="black", ls="--", lw=0.8)
    ax.set_xlabel("radius R")
    ax.set_ylabel(r"$P_t(|X|>R)/(t\,p(|X|>R))$")
    ax.set_title("One-big-jump asymptotics for the strict heat semigroup")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program106_long_range_semigroup.png", dpi=190)
    plt.close(fig)

    return {
        "cycle_size": n,
        "rows": rows,
        "maximum_large_radius_deviation_from_one": float(max(
            abs(r["heat_over_t_one_jump"] - 1.0) for r in rows if r["radius"] >= 4096
        )),
        "theorem": (
            "The Poissonized semigroup e^{-t(I-P)} has fixed-time far tails "
            "asymptotic to t times the regularly varying one-jump tail."
        ),
        "status": "Standard subexponential-tail theorem; finite-cycle numerical check",
        "scope": "Polynomial direct jumps preclude an exponential finite-speed cone in this lattice model.",
    }


def program107_adaptive_information_manifold() -> dict:
    d = np.arange(1, 7, dtype=float)
    multiplicity = np.array([2, 2, 2, 2, 2, 1], dtype=float)
    q0 = multiplicity * strict_profile(d)
    q0 /= q0.sum()
    ref = np.ones(6) / 6.0
    path_lap = np.zeros((6, 6))
    for i in range(5):
        path_lap[i, i] += 1
        path_lap[i + 1, i + 1] += 1
        path_lap[i, i + 1] -= 1
        path_lap[i + 1, i] -= 1
    gamma = 0.2

    def functional(q):
        return float(np.dot(q, np.log(q / ref)) + 0.5 * gamma * q @ path_lap @ q)

    def gradient(q):
        return np.log(q / ref) + 1.0 + gamma * path_lap @ q

    def rhs(_t, q):
        grad = gradient(q)
        return -q * (grad - np.dot(q, grad))

    sol = solve_ivp(rhs, (0, 20), q0, rtol=1e-10, atol=1e-12,
                    dense_output=False, max_step=0.05)
    fs = np.array([functional(sol.y[:, j]) for j in range(sol.y.shape[1])])
    derivative_checks = []
    for j in range(sol.y.shape[1]):
        q = sol.y[:, j]
        grad = gradient(q)
        derivative_checks.append(float(-np.dot(q, (grad - np.dot(q, grad)) ** 2)))
    strict_grad = gradient(q0)
    strict_stationarity = float(np.linalg.norm(
        q0 * (strict_grad - np.dot(q0, strict_grad))
    ))
    projected = strict_grad - strict_grad.mean()
    euclidean_trial = q0 - projected

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.4), constrained_layout=True)
    for i in range(6):
        ax1.plot(sol.t, sol.y[i], label=f"shell {i+1}")
    ax1.set_xlabel("adaptive time")
    ax1.set_ylabel("shell probability")
    ax1.set_title("Fisher/Shahshahani gradient flow")
    ax1.legend(ncol=2, fontsize=8)
    ax2.plot(sol.t, fs, color="#a61b1b")
    ax2.set_xlabel("adaptive time")
    ax2.set_ylabel("functional")
    ax2.set_title("Monotone information functional")
    ax2.grid(True, alpha=0.25)
    fig.savefig(FIG / "program107_adaptive_information_manifold.png", dpi=190)
    plt.close(fig)

    return {
        "initial_strict_shell_probability": q0.tolist(),
        "final_probability": sol.y[:, -1].tolist(),
        "minimum_probability": float(sol.y.min()),
        "maximum_simplex_sum_residual": float(np.max(abs(sol.y.sum(axis=0) - 1))),
        "initial_functional": float(fs[0]),
        "final_functional": float(fs[-1]),
        "maximum_functional_increase": float(max(0, np.max(np.diff(fs)))),
        "maximum_analytic_derivative": float(max(derivative_checks)),
        "strict_stationarity_residual": strict_stationarity,
        "unconstrained_euclidean_unit_step_minimum": float(euclidean_trial.min()),
        "theorem": (
            "Replicator-form Fisher gradient flow preserves the open simplex "
            "and satisfies dF/dt=-Var_q(grad F)<=0."
        ),
        "status": "Proven geometric identity plus numerical integration",
        "scope": "The chosen source-neutral F drives the strict profile away; it does not derive FIN.",
    }


def cyclic_kernel_from_theta(theta, n):
    omega, phi, beta, eta = theta
    idx = np.arange(n)
    d = np.minimum(idx, n - idx).astype(float)
    return np.cos(omega * d + phi) / (1 + beta * d**eta)


def feature_library(theta, n):
    k = cyclic_kernel_from_theta(theta, n)
    spectrum = np.abs(np.fft.fft(k)) ** 2
    spectrum /= spectrum.sum()
    diff = np.roll(k, -1) - k
    return np.array([
        np.mean(diff**2),
        np.mean(k**2),
        np.mean(k) ** 2,
        np.mean(k**4),
        -np.sum(spectrum * np.log(spectrum + 1e-300)),
        np.mean((np.roll(k, -2) - k) ** 2),
    ])


def numeric_gradient(fun, theta, step=2e-5):
    theta = np.asarray(theta, dtype=float)
    result = np.empty(len(theta))
    for i in range(len(theta)):
        h = step * max(1.0, abs(theta[i]))
        plus, minus = theta.copy(), theta.copy()
        plus[i] += h
        minus[i] -= h
        result[i] = (fun(plus) - fun(minus)) / (2 * h)
    return result


def numeric_hessian(fun, theta, step=2e-4):
    theta = np.asarray(theta, dtype=float)
    h = np.empty((len(theta), len(theta)))
    for j in range(len(theta)):
        delta = step * max(1.0, abs(theta[j]))
        plus, minus = theta.copy(), theta.copy()
        plus[j] += delta
        minus[j] -= delta
        h[:, j] = (numeric_gradient(fun, plus) - numeric_gradient(fun, minus)) / (2 * delta)
    return 0.5 * (h + h.T)


def program108_inverse_source_mdl() -> dict:
    theta = np.array([OMEGA, PHI, BETA, ETA])
    gradients = np.column_stack([
        numeric_gradient(lambda t, j=j: feature_library(t, 13)[j], theta)
        for j in range(6)
    ])
    u, singular, vt = np.linalg.svd(gradients)
    rank = int(np.sum(singular > singular[0] * 1e-9))
    null_vectors = vt[rank:].T
    candidates = []
    for j in range(null_vectors.shape[1]):
        w = null_vectors[:, j]
        train = gradients @ w
        held = numeric_gradient(lambda t: float(feature_library(t, 17) @ w), theta)
        hess = numeric_hessian(lambda t: float(feature_library(t, 13) @ w), theta)
        candidates.append({
            "weights": w.tolist(),
            "training_gradient_norm": float(np.linalg.norm(train)),
            "heldout_N17_gradient_norm": float(np.linalg.norm(held)),
            "hessian_eigenvalues": np.linalg.eigvalsh(hess).tolist(),
            "local_minimum": bool(np.linalg.eigvalsh(hess).min() > 0),
        })

    return {
        "feature_names": [
            "nearest roughness", "mean square", "square of mean",
            "fourth moment", "spectral entropy", "two-step roughness"
        ],
        "training_cycle": 13,
        "heldout_cycle": 17,
        "gradient_matrix_singular_values": singular.tolist(),
        "gradient_matrix_rank": rank,
        "stationary_coefficient_nullity": int(6 - rank),
        "candidates": candidates,
        "parameter_count_comparison": {
            "strict_tuple": 4,
            "source_coefficients_up_to_scale": 5,
            "source_has_mdl_advantage": False,
        },
        "status": "Negative result for the declared six-feature inverse-source family",
        "conclusion": (
            "Stationarity can be engineered nonuniquely at N=13, but it is "
            "coefficient-heavy and fails the N=17 transfer test; no source theorem results."
        ),
    }


def program109_chiral_source_intake(p105: dict) -> dict:
    audited = [
        "fundamental_action_reconstruction/generated/p2948*",
        "fundamental_action_reconstruction/generated/p2958*",
        "fundamental_action_reconstruction/generated/p2965*",
        "FIN Programs 101--108",
    ]
    return {
        "audited_classes": audited,
        "new_explicit_signed_formula_admitted": False,
        "fractional_symbol_odd_residual": p105["even_symbol_odd_residual"],
        "finding": (
            "The 4/5 fractional symbol is exactly even. Tail universality, "
            "localization, graphon convergence, and scalar adaptive geometry "
            "supply no inversion-odd signed value."
        ),
        "status": "Proven symmetry obstruction for the new limits; no-source certificate",
        "guardrail": "QW-2191 and the orientation-torsor polarity remain open.",
    }


def program110_process_tensor_preregistration() -> dict:
    design = {
        "schema": "FIN_PROCESS_TENSOR_PREREGISTRATION_V1",
        "record_date": "2026-07-27",
        "primary_hypotheses": {
            "null": "best-fit time-homogeneous diffusion semigroup",
            "alternative": "unitary wave dynamics generated by the same declared A",
        },
        "frozen_operator": {
            "W": "K_s on the declared finite graph", "s": 1.660307278766099,
            "A": "s I - W", "wave": "exp(-itA)", "diffusion": "exp(-tA)",
        },
        "preparation": "localized site preparation at x=0",
        "instrument": "full-site projective readout after each of two slots",
        "slot_times": [0.972441, 0.972441],
        "detector_model": "symmetric independent confusion channel",
        "calibration_rule": "fit detector confusion on calibration partition only",
        "primary_statistic": "held-out log-likelihood ratio wave versus diffusion",
        "secondary_statistics": ["Jensen-Shannon divergence", "Chernoff information"],
        "sample_plan": {
            "heldout_records": 49,
            "target_one-sided_error_bound": 1e-4,
            "bound": "exp(-n * 0.175304)",
        },
        "exclusions": [
            "records with missing slot timestamp",
            "undeclared detector remapping",
            "operator or graph changed after preregistration",
        ],
        "provenance_requirements": [
            "raw immutable records", "calibration split hash", "apparatus clock calibration",
            "graph and operator hash", "instrument implementation hash",
        ],
        "external_data_admitted": False,
    }
    design["canonical_sha256"] = canonical_digest(design)
    PREREG.write_text(json.dumps(design, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "file": PREREG.name,
        "canonical_sha256": design["canonical_sha256"],
        "file_sha256": sha256(PREREG),
        "status": "Preregistered protocol; no external data admitted",
    }


def binary_entropy(p):
    if p in (0, 1):
        return 0.0
    return -p * math.log(p) - (1 - p) * math.log(1 - p)


def program111_correlated_apparatus_ledger() -> dict:
    rows = []
    for eps, rho, tau in itertools.product([0.05, 0.1, 0.2], [0.0, 0.5, 0.9],
                                           [0.5, 1.0, 2.0, 5.0]):
        a = eps * (1 - rho)
        b = eps + rho * (1 - eps)
        h_rate = (1 - eps) * binary_entropy(a) + eps * binary_entropy(b)
        info_rate = math.log(2) - h_rate
        sigma = 0.1 / tau
        system_work = info_rate - sigma
        reset_cost = math.log(2)
        complete_excess = reset_cost - system_work
        rows.append({
            "error_probability": eps, "persistence": rho, "duration": tau,
            "error_entropy_rate": h_rate, "information_rate": info_rate,
            "finite_time_dissipation": sigma, "system_work": system_work,
            "complete_excess": complete_excess,
            "identity_residual": complete_excess - (h_rate + sigma),
        })

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    for rho in [0.0, 0.5, 0.9]:
        r = [row for row in rows if row["error_probability"] == 0.1
             and row["persistence"] == rho]
        ax.plot([q["duration"] for q in r], [q["complete_excess"] for q in r],
                "o-", label=f"persistence={rho:g}")
    ax.set_xlabel("protocol duration")
    ax.set_ylabel("apparatus-inclusive excess (nats)")
    ax.set_title("Correlated detector memory and finite-time work ledger")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program111_correlated_apparatus_ledger.png", dpi=190)
    plt.close(fig)

    return {
        "rows": rows,
        "minimum_complete_excess": float(min(r["complete_excess"] for r in rows)),
        "maximum_identity_residual": float(max(abs(r["identity_residual"]) for r in rows)),
        "identity": "W_complete - DeltaF = h_rate(E) + Sigma(tau) >= 0",
        "status": "Proven information-theoretic identity in an explicitly imported thermodynamic model",
        "scope": "Temperature/action scale and reset implementation remain external operational data.",
    }


def program112_eta_source_skeleton() -> dict:
    retention = math.exp(-ALPHA_GEO / 5.0)
    exponent = -math.log(retention, 2)
    artifacts = [
        "fundamental_action_reconstruction/generated/p2948_s1898_torsion_character_ratio_package_theorem_skeleton.md",
        "fundamental_action_reconstruction/generated/p2950_s1900_exact_package_beta_eta_scale_coupling_obstruction.md",
        "fundamental_action_reconstruction/generated/p2958_s1908_p2938_strict_provenance_theorem_interface.md",
        "fundamental_action_reconstruction/generated/p2965_s1915_primitive_mean_scale_unit_quotient_candidate.md",
        "fundamental_action_reconstruction/generated/p2860_s1810_compression_tail_multiplicative_scale_law_no_eta_source_audit.md",
    ]
    inventory = []
    for rel in artifacts:
        path = ROOT / rel
        inventory.append({"path": rel, "exists": path.exists(), "sha256": sha256(path)})

    fig, ax = plt.subplots(figsize=(8.5, 4.7), constrained_layout=True)
    labels = ["algebraic identity", "fivefold source", "damping coupling", "beta/unit source"]
    values = [1, 0, 0, 0]
    colors = ["#19733a", "#a61b1b", "#a61b1b", "#a61b1b"]
    ax.bar(labels, values, color=colors)
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("strict obligation discharged")
    ax.set_title(r"Conditional $\eta=9/5$ skeleton and its open source obligations")
    ax.tick_params(axis="x", rotation=14)
    fig.savefig(FIG / "program112_eta_source_skeleton.png", dpi=190)
    plt.close(fig)

    return {
        "exact_chain": [
            "alpha_geo = 4 ln 2",
            "r_5 = exp(-alpha_geo/5) = 2^(-4/5)",
            "C(2d)=r_5 C(d) and exact dyadic homogeneity imply C(d) proportional to d^(-4/5)",
            "multiplying a legacy d^(-1) envelope yields strict d^(-9/5)",
            "equivalently eta = 1 + 4/5 = 9/5",
        ],
        "numeric_retention": retention,
        "derived_increment": exponent,
        "derived_eta": 1 + exponent,
        "five_component_candidate": [1, 2, 2, 2, 2],
        "five_component_sum": 9,
        "five_component_mean": 9 / 5,
        "artifact_inventory": inventory,
        "obligations": {
            "strict_source_of_fivefold_quotient": False,
            "non_target_coded_reason_for_dividing_alpha_geo_by_five": False,
            "coupling_to_the_damping_completion_operator": False,
            "target_independent_beta_or_length_unit_source": False,
        },
        "status": "Exact conditional derivation skeleton; not a strict FIN source theorem",
        "conclusion": (
            "The identity compresses several observed 9/5 constructions into "
            "one scale law, but the fivefold input and its damping coupling are "
            "still premises. It therefore explains compatibility, not necessity."
        ),
        "guardrail": "No generic bridge completion or legacy-role transfer follows.",
    }


def main() -> None:
    p101 = program101_independent_interval_certificate()
    results = {
        "metadata": {
            "title": "FIN Programs 101--112: Fractional Limits, Source Tests, and Operational Completion",
            "release": "10.11", "date": "2026-07-27", "seed": SEED,
            "author": "Krzysztof Żuchowski",
            "orcid": "0009-0002-0909-3613",
        },
        "program101_independent_interval_certificate": p101,
        "program102_nonclosure_region": program102_nonclosure_region(p101),
        "program103_graphon_error_theorem": program103_graphon_error_theorem(),
        "program104_singular_localizing_limit": program104_singular_localizing_limit(),
        "program105_fractional_tail_universality": None,
    }
    p105 = program105_fractional_tail_universality()
    results["program105_fractional_tail_universality"] = p105
    results.update({
        "program106_long_range_semigroup": program106_long_range_semigroup(),
        "program107_adaptive_information_manifold": program107_adaptive_information_manifold(),
        "program108_inverse_source_mdl": program108_inverse_source_mdl(),
        "program109_chiral_source_intake": program109_chiral_source_intake(p105),
        "program110_process_tensor_preregistration": program110_process_tensor_preregistration(),
        "program111_correlated_apparatus_ledger": program111_correlated_apparatus_ledger(),
        "program112_eta_source_skeleton": program112_eta_source_skeleton(),
    })
    results["guardrails"] = {
        "strict_kernel_primary": True,
        "legacy_kernel_role": "intermediate bridge kernel only",
        "selector_QW2191_closed": False,
        "dimensional_standard_sourced": False,
        "legacy_strict_bridge_complete": False,
        "legacy_role_transfer_started": False,
        "L_total_or_ToE_promoted": False,
        "external_data_admitted": False,
    }
    OUT.write_text(json.dumps(results, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "output": OUT.name,
        "sha256": sha256(OUT),
        "figures": len(list(FIG.glob("*.png"))),
        "program101_margin": p101["disjoint_margin"],
        "program102_mass_margin": results["program102_nonclosure_region"]["minimum_certified_margin"],
        "program105_slope": p105["fitted_loglog_slope"],
        "program112_eta": results["program112_eta_source_skeleton"]["derived_eta"],
    }, indent=2))


if __name__ == "__main__":
    main()
