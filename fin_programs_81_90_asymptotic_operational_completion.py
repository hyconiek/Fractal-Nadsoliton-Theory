#!/usr/bin/env python3
"""Execute FIN research Programs 81--90.

This suite studies asymptotic Fourier-symbol naturality, Schur fixed points,
dense-kernel discretization universality, approximate locality bounds,
zero-mode quotient actions, robust dimensional calibration, chiral-state
stability, feedback thermodynamics, two-slot process tensors, and the
admission gate for genuinely external data.

All physical standards, baths, preparations, instruments, feedback rules, and
data provenance fields are explicit conditioning data.  No computation here
closes the strict selector, dimensional-scale, legacy-to-strict completion, or
experimental-validation obligations.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_81_90_Asymptotic_Operational_Completion_Results.json"
INTAKE = ROOT / "FIN_Programs_81_90_External_Data_Intake_Template.json"
FIG = ROOT / "FIN_Programs_81_90_Asymptotic_Operational_Completion_Figures"
FIG.mkdir(exist_ok=True)

ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
OMEGA_L = math.pi / 4.0
PHI_L = math.pi / 6.0
OMEGA_S = 0.18575
PHI_S = 0.16250
BETA_S = 1.0
ETA_S = 1.8
MASS2 = 0.25
SEED = 20260727


def k_legacy(d):
    d = np.asarray(d, dtype=float)
    return ALPHA_GEO * np.cos(OMEGA_L * d + PHI_L) / (
        1.0 + BETA_TORS * d
    )


def k_strict(d):
    d = np.asarray(d, dtype=float)
    return np.cos(OMEGA_S * d + PHI_S) / (
        1.0 + BETA_S * d**ETA_S
    )


def cycle_distance(i: int, j: int, n: int) -> int:
    return min((i - j) % n, (j - i) % n)


def radial_matrix(kernel, n: int, *, absolute: bool = False) -> np.ndarray:
    d = np.minimum(
        (np.arange(n)[:, None] - np.arange(n)[None, :]) % n,
        (np.arange(n)[None, :] - np.arange(n)[:, None]) % n,
    )
    w = np.asarray(kernel(d), dtype=float)
    np.fill_diagonal(w, 0.0)
    return np.abs(w) if absolute else w


def laplacian(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(axis=1)) - w


def first_row_profile(profile, n: int, *, mass2: float = MASS2) -> np.ndarray:
    d = np.minimum(np.arange(n), n - np.arange(n)).astype(float)
    x = d / n
    weights = np.abs(np.asarray(profile(x), dtype=float))
    weights[0] = 0.0
    weights /= weights.sum()
    row = -weights
    row[0] = mass2 + 1.0
    return row


def first_row_precision(
    kernel,
    n: int,
    *,
    coordinate_scaled: bool = False,
    mass2: float = MASS2,
) -> np.ndarray:
    d = np.minimum(np.arange(n), n - np.arange(n)).astype(float)
    arg = d / n if coordinate_scaled else d
    weights = np.abs(np.asarray(kernel(arg), dtype=float))
    weights[0] = 0.0
    weights /= weights.sum()
    row = -weights
    row[0] = mass2 + 1.0
    return row


def circulant_from_row(row: np.ndarray) -> np.ndarray:
    return np.array([np.roll(row, i) for i in range(len(row))])


def normalize_row(row: np.ndarray) -> np.ndarray:
    return np.asarray(row, dtype=float) / float(row[0])


def schur_row_symbol(parent_row: np.ndarray) -> np.ndarray:
    even_row = np.asarray(parent_row[::2], dtype=float)
    odd_row = np.asarray(parent_row[1::2], dtype=float)
    a = np.fft.fft(even_row)
    b = np.fft.fft(odd_row)
    s = a - np.abs(b) ** 2 / a
    row = np.fft.ifft(s).real
    return 0.5 * (row + np.r_[row[0], row[:0:-1]])


def harmonic_alias_rg(eigenvalues: np.ndarray) -> np.ndarray:
    """Exact retained precision spectrum after eliminating alternating sites."""
    n2 = len(eigenvalues)
    if n2 % 2:
        raise ValueError("Need an even number of Fourier modes")
    n = n2 // 2
    a = np.asarray(eigenvalues[:n], dtype=float)
    b = np.asarray(eigenvalues[n:], dtype=float)
    return 2.0 * a * b / (a + b)


def relative_l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def naturality_metrics(
    kernel,
    n: int,
    *,
    coordinate_scaled: bool,
) -> dict:
    parent = first_row_precision(
        kernel, 2 * n, coordinate_scaled=coordinate_scaled
    )
    compressed = normalize_row(schur_row_symbol(parent))
    native = normalize_row(
        first_row_precision(kernel, n, coordinate_scaled=coordinate_scaled)
    )
    lc = np.fft.fft(compressed).real
    ln = np.fft.fft(native).real
    return {
        "retained_n": n,
        "frobenius_operator_defect": relative_l2(compressed, native),
        "sup_symbol_defect": float(np.max(np.abs(lc - ln)) / np.max(np.abs(ln))),
        "mode1_relative_defect": float(abs(lc[1] - ln[1]) / abs(ln[1])),
        "green_frobenius_defect": relative_l2(1.0 / lc, 1.0 / ln),
    }


def program81_asymptotic_symbol_audit() -> dict:
    sizes = [192, 768, 3072, 12288, 49152]
    cases = {}
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), constrained_layout=True)
    for name, scaled in [
        ("strict_lattice_distance", False),
        ("strict_coordinate_rescaled", True),
    ]:
        rows = [
            naturality_metrics(k_strict, n, coordinate_scaled=scaled)
            for n in sizes
        ]
        fro = np.array([r["frobenius_operator_defect"] for r in rows])
        sup = np.array([r["sup_symbol_defect"] for r in rows])
        mode1 = np.array([r["mode1_relative_defect"] for r in rows])
        cases[name] = {
            "rows": rows,
            "frobenius_loglog_slope": float(
                np.polyfit(np.log(sizes), np.log(fro), 1)[0]
            ),
            "sup_symbol_loglog_slope": float(
                np.polyfit(np.log(sizes), np.log(sup), 1)[0]
            ),
            "mode1_loglog_slope": float(
                np.polyfit(np.log(sizes), np.log(mode1), 1)[0]
            ),
        }
        axes[0].loglog(sizes, fro, "o-", label=name.replace("_", " "))
        axes[1].semilogx(sizes, sup, "s-", label=name.replace("_", " "))
    axes[0].set_xlabel("retained N")
    axes[0].set_ylabel("relative Frobenius defect")
    axes[0].set_title("Averaged naturality metric")
    axes[1].set_xlabel("retained N")
    axes[1].set_ylabel("normalized sup symbol defect")
    axes[1].set_title("Uniform spectral obstruction")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=8)
    fig.savefig(FIG / "program81_asymptotic_symbol_audit.png", dpi=190)
    plt.close(fig)

    # Verify the exact alias/harmonic-mean form against the block Schur symbol.
    row = first_row_precision(k_strict, 4096, coordinate_scaled=True)
    lam = np.fft.fft(row).real
    harmonic = harmonic_alias_rg(lam)
    direct = np.fft.fft(schur_row_symbol(row)).real
    alias_residual = relative_l2(harmonic, direct)
    return {
        "sizes": sizes,
        "cases": cases,
        "alias_harmonic_mean_formula": (
            "lambda_coarse(k)=2 lambda(k) lambda(k+N) / "
            "(lambda(k)+lambda(k+N))"
        ),
        "alias_formula_relative_residual": alias_residual,
        "dilution_theorem": (
            "If two dense circulant precision rows have the same O(1) "
            "diagonal and uniformly O(1/N) off-diagonal differences, their "
            "relative Frobenius defect is O(N^-1/2). This does not imply "
            "uniform convergence of their Fourier symbols."
        ),
        "status": "Exact alias theorem plus asymptotic numerical obstruction",
    }


def normalized_symbol(symbol: np.ndarray) -> np.ndarray:
    return np.asarray(symbol, dtype=float) / float(np.mean(symbol))


def iterate_rg(symbol: np.ndarray, steps: int) -> list[dict]:
    lam = normalized_symbol(symbol)
    rows = []
    for step in range(steps + 1):
        n = len(lam)
        q = 2.0 * math.pi * np.arange(n) / n
        nn = normalized_symbol(2.0 - 2.0 * np.cos(q))
        rows.append(
            {
                "step": step,
                "modes": n,
                "distance_to_constant": relative_l2(lam, np.ones(n)),
                "distance_to_massless_nearest_neighbour": relative_l2(lam, nn),
                "minimum_eigenvalue": float(lam.min()),
                "maximum_eigenvalue": float(lam.max()),
            }
        )
        if step < steps:
            lam = normalized_symbol(harmonic_alias_rg(lam))
    return rows


def program82_fixed_symbol_classification() -> dict:
    n = 4096
    q = 2.0 * math.pi * np.arange(n) / n
    symbols = {
        "constant": np.ones(n),
        "massless_nearest_neighbour": 2.0 - 2.0 * np.cos(q),
        "massive_nearest_neighbour_r_0_1": 0.1 + 2.0 - 2.0 * np.cos(q),
        "strict_positive_precision": np.fft.fft(
            first_row_precision(k_strict, n)
        ).real,
        "legacy_absolute_precision": np.fft.fft(
            first_row_precision(k_legacy, n)
        ).real,
    }
    results = {name: iterate_rg(lam, 8) for name, lam in symbols.items()}
    fig, ax = plt.subplots(figsize=(8.5, 4.7), constrained_layout=True)
    for name, rows in results.items():
        ax.semilogy(
            [r["step"] for r in rows],
            [max(r["distance_to_constant"], 1e-16) for r in rows],
            "o-",
            label=name.replace("_", " "),
        )
    ax.set_xlabel("normalized harmonic-mean RG step")
    ax.set_ylabel("relative distance to constant symbol")
    ax.set_title("Massive symbols flow to the ultralocal constant boundary")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=7)
    fig.savefig(FIG / "program82_fixed_symbol_classification.png", dpi=190)
    plt.close(fig)
    return {
        "families": results,
        "affine_cosine_theorem": {
            "family": "lambda(q)=m+2c(1-cos q), m>=0, c>0",
            "ratio_map": "r'=r(r+4), r=m/c",
            "admissible_finite_fixed_ratio": 0.0,
            "massless_fixed_symbol_after_mean_normalization": True,
            "constant_ultralocal_boundary": "r=infinity",
            "stability": (
                "r=0 is unstable to positive mass; every r>0 flows toward "
                "the constant-symbol boundary under repeated normalized blocking"
            ),
        },
        "status": "Proven in the affine-cosine cone; numerical wider-family audit",
    }


def program83_dense_kernel_universality() -> dict:
    profiles = {
        "constant": lambda x: np.ones_like(x),
        "exponential": lambda x: np.exp(-3.0 * x),
        "gaussian": lambda x: np.exp(-12.0 * x**2),
        "rational_power": lambda x: 1.0 / (1.0 + 4.0 * x**1.8),
        "strict_coordinate_profile": lambda x: np.abs(k_strict(x)),
        "legacy_coordinate_profile": lambda x: np.abs(k_legacy(x)),
    }
    sizes = [96, 192, 384, 768, 1536, 3072, 6144]
    results = {}
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4), constrained_layout=True)
    for name, profile in profiles.items():
        rows = []
        for n in sizes:
            parent = first_row_profile(profile, 2 * n)
            compressed = normalize_row(schur_row_symbol(parent))
            native = normalize_row(first_row_profile(profile, n))
            lc = np.fft.fft(compressed).real
            ln = np.fft.fft(native).real
            rows.append(
                {
                    "retained_n": n,
                    "frobenius_defect": relative_l2(compressed, native),
                    "sup_symbol_defect": float(
                        np.max(np.abs(lc - ln)) / np.max(np.abs(ln))
                    ),
                }
            )
        fro = np.array([r["frobenius_defect"] for r in rows])
        sup = np.array([r["sup_symbol_defect"] for r in rows])
        results[name] = {
            "rows": rows,
            "frobenius_loglog_slope": float(
                np.polyfit(np.log(sizes), np.log(fro), 1)[0]
            ),
            "final_sup_symbol_defect": float(sup[-1]),
            "sup_symbol_loglog_slope": float(
                np.polyfit(np.log(sizes), np.log(sup), 1)[0]
            ),
        }
        axes[0].loglog(sizes, fro, "o-", label=name.replace("_", " "))
        axes[1].semilogx(sizes, sup, "s-", label=name.replace("_", " "))
    axes[0].set_xlabel("retained N")
    axes[0].set_ylabel("Frobenius naturality defect")
    axes[0].set_title("Generic dense-row dilution")
    axes[1].set_xlabel("retained N")
    axes[1].set_ylabel("sup symbol defect")
    axes[1].set_title("Uniform mismatch survives")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=7)
    fig.savefig(FIG / "program83_dense_kernel_universality.png", dpi=190)
    plt.close(fig)
    return {
        "sizes": sizes,
        "profiles": results,
        "constant_profile_exact_limit": {
            "native_normalized_zero_mode": MASS2 / (MASS2 + 1.0),
            "schur_normalized_zero_mode": 2.0 * MASS2 / (2.0 * MASS2 + 1.0),
            "absolute_sup_symbol_gap": abs(
                2.0 * MASS2 / (2.0 * MASS2 + 1.0)
                - MASS2 / (MASS2 + 1.0)
            ),
            "formula": "|2m/(2m+1)-m/(m+1)|",
        },
        "theorem_scope": (
            "Bounded positive coordinate profiles with normalized dense weights "
            "have O(1/N) row entries. Under a fixed positive mass and uniformly "
            "bounded eliminated inverse, Schur rows retain that scale, so any "
            "same-diagonal row comparison has O(N^-1/2) relative Frobenius defect."
        ),
        "interpretation": (
            "An N^-1/2 Frobenius slope is therefore generic and cannot by itself "
            "identify a FIN-specific continuum fixed point."
        ),
        "status": "Analytic dilution theorem and six-profile numerical test",
    }


def exponential_tail(x: float, start: int) -> float:
    term = x**start / math.factorial(start)
    total = term
    n = start
    while term > 1e-18 * max(1.0, total) and n < 500:
        n += 1
        term *= x / n
        total += term
    return total


def program84_locality_bounds() -> dict:
    n = 64
    source, target = 0, n // 2
    w = radial_matrix(k_strict, n, absolute=True)
    w /= w.sum(axis=1)[0]
    l_full = laplacian(w)
    radii = [1, 2, 4, 8, 16]
    times = [0.1, 0.5, 1.0]
    rows = []
    for radius in radii:
        mask = np.zeros_like(w, dtype=bool)
        for i in range(n):
            for j in range(n):
                mask[i, j] = i != j and cycle_distance(i, j, n) <= radius
        wr = np.where(mask, w, 0.0)
        lr = laplacian(wr)
        norm_lr = float(np.linalg.norm(lr, 2))
        delta_norm = float(np.linalg.norm(l_full - lr, 2))
        hop_order = math.ceil((n // 2) / radius)
        for time in times:
            u = expm(-1j * time * l_full)
            p = expm(-time * l_full)
            local_tail = exponential_tail(time * norm_lr, hop_order)
            amplitude_bound = min(1.0, local_tail + time * delta_norm)
            rows.append(
                {
                    "radius": radius,
                    "time": time,
                    "hop_order": hop_order,
                    "truncated_generator_norm": norm_lr,
                    "truncation_operator_norm": delta_norm,
                    "series_tail_bound": local_tail,
                    "full_wave_amplitude_actual": float(abs(u[target, source])),
                    "full_wave_probability_actual": float(
                        abs(u[target, source]) ** 2
                    ),
                    "full_diffusion_entry_actual": float(p[target, source]),
                    "combined_amplitude_bound": amplitude_bound,
                    "combined_wave_probability_bound": amplitude_bound**2,
                    "combined_diffusion_entry_bound": amplitude_bound,
                }
            )
    best = []
    for time in times:
        subset = [r for r in rows if r["time"] == time]
        best.append(min(subset, key=lambda r: r["combined_amplitude_bound"]))
    fig, ax = plt.subplots(figsize=(8.5, 4.7), constrained_layout=True)
    for time in times:
        subset = [r for r in rows if r["time"] == time]
        ax.semilogy(
            [r["radius"] for r in subset],
            [r["combined_amplitude_bound"] for r in subset],
            "o-",
            label=f"t={time}",
        )
    ax.set_xlabel("truncation radius")
    ax.set_ylabel("rigorous combined amplitude bound")
    ax.set_title("Local-series tail plus Duhamel truncation error")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program84_locality_bounds.png", dpi=190)
    plt.close(fig)
    return {
        "cycle_size": n,
        "source_target_distance": n // 2,
        "rows": rows,
        "best_bound_by_time": best,
        "theorem": (
            "For a radius-R truncation L_R and hop order m, the far matrix "
            "element of exp(z L_R) is bounded by sum_(k>=m) "
            "(|z| ||L_R||)^k/k!. Duhamel adds at most "
            "|t| ||L-L_R|| for either unitary self-adjoint or positive "
            "contractive semigroups."
        ),
        "status": "Rigorous finite-volume bound; no Lorentz interpretation",
    }


def quotient_basis(n: int) -> np.ndarray:
    seed = np.column_stack([np.ones(n) / math.sqrt(n), np.eye(n)[:, 1:]])
    q, _ = np.linalg.qr(seed)
    return q[:, 1:]


def quotient_action_record(l: np.ndarray, rng: np.random.Generator) -> dict:
    n = l.shape[0]
    v = quotient_basis(n)
    a = v.T @ l @ v
    eig = np.linalg.eigvalsh(a)
    j = rng.normal(size=n)
    j -= j.mean()
    f = np.linalg.pinv(l, hermitian=True) @ j
    action = lambda x: 0.5 * float(x @ l @ x) - float(j @ x)
    gauge_differences = [abs(action(f + c * np.ones(n)) - action(f)) for c in [-3, -1, 1, 4]]
    log_pdet = float(np.log(eig).sum())
    log_z_pseudo = float(
        0.5 * (n - 1) * math.log(2.0 * math.pi)
        - 0.5 * log_pdet
        + 0.5 * j @ np.linalg.pinv(l, hermitian=True) @ j
    )
    jr = v.T @ j
    log_z_basis = float(
        0.5 * (n - 1) * math.log(2.0 * math.pi)
        - 0.5 * np.linalg.slogdet(a)[1]
        + 0.5 * jr @ np.linalg.solve(a, jr)
    )
    entropy = float(
        0.5 * ((n - 1) * (1.0 + math.log(2.0 * math.pi)) - log_pdet)
    )
    t0 = 0.4
    ft = expm(-t0 * l) @ f
    analytic_derivative = -float((l @ ft) @ (l @ ft))
    eps = 1e-6
    ep = 0.5 * float((expm(-(t0 + eps) * l) @ f) @ l @ (expm(-(t0 + eps) * l) @ f))
    em = 0.5 * float((expm(-(t0 - eps) * l) @ f) @ l @ (expm(-(t0 - eps) * l) @ f))
    finite_derivative = (ep - em) / (2 * eps)
    return {
        "quotient_dimension": n - 1,
        "minimum_quotient_eigenvalue": float(eig.min()),
        "stationarity_residual": float(np.linalg.norm(l @ f - j)),
        "maximum_gauge_action_difference": max(gauge_differences),
        "minimum_action": action(f),
        "log_partition_pseudodeterminant": log_z_pseudo,
        "log_partition_basis": log_z_basis,
        "partition_formula_residual": abs(log_z_pseudo - log_z_basis),
        "gaussian_differential_entropy": entropy,
        "energy_derivative_analytic": analytic_derivative,
        "energy_derivative_finite_difference": finite_derivative,
        "energy_derivative_residual": abs(analytic_derivative - finite_derivative),
    }


def program85_quotient_action_entropy() -> dict:
    n = 12
    rng = np.random.default_rng(SEED)
    realizations = {}
    for name, kernel in [("strict", k_strict), ("legacy_absolute", k_legacy)]:
        w = radial_matrix(kernel, n, absolute=True)
        w /= w.sum(axis=1)[0]
        realizations[name] = quotient_action_record(laplacian(w), rng)
    return {
        "realizations": realizations,
        "theorems": [
            "For J perpendicular to 1, S[f]=1/2<f,Lf>-<J,f> is gauge invariant under f->f+c1.",
            "The stationary class is f*=L^+J+c1 and is unique on 1^perp.",
            "The quotient Gaussian partition uses pdet(L), not an infrared mass.",
            "For df/dt=-Lf, d(1/2<f,Lf>)/dt=-||Lf||^2<=0.",
        ],
        "status": "Proven quotient variational structure; finite verification",
    }


def composition_counts(total_units: int, parts: int):
    for cuts in itertools.combinations(range(1, total_units), parts - 1):
        points = (0,) + cuts + (total_units,)
        yield tuple(points[i + 1] - points[i] for i in range(parts))


def fisher_matrix(vectors: dict, counts: dict, sigmas: dict) -> np.ndarray:
    f = np.zeros((5, 5))
    for name, v in vectors.items():
        f += counts[name] * np.outer(v, v) / sigmas[name] ** 2
    return f


def program86_robust_calibration() -> dict:
    vectors = {
        "length_primary": np.array([1.0, 0.0, 0.0, 0.0, 1.0]),
        "length_crosscheck": np.array([1.0, 0.0, 0.0, 0.0, 0.0]),
        "clock_primary": np.array([0.0, 1.0, 0.0, 1.0, 0.0]),
        "clock_crosscheck": np.array([0.0, 1.0, 0.0, 0.0, 0.0]),
        "energy_standard": np.array([0.0, -1.0, 1.0, 0.0, 0.0]),
        "length_time_link": np.array([1.0, -1.0, 0.0, 0.0, 0.0]),
    }
    nominal = {
        "length_primary": 0.010,
        "length_crosscheck": 0.015,
        "clock_primary": 0.005,
        "clock_crosscheck": 0.008,
        "energy_standard": 0.020,
        "length_time_link": 0.012,
    }
    scenarios = {"nominal": nominal}
    for group, names in {
        "length_degraded": ["length_primary", "length_crosscheck"],
        "clock_degraded": ["clock_primary", "clock_crosscheck"],
        "energy_degraded": ["energy_standard"],
        "link_degraded": ["length_time_link"],
    }.items():
        s = nominal.copy()
        for name in names:
            s[name] *= 1.8
        scenarios[group] = s

    names = list(vectors)
    best = None
    evaluated = 0
    # Counts are even; 30 units correspond to 60 records with >=2 per class.
    for units in composition_counts(30, len(names)):
        counts = {name: 2 * units[i] for i, name in enumerate(names)}
        logdets = {}
        valid = True
        for scenario, sigmas in scenarios.items():
            sign, logdet = np.linalg.slogdet(
                fisher_matrix(vectors, counts, sigmas)
            )
            if sign <= 0:
                valid = False
                break
            logdets[scenario] = float(logdet)
        if valid:
            score = min(logdets.values())
            if best is None or score > best[0]:
                best = (score, counts, logdets)
        evaluated += 1
    assert best is not None
    score, counts, logdets = best
    worst_name = min(logdets, key=logdets.get)
    worst_f = fisher_matrix(vectors, counts, scenarios[worst_name])
    worst_cov = np.linalg.inv(worst_f)

    group_ranks = {}
    for group, omitted in {
        "omit_both_length_classes": ["length_primary", "length_crosscheck"],
        "omit_both_clock_classes": ["clock_primary", "clock_crosscheck"],
        "omit_energy_class": ["energy_standard"],
        "omit_link_only": ["length_time_link"],
    }.items():
        f = np.zeros((5, 5))
        for name, v in vectors.items():
            if name not in omitted:
                f += np.outer(v, v) / nominal[name] ** 2
        group_ranks[group] = int(np.linalg.matrix_rank(f))
    return {
        "parameters": [
            "log_length",
            "log_time",
            "log_action",
            "clock_drift",
            "distance_bias",
        ],
        "measurement_jacobians": {k: v.tolist() for k, v in vectors.items()},
        "records": 60,
        "allocations_evaluated": evaluated,
        "maximin_d_optimal_even_allocation": counts,
        "scenario_log_determinants": logdets,
        "worst_scenario": worst_name,
        "worst_case_logdet": score,
        "worst_case_one_sigma_bounds": np.sqrt(np.diag(worst_cov)).tolist(),
        "worst_case_condition_number": float(np.linalg.cond(worst_f)),
        "ranks_after_group_omissions": group_ranks,
        "status": "Conditional robust-design computation with explicit nuisances",
    }


def landau_flow(q0: float, bias: float, duration: float = 20.0) -> float:
    sol = solve_ivp(
        # The mobility factor keeps the density-mixture coordinate q in [-1,1].
        lambda _t, y: (1.0 - y**2) * (y + bias),
        (0.0, duration),
        [q0],
        rtol=1e-10,
        atol=1e-12,
    )
    return float(sol.y[0, -1])


def program87_chiral_state_stability() -> dict:
    n = 12
    w = radial_matrix(k_strict, n)
    psi_p = np.exp(2j * math.pi * np.arange(n) / n) / math.sqrt(n)
    rho_p = np.outer(psi_p, psi_p.conj())
    lambda_plus = float(
        sum(
            2.0 * np.imag(rho_p[i, (i + 1) % n] * w[(i + 1) % n, i])
            for i in range(n)
        )
    )
    initial = [-0.2, -1e-3, 0.0, 1e-3, 0.2]
    unbiased = [
        {
            "q_initial": q0,
            "q_final": landau_flow(q0, 0.0),
            "chiral_current_final": lambda_plus * landau_flow(q0, 0.0),
        }
        for q0 in initial
    ]
    rng = np.random.default_rng(SEED)
    ensemble = rng.normal(0.0, 1e-3, size=20000)
    positive_fraction_unbiased = float(np.mean(ensemble > 0.0))
    biased_final = np.array(
        [landau_flow(float(q), 1.2, 12.0) for q in ensemble[:1000]]
    )

    q_grid = np.linspace(-1.5, 1.5, 500)
    fig, ax = plt.subplots(figsize=(8.4, 4.6), constrained_layout=True)
    for bias, label in [
        (0.0, "reflection symmetric"),
        (1.2, "explicit bias h=1.2"),
    ]:
        potential = (
            0.25 * q_grid**4
            + (bias / 3.0) * q_grid**3
            - 0.5 * q_grid**2
            - bias * q_grid
        )
        ax.plot(q_grid, potential, label=label)
    ax.set_xlabel("chiral population imbalance q")
    ax.set_ylabel("Landau potential")
    ax.set_title("Two symmetry-related states versus imported sign bias")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program87_chiral_state_stability.png", dpi=190)
    plt.close(fig)
    return {
        "state_family": (
            "rho(q)=(1+q)/2 rho_(+1)+(1-q)/2 rho_(-1), q in [-1,1]"
        ),
        "strict_lambda_for_k_plus_1": lambda_plus,
        "unbiased_gradient_flow": "dq/dt=q(1-q^2)",
        "unbiased_trajectories": unbiased,
        "symmetric_noise_positive_branch_fraction": positive_fraction_unbiased,
        "explicit_bias": {
            "h": 1.2,
            "sample_size": len(biased_final),
            "positive_final_fraction": float(np.mean(biased_final > 0.0)),
            "mean_final_q": float(np.mean(biased_final)),
        },
        "theorem": (
            "Reflection-equivariant gradient flow has paired stable minima "
            "q=+1 and q=-1; it cannot select one globally. A sufficiently "
            "large explicit bias h>1 makes +1 the unique stable interior-reachable "
            "branch while preserving the physical interval."
        ),
        "status": "Proven symmetry statement and numerical stability test",
        "guardrail": (
            "The Landau state law and h are not strict FIN exports. Spontaneous "
            "branching does not discharge QW-2191."
        ),
    }


def binary_entropy(eps: float) -> float:
    return -(1.0 - eps) * math.log(1.0 - eps) - eps * math.log(eps)


def feedback_row(eps: float) -> dict:
    delta = math.log((1.0 - eps) / eps)
    mutual = math.log(2.0) - binary_entropy(eps)
    delta_f = math.log(2.0) - math.log(1.0 + math.exp(-delta))
    mean_work = eps * delta
    # Two trajectory classes: correct and incorrect binary measurement.
    info_correct = math.log(2.0 * (1.0 - eps))
    info_wrong = math.log(2.0 * eps)
    generalized_jarzynski = (
        (1.0 - eps) * math.exp(delta_f - info_correct)
        + eps * math.exp(-delta + delta_f - info_wrong)
    )
    return {
        "measurement_error": eps,
        "feedback_gap_delta": delta,
        "mutual_information": mutual,
        "free_energy_difference": delta_f,
        "mean_feedback_work": mean_work,
        "generalized_second_law_rhs_deltaF_minus_I": delta_f - mutual,
        "second_law_residual": mean_work - (delta_f - mutual),
        "generalized_jarzynski_average": generalized_jarzynski,
        "generalized_jarzynski_residual": abs(generalized_jarzynski - 1.0),
        "memory_entropy": math.log(2.0),
        "conditional_memory_entropy": binary_entropy(eps),
        "memory_reset_saving_equals_I_residual": abs(
            math.log(2.0) - binary_entropy(eps) - mutual
        ),
    }


def program88_feedback_thermodynamics() -> dict:
    rows = [feedback_row(e) for e in [0.01, 0.05, 0.1, 0.2, 0.4]]
    fig, ax = plt.subplots(figsize=(8.4, 4.6), constrained_layout=True)
    ax.plot(
        [r["measurement_error"] for r in rows],
        [r["free_energy_difference"] for r in rows],
        "o-",
        label="Delta F",
    )
    ax.plot(
        [r["measurement_error"] for r in rows],
        [r["mean_feedback_work"] for r in rows],
        "s-",
        label="mean work",
    )
    ax.plot(
        [r["measurement_error"] for r in rows],
        [r["mutual_information"] for r in rows],
        "^-",
        label="mutual information",
    )
    ax.set_xlabel("binary measurement error")
    ax.set_ylabel("dimensionless information/work")
    ax.set_title("Reversible feedback saturates <W> = Delta F - I")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program88_feedback_thermodynamics.png", dpi=190)
    plt.close(fig)
    return {
        "protocol": (
            "Equiprobable two-state system, binary-symmetric measurement, "
            "outcome-conditioned energy raising with Delta=log((1-eps)/eps)"
        ),
        "rows": rows,
        "status": "Exact conditional feedback identity",
        "scope": (
            "The measurement channel, Hamiltonian feedback, bath beta=1, "
            "apparatus memory, and erasure convention are imported. Including "
            "memory reset prevents mutual information from becoming free energy."
        ),
    }


def dephase(rho: np.ndarray) -> np.ndarray:
    return np.diag(np.diag(rho))


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.asarray(p, dtype=float).ravel()
    q = np.asarray(q, dtype=float).ravel()
    p /= p.sum()
    q /= q.sum()
    m = 0.5 * (p + q)
    keep_p = p > 0
    keep_q = q > 0
    return float(
        0.5 * np.sum(p[keep_p] * np.log(p[keep_p] / m[keep_p]))
        + 0.5 * np.sum(q[keep_q] * np.log(q[keep_q] / m[keep_q]))
    )


def program89_process_tensor_witness() -> dict:
    n = 12
    w = radial_matrix(k_strict, n, absolute=True)
    w /= w.sum(axis=1)[0]
    l = laplacian(w)
    rho0 = np.zeros((n, n), dtype=complex)
    rho0[0, 0] = 1.0
    p0 = np.zeros(n)
    p0[0] = 1.0
    rows = []
    for half_time in [0.1, 0.25, 0.5]:
        u = expm(-1j * half_time * l)
        p = expm(-half_time * l)
        rho_mid = u @ rho0 @ u.conj().T
        wave_no = np.real(np.diag(u @ rho_mid @ u.conj().T))
        wave_dephase = np.real(np.diag(u @ dephase(rho_mid) @ u.conj().T))
        diff_mid = p @ p0
        diff_no = p @ diff_mid
        diff_dephase = p @ diff_mid
        joint_wave = np.zeros((n, n))
        joint_diff = np.zeros((n, n))
        for y in range(n):
            joint_wave[y, :] = abs(u[y, 0]) ** 2 * abs(u[:, y]) ** 2
            joint_diff[y, :] = p[y, 0] * p[:, y]
        rows.append(
            {
                "half_time": half_time,
                "wave_dephasing_intervention_tv": float(
                    0.5 * np.sum(np.abs(wave_no - wave_dephase))
                ),
                "diffusion_dephasing_intervention_tv": float(
                    0.5 * np.sum(np.abs(diff_no - diff_dephase))
                ),
                "joint_intermediate_final_js_divergence": js_divergence(
                    joint_wave, joint_diff
                ),
                "wave_joint_normalization_residual": abs(joint_wave.sum() - 1.0),
                "diffusion_joint_normalization_residual": abs(
                    joint_diff.sum() - 1.0
                ),
            }
        )
    return {
        "two_slot_map": (
            "T[M]=E_(t2) o M o E_(t1)(rho0), with M equal to identity, "
            "nonselective site dephasing, or a recorded site instrument"
        ),
        "rows": rows,
        "status": "Finite operational process-tensor witness",
        "interpretation": (
            "Nonselective dephasing changes the wave branch but is the identity "
            "on the classical diffusion branch. The recorded two-time joint "
            "law provides a stronger operational object than a one-time spectrum."
        ),
    }


def wave_model(l: np.ndarray, times: list[float]) -> dict:
    return {t: np.abs(expm(-1j * t * l)) ** 2 for t in times}


def log_likelihood(counts: np.ndarray, probabilities: np.ndarray) -> float:
    p = np.maximum(probabilities, 1e-300)
    return float(counts @ np.log(p))


def model_selection_dry_run(
    true_models: dict,
    candidates: dict,
    times: list[float],
    preparations: list[int],
    shots: int,
    trials: int,
    rng: np.random.Generator,
) -> dict:
    wins = {truth: 0 for truth in true_models}
    margins = {truth: [] for truth in true_models}
    for truth_name, truth in true_models.items():
        for _ in range(trials):
            scores = {name: 0.0 for name in candidates}
            for time in times:
                for source in preparations:
                    pt = truth[time][:, source]
                    counts = rng.multinomial(shots, pt / pt.sum())
                    for name, model in candidates.items():
                        scores[name] += log_likelihood(
                            counts, model[time][:, source]
                        )
            winner = max(scores, key=scores.get)
            wins[truth_name] += int(winner == truth_name)
            other = "nearest_neighbour_null" if truth_name == "strict" else "strict"
            margins[truth_name].append(scores[truth_name] - scores[other])
    return {
        truth: {
            "correct_selection_rate": wins[truth] / trials,
            "mean_log_likelihood_margin": float(np.mean(margins[truth])),
            "minimum_log_likelihood_margin": float(np.min(margins[truth])),
        }
        for truth in true_models
    }


def program90_external_data_admission() -> dict:
    required = [
        "independent_producer",
        "preparation_map",
        "instrument_or_POVM",
        "clock_calibration",
        "length_calibration",
        "environment_and_boundaries",
        "raw_counts",
        "license",
        "frozen_split_assignment",
    ]
    template_core = {
        "schema": "FIN external-data intake 1.0",
        "required_fields": required,
        "dataset": {
            "identifier": None,
            "independent_producer": None,
            "preparation_map": None,
            "instrument_or_POVM": None,
            "clock_calibration": None,
            "length_calibration": None,
            "environment_and_boundaries": None,
            "raw_counts": None,
            "license": None,
            "frozen_split_assignment": None,
        },
        "admission_rule": "All required fields must be non-null before unblinding.",
        "current_state": "template only; no dataset admitted",
    }
    canonical = json.dumps(
        template_core, sort_keys=True, separators=(",", ":")
    ).encode()
    digest = hashlib.sha256(canonical).hexdigest()
    template = {**template_core, "canonical_core_sha256": digest}
    INTAKE.write_text(json.dumps(template, indent=2) + "\n", encoding="utf-8")

    intake_dir = ROOT / "external_data"
    candidate_manifests = (
        sorted(intake_dir.glob("*.json")) if intake_dir.is_dir() else []
    )
    admitted = []
    rejected = []
    for path in candidate_manifests:
        data = json.loads(path.read_text(encoding="utf-8"))
        missing = [field for field in required if not data.get(field)]
        if missing:
            rejected.append({"file": path.name, "missing": missing})
        else:
            admitted.append(path.name)

    n = 12
    strict_w = radial_matrix(k_strict, n, absolute=True)
    strict_w /= strict_w.sum(axis=1)[0]
    nn_w = np.zeros((n, n))
    for i in range(n):
        nn_w[i, (i - 1) % n] = 0.5
        nn_w[i, (i + 1) % n] = 0.5
    times = [0.2, 0.5, 1.0]
    strict = wave_model(laplacian(strict_w), times)
    null = wave_model(laplacian(nn_w), times)
    rng = np.random.default_rng(SEED)
    dry = model_selection_dry_run(
        {"strict": strict, "nearest_neighbour_null": null},
        {"strict": strict, "nearest_neighbour_null": null},
        times,
        [0, 3, 6],
        100,
        300,
        rng,
    )
    return {
        "intake_template": INTAKE.name,
        "template_sha256": digest,
        "intake_directory_exists": intake_dir.is_dir(),
        "candidate_manifest_count": len(candidate_manifests),
        "admitted_external_datasets": admitted,
        "rejected_candidate_manifests": rejected,
        "synthetic_runner_dry_run": {
            "warning": "Method validation only; not physical evidence",
            "times": times,
            "preparations": [0, 3, 6],
            "shots_per_record": 100,
            "trials_per_truth": 300,
            "results": dry,
        },
        "status": (
            "No external dataset admitted; protocol runner passes a synthetic "
            "strict-versus-null directionality check"
        ),
    }


def main() -> None:
    results = {
        "release": "10.9",
        "title": "Programs 81--90: Asymptotic and Operational Completion Tests",
        "date": "2026-07-27",
        "author": "Krzysztof Żuchowski",
        "orcid": "0009-0002-0909-3613",
        "ontology_guardrail": (
            "The nadsoliton remains primordial fractal information in a "
            "solitonic state; no lower information layer is introduced."
        ),
        "program81": program81_asymptotic_symbol_audit(),
        "program82": program82_fixed_symbol_classification(),
        "program83": program83_dense_kernel_universality(),
        "program84": program84_locality_bounds(),
        "program85": program85_quotient_action_entropy(),
        "program86": program86_robust_calibration(),
        "program87": program87_chiral_state_stability(),
        "program88": program88_feedback_thermodynamics(),
        "program89": program89_process_tensor_witness(),
        "program90": program90_external_data_admission(),
        "global_guardrail": (
            "No result closes QW-2191, derives a dimensional standard, proves "
            "Lorentz invariance, validates a legacy-to-strict completion map, "
            "transfers legacy physical roles, or supplies external physical evidence."
        ),
    }
    OUT.write_text(
        json.dumps(results, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
