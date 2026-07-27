#!/usr/bin/env python3
"""Execute FIN research Programs 71--80.

The suite continues Programs 61--70 without promoting finite computations to
physical closure.  It studies:

71. exact Fourier-symbol Schur elimination and an analytic nearest-neighbour RG;
72. large-N projective/naturality defects under two declared scaling semantics;
73. regulator-free Gaussian information geometry on the zero-mode quotient;
74. the quantitative trade-off between kernel truncation and propagation order;
75. noisy D-optimal calibration of length, time, and action conversion scales;
76. a state-dependent chiral current and its unresolved sign source;
77. a finite-time two-state Landauer/Jarzynski protocol;
78. shot-noise-limited wave-versus-diffusion process discrimination;
79. a fixed, no-fit legacy-to-strict Schur-flow bridge test; and
80. a cryptographically frozen external-data preregistration protocol.

All FIN kernel inputs are dimensionless.  Absolute-value repairs, coordinate
rescalings, baths, clocks, standards, preparations, and instruments are
explicit conditioning data and are not strict FIN exports.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_71_80_Symbol_Operational_Bridge_Results.json"
PREREG = ROOT / "FIN_Programs_71_80_External_Data_Preregistration.json"
FIG = ROOT / "FIN_Programs_71_80_Symbol_Operational_Bridge_Figures"
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


def radial_matrix(kernel, n: int, absolute: bool = False) -> np.ndarray:
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j:
                w[i, j] = float(kernel(cycle_distance(i, j, n)))
    return np.abs(w) if absolute else w


def laplacian(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(axis=1)) - w


def first_row_precision(
    kernel,
    n: int,
    *,
    coordinate_scaled: bool = False,
    mass2: float = MASS2,
) -> np.ndarray:
    """Positive, row-normalized, absolute-value precision on a cycle."""
    d = np.minimum(np.arange(n), n - np.arange(n)).astype(float)
    arg = d / n if coordinate_scaled else d
    weights = np.abs(np.asarray(kernel(arg), dtype=float))
    weights[0] = 0.0
    weights /= weights.sum()
    row = -weights
    row[0] = mass2 + 1.0
    return row


def circulant_from_row(row: np.ndarray) -> np.ndarray:
    n = len(row)
    return np.array([np.roll(row, i) for i in range(n)])


def normalize_row(row: np.ndarray) -> np.ndarray:
    return np.asarray(row, dtype=float) / float(np.real(row[0]))


def schur_dense(a: np.ndarray, keep: np.ndarray) -> np.ndarray:
    keep = np.asarray(keep, dtype=int)
    keep_set = set(map(int, keep))
    drop = np.array([i for i in range(a.shape[0]) if i not in keep_set])
    akk = a[np.ix_(keep, keep)]
    ako = a[np.ix_(keep, drop)]
    aoo = a[np.ix_(drop, drop)]
    return akk - ako @ np.linalg.solve(aoo, ako.T)


def schur_row_symbol(parent_row: np.ndarray) -> np.ndarray:
    """Eliminate odd sites from a real symmetric circulant precision."""
    if len(parent_row) % 2:
        raise ValueError("The parent cycle must have even cardinality")
    even_row = np.asarray(parent_row[::2], dtype=float)
    odd_row = np.asarray(parent_row[1::2], dtype=float)
    a_symbol = np.fft.fft(even_row)
    b_symbol = np.fft.fft(odd_row)
    if np.min(np.abs(a_symbol)) < 1e-14:
        raise np.linalg.LinAlgError("Odd block is singular at a Fourier mode")
    schur_symbol = a_symbol - np.abs(b_symbol) ** 2 / a_symbol
    row = np.fft.ifft(schur_symbol).real
    return 0.5 * (row + np.r_[row[0], row[:0:-1]])


def row_operator_defect(row_a: np.ndarray, row_b: np.ndarray) -> float:
    return float(np.linalg.norm(row_a - row_b) / np.linalg.norm(row_b))


def row_green_defect(row_a: np.ndarray, row_b: np.ndarray) -> float:
    la = np.fft.fft(row_a).real
    lb = np.fft.fft(row_b).real
    return float(np.linalg.norm(1.0 / la - 1.0 / lb) / np.linalg.norm(1.0 / lb))


def program71_fourier_schur_theorem() -> dict:
    parent_row = first_row_precision(k_strict, 24)
    dense = circulant_from_row(parent_row)
    dense_s = schur_dense(dense, np.arange(0, 24, 2))
    symbol_row = schur_row_symbol(parent_row)
    symbol_s = circulant_from_row(symbol_row)
    dense_symbol_relative_residual = float(
        np.linalg.norm(dense_s - symbol_s, "fro") / np.linalg.norm(dense_s, "fro")
    )

    # Exact nearest-neighbour massive-chain renormalization.
    m, c = 0.37, 1.21
    nn_row = np.zeros(24)
    nn_row[0] = m + 2.0 * c
    nn_row[1] = nn_row[-1] = -c
    nn_s = schur_row_symbol(nn_row)
    c_prime = c * c / (m + 2.0 * c)
    m_prime = m * (m + 4.0 * c) / (m + 2.0 * c)
    predicted = np.zeros(12)
    predicted[0] = m_prime + 2.0 * c_prime
    predicted[1] = predicted[-1] = -c_prime
    nn_relative_residual = float(
        np.linalg.norm(nn_s - predicted) / np.linalg.norm(predicted)
    )
    ratio = m / c
    ratio_prime_numeric = m_prime / c_prime
    ratio_prime_formula = ratio * (ratio + 4.0)

    return {
        "fourier_block_formula": (
            "For the even/odd block symbol B(k), "
            "s(k)=B00(k)-B01(k) B11(k)^(-1) B10(k)."
        ),
        "dense_symbol_relative_residual": dense_symbol_relative_residual,
        "nearest_neighbour_exact_map": {
            "m": m,
            "c": c,
            "c_prime": c_prime,
            "m_prime": m_prime,
            "ratio_m_over_c": ratio,
            "ratio_prime_numeric": ratio_prime_numeric,
            "ratio_prime_formula_r_times_r_plus_4": ratio_prime_formula,
            "relative_residual": nn_relative_residual,
            "nonnegative_fixed_points": [0.0],
            "linearized_derivative_at_zero": 4.0,
        },
        "status": "Proven (finite Fourier diagonalization and exact algebra)",
        "interpretation": (
            "Schur elimination is an exact symbol-level operation.  Exactness of "
            "this operation does not imply equality with any native coarse FIN "
            "kernel; that separate naturality question is tested in Program 72."
        ),
    }


def projective_row_defect(kernel, n: int, coordinate_scaled: bool) -> dict:
    parent = first_row_precision(
        kernel, 2 * n, coordinate_scaled=coordinate_scaled
    )
    compressed = normalize_row(schur_row_symbol(parent))
    native = normalize_row(
        first_row_precision(kernel, n, coordinate_scaled=coordinate_scaled)
    )
    return {
        "retained_n": n,
        "operator_defect": row_operator_defect(compressed, native),
        "green_defect": row_green_defect(compressed, native),
    }


def program72_large_n_projective_scaling() -> dict:
    sizes = [48, 96, 192, 384, 768, 1536]
    cases = {
        "strict_lattice_distance": (k_strict, False),
        "legacy_absolute_lattice_distance": (k_legacy, False),
        "strict_coordinate_rescaled": (k_strict, True),
        "legacy_absolute_coordinate_rescaled": (k_legacy, True),
    }
    results = {}
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.5), constrained_layout=True)
    for name, (kernel, scaled) in cases.items():
        rows = [projective_row_defect(kernel, n, scaled) for n in sizes]
        op = np.array([r["operator_defect"] for r in rows])
        green = np.array([r["green_defect"] for r in rows])
        op_slope = float(np.polyfit(np.log(sizes), np.log(op), 1)[0])
        green_slope = float(np.polyfit(np.log(sizes), np.log(green), 1)[0])
        results[name] = {
            "rows": rows,
            "operator_loglog_slope": op_slope,
            "green_loglog_slope": green_slope,
            "operator_monotone_decreasing": bool(np.all(np.diff(op) < 0)),
            "green_monotone_decreasing": bool(np.all(np.diff(green) < 0)),
            "final_operator_defect": float(op[-1]),
            "final_green_defect": float(green[-1]),
        }
        label = name.replace("_", " ")
        axes[0].loglog(sizes, op, "o-", label=label)
        axes[1].loglog(sizes, green, "s-", label=label)
    for ax, title in zip(
        axes,
        ["large-N operator naturality", "large-N Green naturality"],
    ):
        ax.set_xlabel("retained cycle size N")
        ax.set_ylabel("relative defect")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.25)
        ax.legend(fontsize=7)
    fig.savefig(FIG / "program72_large_n_projective_scaling.png", dpi=190)
    plt.close(fig)
    return {
        "sizes": sizes,
        "fixed_conventions": {
            "positive_realization": "absolute-value weights",
            "row_normalization": "weights sum to one",
            "mass_squared": MASS2,
            "post_compression_normalization": "unit diagonal",
        },
        "cases": results,
        "semantic_warning": (
            "K(d) and K(d/N) are different model families.  Coordinate "
            "rescaling is imported lattice semantics, not a FIN-derived length."
        ),
        "status": "Numerical evidence under declared discretizations",
    }


def gaussian_kl_from_covariances(c0: np.ndarray, c1: np.ndarray) -> float:
    """KL[N(0,c0)||N(0,c1)]."""
    k = c0.shape[0]
    sign0, ld0 = np.linalg.slogdet(c0)
    sign1, ld1 = np.linalg.slogdet(c1)
    if sign0 <= 0 or sign1 <= 0:
        raise ValueError("Covariances must be positive definite")
    return float(0.5 * (np.trace(np.linalg.solve(c1, c0)) - k + ld1 - ld0))


def program73_zero_mode_quotient_information() -> dict:
    n = 12
    ls = laplacian(radial_matrix(k_strict, n, absolute=True))
    ll = laplacian(radial_matrix(k_legacy, n, absolute=True))
    q, _ = np.linalg.qr(
        np.column_stack(
            [np.ones(n) / math.sqrt(n), np.eye(n)[:, 1:]]
        )
    )
    v = q[:, 1:]
    a_s = v.T @ ls @ v
    a_l = v.T @ ll @ v
    c_s = np.linalg.inv(a_s)
    c_l = np.linalg.inv(a_l)
    quotient_kl_s_to_l = gaussian_kl_from_covariances(c_s, c_l)
    quotient_kl_l_to_s = gaussian_kl_from_covariances(c_l, c_s)

    rows = []
    for delta in np.logspace(-1, -8, 8):
        csd = np.linalg.inv(ls + delta * np.eye(n))
        cld = np.linalg.inv(ll + delta * np.eye(n))
        kl = gaussian_kl_from_covariances(csd, cld)
        rows.append(
            {
                "delta": float(delta),
                "full_shifted_kl_strict_to_legacy": kl,
                "absolute_error_to_quotient": abs(kl - quotient_kl_s_to_l),
            }
        )

    eig_s = np.linalg.eigvalsh(ls)[1:]
    eig_l = np.linalg.eigvalsh(ll)[1:]
    fig, ax = plt.subplots(figsize=(7.8, 4.5), constrained_layout=True)
    ax.loglog(
        [r["delta"] for r in rows],
        [r["absolute_error_to_quotient"] for r in rows],
        "o-",
    )
    ax.invert_xaxis()
    ax.set_xlabel(r"common shift $\delta\to0^+$")
    ax.set_ylabel("absolute KL error")
    ax.set_title(r"Regulator-free Gaussian KL on $\mathbf{1}^{\perp}$")
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program73_zero_mode_quotient_kl.png", dpi=190)
    plt.close(fig)
    return {
        "quotient_dimension": n - 1,
        "strict_nonzero_eigenvalue_min": float(eig_s.min()),
        "legacy_absolute_nonzero_eigenvalue_min": float(eig_l.min()),
        "kl_strict_to_legacy_on_quotient": quotient_kl_s_to_l,
        "kl_legacy_to_strict_on_quotient": quotient_kl_l_to_s,
        "shifted_limit_rows": rows,
        "last_absolute_error": rows[-1]["absolute_error_to_quotient"],
        "status": "Proven quotient construction; numerically verified limit",
        "scope": (
            "The construction removes the common constant zero mode without a "
            "mass regulator.  It still uses a declared positive absolute-value "
            "realization of the legacy weights."
        ),
    }


def first_nonzero_power(a: np.ndarray, source: int, target: int) -> int | None:
    power = np.eye(a.shape[0])
    scale = max(1.0, np.linalg.norm(a, 2))
    for m in range(1, a.shape[0] + 1):
        power = power @ a
        if abs(power[target, source]) > 1e-11 * scale**m:
            return m
    return None


def tv(p: np.ndarray, q: np.ndarray) -> float:
    return float(0.5 * np.sum(np.abs(p - q)))


def program74_locality_truncation() -> dict:
    n = 12
    source, target = 0, n // 2
    w_full = radial_matrix(k_strict, n, absolute=True)
    w_full /= w_full.sum(axis=1)[0]
    l_full = laplacian(w_full)
    wave_full = np.abs(expm(-1j * l_full)[:, source]) ** 2
    diff_full = expm(-l_full)[:, source]
    rows = []
    for radius in range(1, n // 2 + 1):
        mask = np.zeros_like(w_full, dtype=bool)
        for i in range(n):
            for j in range(n):
                mask[i, j] = i != j and cycle_distance(i, j, n) <= radius
        w = np.where(mask, w_full, 0.0)
        l = laplacian(w)
        m = first_nonzero_power(l, source, target)
        wave = np.abs(expm(-1j * l)[:, source]) ** 2
        diff = expm(-l)[:, source]
        rows.append(
            {
                "truncation_radius": radius,
                "operator_relative_error": float(
                    np.linalg.norm(l - l_full, "fro")
                    / np.linalg.norm(l_full, "fro")
                ),
                "first_nonzero_matrix_power_to_opposite_site": m,
                "leading_wave_probability_power": None if m is None else 2 * m,
                "leading_diffusion_entry_power": m,
                "wave_output_tv_at_t1": tv(wave, wave_full),
                "diffusion_output_tv_at_t1": tv(diff, diff_full),
            }
        )
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    radii = [r["truncation_radius"] for r in rows]
    axes[0].semilogy(
        radii, [r["operator_relative_error"] for r in rows], "o-"
    )
    axes[0].set_xlabel("interaction radius")
    axes[0].set_ylabel("operator relative error")
    axes[0].set_title("Approximation to the full strict operator")
    axes[1].plot(
        radii,
        [r["leading_wave_probability_power"] for r in rows],
        "o-",
        label="wave probability",
    )
    axes[1].plot(
        radii,
        [r["leading_diffusion_entry_power"] for r in rows],
        "s-",
        label="diffusion entry",
    )
    axes[1].set_xlabel("interaction radius")
    axes[1].set_ylabel("leading small-time power")
    axes[1].set_title("Locality order to the opposite site")
    axes[1].legend()
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program74_locality_truncation.png", dpi=190)
    plt.close(fig)
    return {
        "rows": rows,
        "status": "Proven short-time order; numerical finite-time comparison",
        "interpretation": (
            "Finite-range truncation increases the perturbative order required "
            "to reach a distant site, but changes the operator and its t=1 "
            "outputs.  The untruncated radial kernel has immediate graph edges "
            "to every site and therefore does not by itself define a strict "
            "finite propagation cone."
        ),
    }


def program75_noisy_scale_calibration() -> dict:
    # Parameters are log(ell_*), log(tau_*), log(hbar_*).
    vectors = {
        "length_standard": np.array([1.0, 0.0, 0.0]),
        "clock_standard": np.array([0.0, 1.0, 0.0]),
        "energy_standard": np.array([0.0, -1.0, 1.0]),
    }
    sigmas = {
        "length_standard": 0.010,
        "clock_standard": 0.005,
        "energy_standard": 0.020,
    }
    total = 30
    candidates = []
    for nl in range(1, total - 1):
        for nt in range(1, total - nl):
            ne = total - nl - nt
            if ne < 1:
                continue
            counts = {
                "length_standard": nl,
                "clock_standard": nt,
                "energy_standard": ne,
            }
            fisher = np.zeros((3, 3))
            for name, count in counts.items():
                v = vectors[name]
                fisher += count * np.outer(v, v) / sigmas[name] ** 2
            sign, logdet = np.linalg.slogdet(fisher)
            cov = np.linalg.inv(fisher)
            candidates.append(
                (float(logdet if sign > 0 else -np.inf), counts, fisher, cov)
            )
    best = max(candidates, key=lambda x: x[0])
    logdet, counts, fisher, cov = best
    missing_rank = {}
    for omitted in vectors:
        f = np.zeros((3, 3))
        for name, v in vectors.items():
            if name != omitted:
                f += np.outer(v, v) / sigmas[name] ** 2
        missing_rank[omitted] = int(np.linalg.matrix_rank(f))
    return {
        "parameterization": ["log_length", "log_time", "log_action"],
        "measurement_jacobians": {
            k: v.tolist() for k, v in vectors.items()
        },
        "assumed_log_noise_standard_deviations": sigmas,
        "total_measurements": total,
        "d_optimal_allocation": counts,
        "fisher_matrix": fisher.tolist(),
        "covariance_lower_bound": cov.tolist(),
        "one_sigma_log_parameter_bounds": np.sqrt(np.diag(cov)).tolist(),
        "fisher_condition_number": float(np.linalg.cond(fisher)),
        "log_determinant": logdet,
        "rank_after_omitting_each_class": missing_rank,
        "status": "Proven identifiability; conditional numerical design",
        "scope": (
            "The design proves that three independent operational standards are "
            "sufficient and individually necessary for this parameterization. "
            "Their values and noise levels are external conditioning data; FIN "
            "does not generate an SI unit."
        ),
    }


def chiral_current_scalar(w: np.ndarray, rho: np.ndarray) -> float:
    n = w.shape[0]
    return float(
        sum(2.0 * np.imag(rho[i, (i + 1) % n] * w[(i + 1) % n, i]) for i in range(n))
    )


def reflection_matrix(n: int) -> np.ndarray:
    r = np.zeros((n, n))
    for i in range(n):
        r[i, (-i) % n] = 1.0
    return r


def shift_matrix(n: int) -> np.ndarray:
    s = np.zeros((n, n))
    for i in range(n):
        s[(i + 1) % n, i] = 1.0
    return s


def program76_state_dependent_chirality() -> dict:
    n = 12
    r = reflection_matrix(n)
    t = shift_matrix(n)
    kernels = {
        "strict": radial_matrix(k_strict, n),
        "legacy_raw_signed": radial_matrix(k_legacy, n),
    }
    results = {}
    for name, w in kernels.items():
        rows = []
        for k in [0, 1, -1]:
            psi = np.exp(2j * math.pi * k * np.arange(n) / n) / math.sqrt(n)
            rho = np.outer(psi, psi.conj())
            lam = chiral_current_scalar(w, rho)
            reflected = chiral_current_scalar(w, r @ rho @ r.T)
            shifted = chiral_current_scalar(w, t @ rho @ t.T)
            rows.append(
                {
                    "fourier_mode_k": k,
                    "lambda": lam,
                    "reflection_sum_residual": abs(reflected + lam),
                    "translation_invariance_residual": abs(shifted - lam),
                }
            )
        results[name] = rows
    return {
        "formula": (
            "J_ij=2 Im(rho_ij W_ji), "
            "Lambda(rho,W)=sum_i J_{i,i+1}."
        ),
        "kernels": results,
        "status": "Proven covariance properties; source remains open",
        "interpretation": (
            "A nonradial state can supply a signed pseudoscalar current even "
            "when the kernel is real and radial.  Reflection reverses the sign "
            "and translation preserves it.  However +k and -k remain paired: "
            "the selector obstruction is transferred to the preparation/state "
            "law unless a strict chiral nadsoliton state is independently "
            "derived.  QW-2191 is therefore not discharged."
        ),
    }


def two_state_protocol(duration: float, steps: int = 600, gamma: float = 3.0):
    beta = 1.0
    delta_e = 5.0
    dt = duration / steps
    p = np.array([0.5, 0.5])
    first = np.zeros(2)
    second = np.zeros(2)
    jarzynski_weight = p.copy()
    previous_e = 0.0
    for step in range(1, steps + 1):
        energy = delta_e * step / steps
        work = energy - previous_e
        first[1] += work * p[1]
        second[1] += 2.0 * work * (first[1] - work * p[1]) + work**2 * p[1]
        jarzynski_weight[1] *= math.exp(-beta * work)

        gibbs = np.array([1.0, math.exp(-beta * energy)])
        gibbs /= gibbs.sum()
        relax = 1.0 - math.exp(-gamma * dt)
        transition = (1.0 - relax) * np.eye(2) + relax * gibbs[:, None]
        p = transition @ p
        first = transition @ first
        second = transition @ second
        jarzynski_weight = transition @ jarzynski_weight
        previous_e = energy
    mean = float(first.sum())
    second_moment = float(second.sum())
    delta_f = math.log(2.0) - math.log(1.0 + math.exp(-delta_e))
    return {
        "duration": duration,
        "mean_work": mean,
        "work_variance": max(0.0, second_moment - mean**2),
        "free_energy_difference": delta_f,
        "dissipated_work": mean - delta_f,
        "jarzynski_average_exp_minus_work": float(jarzynski_weight.sum()),
        "exp_minus_delta_free_energy": math.exp(-delta_f),
        "jarzynski_absolute_residual": abs(
            float(jarzynski_weight.sum()) - math.exp(-delta_f)
        ),
        "final_excited_probability": float(p[1]),
    }


def program77_finite_time_landauer_jarzynski() -> dict:
    durations = [0.1, 0.3, 1.0, 3.0, 10.0]
    rows = [two_state_protocol(t) for t in durations]
    fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.2), constrained_layout=True)
    axes[0].semilogx(
        durations, [r["mean_work"] for r in rows], "o-", label="mean work"
    )
    axes[0].axhline(
        rows[0]["free_energy_difference"],
        color="black",
        linestyle="--",
        label="free-energy difference",
    )
    axes[0].set_xlabel("protocol duration")
    axes[0].set_ylabel("dimensionless work")
    axes[0].set_title("Finite-time erasure cost")
    axes[0].legend()
    axes[1].loglog(
        durations,
        [r["jarzynski_absolute_residual"] for r in rows],
        "s-",
    )
    axes[1].set_xlabel("protocol duration")
    axes[1].set_ylabel("absolute Jarzynski residual")
    axes[1].set_title("Numerical fluctuation-identity check")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program77_landauer_jarzynski.png", dpi=190)
    plt.close(fig)
    return {
        "protocol": (
            "Two-state memory, beta=1, E1:0->5, Markov relaxation toward the "
            "instantaneous Gibbs distribution, 600 quenches."
        ),
        "rows": rows,
        "mean_work_lower_bound_satisfied": bool(
            all(r["dissipated_work"] >= -1e-11 for r in rows)
        ),
        "dissipation_monotone_with_duration": bool(
            np.all(np.diff([r["dissipated_work"] for r in rows]) < 0)
        ),
        "status": "Conditional theorem instance and numerical protocol",
        "scope": (
            "The bath, inverse temperature, energy schedule, Markov rule, and "
            "work convention are imported operational structure.  The result "
            "does not derive thermodynamic temperature or energy from FIN."
        ),
    }


def aggregate_probabilities(p: np.ndarray, instrument: str) -> np.ndarray:
    if instrument == "full":
        return p
    if instrument == "adjacent_pairs":
        return p.reshape(-1, 2).sum(axis=1)
    raise ValueError(instrument)


def process_models(n: int, times: list[float]) -> dict:
    w = radial_matrix(k_strict, n, absolute=True)
    w /= w.sum(axis=1)[0]
    l = laplacian(w)
    return {
        "wave": {t: np.abs(expm(-1j * t * l)) ** 2 for t in times},
        "diffusion": {t: expm(-t * l) for t in times},
    }


def process_classification_accuracy(
    models: dict,
    preps: list[int],
    times: list[float],
    instrument: str,
    shots: int,
    rng: np.random.Generator,
    trials: int = 300,
) -> float:
    correct = 0
    labels = ["wave", "diffusion"]
    for trial in range(trials):
        truth = labels[trial % 2]
        loglike = {label: 0.0 for label in labels}
        for source in preps:
            for time in times:
                p_true = aggregate_probabilities(
                    models[truth][time][:, source], instrument
                )
                counts = rng.multinomial(shots, p_true / p_true.sum())
                for label in labels:
                    pred = aggregate_probabilities(
                        models[label][time][:, source], instrument
                    )
                    pred = np.maximum(pred / pred.sum(), 1e-300)
                    loglike[label] += float(counts @ np.log(pred))
        chosen = max(labels, key=lambda x: loglike[x])
        correct += int(chosen == truth)
    return correct / trials


def program78_noisy_process_tomography() -> dict:
    n = 12
    # The short-time window avoids a trivial perfect-separation regime and
    # exposes the O(t) diffusion versus O(t^2) wave-probability distinction.
    times = [0.03, 0.06, 0.12]
    models = process_models(n, times)
    rng = np.random.default_rng(SEED)
    configs = {
        "one_prep_full": ([0], "full"),
        "one_prep_pair_coarse": ([0], "adjacent_pairs"),
        "three_preps_full": ([0, 3, 6], "full"),
        "three_preps_pair_coarse": ([0, 3, 6], "adjacent_pairs"),
    }
    shots_grid = [2, 10, 50]
    rows = []
    for name, (preps, instrument) in configs.items():
        for shots in shots_grid:
            rows.append(
                {
                    "configuration": name,
                    "preparations": preps,
                    "instrument": instrument,
                    "shots_per_preparation_time": shots,
                    "classification_accuracy": process_classification_accuracy(
                        models,
                        preps,
                        times,
                        instrument,
                        shots,
                        rng,
                        trials=500,
                    ),
                }
            )
    fig, ax = plt.subplots(figsize=(8.8, 4.8), constrained_layout=True)
    for name in configs:
        subset = [r for r in rows if r["configuration"] == name]
        ax.semilogx(
            [r["shots_per_preparation_time"] for r in subset],
            [r["classification_accuracy"] for r in subset],
            "o-",
            label=name.replace("_", " "),
        )
    ax.axhline(0.5, color="black", linestyle="--", linewidth=1)
    ax.set_xlabel("shots per preparation and time")
    ax.set_ylabel("wave/diffusion classification accuracy")
    ax.set_ylim(0.45, 1.02)
    ax.set_title("Instrument- and shot-noise-limited process discrimination")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program78_process_tomography.png", dpi=190)
    plt.close(fig)
    return {
        "times": times,
        "rows": rows,
        "status": "Synthetic operational numerical evidence",
        "scope": (
            "The computation shows which preparations, instruments, and shot "
            "counts distinguish two dynamics generated from the same operator. "
            "It is not external physical evidence and does not select which "
            "dynamics nature realizes."
        ),
    }


def program79_legacy_to_strict_schur_flow() -> dict:
    start_n = 1536
    legacy_row = normalize_row(first_row_precision(k_legacy, start_n))
    rows = []
    current = legacy_row
    current_n = start_n
    while current_n >= 96:
        current = normalize_row(schur_row_symbol(current))
        current_n //= 2
        strict_native = normalize_row(first_row_precision(k_strict, current_n))
        rows.append(
            {
                "retained_n": current_n,
                "schur_steps": int(round(math.log2(start_n / current_n))),
                "operator_defect_to_native_strict": row_operator_defect(
                    current, strict_native
                ),
                "green_defect_to_native_strict": row_green_defect(
                    current, strict_native
                ),
            }
        )
    op = np.array([r["operator_defect_to_native_strict"] for r in rows])
    green = np.array([r["green_defect_to_native_strict"] for r in rows])
    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.semilogx(
        [r["retained_n"] for r in rows],
        op,
        "o-",
        label="operator defect",
    )
    ax.semilogx(
        [r["retained_n"] for r in rows],
        green,
        "s-",
        label="Green defect",
    )
    ax.invert_xaxis()
    ax.set_xlabel("retained N along repeated Schur flow")
    ax.set_ylabel("relative defect to native strict kernel")
    ax.set_title("No-fit legacy-absolute to strict-absolute bridge test")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program79_legacy_strict_schur_flow.png", dpi=190)
    plt.close(fig)
    return {
        "start_n": start_n,
        "map": (
            "absolute-value row-normalized legacy precision; repeated even-site "
            "Schur elimination; unit-diagonal normalization; no fitted parameter"
        ),
        "rows": rows,
        "operator_defect_monotone_along_coarsening": bool(np.all(np.diff(op) < 0)),
        "green_defect_monotone_along_coarsening": bool(
            np.all(np.diff(green) < 0)
        ),
        "status": "Numerical falsification test under one fixed positive map",
        "guardrail": (
            "Even a decreasing defect would not establish the repository's "
            "legacy-to-strict completion map, phase/topological completion, "
            "selector closure, or physical-role transfer.  This experiment "
            "tests only the declared repaired operator flow."
        ),
    }


def program80_external_preregistration() -> dict:
    core = {
        "title": "FIN operator external-data falsification protocol",
        "version": "1.0.0-frozen",
        "created": "2026-07-27",
        "scope": (
            "A future external-data test of dimensionless operator-shape "
            "predictions after independently calibrated conversion maps."
        ),
        "candidate_models": {
            "strict": (
                "K(d)=cos(omega*d+phi)/(1+beta*d^eta), parameters frozen "
                "from repository provenance before unblinding"
            ),
            "legacy_intermediate": (
                "K(d)=alpha_geo*cos(omega*d+phi)/(1+beta_tors*d), retained "
                "only as an intermediate comparator"
            ),
            "nulls": [
                "nearest-neighbour Laplacian",
                "screened exponential radial kernel",
                "power-law radial kernel",
                "matched-spectrum phase-randomized operator",
            ],
        },
        "operational_requirements": [
            "record preparation map",
            "record instrument/POVM or coarse-graining map",
            "record clock and length calibration with uncertainty",
            "record environment and boundary conditions",
            "record raw counts before normalization",
        ],
        "primary_metrics": [
            "held-out predictive log likelihood",
            "calibration-marginalized KL divergence",
            "posterior predictive residual",
            "complexity-penalized model score",
        ],
        "split": {
            "training_fraction": 0.60,
            "validation_fraction": 0.20,
            "hidden_test_fraction": 0.20,
            "stratification": "by preparation, instrument, and time window",
        },
        "decision_rules": {
            "support": (
                "strict model beats every frozen null on hidden-test mean log "
                "likelihood with a preregistered uncertainty interval excluding 0"
            ),
            "falsification": (
                "a frozen null dominates strict, or strict residuals fail the "
                "preregistered calibration check"
            ),
            "inconclusive": "all other outcomes",
        },
        "forbidden_after_unblinding": [
            "changing kernel parameters",
            "changing distance or coordinate semantics",
            "adding observables",
            "changing binning or instrument maps",
            "dropping failed preparations",
            "using synthetic FIN-generated data as external evidence",
        ],
        "current_data_state": "No external dataset selected or inspected",
    }
    canonical = json.dumps(
        core, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    digest = hashlib.sha256(canonical).hexdigest()
    record = {**core, "canonical_core_sha256": digest}
    PREREG.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return {
        "preregistration_file": PREREG.name,
        "canonical_core_sha256": digest,
        "status": "Protocol frozen; no physical evidence claimed",
        "next_gate": (
            "Select an independently produced dataset whose preparation, "
            "instrument, calibration, environment, and raw records satisfy the "
            "frozen operational requirements."
        ),
    }


def main() -> None:
    results = {
        "release": "10.8",
        "title": "Programs 71--80: Symbol, Quotient, and Operational Bridge",
        "date": "2026-07-27",
        "author": "Krzysztof Żuchowski",
        "orcid": "0009-0002-0909-3613",
        "ontology_guardrail": (
            "The nadsoliton is treated as primordial fractal information in a "
            "solitonic state; no lower informational substrate is introduced."
        ),
        "program71": program71_fourier_schur_theorem(),
        "program72": program72_large_n_projective_scaling(),
        "program73": program73_zero_mode_quotient_information(),
        "program74": program74_locality_truncation(),
        "program75": program75_noisy_scale_calibration(),
        "program76": program76_state_dependent_chirality(),
        "program77": program77_finite_time_landauer_jarzynski(),
        "program78": program78_noisy_process_tomography(),
        "program79": program79_legacy_to_strict_schur_flow(),
        "program80": program80_external_preregistration(),
        "global_guardrail": (
            "No result in this suite closes QW-2191, derives a dimensional "
            "standard, establishes Lorentz invariance, transfers legacy physical "
            "roles, proves the legacy-to-strict completion bridge, or converts "
            "FIN into an experimentally validated physical theory."
        ),
    }
    OUT.write_text(
        json.dumps(results, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(results, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
