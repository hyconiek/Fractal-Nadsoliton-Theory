#!/usr/bin/env python3
"""Execute FIN Programs 51--60: Green geometry, fractal compression, physics bridge.

All core calculations are finite and dimensionless.  Any mass shift, unit
triple, orientation, state, or instrument is explicitly declared as a
regularizer or conditioning datum.  The script exports no strict selector,
physical unit, legacy-role transfer, completed legacy-to-strict bridge, or ToE.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_51_60_Green_Fractal_Physics_Results.json"
FIG = ROOT / "FIN_Programs_51_60_Green_Fractal_Physics_Figures"
FIG.mkdir(exist_ok=True)

ALPHA_GEO = 4.0 * math.log(2.0)
BETA_TORS = 0.01
OMEGA_L = math.pi / 4.0
PHI_L = math.pi / 6.0
OMEGA_S = 0.18575
PHI_S = 0.16250
BETA_S = 1.0
ETA_S = 1.8
S_REPORTED = 1.660307278766099
MASS2 = 0.25  # dimensionless infrared regularizer, not a physical mass


def k_legacy(d):
    d = np.asarray(d, dtype=float)
    return ALPHA_GEO * np.cos(OMEGA_L * d + PHI_L) / (1.0 + BETA_TORS * d)


def k_strict(d):
    d = np.asarray(d, dtype=float)
    return np.cos(OMEGA_S * d + PHI_S) / (1.0 + BETA_S * d**ETA_S)


def cycle_distance(i: int, j: int, n: int) -> int:
    return min((i - j) % n, (j - i) % n)


def radial_matrix(kernel, n: int, diagonal: bool = False) -> np.ndarray:
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j and not diagonal:
                continue
            w[i, j] = float(kernel(cycle_distance(i, j, n)))
    return w


def directed_matrix(kernel, n: int) -> np.ndarray:
    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i != j:
                w[i, j] = float(kernel((j - i) % n))
    return w


def laplacian(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(axis=1)) - w


def reflection(n: int) -> np.ndarray:
    r = np.zeros((n, n))
    for i in range(n):
        r[i, (-i) % n] = 1.0
    return r


def make_spd(a: np.ndarray, floor: float = MASS2) -> tuple[np.ndarray, float]:
    mineig = float(np.linalg.eigvalsh(a).min())
    shift = max(0.0, floor - mineig)
    return a + shift * np.eye(a.shape[0]), shift


def schur_complement(a: np.ndarray, keep: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    drop = np.array([i for i in range(n) if i not in set(map(int, keep))])
    akk = a[np.ix_(keep, keep)]
    akd = a[np.ix_(keep, drop)]
    add = a[np.ix_(drop, drop)]
    return akk - akd @ np.linalg.solve(add, akd.conj().T)


def metric_audit(d: np.ndarray, tol: float = 1e-10) -> dict:
    n = len(d)
    violations = []
    for i in range(n):
        for j in range(n):
            for k in range(n):
                excess = float(d[i, k] - d[i, j] - d[j, k])
                if excess > tol:
                    violations.append(excess)
    diagonal = float(np.max(np.abs(np.diag(d))))
    symmetry = float(np.max(np.abs(d - d.T)))
    nonnegative = float(d.min())
    separation = min(
        float(d[i, j]) for i in range(n) for j in range(n) if i != j
    )
    return {
        "diagonal_max_abs": diagonal,
        "symmetry_max_abs": symmetry,
        "minimum_entry": nonnegative,
        "minimum_offdiagonal": separation,
        "triangle_violations": len(violations),
        "maximum_triangle_excess": max(violations) if violations else 0.0,
        "is_metric": bool(
            diagonal <= tol
            and symmetry <= tol
            and nonnegative >= -tol
            and separation > tol
            and not violations
        ),
    }


def green_squared_distance(a: np.ndarray) -> np.ndarray:
    g = np.linalg.inv(a)
    diag = np.diag(g).real
    return np.maximum(diag[:, None] + diag[None, :] - 2.0 * g.real, 0.0)


def effective_resistance(l: np.ndarray) -> np.ndarray:
    g = np.linalg.pinv(l, hermitian=True)
    diag = np.diag(g)
    return np.maximum(diag[:, None] + diag[None, :] - 2.0 * g, 0.0)


def program51_resolvent_domain() -> dict:
    n = 12
    ws = radial_matrix(k_strict, n)
    wl = radial_matrix(k_legacy, n)
    ss = float(ws.sum(axis=1)[0])
    sl = float(wl.sum(axis=1)[0])
    fixed_s_strict = S_REPORTED * np.eye(n) - ws
    fixed_s_legacy = S_REPORTED * np.eye(n) - wl
    ls = laplacian(ws)
    ll = laplacian(wl)
    als, shift_s = make_spd(ls)
    all_, shift_l = make_spd(ll)
    return {
        "strict_row_sum": ss,
        "reported_s": S_REPORTED,
        "reported_s_minus_strict_row_sum": S_REPORTED - ss,
        "legacy_row_sum": sl,
        "fixed_s_strict_min_eigenvalue": float(np.linalg.eigvalsh(fixed_s_strict).min()),
        "fixed_s_legacy_min_eigenvalue": float(np.linalg.eigvalsh(fixed_s_legacy).min()),
        "strict_laplacian_min_eigenvalue": float(np.linalg.eigvalsh(ls).min()),
        "legacy_signed_laplacian_min_eigenvalue": float(np.linalg.eigvalsh(ll).min()),
        "strict_spd_regularizer": shift_s,
        "legacy_spd_regularizer": shift_l,
        "static_green_inverse_residual_strict": float(
            np.linalg.norm(als @ np.linalg.inv(als) - np.eye(n), "fro")
        ),
        "static_green_inverse_residual_legacy": float(
            np.linalg.norm(all_ @ np.linalg.inv(all_) - np.eye(n), "fro")
        ),
        "theorem": (
            "The resolvent (zI-A)^-1 exists exactly for z outside the spectrum. "
            "Positive definiteness is sufficient for G(0) only when zero is "
            "excluded; a graph Laplacian itself has a zero mode."
        ),
        "correction_to_source_text": (
            "The beta Gram threshold does not imply positivity of an unrelated "
            "fixed-s generator, and the reported strict row sum must not be "
            "silently reused for canonical legacy."
        ),
    }


def program52_green_geometry() -> dict:
    n = 12
    ws = radial_matrix(k_strict, n)
    ls = laplacian(ws)
    a = ls + MASS2 * np.eye(n)
    d2 = green_squared_distance(a)
    droot = np.sqrt(d2)

    # Explicit SPD counterexample: squared Hilbert distances need not be metrics.
    eps = 0.1
    b = np.array([[0.0, 1.0, 2.0], [eps, 0.0, 0.0], [0.0, eps, 0.0]])
    g_counter = b.T @ b
    a_counter = np.linalg.inv(g_counter)
    d2_counter = green_squared_distance(a_counter)

    wlp = np.maximum(radial_matrix(k_legacy, n), 0.0)
    wla = np.abs(radial_matrix(k_legacy, n))
    audits = {
        "strict_green_R_squared": metric_audit(d2),
        "strict_green_sqrt_R": metric_audit(droot),
        "arbitrary_SPD_R_counterexample": metric_audit(d2_counter),
        "strict_effective_resistance": metric_audit(effective_resistance(ls)),
        "legacy_positive_part_effective_resistance": metric_audit(
            effective_resistance(laplacian(wlp))
        ),
        "legacy_absolute_effective_resistance": metric_audit(
            effective_resistance(laplacian(wla))
        ),
    }

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1), constrained_layout=True)
    im0 = axes[0].imshow(d2, cmap="magma")
    axes[0].set_title(r"$R_{ij}$: squared Green embedding")
    axes[0].set_xlabel("state j")
    axes[0].set_ylabel("state i")
    fig.colorbar(im0, ax=axes[0], fraction=0.046)
    im1 = axes[1].imshow(droot, cmap="viridis")
    axes[1].set_title(r"$\sqrt{R_{ij}}$: guaranteed Hilbert metric")
    axes[1].set_xlabel("state j")
    axes[1].set_ylabel("state i")
    fig.colorbar(im1, ax=axes[1], fraction=0.046)
    fig.savefig(FIG / "program52_green_geometry.png", dpi=190)
    plt.close(fig)
    return {
        "audits": audits,
        "counterexample_R_matrix": d2_counter.tolist(),
        "corrected_theorem": (
            "For arbitrary positive-definite A, sqrt(Gii+Gjj-2Re Gij) is a "
            "metric. The unsquared expression is not generally a metric. "
            "Effective resistance itself is a metric only with the additional "
            "connected graph-Laplacian structure."
        ),
    }


def contour_function(a: np.ndarray, function, points: int = 1024) -> np.ndarray:
    eig = np.linalg.eigvalsh(a)
    lo, hi = float(eig.min()), float(eig.max())
    center = 0.5 * (lo + hi)
    radius = 0.5 * (hi - lo) + 0.1 * max(lo, 1e-3)
    out = np.zeros_like(a, dtype=complex)
    eye = np.eye(a.shape[0])
    for theta in np.linspace(0.0, 2.0 * np.pi, points, endpoint=False):
        phase = np.exp(1j * theta)
        z = center + radius * phase
        out += function(z) * np.linalg.inv(z * eye - a) * radius * phase / points
    return out


def program53_resolvent_dual_calculus() -> dict:
    n = 12
    a = laplacian(radial_matrix(k_strict, n)) + MASS2 * np.eye(n)
    t = 0.7
    u_direct = expm(-1j * t * a)
    p_direct = expm(-t * a)
    g_direct = np.linalg.inv(a)
    u_green = contour_function(a, lambda z: np.exp(-1j * t * z))
    p_green = contour_function(a, lambda z: np.exp(-t * z))
    g_green = contour_function(a, lambda z: 1.0 / z)
    return {
        "contour_points": 1024,
        "unitary_relative_error": float(
            np.linalg.norm(u_green - u_direct, "fro") / np.linalg.norm(u_direct, "fro")
        ),
        "diffusive_relative_error": float(
            np.linalg.norm(p_green - p_direct, "fro") / np.linalg.norm(p_direct, "fro")
        ),
        "static_green_relative_error": float(
            np.linalg.norm(g_green - g_direct, "fro") / np.linalg.norm(g_direct, "fro")
        ),
        "unitarity_residual": float(
            np.linalg.norm(u_direct.conj().T @ u_direct - np.eye(n), "fro")
        ),
        "verdict": (
            "The same resolvent reconstructs both semigroups by holomorphic "
            "functional calculus. Their operational outcome difference is not "
            "merely a difference of contours: it also depends on state, "
            "preparation, instrument, and record map."
        ),
    }


def circulant_eigenvalues(kernel, n: int) -> np.ndarray:
    first = radial_matrix(kernel, n)[0]
    return np.fft.fft(first).real


def functional_fit(source: np.ndarray, target: np.ndarray, tol: float = 1e-9):
    order = np.argsort(source)
    groups = []
    for idx in order:
        if not groups or abs(source[idx] - source[groups[-1][0]]) > tol:
            groups.append([int(idx)])
        else:
            groups[-1].append(int(idx))
    x = np.array([np.mean(source[g]) for g in groups])
    y = np.array([np.mean(target[g]) for g in groups])
    spread = max(float(np.ptp(target[g])) for g in groups)
    poly = np.polynomial.Polynomial.fit(x, y, len(x) - 1).convert()
    pred = poly(source)
    return groups, poly, spread, pred


def program54_finite_resolvent_bridge() -> dict:
    mu_l12 = circulant_eigenvalues(k_legacy, 12)
    mu_s12 = circulant_eigenvalues(k_strict, 12)
    groups, poly, spread, pred12 = functional_fit(mu_l12, mu_s12)
    mu_l16 = circulant_eigenvalues(k_legacy, 16)
    mu_s16 = circulant_eigenvalues(k_strict, 16)
    pred16 = poly(mu_l16)
    return {
        "C12_distinct_legacy_eigenvalues": len(groups),
        "C12_degeneracy_obstruction_max_strict_spread": spread,
        "C12_polynomial_degree": int(poly.degree()),
        "C12_relative_functional_residual": float(
            np.linalg.norm(pred12 - mu_s12) / np.linalg.norm(mu_s12)
        ),
        "C16_transfer_relative_residual_using_C12_map": float(
            np.linalg.norm(pred16 - mu_s16) / np.linalg.norm(mu_s16)
        ),
        "C12_polynomial_coefficients_in_power_basis": poly.coef.tolist(),
        "verdict": (
            "A support-specific functional calculus map exists on C12 because "
            "the two radial circulants share Fourier eigenspaces and legacy "
            "degeneracies do not split strict eigenvalues. Transfer to C16 "
            "tests whether this is more than finite interpolation; failure is "
            "not a continuum bridge theorem."
        ),
    }


def radial_offdiagonal_profile(a: np.ndarray) -> np.ndarray:
    n = a.shape[0]
    out = []
    for d in range(1, n // 2 + 1):
        vals = [-a[i, j].real for i in range(n) for j in range(n) if cycle_distance(i, j, n) == d]
        out.append(float(np.mean(vals)))
    return np.array(out)


def amplitude_fit_error(profile: np.ndarray, target: np.ndarray) -> tuple[float, float]:
    c = float(np.dot(profile, target) / max(np.dot(profile, profile), 1e-30))
    err = float(np.linalg.norm(c * profile - target) / np.linalg.norm(target))
    return c, err


def strict_tail_fit(profile: np.ndarray) -> dict:
    d = np.arange(1.0, len(profile) + 1.0)

    def residual(x):
        c, beta, eta = np.exp(x[0]), np.exp(x[1]), np.exp(x[2])
        pred = c * np.cos(OMEGA_S * d + PHI_S) / (1.0 + beta * d**eta)
        return pred - profile

    fit = least_squares(
        residual,
        np.log([max(abs(profile[0]), 1e-3), 1.0, ETA_S]),
        max_nfev=20000,
    )
    c, beta, eta = np.exp(fit.x)
    rel = float(np.linalg.norm(fit.fun) / max(np.linalg.norm(profile), 1e-30))
    return {"c": float(c), "beta": float(beta), "eta": float(eta), "relative_error": rel}


def compressed_legacy(n: int) -> tuple[np.ndarray, np.ndarray, float, float]:
    l = laplacian(radial_matrix(k_legacy, n))
    a, shift = make_spd(l)
    keep = np.arange(0, n, 2)
    s = schur_complement(a, keep)
    identity = float(
        np.linalg.norm(
            np.linalg.inv(s) - np.linalg.inv(a)[np.ix_(keep, keep)], "fro"
        )
    )
    return a, s, shift, identity


def program55_fractal_schur_compression() -> dict:
    target = k_strict(np.arange(1.0, 7.0))
    raw = k_legacy(np.arange(1.0, 7.0))
    _, s12, shift24, identity24 = compressed_legacy(24)
    eff12 = radial_offdiagonal_profile(s12)
    c_raw, e_raw = amplitude_fit_error(raw, target)
    c_eff, e_eff = amplitude_fit_error(eff12, target)

    # Two exact binary eliminations, 48 -> 24 -> 12.
    l48 = laplacian(radial_matrix(k_legacy, 48))
    a48, shift48 = make_spd(l48)
    keep48 = np.arange(0, 48, 2)
    s24 = schur_complement(a48, keep48)
    keep24 = np.arange(0, 24, 2)
    s12_twostep = schur_complement(s24, keep24)
    eff12_twostep = radial_offdiagonal_profile(s12_twostep)
    c_eff2, e_eff2 = amplitude_fit_error(eff12_twostep, target)

    fit1 = strict_tail_fit(eff12)
    fit2 = strict_tail_fit(eff12_twostep)
    d = np.arange(1, 7)
    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    ax.plot(d, target / np.linalg.norm(target), "o-", label="strict target")
    ax.plot(d, raw / np.linalg.norm(raw), "s--", label="canonical legacy")
    ax.plot(d, eff12 / np.linalg.norm(eff12), "^-", label="Green-Schur 24→12")
    ax.plot(
        d,
        eff12_twostep / np.linalg.norm(eff12_twostep),
        "d-",
        label="Green-Schur 48→24→12",
    )
    ax.axhline(0.0, color="black", linewidth=0.6)
    ax.set_xlabel("cycle distance d")
    ax.set_ylabel("unit-L2 radial coupling")
    ax.set_title("Fractal binary compression preserves retained Green data, not strict shape")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "program55_fractal_compression.png", dpi=190)
    plt.close(fig)
    return {
        "compression_rule": "exact Schur complement after retaining even sites",
        "scale_ratios_tested": [2, 4],
        "dimensionless_regularizer_24": shift24,
        "dimensionless_regularizer_48": shift48,
        "retained_green_identity_residual_24_to_12": identity24,
        "raw_legacy_to_strict_amplitude_fit": {"amplitude": c_raw, "relative_error": e_raw},
        "one_step_to_strict_amplitude_fit": {"amplitude": c_eff, "relative_error": e_eff},
        "two_step_to_strict_amplitude_fit": {"amplitude": c_eff2, "relative_error": e_eff2},
        "one_step_strict_phase_tail_fit": fit1,
        "two_step_strict_phase_tail_fit": fit2,
        "verdict": (
            "Schur compression is an exact information-preserving marginalization "
            "for retained Green responses. Whether it approaches strict is an "
            "empirical finite test; a fitted exponent is not a source theorem. "
            "The binary ratio supplies relative scale only, never an absolute length."
        ),
    }


def spectral_dimension_curve(l: np.ndarray, times: np.ndarray):
    eig = np.linalg.eigvalsh(l)
    positive = eig[eig > 1e-10]
    z = np.array([np.exp(-t * positive).sum() for t in times])
    ds = -2.0 * np.gradient(np.log(np.maximum(z, 1e-300)), np.log(times))
    reliable = z > 1e-8
    return eig, z, ds, reliable


def plateau_width(times: np.ndarray, ds: np.ndarray, center=1.0, tol=0.15) -> float:
    mask = np.abs(ds - center) <= tol
    best = 0.0
    start = None
    for i, ok in enumerate(mask):
        if ok and start is None:
            start = i
        if start is not None and (not ok or i == len(mask) - 1):
            end = i if ok and i == len(mask) - 1 else i - 1
            best = max(best, math.log10(times[end] / times[start]) if end > start else 0.0)
            start = None
    return float(best)


def program56_spectral_dimension() -> dict:
    times = np.logspace(-3, 3, 700)
    operators = {
        "strict_C12_raw": (12, laplacian(radial_matrix(k_strict, 12))),
        "strict_C48_positive_part": (
            48,
            laplacian(np.maximum(radial_matrix(k_strict, 48), 0.0)),
        ),
        "strict_C48_absolute": (
            48,
            laplacian(np.abs(radial_matrix(k_strict, 48))),
        ),
        "legacy_C48_positive_part": (
            48,
            laplacian(np.maximum(radial_matrix(k_legacy, 48), 0.0)),
        ),
        "legacy_C48_absolute": (
            48,
            laplacian(np.abs(radial_matrix(k_legacy, 48))),
        ),
    }
    results = {}
    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    for name, (n, l) in operators.items():
        eig, z, ds, reliable = spectral_dimension_curve(l, times)
        positive = eig[eig > 1e-10]
        ds_plot = np.where(reliable, ds, np.nan)
        reliable_indices = np.flatnonzero(reliable)
        closest = reliable_indices[np.argmin(abs(ds[reliable] - 1.0))]
        results[name] = {
            "cycle_size": n,
            "minimum_eigenvalue": float(eig.min()),
            "spectral_gap": float(positive.min()),
            "reliability_cut": "nonzero-mode heat trace > 1e-8",
            "time_at_dimension_closest_to_one": float(times[closest]),
            "plateau_width_decades_within_1_plusminus_0p15": plateau_width(
                times, ds_plot
            ),
        }
        ax.semilogx(times, ds_plot, label=name.replace("_", " "))
    ax.axhline(1.0, color="black", linestyle="--", linewidth=0.8, label="d_s=1")
    ax.set_ylim(0, 3)
    ax.set_xlabel("dimensionless diffusion time")
    ax.set_ylabel(r"local spectral dimension $d_s(t)$")
    ax.set_title("Finite spectral dimension is window- and repair-dependent")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "program56_spectral_dimension.png", dpi=190)
    plt.close(fig)
    return {
        "results": results,
        "verdict": (
            "Finite heat traces define a scale-dependent spectral-dimension "
            "diagnostic, but no robust plateau by itself exports spacetime "
            "dimension. Beyond the positive C12 strict support, both strict and "
            "legacy extrapolations depend on a declared positivity repair."
        ),
    }


def gaussian_kl(cov_p: np.ndarray, cov_q: np.ndarray) -> float:
    n = cov_p.shape[0]
    signp, logdetp = np.linalg.slogdet(cov_p)
    signq, logdetq = np.linalg.slogdet(cov_q)
    assert signp > 0 and signq > 0
    return float(
        0.5
        * (
            np.trace(np.linalg.solve(cov_q, cov_p)).real
            - n
            + logdetq
            - logdetp
        )
    )


def program57_information_action_compression() -> dict:
    n = 24
    al, shift_l = make_spd(laplacian(radial_matrix(k_legacy, n)))
    ass, shift_s = make_spd(laplacian(radial_matrix(k_strict, n)))
    gl = np.linalg.inv(al)
    gs = np.linalg.inv(ass)
    keep = np.arange(0, n, 2)
    glk = gl[np.ix_(keep, keep)]
    gsk = gs[np.ix_(keep, keep)]
    kl_full = gaussian_kl(gl, gs)
    kl_coarse = gaussian_kl(glk, gsk)
    sl = schur_complement(al, keep)
    ss = schur_complement(ass, keep)
    return {
        "legacy_regularizer": shift_l,
        "strict_regularizer": shift_s,
        "KL_legacy_to_strict_full": kl_full,
        "KL_legacy_to_strict_retained_even_states": kl_coarse,
        "data_processing_margin": kl_full - kl_coarse,
        "legacy_marginal_precision_identity_residual": float(
            np.linalg.norm(np.linalg.inv(glk) - sl, "fro")
        ),
        "strict_marginal_precision_identity_residual": float(
            np.linalg.norm(np.linalg.inv(gsk) - ss, "fro")
        ),
        "verdict": (
            "Green-Schur compression is a Gaussian information marginal. It "
            "reduces model distinguishability as required by data processing. "
            "The resulting log-determinant/KL actions are dimensionless and "
            "coordinate/regularizer conditioned, not physical action or entropy."
        ),
    }


def program58_chiral_green_compression() -> dict:
    n = 24
    r = reflection(n)
    base, shift = make_spd(laplacian(radial_matrix(k_legacy, n)), floor=0.4)
    directed = directed_matrix(k_legacy, n)
    c = 0.5j * (directed - directed.T)
    c = c / np.linalg.norm(c, "fro")
    lam = 0.2
    ap = base + lam * c
    am = base - lam * c
    keep = np.arange(0, n, 2)
    sp = schur_complement(ap, keep)
    sm = schur_complement(am, keep)
    rr = reflection(n // 2)
    derivative = (sp - sm) / (2.0 * lam)
    return {
        "legacy_regularizer": shift,
        "full_odd_carrier_residual": float(np.linalg.norm(r @ c @ r.T + c, "fro")),
        "full_branch_covariance_residual": float(np.linalg.norm(r @ ap @ r.T - am, "fro")),
        "compressed_branch_covariance_residual": float(
            np.linalg.norm(rr @ sp @ rr.T - sm, "fro")
        ),
        "compressed_linear_response_odd_residual": float(
            np.linalg.norm(rr @ derivative @ rr.T + derivative, "fro")
        ),
        "compressed_chiral_response_fro_norm": float(np.linalg.norm(derivative, "fro")),
        "plus_minus_isospectral_residual_after_compression": float(
            np.linalg.norm(np.linalg.eigvalsh(sp) - np.linalg.eigvalsh(sm))
        ),
        "verdict": (
            "Exact Green-Schur compression transports the inversion-odd receiver "
            "and mirror-paired branches. It does not generate the directed chart "
            "or choose the sign of lambda; QW-2191 remains open."
        ),
    }


def js_divergence(p: np.ndarray, q: np.ndarray) -> float:
    p = np.maximum(np.asarray(p, float), 1e-300)
    q = np.maximum(np.asarray(q, float), 1e-300)
    p /= p.sum()
    q /= q.sum()
    m = 0.5 * (p + q)
    return float(0.5 * np.sum(p * np.log(p / m)) + 0.5 * np.sum(q * np.log(q / m)))


def program59_operational_records() -> dict:
    n = 12
    l = laplacian(radial_matrix(k_strict, n))
    a, b = 2, 10
    psi = np.zeros(n, complex)
    psi[a] = 1.0 / math.sqrt(2.0)
    psi[b] = 1.0 / math.sqrt(2.0)
    p0 = np.abs(psi) ** 2
    t = 1.0
    u = expm(-1j * t * l)
    heat = expm(-t * l)
    p_wave = np.abs(u @ psi) ** 2
    p_diff = heat @ p0
    p_incoherent = 0.5 * (np.abs(u[:, a]) ** 2 + np.abs(u[:, b]) ** 2)

    coarse = np.zeros((n // 2, n))
    for j in range(n):
        coarse[j // 2, j] = 1.0
    wave_c = coarse @ p_wave
    diff_c = coarse @ p_diff
    js_fine = js_divergence(p_wave, p_diff)
    js_coarse = js_divergence(wave_c, diff_c)

    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    x = np.arange(n)
    ax.plot(x, p_wave, "o-", label="unitary two-slit Born record")
    ax.plot(x, p_diff, "s--", label="diffusive two-source record")
    ax.plot(x, p_incoherent, "^-.", label="unitary incoherent control")
    ax.set_xlabel("detector site")
    ax.set_ylabel("normalized record probability")
    ax.set_title("One strict generator, different operational records")
    ax.legend(fontsize=8)
    fig.savefig(FIG / "program59_operational_records.png", dpi=190)
    plt.close(fig)
    return {
        "dimensionless_time": t,
        "wave_probability_sum": float(p_wave.sum()),
        "diffusion_probability_sum": float(p_diff.sum()),
        "heat_minimum_entry": float(heat.min()),
        "coherent_vs_incoherent_L1": float(np.sum(np.abs(p_wave - p_incoherent))),
        "wave_vs_diffusion_JSD_fine": js_fine,
        "wave_vs_diffusion_JSD_after_pair_coarse_graining": js_coarse,
        "data_processing_margin": js_fine - js_coarse,
        "verdict": (
            "Preparation plus dynamical calculus plus instrument coarse-graining "
            "produces distinct records without invoking observer consciousness. "
            "The kernel does not itself choose the preparation, clock, Born "
            "instrument, environment, apparatus, or physical calibration."
        ),
    }


def program60_physical_scale_obstruction() -> dict:
    n = 12
    eig = np.linalg.eigvalsh(laplacian(radial_matrix(k_strict, n)))
    gap = float(eig[1])
    # Purely illustrative conversion triples; none is sourced by W0.
    triples = [
        {"name": "gauge_A", "ell_star": 1.0, "tau_star": 1.0, "hbar_star": 1.0},
        {"name": "gauge_B", "ell_star": 1e-6, "tau_star": 1e2, "hbar_star": 1e-3},
        {"name": "gauge_C", "ell_star": 1e9, "tau_star": 1e-4, "hbar_star": 1e5},
    ]
    rows = []
    for row in triples:
        rows.append(
            {
                **row,
                "derived_speed_unit": row["ell_star"] / row["tau_star"],
                "derived_energy_gap": row["hbar_star"] * gap / row["tau_star"],
                "physical_length_of_one_binary_compression_cell": 2.0 * row["ell_star"],
            }
        )
    dimension_matrix = np.eye(3)
    return {
        "dimensionless_strict_gap": gap,
        "conversion_examples": rows,
        "minimal_independent_conversion_rank": int(np.linalg.matrix_rank(dimension_matrix)),
        "fractal_compression_exports": "relative scale ratio 2^k and dimensionless flow",
        "fractal_compression_does_not_export": "ell_star, tau_star, hbar_star, SI calibration",
        "minimal_conditioned_bridge": {
            "W0": "informational nadsoliton, canonical kernels, Green/Schur relations",
            "CA": "ell_star, tau_star, hbar_star or an equivalent calibrated unit basis",
            "OA": "state, preparation, clock, dynamics, instrument, environment, apparatus, record",
            "SA_when_orientation_is_measured": "origin/polarity representative and signed coupling law",
            "W3": "experimentally testable predictions conditioned on CA+OA(+SA)",
        },
        "theorem": (
            "Dimensionless kernel, Green, entropy, and fractal scale-ratio data "
            "are invariant under arbitrary changes of the conversion triple. "
            "Consequently they cannot determine dimensional experimental values "
            "without at least one calibrated scale source and an operational map."
        ),
    }


def main() -> None:
    programs = {
        "program_51_resolvent_domain": program51_resolvent_domain(),
        "program_52_green_geometry": program52_green_geometry(),
        "program_53_resolvent_dual_calculus": program53_resolvent_dual_calculus(),
        "program_54_finite_resolvent_bridge": program54_finite_resolvent_bridge(),
        "program_55_fractal_schur_compression": program55_fractal_schur_compression(),
        "program_56_spectral_dimension": program56_spectral_dimension(),
        "program_57_information_action_compression": program57_information_action_compression(),
        "program_58_chiral_green_compression": program58_chiral_green_compression(),
        "program_59_operational_records": program59_operational_records(),
        "program_60_physical_scale_obstruction": program60_physical_scale_obstruction(),
    }
    payload = {
        "release": "10.6",
        "programs": "51-60",
        "seed": "deterministic/no stochastic fitting",
        "scope": (
            "finite dimensionless Green/operator/information calculations; "
            "conditioned physics bridge only"
        ),
        "source_files_audited": [
            "Funkcja Greena w kontekście FIN qwen.md",
            "geometria przestrzeni stanow.md",
        ],
        "programs_data": programs,
        "global_verdict": (
            "Green resolvents and exact Schur complements unify finite response, "
            "geometry, marginalization, dual dynamics, and chiral receiver "
            "transport. Fractal compression generates relative scale flow but "
            "neither strict eta=1.8, an absolute unit, nor a selector. The shortest "
            "honest physics bridge is W0 + conversion axioms + operational map, "
            "with a sector axiom only for orientation-sensitive records."
        ),
    }
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps({"output": str(OUT), "figures": sorted(p.name for p in FIG.iterdir())}, indent=2))


if __name__ == "__main__":
    main()
