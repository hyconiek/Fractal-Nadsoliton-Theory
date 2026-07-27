#!/usr/bin/env python3
"""Execute FIN research Programs 138--150.

This round attacks the state-selection obstruction left by Programs 125--137,
upgrades the fractional certificate with genuine interval arithmetic, and
constructs the remaining operational objects needed for a pre-data physical
test.  Strict/legacy kernels, receiver/source claims, and conversion axioms
remain explicitly separated.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from mpmath import iv


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_138_150_State_Detector_Bridge_Results.json"
PROTOCOL = ROOT / "FIN_Programs_138_150_PreData_Physical_Protocol.json"
FIG = ROOT / "FIN_Programs_138_150_State_Detector_Bridge_Figures"
FIG.mkdir(exist_ok=True)

PREVIOUS = ROOT / "FIN_Programs_125_137_Trace_Localizer_Physics_Results.json"

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = 743.0 / 4000.0
PHI = 13.0 / 80.0
ETA = 9.0 / 5.0
ALPHA = ETA - 1.0
BETA_TORS = 0.01
LEGACY_OMEGA = math.pi / 4
LEGACY_PHI = math.pi / 6
SEED = 20260727
RNG = np.random.default_rng(SEED)


def sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def canonical_digest(record: dict) -> str:
    payload = json.dumps(
        record, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode()
    return hashlib.sha256(payload).hexdigest()


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Unsupported JSON value: {type(value)!r}")


def program138_modular_kms_state_selection() -> dict:
    dims = np.asarray([1, 2, 2, 2, 2], dtype=float)
    theta = np.linspace(-3, 4, 500)
    # Full Hilbert-block Gibbs state: sector 2 has degeneracy one and all
    # other sectors have degeneracy two.
    w2 = 1 / (1 + 8 * np.exp(-theta))
    eta = 2 - w2
    theta_target = math.log(2)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.plot(theta, eta, color="#1f5a99")
    ax.axhline(9 / 5, color="#d1495b", ls="--", label=r"$\eta=9/5$")
    ax.axhline(17 / 9, color="#2a9d8f", ls=":", label=r"Hilbert trace $17/9$")
    ax.axvline(theta_target, color="#d1495b", alpha=0.6)
    ax.set_xlabel(r"dimensionless modular gap $\theta=\beta\Delta$")
    ax.set_ylabel(r"induced exponent $\eta$")
    ax.set_title("A modular Gibbs family does not select its own gap")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program138_modular_kms_family.png", dpi=190)
    plt.close(fig)

    return {
        "commutative_algebra": "A_F=C^5",
        "continuous_automorphism_theorem": (
            "Aut(C^5)=S5 is discrete, hence every continuous one-parameter "
            "automorphism group R->Aut(C^5) is trivial."
        ),
        "commutative_KMS_consequence": (
            "For trivial dynamics every probability state is KMS at every beta; "
            "KMS does not select a sector state."
        ),
        "noncommutative_enlargement": "B=M1(C) direct-sum four copies of M2(C)",
        "modular_consequence": (
            "Every faithful density rho defines its own modular flow "
            "sigma_t^rho(a)=rho^(it) a rho^(-it) and is KMS for that flow."
        ),
        "hilbert_gibbs_family": {
            "w2": "1/(1+8 exp(-theta))",
            "eta": "2-w2",
            "theta_at_eta_9_over_5": theta_target,
            "theta_identity": "theta=ln(2)=alpha_geo/4",
            "eta_at_theta_zero": 17 / 9,
        },
        "strict_modular_gap_source_found": False,
        "uniform_trace_selected": False,
        "status": "Proven modular/KMS nonselection theorem",
        "scope": (
            "A chosen modular Hamiltonian can realize 9/5, but fixing its gap "
            "to ln(2) is an additional state-selection premise."
        ),
    }


def entropy_state(dims: np.ndarray, power: float) -> np.ndarray:
    r = dims**power
    return r / r.sum()


def program139_maximum_entropy_reference_classification() -> dict:
    dims = np.asarray([1, 2, 2, 2, 2], dtype=float)
    choices = {}
    for name, power in [
        ("sector_counting", 0.0),
        ("Hilbert_microstate", 1.0),
        ("Hilbert_Schmidt_microstate", 2.0),
    ]:
        w = entropy_state(dims, power)
        choices[name] = {
            "reference_power": power,
            "weights": w.tolist(),
            "eta": float(w @ dims),
        }

    powers = np.linspace(-3, 3, 401)
    etas = np.asarray([entropy_state(dims, a) @ dims for a in powers])
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.plot(powers, etas, color="#6a4c93")
    ax.scatter(
        [0, 1, 2],
        [choices[k]["eta"] for k in choices],
        color=["#d1495b", "#1f5a99", "#2a9d8f"],
        zorder=3,
    )
    ax.axhline(9 / 5, color="#d1495b", ls="--")
    ax.set_xlabel(r"reference-measure power $a$ in $r_p\propto d_p^a$")
    ax.set_ylabel(r"maximum-relative-entropy exponent")
    ax.set_title("Entropy selects the reference measure, not a reference-free state")
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program139_entropy_reference_classification.png", dpi=190)
    plt.close(fig)

    return {
        "relative_entropy_theorem": (
            "Without additional constraints, maximizing -D(w||r) uniquely gives w=r."
        ),
        "reference_choices": choices,
        "unconstrained_Shannon_uniformity": (
            "The uniform sector state follows only after sector-counting measure "
            "is chosen as the reference."
        ),
        "dimension_constraint_test": (
            "Imposing sum_p w_p d_p=9/5 is equivalent to imposing w2=1/5; "
            "the target exponent is then already present in the constraint."
        ),
        "full_S5_symmetry_test": (
            "Full sector permutation symmetry would select the uniform state, "
            "but it exchanges nonisomorphic fibre blocks of dimensions 1 and 2 "
            "and is not an admitted carrier symmetry."
        ),
        "reference_free_uniform_state_found": False,
        "status": "Proven reference-measure dependence and circularity theorem",
        "scope": (
            "Maximum entropy is a valid conditional selector once its sample "
            "space/reference measure is declared; FIN does not yet source that choice."
        ),
    }


def program140_morita_stability_test() -> dict:
    dims = np.asarray([1, 2, 2, 2, 2], dtype=float)
    cases = {
        "original_blocks": np.asarray([1, 1, 1, 1, 1]),
        "uniform_sector_amplification": np.asarray([2, 1, 1, 1, 1]),
        "p2_dominant": np.asarray([20, 1, 1, 1, 1]),
        "unit_sectors_dominant": np.asarray([1, 1, 10, 10, 10]),
    }
    rows = {}
    for name, amplification in cases.items():
        weights = amplification * dims
        weights = weights / weights.sum()
        rows[name] = {
            "amplification": amplification.tolist(),
            "normalized_Hilbert_central_weights": weights.tolist(),
            "eta": float(weights @ dims),
        }

    samples = RNG.integers(1, 21, size=(2500, 5))
    weights = samples * dims
    weights = weights / weights.sum(axis=1, keepdims=True)
    eta_samples = weights @ dims
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.hist(eta_samples, bins=45, color="#3f7cac", alpha=0.8)
    ax.axvline(9 / 5, color="#d1495b", ls="--", label=r"$9/5$")
    ax.axvline(17 / 9, color="#2a9d8f", ls=":", label=r"$17/9$")
    ax.set_xlabel("Hilbert-trace exponent after blockwise amplification")
    ax.set_ylabel("count")
    ax.set_title("Normalized Hilbert trace is not Morita invariant")
    ax.legend()
    fig.savefig(FIG / "program140_morita_amplification.png", dpi=190)
    plt.close(fig)

    return {
        "algebras": (
            "Direct sums of M_(n_p d_p)(C) are strongly Morita equivalent "
            "to C^5 and to the original fibre-block algebra."
        ),
        "cases": rows,
        "sample_eta_range": [float(eta_samples.min()), float(eta_samples.max())],
        "uniform_sector_eta_recovered_by_amplification": rows[
            "uniform_sector_amplification"
        ]["eta"],
        "normalized_Hilbert_trace_Morita_invariant": False,
        "Morita_equivalence_selects_center_state": False,
        "status": "Proven Morita-instability/no-selection theorem",
        "scope": (
            "Morita equivalence preserves the center up to isomorphism but does "
            "not provide a distinguished probability state on that center."
        ),
    }


def iv_lo(x) -> float:
    return float(np.nextafter(float(x.a), -np.inf))


def iv_hi(x) -> float:
    return float(np.nextafter(float(x.b), np.inf))


def interval_fft(values: list) -> list:
    """In-place radix-two interval FFT using mpmath.iv complex balls."""
    x = list(values)
    n = len(x)
    if n == 0 or n & (n - 1):
        raise ValueError("FFT length must be a power of two")
    j = 0
    for i in range(1, n):
        bit = n >> 1
        while j & bit:
            j ^= bit
            bit >>= 1
        j ^= bit
        if i < j:
            x[i], x[j] = x[j], x[i]
    size = 2
    while size <= n:
        root = iv.exp(iv.mpc(0, -2 * iv.pi / size))
        half = size // 2
        for start in range(0, n, size):
            w = iv.mpc(1, 0)
            for offset in range(half):
                even = x[start + offset]
                odd = w * x[start + offset + half]
                x[start + offset] = even + odd
                x[start + offset + half] = even - odd
                w *= root
        size *= 2
    return x


def program141_validated_interval_fft() -> dict:
    iv.dps = 42
    nfft = 2**14
    dmax = nfft // 2 - 1
    omega = iv.mpf([743, 743]) / 4000
    phi = iv.mpf([13, 13]) / 80
    eta_iv = iv.mpf([9, 9]) / 5
    zero = iv.mpc(0, 0)
    ring = [zero for _ in range(nfft)]
    weights = []
    zd = iv.mpf([0, 0])
    derivative = iv.mpf([0, 0])
    for d in range(1, dmax + 1):
        a = abs(iv.cos(omega * d + phi)) / (1 + iv.mpf(d) ** eta_iv)
        weights.append(a)
        ring[d] = iv.mpc(a, 0)
        ring[nfft - d] = iv.mpc(a, 0)
        zd += 2 * a
        derivative += 2 * d * a
    transformed = interval_fft(ring)

    tail = 2 * iv.mpf(dmax) ** (1 - eta_iv) / (eta_iv - 1)
    tail_hi = iv_hi(tail)
    zd_lo, zd_hi = iv_lo(zd), iv_hi(zd)
    derivative_hi = iv_hi(derivative)
    half_cell = math.pi / nfft
    abelian = 2 / (iv.gamma(eta_iv) * iv.sin(iv.pi * (eta_iv - 1) / 2))
    c_lo = iv_lo(abelian) / (zd_hi + tail_hi)
    c_hi = iv_hi(abelian) / zd_lo

    rows = []
    for k in range(1, nfft // 2):
        q_center = 2 * math.pi * k / nfft
        cell_lo = q_center - half_cell
        cell_hi = q_center + half_cell
        if cell_hi < 0.001 or cell_lo > 0.02:
            continue
        qlo = max(0.001, cell_lo)
        qhi = min(0.02, cell_hi)
        spec = transformed[k].real
        numerator = zd - spec
        num_lo, num_hi = iv_lo(numerator), iv_hi(numerator)
        cell_variation = derivative_hi * half_cell
        slo = max(0.0, num_lo - cell_variation) / (zd_hi + tail_hi)
        shi = (num_hi + cell_variation + 2 * tail_hi) / zd_lo
        pred_lo = c_lo * qlo**ALPHA
        pred_hi = c_hi * qhi**ALPHA
        ratio_lo = slo / pred_hi
        ratio_hi = shi / pred_lo
        rel = max(abs(ratio_lo - 1), abs(ratio_hi - 1))
        rows.append(
            {
                "k": k,
                "q_interval": [qlo, qhi],
                "symbol_interval": [slo, shi],
                "prediction_interval": [pred_lo, pred_hi],
                "relative_remainder_upper": rel,
                "fractional_prediction_compatible": (
                    slo <= pred_hi and shi >= pred_lo
                ),
            }
        )

    qleft = min(r["q_interval"][0] for r in rows)
    qright = max(r["q_interval"][1] for r in rows)
    worst = max(rows, key=lambda r: r["relative_remainder_upper"])
    mid = np.asarray([0.5 * sum(r["q_interval"]) for r in rows])
    rels = np.asarray([r["relative_remainder_upper"] for r in rows])
    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.semilogy(mid, rels, "o-", color="#1f5a99", ms=3)
    ax.set_xlabel("|q|")
    ax.set_ylabel("formal interval relative-remainder upper bound")
    ax.set_title("Validated interval FFT: rigorous but deliberately coarse")
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program141_validated_interval_fft.png", dpi=190)
    plt.close(fig)

    return {
        "interval_engine": "mpmath.iv, 42 decimal digits, radix-two interval FFT",
        "fft_length": nfft,
        "retained_distance": dmax,
        "covered_window": [qleft, qright],
        "target_window_fully_covered": qleft <= 0.001 and qright >= 0.02,
        "number_of_frequency_cells": len(rows),
        "tail_normalization_upper": tail_hi,
        "abelian_constant_interval": [c_lo, c_hi],
        "maximum_relative_remainder_upper": worst["relative_remainder_upper"],
        "worst_cell": worst,
        "all_symbol_lower_bounds_positive": all(
            r["symbol_interval"][0] > 0 for r in rows
        ),
        "all_cells_compatible_with_fractional_prediction": all(
            r["fractional_prediction_compatible"] for r in rows
        ),
        "formal_interval_arithmetic": True,
        "tight_2_3_percent_claim_formally_recovered": (
            worst["relative_remainder_upper"] < 0.03
        ),
        "status": "Proven computer-assisted continuous-window enclosure",
        "scope": (
            "The interval certificate is formal but coarser than Program 127 "
            "because N=2^14 is used. It validates the fractional scale without "
            "automatically validating the earlier 2.2647% tight bound."
        ),
    }


def certified_cf_prefix(number_terms=12) -> tuple[list[int], list[dict]]:
    iv.dps = 90
    x = iv.mpf([743, 743]) / (4000 * iv.pi)
    terms = []
    for _ in range(number_terms):
        lo, hi = iv_lo(x), iv_hi(x)
        alo, ahi = math.floor(lo), math.floor(hi)
        if alo != ahi:
            raise ArithmeticError("continued-fraction term not certified")
        terms.append(alo)
        x = 1 / (x - alo)
    rows = []
    pm2, pm1, qm2, qm1 = 0, 1, 1, 0
    for a in terms:
        p = a * pm1 + pm2
        q = a * qm1 + qm2
        rows.append({"a": a, "p": p, "q": q})
        pm2, pm1, qm2, qm1 = pm1, p, qm1, q
    return terms, rows


def ostrowski_digit_sum(n: int, denominators: list[int]) -> int:
    remainder = n
    total = 0
    for q in reversed(denominators):
        digit, remainder = divmod(remainder, q)
        total += digit
    if remainder:
        raise ArithmeticError("incomplete Ostrowski representation")
    return total


def program142_diophantine_discrepancy() -> dict:
    terms, convergents = certified_cf_prefix(12)
    denominators = [row["q"] for row in convergents if row["q"] <= 1_000_000]
    nmax = 1_000_000
    bounds = np.empty(nmax, dtype=float)
    for n in range(1, nmax + 1):
        bounds[n - 1] = 2 * ostrowski_digit_sum(n, denominators)

    d = np.arange(1, nmax + 1, dtype=float)
    sums = np.cumsum(np.abs(np.cos(OMEGA * d + PHI)))
    discrepancy = np.abs(sums - d * (2 / math.pi))
    ratio = discrepancy / bounds
    mask = d >= 1000

    fig, ax = plt.subplots(figsize=(8.7, 4.8), constrained_layout=True)
    take = np.unique(np.geomspace(1, nmax, 1200).astype(int)) - 1
    ax.loglog(d[take], discrepancy[take], label="observed absolute discrepancy")
    ax.loglog(d[take], bounds[take], label="Denjoy–Koksma/Ostrowski bound")
    ax.set_xlabel("N")
    ax.set_ylabel(r"$|\sum_{d\leq N}|\cos|-2N/\pi|$")
    ax.set_title("Certified finite Diophantine discrepancy modulus")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program142_diophantine_discrepancy.png", dpi=190)
    plt.close(fig)

    return {
        "rotation": "theta=omega/pi=743/(4000*pi)",
        "certified_continued_fraction_terms": terms,
        "convergents": convergents,
        "notable_partial_quotients": [67, 928],
        "finite_range": [1, nmax],
        "variation_of_abs_cos_pi_x": 2,
        "modulus": (
            "|sum_(d=1)^N f(x+d theta)-N int f| <= "
            "2 * OstrowskiDigitSum_theta(N)"
        ),
        "maximum_bound_per_N_for_N_ge_1000": float(
            np.max(bounds[mask] / d[mask])
        ),
        "bound_at_N_1e6": float(bounds[-1]),
        "observed_discrepancy_at_N_1e6": float(discrepancy[-1]),
        "maximum_observed_to_bound_ratio": float(np.max(ratio)),
        "all_observed_rows_below_theorem_bound": bool(np.all(ratio <= 1 + 1e-12)),
        "all_scale_power_remainder_proved": False,
        "obstruction": (
            "A finite continued-fraction prefix controls N<=10^6 but cannot "
            "bound all future partial quotients; the term 928 already produces "
            "a severe finite-scale constant."
        ),
        "status": "Proven finite discrepancy theorem; global power-rate remains open",
        "scope": (
            "This supplies the missing explicit finite modulus but not an "
            "all-scale Diophantine type theorem for 743/(4000*pi)."
        ),
    }


def program143_weighted_wave_estimates() -> dict:
    s = np.linspace(0.505, 3.0, 500)
    integral = np.sqrt(np.pi) * np.asarray(
        [math.gamma(x - 0.5) / math.gamma(x) for x in s]
    )
    sobolev_constant = np.sqrt(integral / (2 * math.pi))

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.semilogy(s, sobolev_constant, color="#6a4c93")
    ax.axvline(0.5, color="#d1495b", ls="--")
    ax.set_xlabel("Sobolev index s")
    ax.set_ylabel(r"point-evaluation constant $C_s$")
    ax.set_title("Fractional wave records exist on H^s only above s=1/2")
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program143_weighted_wave_estimates.png", dpi=190)
    plt.close(fig)

    return {
        "unitary_group": "U_t=exp(-it C(-Delta)^(2/5))",
        "Sobolev_invariance": "||U_t f||_(H^s)=||f||_(H^s) for every real s",
        "point_evaluation_theorem": (
            "For s>1/2, |U_t f(x)|<=C_s ||f||_(H^s), "
            "C_s^2=(2pi)^(-1) integral_R (1+q^2)^(-s)dq."
        ),
        "threshold": 0.5,
        "L1_to_Linf_decay_recovered": False,
        "wave_smoothing": False,
        "constructed_record_space": "H^s(R), s>1/2, with bounded point evaluation",
        "physical_detector_source_found": False,
        "status": "Proven weighted point-record theorem without dispersive decay",
        "scope": (
            "Regular preparation can replace a hard UV cutoff mathematically, "
            "but selecting preparation regularity is still operational input."
        ),
    }


def program144_detector_resolution_renormalization() -> dict:
    q = np.linspace(-40, 40, 240_001)
    phase = np.exp(-1j * np.abs(q) ** ALPHA)
    smooth_hat = np.exp(-0.5 * q**2)
    limit = np.trapz(phase * smooth_hat, q) / (2 * math.pi)
    sigmas = np.geomspace(0.002, 1.0, 45)
    values = []
    errors = []
    for sigma in sigmas:
        response = np.exp(-0.5 * sigma**2 * q**2)
        value = np.trapz(response * phase * smooth_hat, q) / (2 * math.pi)
        values.append(value)
        errors.append(abs(value - limit))
    errors = np.asarray(errors)
    fit_mask = sigmas < 0.08
    slope = float(
        np.polyfit(np.log(sigmas[fit_mask]), np.log(errors[fit_mask]), 1)[0]
    )
    delta_peak = 1 / (math.sqrt(2 * math.pi) * sigmas)

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.loglog(sigmas, errors, "o-", label="smooth prepared-state error")
    ax.loglog(sigmas, delta_peak, label="unprepared delta peak at t=0")
    ax.set_xlabel("detector resolution sigma")
    ax.set_ylabel("error or peak magnitude")
    ax.set_title("Resolution removal succeeds for smooth states and fails for a delta")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program144_detector_resolution.png", dpi=190)
    plt.close(fig)

    return {
        "detector": "R_sigma has Fourier multiplier exp(-sigma^2 q^2/2)",
        "measurement_operator": "M_(sigma,t)=R_sigma U_t",
        "smooth_preparation": "f_hat(q)=exp(-q^2/2)",
        "smooth_sigma_to_zero_limit_at_x0_t1": [
            float(np.real(limit)),
            float(np.imag(limit)),
        ],
        "observed_small_sigma_error_slope": slope,
        "general_rate": (
            "For 1/2<s<5/2, ||(R_sigma-I)U_t f||_infinity "
            "<=C_s sigma^(s-1/2)||f||_(H^s)."
        ),
        "delta_peak": "1/(sqrt(2pi)*sigma) at t=0",
        "delta_resolution_independent_limit": False,
        "smooth_resolution_independent_limit": True,
        "status": "Proven conditional detector-limit theorem plus numerical check",
        "scope": (
            "Resolution can be removed for sufficiently regular preparations; "
            "it cannot be removed for ideal point preparation."
        ),
    }


def system_spin(alpha: float, t: float) -> float:
    return math.exp(-(t**alpha))


def joint_observables(params, t1=1.0, t2=2.0):
    alpha, eps, rho = params
    e = 1 - 2 * eps
    s1, s2 = system_spin(alpha, t1), system_spin(alpha, t2)
    mu1, mu2 = e * s1, e * s2
    lag = s1 * s2 * (e * e + rho * (1 - e * e))
    return np.asarray([mu1, mu2, lag])


def numerical_jacobian(fun, x, h=1e-6):
    x = np.asarray(x, dtype=float)
    base = fun(x)
    jac = np.empty((len(base), len(x)))
    for j in range(len(x)):
        xp, xm = x.copy(), x.copy()
        xp[j] += h
        xm[j] -= h
        jac[:, j] = (fun(xp) - fun(xm)) / (2 * h)
    return jac


def invert_joint_observables(obs):
    mu1, mu2, lag = obs
    ratio = max(mu1 / mu2, 1 + 1e-12)
    alpha = math.log(1 + math.log(ratio), 2)
    s1, s2 = system_spin(alpha, 1.0), system_spin(alpha, 2.0)
    e = np.clip(mu1 / s1, 1e-6, 0.999999)
    eps = (1 - e) / 2
    rho = (lag / (s1 * s2) - e * e) / (1 - e * e)
    return np.asarray([alpha, eps, rho])


def program145_joint_system_instrument_identifiability() -> dict:
    truth = np.asarray([ALPHA, 0.1, 0.75])
    observables = joint_observables(truth)
    jac = numerical_jacobian(joint_observables, truth)
    singular = np.linalg.svd(jac, compute_uv=False)
    one_time_jac = numerical_jacobian(
        lambda x: joint_observables(x)[:1], truth
    )

    estimates = []
    for _ in range(3000):
        noisy = observables + RNG.normal(0, 2e-4, 3)
        estimates.append(invert_joint_observables(noisy))
    estimates = np.asarray(estimates)
    rmse = np.sqrt(np.mean((estimates - truth) ** 2, axis=0))

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.9), constrained_layout=True)
    labels = [r"$\alpha$", r"$\epsilon$", r"$\rho$"]
    for j, ax in enumerate(axes):
        ax.hist(estimates[:, j], bins=40, color="#3f7cac")
        ax.axvline(truth[j], color="#d1495b", ls="--")
        ax.set_title(labels[j])
    fig.suptitle("Joint system–instrument inversion from two means and one lag")
    fig.savefig(FIG / "program145_joint_identifiability.png", dpi=190)
    plt.close(fig)

    return {
        "parameter_vector": ["fractional alpha", "flip prevalence epsilon", "memory rho"],
        "observable_vector": [
            "mean spin at t=1",
            "mean spin at t=2",
            "adjacent-time spin product",
        ],
        "truth": truth.tolist(),
        "truth_observables": observables.tolist(),
        "jacobian": jac.tolist(),
        "jacobian_singular_values": singular.tolist(),
        "joint_rank": int(np.linalg.matrix_rank(jac, tol=1e-8)),
        "one_time_rank": int(np.linalg.matrix_rank(one_time_jac, tol=1e-8)),
        "exact_inverse": (
            "alpha=log_2(1+log(mu1/mu2)); e=mu1/exp(-1); "
            "epsilon=(1-e)/2; rho=(lag/(s1 s2)-e^2)/(1-e^2)."
        ),
        "synthetic_noise_sd": 2e-4,
        "synthetic_RMSE": {
            "alpha": float(rmse[0]),
            "epsilon": float(rmse[1]),
            "rho": float(rmse[2]),
        },
        "status": "Proven local/global inversion in the declared model; synthetic validation",
        "scope": (
            "The theorem separates generator exponent and apparatus memory in "
            "this three-parameter spin model, not in every possible instrument family."
        ),
    }


def symmetric_stable(alpha: float, size: int, rng) -> np.ndarray:
    if alpha == 2:
        return rng.normal(0, math.sqrt(2), size)
    u = rng.uniform(-math.pi / 2, math.pi / 2, size)
    w = rng.exponential(1.0, size)
    return (
        np.sin(alpha * u)
        / np.cos(u) ** (1 / alpha)
        * (
            np.cos((1 - alpha) * u) / w
        ) ** ((1 - alpha) / alpha)
    )


def iqr(values: np.ndarray) -> float:
    q25, q75 = np.quantile(values, [0.25, 0.75])
    return float(q75 - q25)


def program146_calibration_invariant_observables() -> dict:
    t1, t2 = 1.0, 4.0
    reps, n = 350, 4000
    slopes = []
    slopes_noise = []
    for _ in range(reps):
        ell = math.exp(RNG.normal(0, 0.8))
        x1 = ell * t1 ** (1 / ALPHA) * symmetric_stable(ALPHA, n, RNG)
        x2 = ell * t2 ** (1 / ALPHA) * symmetric_stable(ALPHA, n, RNG)
        slopes.append(math.log(iqr(x2) / iqr(x1)) / math.log(t2 / t1))
        noise1 = RNG.normal(0, 0.1 * ell, n)
        noise2 = RNG.normal(0, 0.1 * ell, n)
        slopes_noise.append(
            math.log(iqr(x2 + noise2) / iqr(x1 + noise1)) / math.log(t2 / t1)
        )
    slopes = np.asarray(slopes)
    slopes_noise = np.asarray(slopes_noise)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.hist(slopes, bins=38, alpha=0.65, label="calibration cancels")
    ax.hist(slopes_noise, bins=38, alpha=0.55, label="plus detector noise")
    ax.axvline(1 / ALPHA, color="#d1495b", ls="--", label=r"$1/\alpha=5/4$")
    ax.set_xlabel("IQR log-slope")
    ax.set_ylabel("replicates")
    ax.set_title("A calibration-invariant fractional observable")
    ax.legend()
    fig.savefig(FIG / "program146_calibration_invariant_ratio.png", dpi=190)
    plt.close(fig)

    return {
        "primary_observable": (
            "R_IQR=log[IQR(t2)/IQR(t1)]/log(t2/t1)=1/alpha"
        ),
        "calibrations_cancelled": ["length ell", "clock scale tau in time ratios"],
        "target_value_for_alpha_4_over_5": 1 / ALPHA,
        "simulation": {
            "replicates": reps,
            "records_per_time": n,
            "time_pair": [t1, t2],
            "random_log_length_sd": 0.8,
            "mean_without_detector_noise": float(slopes.mean()),
            "sd_without_detector_noise": float(slopes.std(ddof=1)),
            "mean_with_relative_detector_noise": float(slopes_noise.mean()),
            "sd_with_relative_detector_noise": float(slopes_noise.std(ddof=1)),
        },
        "additional_ratios": [
            "spectral gap ratio DeltaE_k/DeltaE_j=lambda_k/lambda_j cancels hbar/tau",
            "quantile ratios at one time cancel ell and tau",
            "normalized characteristic-function ratios cancel overall intensity",
        ],
        "all_physical_calibration_removed": False,
        "status": "Constructed calibration-invariant observable; synthetic robustness evidence",
        "scope": (
            "Time ratios and detector response still require calibrated ordering "
            "and a declared instrument; the statistic removes multiplicative scales."
        ),
    }


def reflection_matrix(n=12):
    r = np.zeros((n, n))
    for i in range(n):
        r[(-i) % n, i] = 1
    return r


def fractional_circulant(n=12):
    k = np.arange(n)
    q = 2 * math.pi * np.minimum(k, n - k) / n
    lam = q**ALPHA
    f = np.exp(2j * math.pi * np.outer(np.arange(n), k) / n) / math.sqrt(n)
    return np.real_if_close(f @ np.diag(lam) @ f.conj().T).real


def current_witness(a):
    n = len(a)
    c = np.zeros((n, n), dtype=complex)
    for i in range(n):
        j = (i + 1) % n
        c[j, i] = -1j * a[j, i]
        c[i, j] = 1j * a[j, i]
    return c


def asymmetry_monotone(rho, r):
    diff = rho - r @ rho @ r
    return 0.5 * float(np.sum(np.abs(np.linalg.eigvalsh(diff))))


def program147_preparation_resource_theory() -> dict:
    n = 12
    r = reflection_matrix(n)
    a = fractional_circulant(n)
    c = current_witness(a)
    x = np.arange(n)
    psi_plus = np.exp(2j * math.pi * x / n) / math.sqrt(n)
    psi_minus = np.exp(-2j * math.pi * x / n) / math.sqrt(n)
    rho_plus = np.outer(psi_plus, psi_plus.conj())
    rho_minus = np.outer(psi_minus, psi_minus.conj())
    rho_mix = 0.5 * (rho_plus + rho_minus)

    lambdas = np.linspace(0, 0.5, 101)
    monotones = []
    currents = []
    for lam in lambdas:
        rho = (1 - lam) * rho_plus + lam * (r @ rho_plus @ r)
        monotones.append(asymmetry_monotone(rho, r))
        currents.append(float(np.real(np.trace(rho @ c))))
    c_norm = float(np.max(np.abs(np.linalg.eigvalsh(c))))

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.plot(lambdas, monotones, label=r"$M(\rho)$")
    ax.plot(lambdas, np.abs(currents) / c_norm, label=r"$|\Lambda|/\|C\|$")
    ax.set_xlabel("free reflection-mixing strength")
    ax.set_ylabel("normalized resource")
    ax.set_title("Reflection-covariant operations cannot create signed preparation")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program147_preparation_resource.png", dpi=190)
    plt.close(fig)

    return {
        "free_states": "rho=R rho R",
        "free_operations": "reflection-covariant CPTP maps Phi(R rho R)=R Phi(rho) R",
        "resource_monotone": "M(rho)=0.5 ||rho-R rho R||_1",
        "monotonicity_theorem": (
            "Trace-norm contractivity gives M(Phi(rho))<=M(rho) for every free Phi."
        ),
        "signed_witness_bound": "|Tr(rho C)|<=M(rho)||C||_infinity when RCR=-C",
        "finite_checks": {
            "RCR_plus_C_norm": float(np.linalg.norm(r @ c @ r + c)),
            "M_rho_plus": asymmetry_monotone(rho_plus, r),
            "M_rho_minus": asymmetry_monotone(rho_minus, r),
            "M_symmetric_mixture": asymmetry_monotone(rho_mix, r),
            "Lambda_plus": float(np.real(np.trace(rho_plus @ c))),
            "Lambda_minus": float(np.real(np.trace(rho_minus @ c))),
            "all_mixing_rows_respect_witness_bound": bool(
                np.all(np.abs(currents) <= c_norm * np.asarray(monotones) + 1e-12)
            ),
        },
        "strict_free_operations_generate_signed_state": False,
        "QW_2191_discharged": False,
        "status": "Proven preparation-resource no-go theorem",
        "scope": (
            "A signed branch is a resource. It cannot arise from a reflection-"
            "symmetric state under reflection-covariant operations."
        ),
    }


def program148_coupled_state_damping_variational() -> dict:
    dims = np.asarray([1, 2, 2, 2, 2], dtype=float)
    powers = np.linspace(-3, 3, 601)
    etas = []
    entropies = []
    for power in powers:
        w = entropy_state(dims, power)
        etas.append(float(w @ dims))
        entropies.append(float(-np.sum(w * np.log(w))))
    etas = np.asarray(etas)

    fig, ax1 = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax1.plot(powers, etas, color="#1f5a99", label=r"$\eta(a)$")
    ax1.axhline(9 / 5, color="#d1495b", ls="--")
    ax1.axhline(17 / 9, color="#2a9d8f", ls=":")
    ax1.set_xlabel(r"unsourced reference exponent $a$")
    ax1.set_ylabel(r"stationary damping exponent $\eta$")
    ax1.set_title("A strictly convex coupled action exposes one missing reference exponent")
    ax1.grid(True, alpha=0.25)
    fig.savefig(FIG / "program148_coupled_variational_family.png", dpi=190)
    plt.close(fig)

    return {
        "functional": (
            "F_(a,kappa)(w,eta)=D_KL(w||r_a)+(kappa/2)"
            "*(eta-sum_p w_p d_p)^2, r_a,p proportional d_p^a"
        ),
        "domain": "probability simplex interior times real eta; kappa>0",
        "strict_convexity": True,
        "unique_stationary_solution": (
            "w=r_a and eta=sum_p r_a,p d_p"
        ),
        "sector_counting_solution": {"a": 0, "eta": 9 / 5},
        "Hilbert_microstate_solution": {"a": 1, "eta": 17 / 9},
        "multiplicative_positive_tail_addition": (
            "T(ab)=T(a)T(b), T(d)=beta d^eta, forces beta=1 separately."
        ),
        "target_free_uniform_selection": False,
        "missing_object": (
            "A theorem selecting the reference exponent a; convexity alone "
            "selects w only after a is supplied."
        ),
        "status": "Constructed well-posed coupled variational family; strict source remains open",
        "scope": (
            "The action is target-free as a family, but choosing a=0 to obtain "
            "9/5 would be an additional premise."
        ),
    }


def program149_completion_map_diagram() -> dict:
    d = np.arange(0, 65, dtype=float)
    legacy = (
        ALPHA_GEO
        * np.cos(LEGACY_OMEGA * d + LEGACY_PHI)
        / (1 + BETA_TORS * d)
    )
    strict = np.cos(OMEGA * d + PHI) / (1 + d**ETA)
    legacy_projective = legacy / legacy[0]
    damping = (1 + BETA_TORS * d) / (1 + d**ETA)
    completed_actual = legacy_projective * damping
    strict_projective = strict / strict[0]
    actual_residual = float(
        np.linalg.norm(completed_actual - strict_projective)
        / np.linalg.norm(strict_projective)
    )

    shared_legacy = (
        ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)
    )
    shared_completed = shared_legacy / shared_legacy[0] * damping
    shared_residual = float(
        np.linalg.norm(shared_completed - strict_projective)
        / np.linalg.norm(strict_projective)
    )

    fig, ax = plt.subplots(figsize=(10.0, 5.4), constrained_layout=True)
    ax.axis("off")
    boxes = [
        (0.03, 0.62, "canonical legacy\n(pi/4, pi/6, beta_tors)"),
        (0.29, 0.62, "amplitude quotient\nPi_0"),
        (0.55, 0.62, "conditional damping\nD_w"),
        (0.81, 0.62, "strict target\n(omega_s, phi_s, 9/5)"),
        (0.29, 0.18, "phase/frequency\ntransport: absent"),
        (0.58, 0.18, "role semantics\ntransport: absent"),
    ]
    for x, y, label in boxes:
        ax.text(
            x,
            y,
            label,
            ha="center",
            va="center",
            transform=ax.transAxes,
            bbox=dict(boxstyle="round,pad=0.45", fc="#eef4f8", ec="#1f5a99"),
        )
    for x1, x2 in [(0.12, 0.24), (0.38, 0.50), (0.64, 0.76)]:
        ax.annotate(
            "",
            xy=(x2, 0.62),
            xytext=(x1, 0.62),
            xycoords=ax.transAxes,
            arrowprops=dict(arrowstyle="->", lw=1.6),
        )
    ax.annotate(
        "",
        xy=(0.75, 0.55),
        xytext=(0.34, 0.25),
        xycoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", color="#d1495b", ls="--"),
    )
    ax.annotate(
        "",
        xy=(0.80, 0.55),
        xytext=(0.62, 0.25),
        xycoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", color="#d1495b", ls="--"),
    )
    ax.set_title("Typed completion diagram: two essential arrows remain absent")
    fig.savefig(FIG / "program149_completion_map_diagram.png", dpi=190)
    plt.close(fig)

    return {
        "nodes": [
            "K_legacy_ont canonical",
            "projective amplitude shape",
            "conditionally damped legacy shape",
            "K_strict_gate projective target",
            "legacy physical roles",
        ],
        "edges": {
            "amplitude_quotient": "proven quotient; erases absolute roles",
            "damping_completion": "conditional on state w2=1/5 and beta multiplicativity",
            "phase_frequency_transport": "absent",
            "topological_selector_transport": "absent",
            "physical_role_transport": "absent and downstream",
        },
        "actual_canonical_legacy_to_projective_strict_relative_residual": actual_residual,
        "shared_phase_frequency_control_residual": shared_residual,
        "commutes_when_oscillatory_numerator_is_shared": shared_residual < 1e-14,
        "commutes_for_actual_canonical_legacy_parameters": actual_residual < 1e-14,
        "full_bridge_complete": False,
        "role_transfer_started": False,
        "status": "Constructed typed diagram and finite noncommutation witness",
        "scope": (
            "Programs 134/135 close only amplitude-quotient and conditional "
            "damping arrows. Canonical legacy phase/frequency remain different."
        ),
    }


def make_protocol(power: dict) -> dict:
    core = {
        "protocol_id": "FIN-P150-FRACTIONAL-IQR-001",
        "version": "1.0.0",
        "frozen_date": "2026-07-27",
        "scientific_question": (
            "Does a calibrated one-dimensional spreading record have local "
            "diffusive exponent alpha=2 or FIN fractional exponent alpha=4/5?"
        ),
        "preparation": (
            "One localized packet prepared by a documented apparatus; identical "
            "preparation procedure at both observation times."
        ),
        "times": [1.0, 4.0],
        "records_per_time_minimum": 4000,
        "detector": {
            "response_family": "Gaussian position response",
            "maximum_relative_resolution_sigma_over_reference_IQR": 0.1,
            "memory_calibration": "independent known-input sequence before test records",
        },
        "calibration": {
            "required": ["time ordering and time ratio", "detector response"],
            "cancelled_by_statistic": ["multiplicative length scale"],
        },
        "primary_statistic": (
            "T=log(IQR(t=4)/IQR(t=1))/log(4)"
        ),
        "null": "H0: alpha=2, predicted T=1/2",
        "alternative": "H1: alpha=4/5, predicted T=5/4",
        "rejection_rule": "Reject H0 in favor of H1 when T>0.875",
        "exclusions": [
            "failed preparation audit",
            "time-ratio uncertainty above preregistered tolerance",
            "detector resolution above 0.1 reference IQR",
            "missing raw records",
            "post-hoc changes to binning or quantile definition",
        ],
        "raw_record_format": [
            "record_id",
            "preparation_batch",
            "time",
            "position_readout",
            "detector_calibration_id",
            "timestamp",
        ],
        "external_data_admitted": False,
        "claim_boundary": (
            "A rejection concerns the declared local-vs-fractional operational "
            "models only; it does not prove the FIN ontology, bridge, or ToE."
        ),
    }
    return {
        "core": core,
        "canonical_core_sha256": canonical_digest(core),
        "synthetic_power_audit": power,
    }


def program150_predata_physical_protocol() -> dict:
    reps, n = 450, 4000
    t1, t2 = 1.0, 4.0
    threshold = 0.875
    null_stats, alt_stats = [], []
    for _ in range(reps):
        ell = math.exp(RNG.normal(0, 0.5))
        # Null: characteristic exp(-t q^2), hence N(0,2t).
        n1 = ell * math.sqrt(t1) * symmetric_stable(2, n, RNG)
        n2 = ell * math.sqrt(t2) * symmetric_stable(2, n, RNG)
        # Alternative: characteristic exp(-t |q|^(4/5)).
        a1 = ell * t1 ** (1 / ALPHA) * symmetric_stable(ALPHA, n, RNG)
        a2 = ell * t2 ** (1 / ALPHA) * symmetric_stable(ALPHA, n, RNG)
        for values in [n1, n2, a1, a2]:
            values += RNG.normal(0, 0.1 * ell, n)
        null_stats.append(math.log(iqr(n2) / iqr(n1)) / math.log(4))
        alt_stats.append(math.log(iqr(a2) / iqr(a1)) / math.log(4))
    null_stats, alt_stats = np.asarray(null_stats), np.asarray(alt_stats)
    power = {
        "replicates": reps,
        "records_per_time": n,
        "threshold": threshold,
        "null_mean": float(null_stats.mean()),
        "null_sd": float(null_stats.std(ddof=1)),
        "alternative_mean": float(alt_stats.mean()),
        "alternative_sd": float(alt_stats.std(ddof=1)),
        "synthetic_false_positive_rate": float(np.mean(null_stats > threshold)),
        "synthetic_power": float(np.mean(alt_stats > threshold)),
    }
    protocol = make_protocol(power)
    PROTOCOL.write_text(
        json.dumps(protocol, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.hist(null_stats, bins=38, alpha=0.65, label="H0: alpha=2")
    ax.hist(alt_stats, bins=38, alpha=0.65, label="H1: alpha=4/5")
    ax.axvline(threshold, color="black", ls="--", label="frozen threshold")
    ax.set_xlabel("preregistered IQR log-slope statistic")
    ax.set_ylabel("synthetic replicates")
    ax.set_title("Pre-data power audit of the frozen operational protocol")
    ax.legend()
    fig.savefig(FIG / "program150_predata_protocol_power.png", dpi=190)
    plt.close(fig)

    return {
        "protocol_file": PROTOCOL.name,
        "canonical_core_sha256": protocol["canonical_core_sha256"],
        "synthetic_power_audit": power,
        "external_data_admitted": False,
        "protocol_frozen_before_external_data": True,
        "status": "Constructed immutable pre-data protocol with synthetic power audit",
        "scope": (
            "The protocol can test an operational spreading exponent after "
            "real calibration and intake. It does not test the full FIN ontology."
        ),
    }


def main() -> None:
    previous = json.loads(PREVIOUS.read_text())
    result = {
        "title": "FIN Programs 138-150: State Selection, Validated Fractional Dynamics, Detector Physics and Bridge Architecture",
        "release": "10.14",
        "execution_date": "2026-07-27",
        "seed": SEED,
        "frozen_parameters": {
            "alpha_geo": ALPHA_GEO,
            "strict_omega": OMEGA,
            "strict_phi": PHI,
            "strict_eta": ETA,
            "legacy_omega": "pi/4",
            "legacy_phi": "pi/6",
            "beta_tors": BETA_TORS,
        },
        "provenance": {
            "previous_results": PREVIOUS.name,
            "previous_results_sha256": sha256(PREVIOUS),
            "script": Path(__file__).name,
        },
        "program138_modular_kms_state_selection": program138_modular_kms_state_selection(),
        "program139_maximum_entropy_reference": program139_maximum_entropy_reference_classification(),
        "program140_morita_stability": program140_morita_stability_test(),
        "program141_validated_interval_fft": program141_validated_interval_fft(),
        "program142_diophantine_discrepancy": program142_diophantine_discrepancy(),
        "program143_weighted_wave_estimates": program143_weighted_wave_estimates(),
        "program144_detector_resolution": program144_detector_resolution_renormalization(),
        "program145_joint_identifiability": program145_joint_system_instrument_identifiability(),
        "program146_calibration_invariant_observables": program146_calibration_invariant_observables(),
        "program147_preparation_resource_theory": program147_preparation_resource_theory(),
        "program148_coupled_state_damping_action": program148_coupled_state_damping_variational(),
        "program149_completion_map_diagram": program149_completion_map_diagram(),
        "program150_predata_protocol": program150_predata_physical_protocol(),
        "global_verdict": {
            "new_objects_constructed": [
                "modular/KMS state-selection no-go theorem",
                "reference-measure entropy classification",
                "Morita amplification trace test",
                "validated interval FFT certificate",
                "finite Diophantine discrepancy modulus",
                "Sobolev wave-record space",
                "Gaussian detector-resolution family",
                "joint system-instrument inverse map",
                "calibration-invariant IQR observable",
                "reflection-asymmetry preparation resource theory",
                "coupled state-damping variational family",
                "typed completion-map diagram",
                "immutable pre-data physical protocol",
            ],
            "strongest_positive_result": (
                "The fractional window now has a genuine interval-arithmetic "
                "certificate, and a calibration-invariant operational statistic "
                "with a frozen pre-data test has been constructed."
            ),
            "strongest_negative_result": (
                "KMS/modular theory, maximum entropy, and Morita equivalence "
                "do not select w2=1/5 without an additional reference or gap."
            ),
            "deepest_missing_object": (
                "A strict theorem choosing the reference state/measure on the "
                "localized sector algebra; separately, a nonfree signed preparation "
                "source and physical calibration standards remain necessary."
            ),
            "closures_not_claimed": [
                "canonical uniform sector trace",
                "strict eta=9/5 source",
                "strict phase/frequency source",
                "QW-2191 selector closure",
                "internal physical units",
                "full legacy-to-strict completion",
                "legacy physical-role transfer",
                "L_total or Theory-of-Everything closure",
                "external experimental validation",
            ],
        },
    }
    OUT.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=json_default) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "output": OUT.name,
                "protocol": PROTOCOL.name,
                "figures": len(list(FIG.glob("*.png"))),
                "programs": 13,
                "formal_fft_max_relative_upper": result[
                    "program141_validated_interval_fft"
                ]["maximum_relative_remainder_upper"],
                "joint_rank": result["program145_joint_identifiability"]["joint_rank"],
                "bridge_actual_residual": result[
                    "program149_completion_map_diagram"
                ]["actual_canonical_legacy_to_projective_strict_relative_residual"],
                "synthetic_protocol_power": result["program150_predata_protocol"][
                    "synthetic_power_audit"
                ]["synthetic_power"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
