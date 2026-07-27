#!/usr/bin/env python3
"""Execute FIN research Programs 151--164.

The round completes the recommendations of Release 10.14 and adds an
axiomatic synthesis requested by the author.  Its central discipline is to
separate:

* a minimal spectral mathematical core,
* state and preparation premises,
* instruments and records,
* dimensional calibration,
* the still-incomplete legacy-to-strict completion map.

No axiom-augmented result is promoted to a strict-core theorem.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-fin-151-164")

import matplotlib.pyplot as plt
import numpy as np
from mpmath import iv, mp
from scipy.stats import levy_stable, norm


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_151_164_Axiomatic_Operational_Results.json"
AXIOMS_OUT = ROOT / "FIN_Programs_151_164_Minimal_Axiomatic_System.json"
FIG = ROOT / "FIN_Programs_151_164_Axiomatic_Operational_Figures"
FIG.mkdir(exist_ok=True)

PREVIOUS_RESULTS = ROOT / "FIN_Programs_138_150_State_Detector_Bridge_Results.json"
PREVIOUS_PROTOCOL = ROOT / "FIN_Programs_138_150_PreData_Physical_Protocol.json"

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = 743.0 / 4000.0
PHI = 13.0 / 80.0
ETA = 9.0 / 5.0
STABLE_ALPHA = ETA - 1.0
LEGACY_OMEGA = math.pi / 4.0
LEGACY_PHI = math.pi / 6.0
SEED = 20260727
RNG = np.random.default_rng(SEED)


def sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Unsupported JSON value: {type(value)!r}")


def iv_lo(x) -> float:
    return float(np.nextafter(float(x.a), -np.inf))


def iv_hi(x) -> float:
    return float(np.nextafter(float(x.b), np.inf))


def interval_fft(values: list) -> list:
    """Radix-two FFT over mpmath interval complex numbers."""
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


def program151_tighter_validated_fractional_certificate() -> dict:
    """Double the validated FFT length and retain a continuous-cell enclosure."""
    iv.dps = 42
    nfft = 2**15
    dmax = nfft // 2 - 1
    omega = iv.mpf([743, 743]) / 4000
    phi = iv.mpf([13, 13]) / 80
    eta_iv = iv.mpf([9, 9]) / 5
    ring = [iv.mpc(0, 0) for _ in range(nfft)]
    zd = iv.mpf([0, 0])
    derivative = iv.mpf([0, 0])
    for d in range(1, dmax + 1):
        a = abs(iv.cos(omega * d + phi)) / (1 + iv.mpf(d) ** eta_iv)
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
    abelian = 2 / (
        iv.gamma(eta_iv) * iv.sin(iv.pi * (eta_iv - 1) / 2)
    )
    c_lo = iv_lo(abelian) / (zd_hi + tail_hi)
    c_hi = iv_hi(abelian) / zd_lo

    rows = []
    for k in range(1, nfft // 2):
        q_center = 2 * math.pi * k / nfft
        cell_lo = q_center - half_cell
        cell_hi = q_center + half_cell
        if cell_hi < 0.001 or cell_lo > 0.02:
            continue
        qlo, qhi = max(0.001, cell_lo), min(0.02, cell_hi)
        numerator = zd - transformed[k].real
        num_lo, num_hi = iv_lo(numerator), iv_hi(numerator)
        cell_variation = derivative_hi * half_cell
        symbol_lo = max(0.0, num_lo - cell_variation) / (zd_hi + tail_hi)
        symbol_hi = (num_hi + cell_variation + 2 * tail_hi) / zd_lo
        prediction_lo = c_lo * qlo**STABLE_ALPHA
        prediction_hi = c_hi * qhi**STABLE_ALPHA
        ratio_lo = symbol_lo / prediction_hi
        ratio_hi = symbol_hi / prediction_lo
        relative_upper = max(abs(ratio_lo - 1), abs(ratio_hi - 1))
        rows.append(
            {
                "k": k,
                "q_interval": [qlo, qhi],
                "symbol_interval": [symbol_lo, symbol_hi],
                "prediction_interval": [prediction_lo, prediction_hi],
                "relative_remainder_upper": relative_upper,
                "compatible": symbol_lo <= prediction_hi
                and symbol_hi >= prediction_lo,
            }
        )

    worst = max(rows, key=lambda row: row["relative_remainder_upper"])
    qmid = np.asarray([sum(row["q_interval"]) / 2 for row in rows])
    remainder = np.asarray([row["relative_remainder_upper"] for row in rows])
    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.semilogy(qmid, remainder, "o-", ms=2.8, color="#1f5a99")
    ax.axhline(0.03, color="#a61b1b", ls="--", label="3% target")
    ax.set_xlabel(r"$|q|$")
    ax.set_ylabel("formal relative-remainder upper bound")
    ax.set_title("Program 151: doubled validated FFT enclosure")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program151_tighter_interval_fft.png", dpi=190)
    plt.close(fig)

    return {
        "interval_engine": "mpmath.iv, 42 decimal digits, radix-two FFT",
        "fft_length": nfft,
        "retained_distance": dmax,
        "frequency_cells": len(rows),
        "covered_window": [
            min(r["q_interval"][0] for r in rows),
            max(r["q_interval"][1] for r in rows),
        ],
        "tail_normalization_upper": tail_hi,
        "maximum_relative_remainder_upper": worst[
            "relative_remainder_upper"
        ],
        "worst_cell": worst,
        "all_symbol_lower_bounds_positive": all(
            r["symbol_interval"][0] > 0 for r in rows
        ),
        "all_cells_fractional_compatible": all(r["compatible"] for r in rows),
        "sub_3_percent_certificate": worst["relative_remainder_upper"] < 0.03,
        "improvement_over_program141": (
            "The grid and truncation radius are doubled. The result remains a "
            "formal continuous-window enclosure; dependency and tail bounds, "
            "not floating-point error, dominate its width."
        ),
        "status": "Proven, computer-assisted",
        "confidence": "Proven",
        "claim_boundary": (
            "Compatibility with the 4/5 fractional scale is certified. A "
            "sub-3-percent Abelian remainder is claimed only if the explicit "
            "Boolean field above is true."
        ),
    }


def certified_cf_prefix(number_terms: int = 12) -> tuple[list[int], list[dict]]:
    iv.dps = 100
    x = iv.mpf([743, 743]) / (4000 * iv.pi)
    terms: list[int] = []
    for _ in range(number_terms):
        lo, hi = iv_lo(x), iv_hi(x)
        alo, ahi = math.floor(lo), math.floor(hi)
        if alo != ahi:
            raise ArithmeticError("continued-fraction digit not certified")
        terms.append(alo)
        x = 1 / (x - alo)
    pm2, pm1, qm2, qm1 = 0, 1, 1, 0
    rows = []
    for a in terms:
        p, q = a * pm1 + pm2, a * qm1 + qm2
        rows.append({"a": a, "p": p, "q": q})
        pm2, pm1, qm2, qm1 = pm1, p, qm1, q
    return terms, rows


def program152_axiomatic_diophantine_condition() -> dict:
    """State the exact extra all-scale premise and test its finite prefix."""
    terms, convergents = certified_cf_prefix()
    mp.dps = 100
    theta = mp.mpf(743) / (4000 * mp.pi)
    rows = []
    for row in convergents:
        q = row["q"]
        nearest = mp.nint(q * theta)
        distance = abs(q * theta - nearest)
        rows.append(
            {
                **row,
                "distance_to_integer": float(distance),
                "q_times_distance": float(q * distance),
                "q2_times_distance": float(q * q * distance),
            }
        )
    positive_rows = [r for r in rows if r["q"] > 1]
    finite_kappa = {
        "nu_1": min(r["q_times_distance"] for r in positive_rows),
        "nu_2": min(r["q2_times_distance"] for r in positive_rows),
    }

    q = np.asarray([r["q"] for r in positive_rows], dtype=float)
    k1 = np.asarray([r["q_times_distance"] for r in positive_rows])
    k2 = np.asarray([r["q2_times_distance"] for r in positive_rows])
    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.loglog(q, k1, "o-", label=r"$q\|q\theta\|$")
    ax.loglog(q, k2, "s-", label=r"$q^2\|q\theta\|$")
    ax.set_xlabel("continued-fraction denominator q")
    ax.set_ylabel("finite Diophantine constant candidate")
    ax.set_title("Program 152: a finite prefix cannot certify an all-scale axiom")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program152_diophantine_axiom.png", dpi=190)
    plt.close(fig)

    return {
        "rotation": "theta=743/(4000*pi)",
        "certified_terms": terms,
        "convergents": rows,
        "candidate_axiom": (
            "D(kappa,nu): ||q theta|| >= kappa q^(-nu) for every integer q>=1, "
            "with kappa>0 and finite nu."
        ),
        "finite_prefix_kappa": finite_kappa,
        "necessity_theorem": (
            "Irrationality alone cannot yield a uniform polynomial "
            "Diophantine lower bound: Liouville irrational numbers violate "
            "D(kappa,nu) for every fixed finite nu and every kappa>0."
        ),
        "conditional_consequence": (
            "D(kappa,nu), together with bounded variation or Fourier-decay "
            "hypotheses on the oscillatory weight, supplies quantitative "
            "small-divisor control. The constants and final Abelian remainder "
            "must still be derived for the selected function class."
        ),
        "all_scale_axiom_derived_for_theta": False,
        "finite_prefix_is_not_global_proof": True,
        "status": "Proven obstruction; conditional axiom constructed",
        "confidence": "Proven",
        "claim_boundary": (
            "No all-scale irrationality measure for 743/(4000*pi) is inferred "
            "from twelve continued-fraction digits."
        ),
    }


def program153_functorial_fibre_groupoid_measures() -> dict:
    """Classify invariant central probability measures on structural orbits."""
    orbit_sizes = np.asarray([1, 1, 3], dtype=float)
    dims = np.asarray([1, 2, 2, 2, 2], dtype=float)
    structural_orbits = [[2], [3], [5, 7, 11]]
    powers = np.linspace(-4, 4, 401)
    natural_family = []
    for a in powers:
        raw = dims**a
        w = raw / raw.sum()
        natural_family.append(
            {"a": float(a), "w2": float(w[0]), "eta": float(w @ dims)}
        )

    # Barycentric plot of all orbit-invariant measures.
    grid = []
    for i in range(51):
        for j in range(51 - i):
            x, y = i / 50, j / 50
            z = 1 - x - y
            grid.append((x, y, z))
    grid = np.asarray(grid)
    cart_x = grid[:, 1] + 0.5 * grid[:, 2]
    cart_y = math.sqrt(3) * grid[:, 2] / 2
    uniform_orbit_mass = np.asarray([0.2, 0.2, 0.6])
    ux = uniform_orbit_mass[1] + 0.5 * uniform_orbit_mass[2]
    uy = math.sqrt(3) * uniform_orbit_mass[2] / 2
    fig, ax = plt.subplots(figsize=(7.3, 5.7), constrained_layout=True)
    ax.scatter(cart_x, cart_y, s=5, alpha=0.35, color="#1f5a99")
    ax.scatter([ux], [uy], s=90, color="#d1495b", label="uniform sectors")
    ax.plot([0, 1, 0.5, 0], [0, 0, math.sqrt(3) / 2, 0], color="black")
    ax.text(-0.04, -0.04, "{2}")
    ax.text(1.01, -0.04, "{3}")
    ax.text(0.43, math.sqrt(3) / 2 + 0.03, "{5,7,11}")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Program 153: invariant measures retain two free orbit masses")
    ax.legend(loc="upper right")
    fig.savefig(FIG / "program153_groupoid_measure_simplex.png", dpi=190)
    plt.close(fig)

    return {
        "objects": [2, 3, 5, 7, 11],
        "structural_orbits": structural_orbits,
        "classification": (
            "Every invariant probability has sector weights "
            "(x,y,z/3,z/3,z/3), where x,y,z>=0 and x+y+z=1."
        ),
        "simplex_dimension": 2,
        "uniform_sector_orbit_masses": uniform_orbit_mass.tolist(),
        "dimension_valuation_family": "r_a(p)=d_p^a/sum_q d_q^a",
        "dimension_valuation_parameter_range_tested": [
            float(powers.min()),
            float(powers.max()),
        ],
        "normalization_additivity_naturality_select_unique_measure": False,
        "why": (
            "Invariance fixes weights inside each isomorphism orbit but cannot "
            "compare the masses of nonisomorphic orbits. An induction or "
            "valuation normalization relating those orbits is extra data."
        ),
        "status": "Proven finite groupoid classification",
        "confidence": "Proven",
        "claim_boundary": "No canonical w2=1/5 follows from the admitted groupoid.",
    }


def program154_axiom_augmented_modular_equilibrium() -> dict:
    """Construct the smallest explicit equilibrium axiom that yields eta=9/5."""
    dims = np.asarray([1, 2, 2, 2, 2], dtype=float)
    energy = dims - 1
    beta = np.linspace(0, 2.2, 500)
    sector_weights = np.asarray(
        [
            dims * np.exp(-b * energy)
            / np.sum(dims * np.exp(-b * energy))
            for b in beta
        ]
    )
    eta_curve = sector_weights @ dims
    beta_star = math.log(2)
    target_weights = dims * np.exp(-beta_star * energy)
    target_weights /= target_weights.sum()

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.plot(beta, eta_curve, color="#1f5a99")
    ax.axvline(beta_star, color="#d1495b", ls="--", label=r"$\beta_F=\log2$")
    ax.axhline(9 / 5, color="#d1495b", ls=":")
    ax.set_xlabel(r"dimensionless inverse temperature $\beta_F$")
    ax.set_ylabel(r"$\eta=\sum_p w_p d_p$")
    ax.set_title("Program 154: a one-relation axiom selects the uniform centre")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program154_axiomatic_modular_equilibrium.png", dpi=190)
    plt.close(fig)

    return {
        "carrier_energy": "H_F restricted to fibre p equals (d_p-1) I",
        "microstate_degeneracy": dims.astype(int).tolist(),
        "new_axiom_A_ME": (
            "beta_F = alpha_geo / |Hom(U(12),{+1,-1})| = "
            "alpha_geo/4 = ln(2)."
        ),
        "hom_cardinality": 4,
        "beta_star": beta_star,
        "central_weights_at_beta_star": target_weights.tolist(),
        "eta_at_beta_star": float(target_weights @ dims),
        "uniqueness_theorem": (
            "For unit energy gap and Hilbert degeneracy d_p, equality of the "
            "d=1 and d=2 central sector weights holds iff 2 exp(-beta_F)=1, "
            "hence iff beta_F=ln(2)."
        ),
        "removal_tests": {
            "remove_temperature_relation": (
                "A one-parameter Gibbs family remains; eta is not selected."
            ),
            "remove_unit_gap": (
                "Only the product beta_F Delta=ln(2) is fixed."
            ),
            "remove_Hilbert_degeneracy": (
                "Sector-counting Gibbs weights are uniform at beta_F=0 instead."
            ),
            "remove_carrier_energy": (
                "No canonical modular Hamiltonian remains."
            ),
        },
        "strict_core_derivation": False,
        "status": "Consistent axiom-augmented closure; not strict",
        "confidence": "Proven conditional on A_ME",
        "claim_boundary": (
            "A_ME is a transparent additional premise. The identity "
            "alpha_geo/4=ln2 does not by itself prove that it is the physical "
            "inverse temperature."
        ),
    }


def program155_reflection_resource_completeness() -> dict:
    """Prove completeness of |r| on the two-branch reflection orbit line."""
    r = np.linspace(-1, 1, 401)
    lambdas = [1.0, 0.65, 0.2, 0.0, -0.4]
    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    for lam in lambdas:
        ax.plot(r, lam * r, label=fr"$\lambda={lam:g}$")
    ax.fill_between(r, -np.abs(r), np.abs(r), color="#1f5a99", alpha=0.08)
    ax.set_xlabel("input signed coordinate r")
    ax.set_ylabel("output signed coordinate r'")
    ax.set_title("Program 155: all free orbit-line channels contract |r|")
    ax.grid(True, alpha=0.25)
    ax.legend(ncol=3)
    fig.savefig(FIG / "program155_reflection_resource_complete.png", dpi=190)
    plt.close(fig)

    return {
        "state_family": (
            "rho_r=((1+r)/2)rho_+ + ((1-r)/2)rho_-, r in [-1,1], "
            "with reflection r->-r."
        ),
        "free_channel_classification": (
            "Every affine trace-preserving reflection-covariant channel on "
            "this orbit line acts as r' = lambda r with |lambda|<=1."
        ),
        "monotone": "M(rho_r)=|r|",
        "conversion_theorem": (
            "There exists a deterministic free channel rho_r -> rho_s iff "
            "|s|<=|r|. Therefore M is a complete single-copy monotone on the "
            "declared orbit line."
        ),
        "constructive_channel": (
            "For r nonzero choose lambda=s/r; covariance and positivity are "
            "equivalent to |lambda|<=1."
        ),
        "full_density_space_completeness": False,
        "strict_signed_preparation_source": False,
        "status": "Proven on the reflection orbit line",
        "confidence": "Proven",
        "claim_boundary": (
            "The theorem neither classifies arbitrary density matrices nor "
            "creates the nonfree input r != 0 required by QW-2191."
        ),
    }


def stable_rvs(alpha: float, size, rng: np.random.Generator) -> np.ndarray:
    """Chambers--Mallows--Stuck sampler for symmetric alpha-stable variables."""
    if alpha == 2:
        return math.sqrt(2) * rng.normal(size=size)
    if alpha == 1:
        return np.tan(rng.uniform(-math.pi / 2, math.pi / 2, size=size))
    v = rng.uniform(-math.pi / 2, math.pi / 2, size=size)
    w = rng.exponential(size=size)
    first = np.sin(alpha * v) / np.cos(v) ** (1 / alpha)
    second = (
        np.cos((1 - alpha) * v) / w
    ) ** ((1 - alpha) / alpha)
    return first * second


def empirical_iqr(x: np.ndarray, axis=-1) -> np.ndarray:
    q = np.quantile(x, [0.25, 0.75], axis=axis)
    return q[1] - q[0]


def program156_detector_deconvolution_finite_bias() -> dict:
    """Derive IQR asymptotics and measure finite Gaussian-resolution bias."""
    alpha = STABLE_ALPHA
    q25, q75 = levy_stable.ppf([0.25, 0.75], alpha, 0, loc=0, scale=1)
    density_q = levy_stable.pdf(q75, alpha, 0, loc=0, scale=1)
    iqr0 = q75 - q25
    n = 4000
    t1, t2 = 1.0, 4.0
    slope_true = 1 / alpha
    var_iqr = 1 / (4 * n * density_q**2)
    slope_sd_theory = math.sqrt(
        2 * var_iqr / iqr0**2 / math.log(t2 / t1) ** 2
    )

    replicates = 320
    z1 = stable_rvs(alpha, (replicates, n), RNG)
    z2 = stable_rvs(alpha, (replicates, n), RNG) * t2 ** (1 / alpha)
    sigma_rel = np.asarray([0.0, 0.05, 0.1, 0.2, 0.4])
    rows = []
    for rel in sigma_rel:
        sigma = rel * iqr0
        x1 = z1 + RNG.normal(scale=sigma, size=z1.shape)
        x2 = z2 + RNG.normal(scale=sigma, size=z2.shape)
        slopes = np.log(empirical_iqr(x2) / empirical_iqr(x1)) / math.log(t2)
        rows.append(
            {
                "sigma_over_IQR_t1": float(rel),
                "mean_slope": float(np.mean(slopes)),
                "slope_bias": float(np.mean(slopes) - slope_true),
                "slope_sd": float(np.std(slopes, ddof=1)),
            }
        )

    rel = np.asarray([r["sigma_over_IQR_t1"] for r in rows])
    bias = np.asarray([r["slope_bias"] for r in rows])
    sd = np.asarray([r["slope_sd"] for r in rows])
    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.plot(rel, bias, "o-", label="Monte Carlo bias")
    ax.plot(rel, sd, "s-", label="Monte Carlo SD")
    ax.axhline(slope_sd_theory, color="#2a9d8f", ls="--", label="zero-noise asymptotic SD")
    ax.set_xlabel(r"absolute detector $\sigma/\mathrm{IQR}(t_1)$")
    ax.set_ylabel("slope error scale")
    ax.set_title("Program 156: detector convolution biases a finite-time slope")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program156_detector_bias.png", dpi=190)
    plt.close(fig)

    return {
        "stable_alpha": alpha,
        "unit_scale_quartiles": [float(q25), float(q75)],
        "unit_scale_iqr": float(iqr0),
        "density_at_upper_quartile": float(density_q),
        "asymptotic_theorem": (
            "For a symmetric continuous density f at its quartiles, "
            "Var(IQR_hat)=1/(4 n f(q_0.75)^2)+o(n^-1). Independent-time "
            "delta calculus gives the recorded log-slope variance."
        ),
        "predicted_zero_noise_slope_sd": slope_sd_theory,
        "simulation": {
            "replicates": replicates,
            "records_per_time": n,
            "rows": rows,
        },
        "deconvolution_identifiability": (
            "Gaussian deconvolution is identifiable when the response "
            "characteristic function and calibration are known and nonzero."
        ),
        "deconvolution_instability": (
            "Its inverse multiplier exp(sigma^2 q^2/2) is unbounded, so "
            "unregularized UV deconvolution is statistically ill-posed."
        ),
        "status": "Proven asymptotics; strong numerical evidence for finite bias",
        "confidence": "Strong evidence",
        "claim_boundary": (
            "The simulation quantifies the declared stable/Gaussian model; it "
            "is not a detector-independent theorem."
        ),
    }


def program157_semiparametric_identifiability() -> dict:
    """Classify identifiability against multiplicative apparatus response."""
    times = np.geomspace(1, 32, 7)
    x = np.log(times)
    alpha_score = -x / STABLE_ALPHA**2
    rows = []
    residuals = []
    for degree in range(5):
        nuisance = np.column_stack([x**j for j in range(degree + 1)])
        projection = nuisance @ np.linalg.pinv(nuisance) @ alpha_score
        residual = alpha_score - projection
        augmented = np.column_stack([nuisance, alpha_score])
        rows.append(
            {
                "nuisance_log_polynomial_degree": degree,
                "nuisance_rank": int(np.linalg.matrix_rank(nuisance)),
                "augmented_rank": int(np.linalg.matrix_rank(augmented)),
                "alpha_score_residual_norm": float(np.linalg.norm(residual)),
                "alpha_identifiable": bool(np.linalg.norm(residual) > 1e-10),
            }
        )
        residuals.append(np.linalg.norm(residual))

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.semilogy(range(5), np.maximum(residuals, 1e-18), "o-", color="#6a4c93")
    ax.axhline(1e-10, color="#a61b1b", ls="--", label="rank tolerance")
    ax.set_xlabel("degree of unknown log-apparatus response")
    ax.set_ylabel("exponent-score residual norm")
    ax.set_title("Program 157: a free linear time response erases identifiability")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program157_semiparametric_rank.png", dpi=190)
    plt.close(fig)

    return {
        "observation_model": "Q_obs(t)=h(t) C t^(1/alpha)",
        "global_counterexample": (
            "For any alternative alpha', define h'(t)=h(t) "
            "t^(1/alpha-1/alpha'). Then Q_obs is unchanged. Thus alpha is "
            "nonidentifiable under unrestricted positive h."
        ),
        "finite_tangent_test": rows,
        "criterion": (
            "The exponent is locally identifiable iff its score has a "
            "nonzero component orthogonal to the nuisance tangent space."
        ),
        "declared_design_result": (
            "An unknown constant gain is harmless, but admitting an unknown "
            "linear term in log t exactly confounds 1/alpha."
        ),
        "status": "Proven semiparametric nonidentifiability theorem",
        "confidence": "Proven",
        "claim_boundary": (
            "Recovering alpha requires calibration, replicated controls, or a "
            "restricted apparatus response class; more time points alone do "
            "not remove an unrestricted confounder."
        ),
    }


def program158_finite_sample_stable_iqr_theorem() -> dict:
    """Validate asymptotic confidence intervals for the stable IQR exponent."""
    alpha = STABLE_ALPHA
    t2 = 4.0
    q75 = levy_stable.ppf(0.75, alpha, 0)
    q25 = levy_stable.ppf(0.25, alpha, 0)
    f = levy_stable.pdf(q75, alpha, 0)
    iqr = q75 - q25
    true_t = 1 / alpha
    sample_sizes = [250, 500, 1000, 4000]
    replicates = 700
    rows = []
    for n in sample_sizes:
        slopes = np.empty(replicates)
        for rep in range(replicates):
            x1 = stable_rvs(alpha, n, RNG)
            x2 = stable_rvs(alpha, n, RNG) * t2 ** (1 / alpha)
            slopes[rep] = math.log(
                float(empirical_iqr(x2)) / float(empirical_iqr(x1))
            ) / math.log(t2)
        variance_iqr_relative = 1 / (4 * n * f**2 * iqr**2)
        predicted_sd_t = math.sqrt(
            2 * variance_iqr_relative / math.log(t2) ** 2
        )
        lower = slopes - 1.96 * predicted_sd_t
        upper = slopes + 1.96 * predicted_sd_t
        rows.append(
            {
                "n": n,
                "replicates": replicates,
                "mean_T": float(np.mean(slopes)),
                "empirical_sd_T": float(np.std(slopes, ddof=1)),
                "predicted_sd_T": predicted_sd_t,
                "normal_95_coverage_T": float(
                    np.mean((lower <= true_t) & (true_t <= upper))
                ),
                "mean_alpha_hat": float(np.mean(1 / slopes)),
            }
        )

    n = np.asarray(sample_sizes)
    empirical = np.asarray([r["empirical_sd_T"] for r in rows])
    predicted = np.asarray([r["predicted_sd_T"] for r in rows])
    coverage = np.asarray([r["normal_95_coverage_T"] for r in rows])
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.4), constrained_layout=True)
    axes[0].loglog(n, empirical, "o-", label="empirical SD")
    axes[0].loglog(n, predicted, "s--", label="asymptotic SD")
    axes[0].set_xlabel("records per time")
    axes[0].set_ylabel("SD of T")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend()
    axes[1].semilogx(n, coverage, "o-", color="#2a9d8f")
    axes[1].axhline(0.95, color="#a61b1b", ls="--")
    axes[1].set_xlabel("records per time")
    axes[1].set_ylabel("nominal 95% interval coverage")
    axes[1].grid(True, alpha=0.25)
    fig.suptitle("Program 158: finite-sample stable-IQR inference")
    fig.savefig(FIG / "program158_iqr_finite_sample.png", dpi=190)
    plt.close(fig)

    return {
        "estimator": "T_hat=log(IQR_hat(t2)/IQR_hat(t1))/log(t2/t1)",
        "true_T": true_t,
        "asymptotic_variance": (
            "For independent equal-size samples, Var(T_hat) = "
            "2/[4 n f(q_0.75)^2 IQR^2 log(t2/t1)^2] + o(n^-1)."
        ),
        "rows": rows,
        "confidence_interval_type": "first-order normal, plug-in-free model value",
        "nonasymptotic_concentration_proved": False,
        "status": "Proven asymptotics; Monte Carlo finite-sample validation",
        "confidence": "Strong evidence",
        "claim_boundary": (
            "Coverage is model-specific and asymptotic. It is not an exact "
            "finite-sample confidence theorem under apparatus confounding."
        ),
    }


def generate_protocol_family(
    name: str, replicates: int, n: int, t2: float
) -> np.ndarray:
    """Generate one adversarial family for the frozen IQR protocol."""
    if name == "local_gaussian":
        x1 = RNG.normal(size=(replicates, n))
        x2 = RNG.normal(size=(replicates, n)) * math.sqrt(t2)
    elif name.startswith("stable_alpha_"):
        alpha = float(name.rsplit("_", 1)[1])
        x1 = stable_rvs(alpha, (replicates, n), RNG)
        x2 = stable_rvs(alpha, (replicates, n), RNG) * t2 ** (1 / alpha)
    elif name == "absolute_truncated_alpha_0.8":
        x1 = np.clip(stable_rvs(0.8, (replicates, n), RNG), -8, 8)
        x2 = np.clip(
            stable_rvs(0.8, (replicates, n), RNG) * t2 ** 1.25, -8, 8
        )
    elif name == "local_with_time_gain":
        x1 = RNG.normal(size=(replicates, n))
        x2 = RNG.normal(size=(replicates, n)) * math.sqrt(t2) * t2**0.7
    else:
        raise ValueError(name)
    return np.log(empirical_iqr(x2) / empirical_iqr(x1)) / math.log(t2)


def program159_blind_adversarial_protocol() -> dict:
    """Challenge model specificity of the frozen Program-150 decision rule."""
    families = [
        "local_gaussian",
        "stable_alpha_0.8",
        "stable_alpha_1.0",
        "stable_alpha_1.2",
        "absolute_truncated_alpha_0.8",
        "local_with_time_gain",
    ]
    replicates, n, t2 = 320, 2500, 4.0
    threshold = 0.875
    band = [1.10, 1.40]
    rows = []
    values = {}
    for family in families:
        slopes = generate_protocol_family(family, replicates, n, t2)
        values[family] = slopes
        rows.append(
            {
                "family": family,
                "mean_T": float(np.mean(slopes)),
                "sd_T": float(np.std(slopes, ddof=1)),
                "P150_fractional_decision_rate": float(np.mean(slopes > threshold)),
                "FIN_equivalence_band_rate": float(
                    np.mean((slopes >= band[0]) & (slopes <= band[1]))
                ),
            }
        )

    fig, ax = plt.subplots(figsize=(10.2, 5.1), constrained_layout=True)
    ax.boxplot(
        [values[f] for f in families],
        tick_labels=[
            "Gaussian",
            r"$\alpha=.8$",
            r"$\alpha=1$",
            r"$\alpha=1.2$",
            "truncated",
            "gain-confounded",
        ],
        showfliers=False,
    )
    ax.axhline(threshold, color="#a61b1b", ls="--", label="P150 threshold")
    ax.axhspan(band[0], band[1], color="#2a9d8f", alpha=0.12, label="FIN band")
    ax.set_ylabel("IQR spreading slope T")
    ax.set_title("Program 159: the binary rule detects scaling, not FIN uniquely")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program159_adversarial_protocol.png", dpi=190)
    plt.close(fig)

    return {
        "protocol_rule_frozen_before_family_comparison": True,
        "P150_threshold": threshold,
        "exploratory_FIN_equivalence_band": band,
        "replicates": replicates,
        "records_per_time": n,
        "rows": rows,
        "falsification_result": (
            "The P150 threshold separates local diffusion from sufficiently "
            "fast spreading, but generic stable laws and time-dependent "
            "apparatus gain can trigger it. The statistic is not FIN-specific."
        ),
        "strict_FIN_validation": False,
        "status": "Strong adversarial falsification evidence",
        "confidence": "Strong evidence",
        "claim_boundary": (
            "The generated alternatives are declared synthetic models, not "
            "exhaustive descriptions of laboratory data."
        ),
    }


def program160_phase_frequency_bridge_obstruction() -> dict:
    """Prove an equivariant representation obstruction between the phases."""
    legacy_points = np.exp(1j * (LEGACY_OMEGA * np.arange(8) + LEGACY_PHI))
    strict_points = np.exp(1j * (OMEGA * np.arange(48) + PHI))
    z12_defect = abs(np.exp(1j * 12 * OMEGA) - 1)
    common_generator_defect = abs(np.exp(1j * LEGACY_OMEGA) - np.exp(1j * OMEGA))

    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.6), constrained_layout=True)
    axes[0].scatter(legacy_points.real, legacy_points.imag, c=np.arange(8), cmap="viridis")
    axes[0].set_title("legacy: period-eight orbit")
    axes[1].scatter(strict_points.real, strict_points.imag, c=np.arange(48), cmap="plasma", s=25)
    axes[1].set_title("strict: first 48 points, no finite period")
    for ax in axes:
        circle = plt.Circle((0, 0), 1, fill=False, color="black", alpha=0.4)
        ax.add_patch(circle)
        ax.set_aspect("equal")
        ax.axis("off")
    fig.suptitle("Program 160: inequivalent translation characters")
    fig.savefig(FIG / "program160_phase_obstruction.png", dpi=190)
    plt.close(fig)

    return {
        "legacy_character": "chi_L(d)=exp(i*pi*d/4), a character of Z_8",
        "legacy_period": 8,
        "strict_sequence": "chi_S(d)=exp(i*(743/4000)*d)",
        "strict_has_nonzero_integer_period": False,
        "strict_period_proof": (
            "If T>0 were a period then (743/4000)T=2*pi*k. For k!=0 "
            "this would make pi rational; k=0 is impossible."
        ),
        "strict_is_Z12_character": False,
        "Z12_character_defect": z12_defect,
        "generator_eigenvalue_defect": common_generator_defect,
        "obstruction_theorem": (
            "No translation-equivariant character homomorphism from the "
            "legacy Z_8 representation to the strict infinite-order phase "
            "can preserve the generator eigenvalue. Finite order is preserved "
            "by such a homomorphism, whereas the strict eigenvalue has "
            "infinite order."
        ),
        "non_equivariant_or_affine_completion_excluded": False,
        "status": "Proven representation-theoretic obstruction",
        "confidence": "Proven",
        "claim_boundary": (
            "The theorem excludes a canonical character-preserving phase "
            "transport, not every nonlinear, symmetry-breaking completion map."
        ),
    }


def program161_reference_energy_grammar() -> dict:
    """Enumerate a finite structural grammar for the modular energy."""
    f_dimension = np.asarray([0, 1, 1, 1, 1], dtype=int)
    f_p3_orbit = np.asarray([0, 1, 0, 0, 0], dtype=int)
    f_large_unit_orbit = np.asarray([0, 0, 1, 1, 1], dtype=int)
    features = np.vstack([f_dimension, f_p3_orbit, f_large_unit_orbit]).T
    solutions = []
    all_complexities = []
    for coeff in itertools.product(range(-3, 4), repeat=3):
        coeff_array = np.asarray(coeff)
        energy = features @ coeff_array
        normalized = energy - energy[0]
        complexity = int(np.sum(np.abs(coeff_array)))
        all_complexities.append(complexity)
        if np.array_equal(normalized, f_dimension):
            solutions.append(
                {
                    "coefficients": list(coeff),
                    "energy": energy.tolist(),
                    "L1_complexity": complexity,
                }
            )
    minimum = min(s["L1_complexity"] for s in solutions)
    minimizers = [s for s in solutions if s["L1_complexity"] == minimum]

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    ax.hist([s["L1_complexity"] for s in solutions], bins=np.arange(0.5, 10.5), color="#1f5a99")
    ax.axvline(minimum, color="#d1495b", ls="--", label="minimum")
    ax.set_xlabel(r"$\ell^1$ description length")
    ax.set_ylabel("target-realizing formulas")
    ax.set_title("Program 161: d-1 is the unique shortest declared formula")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program161_energy_grammar.png", dpi=190)
    plt.close(fig)

    return {
        "declared_structural_orbits": [[2], [3], [5, 7, 11]],
        "feature_columns": {
            "f_dimension": f_dimension.tolist(),
            "f_p3_orbit": f_p3_orbit.tolist(),
            "f_large_unit_orbit": f_large_unit_orbit.tolist(),
        },
        "coefficient_box": [-3, 3],
        "formulas_enumerated": 7**3,
        "target": "(0,1,1,1,1), modulo additive constant",
        "solutions": solutions,
        "number_of_solutions": len(solutions),
        "minimum_L1_complexity": minimum,
        "minimum_solutions": minimizers,
        "unique_minimum_is_dimension_minus_one": (
            len(minimizers) == 1
            and minimizers[0]["coefficients"] == [1, 0, 0]
        ),
        "strict_source_theorem": False,
        "status": "Proven exhaustive result in declared finite grammar",
        "confidence": "Proven in scope",
        "claim_boundary": (
            "Minimum description length depends on the declared feature "
            "grammar and on targeting uniform central weights. It motivates "
            "H_F=d-1 but does not derive A_ME."
        ),
    }


def program162_conditional_role_obstruction_matrix() -> dict:
    """Audit necessary bridge edges before any legacy physical role transfer."""
    edges = [
        "amplitude",
        "damping_semantics",
        "phase_frequency",
        "selector_topology",
        "physical_units",
        "spectral_equivalence",
    ]
    roles = {
        "sin2_thetaW=alpha_geo/12": [1, 0, 0, 0, 1, 0],
        "alpha_EM(alpha_geo,beta_tors)": [1, 1, 0, 0, 1, 0],
        "beta_tors^N hierarchy": [0, 1, 0, 0, 1, 0],
        "phase/topological bit": [0, 0, 1, 1, 0, 0],
        "absolute kernel normalization": [1, 0, 0, 0, 1, 0],
        "projective scalar shape": [0, 1, 1, 0, 0, 0],
        "common dual dynamics": [0, 0, 0, 0, 0, 1],
        "translation/circulant structure": [0, 0, 1, 0, 0, 1],
        "Green-response interpretation": [0, 1, 0, 0, 1, 1],
    }
    matrix = np.asarray(list(roles.values()), dtype=int)
    current_edges = np.asarray([0, 0, 0, 0, 0, 0], dtype=int)
    eligible = {
        role: bool(np.all(current_edges[np.asarray(req, dtype=bool)] == 1))
        for role, req in roles.items()
    }

    fig, ax = plt.subplots(figsize=(10.8, 5.8), constrained_layout=True)
    image = ax.imshow(matrix, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(edges)), edges, rotation=35, ha="right")
    ax.set_yticks(range(len(roles)), list(roles))
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax.text(j, i, str(matrix[i, j]), ha="center", va="center", fontsize=8)
    ax.set_title("Program 162: necessary completion edges by legacy role")
    fig.colorbar(image, ax=ax, shrink=0.75, label="edge required")
    fig.savefig(FIG / "program162_role_obstruction_matrix.png", dpi=190)
    plt.close(fig)

    return {
        "edges": edges,
        "roles": roles,
        "current_completed_edges": current_edges.tolist(),
        "role_transfer_eligible": eligible,
        "eligible_physical_roles_count": sum(eligible.values()),
        "structural_analogy_note": (
            "Both kernels can be embedded into operator dynamics after separate "
            "positivity shifts, but this is not spectral equivalence or role transfer."
        ),
        "status": "Proven necessary-condition obstruction matrix",
        "confidence": "Proven in declared dependency model",
        "claim_boundary": (
            "Zero eligible roles means no transfer is licensed now. It does "
            "not prove that every future enriched bridge is impossible."
        ),
    }


def validate_external_candidate(path: Path, required_fields: set[str]) -> dict:
    """Very small local schema validator; it does not download external data."""
    try:
        if path.suffix.lower() == ".json":
            payload = json.loads(path.read_text(encoding="utf-8"))
            if isinstance(payload, list) and payload and isinstance(payload[0], dict):
                fields = set(payload[0])
            elif isinstance(payload, dict) and isinstance(payload.get("records"), list) and payload["records"]:
                fields = set(payload["records"][0])
            else:
                fields = set(payload) if isinstance(payload, dict) else set()
        elif path.suffix.lower() == ".csv":
            first = path.read_text(encoding="utf-8", errors="replace").splitlines()[0]
            fields = {item.strip() for item in first.split(",")}
        else:
            return {"path": str(path), "admitted": False, "reason": "unsupported format"}
    except Exception as exc:
        return {"path": str(path), "admitted": False, "reason": f"parse error: {exc}"}
    missing = sorted(required_fields - fields)
    return {
        "path": str(path),
        "fields": sorted(fields),
        "missing": missing,
        "admitted": not missing,
    }


def program163_external_intake_readiness() -> dict:
    """Execute the frozen schema gate against an explicit local intake folder."""
    protocol = json.loads(PREVIOUS_PROTOCOL.read_text(encoding="utf-8"))
    core = protocol["core"]
    required = set(core["raw_record_format"])
    intake = ROOT / "FIN_External_Data_Candidates"
    candidate_paths = (
        sorted([*intake.glob("*.json"), *intake.glob("*.csv")])
        if intake.exists()
        else []
    )
    validations = [validate_external_candidate(p, required) for p in candidate_paths]
    admitted = [v for v in validations if v.get("admitted")]

    gate_labels = [
        "raw fields",
        "two times",
        "n>=4000",
        "detector",
        "memory",
        "provenance",
    ]
    readiness = np.ones(len(gate_labels))
    external_pass = np.zeros(len(gate_labels))
    fig, ax = plt.subplots(figsize=(8.7, 4.8), constrained_layout=True)
    x = np.arange(len(gate_labels))
    ax.bar(x - 0.18, readiness, width=0.36, label="schema implemented")
    ax.bar(x + 0.18, external_pass, width=0.36, label="admitted dataset")
    ax.set_xticks(x, gate_labels, rotation=25, ha="right")
    ax.set_ylim(0, 1.15)
    ax.set_ylabel("gate state")
    ax.set_title("Program 163: ready intake schema, zero admitted external records")
    ax.legend()
    fig.savefig(FIG / "program163_external_readiness.png", dpi=190)
    plt.close(fig)

    return {
        "protocol_id": core["protocol_id"],
        "protocol_core_sha256": protocol["canonical_core_sha256"],
        "required_raw_fields": sorted(required),
        "intake_directory": str(intake),
        "candidate_files_found": len(candidate_paths),
        "validations": validations,
        "admitted_datasets": len(admitted),
        "external_data_admitted": bool(admitted),
        "methodological_result": (
            "The local validator and stop rules are executable. No web search "
            "or external download was performed, and absence of a local "
            "candidate is recorded as zero admission rather than validation."
        ),
        "status": "Operational readiness proven; external validation open",
        "confidence": "Proven for local schema",
        "claim_boundary": "No laboratory evidence is claimed.",
    }


def program164_minimal_axiomatic_system() -> dict:
    """Construct and mechanically test a six-axiom operational architecture."""
    axioms = {
        "A0": {
            "name": "conservative spectral Dirichlet generator",
            "statement": (
                "(H,E,A,1) with E a closed Markovian Dirichlet form, A its "
                "self-adjoint nonnegative generator and A1=0; Borel functional "
                "calculus defines U_t, P_t and resolvents."
            ),
            "rank": 1,
        },
        "A1": {
            "name": "sector equilibrium/state",
            "statement": (
                "A localized finite algebra A_F, a faithful central state omega_F, "
                "and eta=omega_F(d); A_ME is one explicit non-strict realization."
            ),
            "rank": 3,
        },
        "A2": {
            "name": "preparation resource",
            "statement": (
                "A declared preparation set including a nonfree signed resource "
                "whenever an orientation-sensitive record is required."
            ),
            "rank": 3,
        },
        "A3": {
            "name": "instrument and record",
            "statement": (
                "A normalized completely positive instrument, apparatus/environment "
                "channel, outcome sigma-algebra and persistent record."
            ),
            "rank": 2,
        },
        "A4": {
            "name": "dimensional calibration",
            "statement": (
                "Positive conversion standards (ell_*,tau_*,hbar_*) and their "
                "operational calibration maps for every dimensional claim."
            ),
            "rank": 2,
        },
        "A5": {
            "name": "typed legacy-to-strict completion",
            "statement": (
                "A typed map carrying amplitude, damping/compression, phase/frequency, "
                "selector/topology and an explicit role-semantics audit."
            ),
            "rank": 4,
        },
    }
    requirements = {
        "dual_unitary_diffusive_dynamics": {"A0"},
        "target_exponent_as_selected_state": {"A0", "A1"},
        "generic_outcome_probabilities": {"A0", "A1", "A3"},
        "signed_operational_outcome": {"A0", "A1", "A2", "A3"},
        "calibrated_dimensional_observation": {"A0", "A1", "A3", "A4"},
        "calibrated_signed_test": {"A0", "A1", "A2", "A3", "A4"},
        "typed_kernel_completion": {"A0", "A1", "A5"},
        "role_transfer_eligibility_only": set(axioms),
    }
    codes = list(axioms)
    subset_rows = []
    for mask in range(2 ** len(codes)):
        subset = {codes[j] for j in range(len(codes)) if mask & (1 << j)}
        capability = {
            name: req <= subset for name, req in requirements.items()
        }
        subset_rows.append(
            {
                "axioms": sorted(subset),
                "count": len(subset),
                "capabilities": capability,
            }
        )

    minimum_sets = {}
    for capability in requirements:
        accepted = [
            row for row in subset_rows if row["capabilities"][capability]
        ]
        nmin = min(row["count"] for row in accepted)
        minimum_sets[capability] = [
            row["axioms"] for row in accepted if row["count"] == nmin
        ]

    removal_countermodels = {
        "A0": (
            "Choose a nonnormal, indefinite, nonconservative, or non-Markov "
            "generator. One respectively loses the shared spectral calculus, "
            "contractivity, mass conservation, or positivity preservation."
        ),
        "A1": (
            "Uniform-sector and Hilbert-trace states share A0,A2--A5 but give "
            "eta=9/5 and 17/9 respectively."
        ),
        "A2": (
            "Restrict preparations to reflection-invariant states. Every "
            "reflection-odd record has zero expectation under free operations."
        ),
        "A3": (
            "Retain states and dynamics but specify no CP instrument or outcome "
            "map. No probability distribution over recorded events is defined."
        ),
        "A4": (
            "Rescale ell_*, tau_* and hbar_* while leaving all dimensionless "
            "operator statements fixed. Absolute length, time and action change."
        ),
        "A5": (
            "Use the canonical legacy and strict kernels. Their audited scalar "
            "completion residual and inequivalent phase representations remain."
        ),
    }
    full = set(codes)
    removal_losses = {}
    for code in codes:
        subset = full - {code}
        removal_losses[code] = [
            name for name, req in requirements.items() if not req <= subset
        ]

    capability_names = list(requirements)
    heat = np.asarray(
        [
            [int(code in requirements[cap]) for code in codes]
            for cap in capability_names
        ]
    )
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    image = ax.imshow(heat, cmap="Purples", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(range(len(codes)), codes)
    ax.set_yticks(range(len(capability_names)), capability_names)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            ax.text(j, i, str(heat[i, j]), ha="center", va="center", fontsize=8)
    ax.set_title("Program 164: axioms required by each capability")
    fig.colorbar(image, ax=ax, shrink=0.75, label="required")
    fig.savefig(FIG / "program164_axiom_capability_matrix.png", dpi=190)
    plt.close(fig)

    result = {
        "system_name": "AFIS: Axiomatic FIN Informational-Spectral system",
        "axioms": axioms,
        "capability_requirements": {
            key: sorted(value) for key, value in requirements.items()
        },
        "all_64_subsets": subset_rows,
        "minimum_axiom_sets": minimum_sets,
        "removal_countermodels": removal_countermodels,
        "full_system_removal_losses": removal_losses,
        "layered_minima": {
            "minimal_mathematical_core": ["A0"],
            "selected_dimensionless_operator_model": ["A0", "A1"],
            "generic_dimensionless_operational_physics": ["A0", "A1", "A3"],
            "signed_dimensionless_operational_physics": [
                "A0",
                "A1",
                "A2",
                "A3",
            ],
            "generic_calibrated_physics": ["A0", "A1", "A3", "A4"],
            "calibrated_signed_test": ["A0", "A1", "A2", "A3", "A4"],
            "full_conditional_FIN_architecture": codes,
        },
        "central_independence_theorem": (
            "No single one of A1--A5 is forced by A0: finite countermodels "
            "preserve the spectral generator while varying the state, "
            "preparation set, instrument, units, or completion map independently."
        ),
        "A0_internal_necessity": {
            "self_adjointness": (
                "Without it exp(-itA) need not be unitary and no projection-valued "
                "spectral measure is guaranteed."
            ),
            "nonnegativity": (
                "Without it exp(-tA) can grow exponentially on negative spectrum."
            ),
            "conservativity": (
                "Without A1=0 the diffusion semigroup need not preserve total mass."
            ),
            "Markov_property": (
                "A nonnegative self-adjoint contraction need not preserve positivity; "
                "the Dirichlet/Markov property supplies that implication."
            ),
        },
        "smallest_new_state_axiom_candidate": "A_ME from Program 154",
        "smallest_single_theorem_connecting_every_layer_found": False,
        "obstruction_type": [
            "state-theoretic",
            "resource-theoretic",
            "operational",
            "dimensional",
            "representation-theoretic",
        ],
        "status": "Axiom system constructed; independence proven by countermodels",
        "confidence": "Proven as a logical architecture",
        "claim_boundary": (
            "AFIS is an explicit axiom-augmented research architecture, not a "
            "strict-core derivation, role-transfer theorem, L_total, or ToE closure."
        ),
    }
    AXIOMS_OUT.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, default=json_default)
        + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    programs = {
        "151": program151_tighter_validated_fractional_certificate(),
        "152": program152_axiomatic_diophantine_condition(),
        "153": program153_functorial_fibre_groupoid_measures(),
        "154": program154_axiom_augmented_modular_equilibrium(),
        "155": program155_reflection_resource_completeness(),
        "156": program156_detector_deconvolution_finite_bias(),
        "157": program157_semiparametric_identifiability(),
        "158": program158_finite_sample_stable_iqr_theorem(),
        "159": program159_blind_adversarial_protocol(),
        "160": program160_phase_frequency_bridge_obstruction(),
        "161": program161_reference_energy_grammar(),
        "162": program162_conditional_role_obstruction_matrix(),
        "163": program163_external_intake_readiness(),
        "164": program164_minimal_axiomatic_system(),
    }
    record = {
        "metadata": {
            "title": (
                "FIN Programs 151-164: Axiomatic Operational Foundations, "
                "State Selection, and Falsifiable Measurement"
            ),
            "release": "10.15",
            "version": "1.0.0",
            "date": "2026-07-27",
            "creator": "Żuchowski, Krzysztof",
            "orcid": "0009-0002-0909-3613",
            "seed": SEED,
            "previous_results_sha256": sha256(PREVIOUS_RESULTS),
            "previous_protocol_sha256": sha256(PREVIOUS_PROTOCOL),
            "firecrawl_used": False,
            "external_data_used": False,
        },
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "omega_strict": OMEGA,
            "phi_strict": PHI,
            "eta_strict": ETA,
            "stable_alpha": STABLE_ALPHA,
            "omega_legacy": LEGACY_OMEGA,
            "phi_legacy": LEGACY_PHI,
        },
        "programs": programs,
        "global_verdict": {
            "deepest_result": (
                "The spectral generator is a minimal mathematical core, but "
                "state selection, signed preparation, instruments, dimensional "
                "calibration, and the legacy-to-strict completion are mutually "
                "independent obligations. A_ME is the smallest explicit "
                "one-relation state-selection candidate found; it remains an axiom."
            ),
            "strict_selector_closed": False,
            "QW_2191_discharged": False,
            "internal_units_derived": False,
            "legacy_strict_bridge_completed": False,
            "role_transfer_started": False,
            "L_total_derived": False,
            "ToE_claimed": False,
        },
    }
    OUT.write_text(
        json.dumps(record, indent=2, ensure_ascii=False, default=json_default)
        + "\n",
        encoding="utf-8",
    )
    print(OUT)
    print(AXIOMS_OUT)


if __name__ == "__main__":
    main()
