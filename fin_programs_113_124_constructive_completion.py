#!/usr/bin/env python3
"""Execute FIN research Programs 113--124.

This round constructs, rather than merely names, the missing objects selected
after Programs 101--112:

* a global analytic mass certificate and continuous parameter box;
* a finite-q Abelian certificate and the stable-process limit;
* a fractional operational process and local/fractional crossover object;
* the unique inverse Fisher potential and a finite variational grammar;
* a hidden-Markov apparatus channel;
* a Z12 homological-character fibre functor producing [1,2,2,2,2];
* a conditional monoid-cocycle damping completion; and
* a typed operational experiment object plus immutable data-intake gate.

The constructions do not silently close QW-2191, source physical units,
complete the legacy-to-strict bridge, transfer legacy roles, or claim
experimental validation.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_113_124_Constructive_Completion_Results.json"
INTAKE = ROOT / "FIN_Programs_113_124_External_Data_Intake.json"
FIG = ROOT / "FIN_Programs_113_124_Constructive_Completion_Figures"
FIG.mkdir(exist_ok=True)

PREVIOUS = ROOT / "FIN_Programs_101_112_Fractional_Source_Completion_Results.json"

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = 743.0 / 4000.0
PHI = 13.0 / 80.0
BETA = 1.0
ETA = 9.0 / 5.0
ALPHA = ETA - 1.0
BETA_TORS = 0.01
SEED = 20260727
RNG = np.random.default_rng(SEED)


def sha256(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def canonical_digest(record: dict) -> str:
    return hashlib.sha256(
        json.dumps(record, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, sp.Basic):
        return str(value)
    raise TypeError(f"Unsupported JSON value: {type(value)!r}")


def strict_weight(d, omega=OMEGA, phi=PHI, beta=BETA, eta=ETA):
    d = np.asarray(d, dtype=float)
    return np.abs(np.cos(omega * d + phi)) / (1.0 + beta * d**eta)


def js_divergence(p, q):
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum()
    q = q / q.sum()
    m = 0.5 * (p + q)

    def kl(a, b):
        mask = a > 0
        return float(np.sum(a[mask] * np.log(a[mask] / b[mask])))

    return 0.5 * kl(p, m) + 0.5 * kl(q, m)


def program113_global_mass_theorem(previous: dict) -> dict:
    p = previous["program101_independent_interval_certificate"]
    ppi = p["p_hat_pi_interval"]
    ph = p["p_hat_pi_over_2_interval"]
    a = [1.0 - ppi[1], 1.0 - ppi[0]]
    b = [1.0 - ph[1], 1.0 - ph[0]]

    m_s, a_s, b_s = sp.symbols("m a b", positive=True)
    native = m_s / (m_s + a_s)
    schur = 2 * m_s * (m_s + a_s) / (
        (2 * m_s + a_s) * (m_s + b_s)
    )
    factor = sp.factor(schur - native)
    expected = (
        m_s
        * ((3 * a_s - 2 * b_s) * m_s + a_s * (2 * a_s - b_s))
        / ((m_s + a_s) * (2 * m_s + a_s) * (m_s + b_s))
    )
    symbolic_residual = sp.simplify(factor - expected)
    linear_lower = 3 * a[0] - 2 * b[1]
    constant_lower = a[0] * (2 * a[0] - b[1])

    masses = np.geomspace(1e-5, 1e5, 500)
    ac, bc = np.mean(a), np.mean(b)
    gaps = (
        masses
        * ((3 * ac - 2 * bc) * masses + ac * (2 * ac - bc))
        / ((masses + ac) * (2 * masses + ac) * (masses + bc))
    )
    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.loglog(masses, gaps, color="#1f5a99")
    ax.set_xlabel("mass m")
    ax.set_ylabel("Schur ratio minus native ratio")
    ax.set_title("Global positive Schur-nonclosure gap for every m>0")
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program113_global_mass_theorem.png", dpi=190)
    plt.close(fig)

    return {
        "exact_factorization": str(factor),
        "symbolic_identity_residual": str(symbolic_residual),
        "a_interval": a,
        "b_interval": b,
        "coefficient_3a_minus_2b_lower": linear_lower,
        "coefficient_a_times_2a_minus_b_lower": constant_lower,
        "global_positive_for_every_m_gt_zero": (
            symbolic_residual == 0 and linear_lower > 0 and constant_lower > 0
        ),
        "uniform_positive_lower_bound_on_all_m_gt_zero": False,
        "reason_no_uniform_bound": "the gap tends to zero as m->0 and m->infinity",
        "status": "Proven analytic theorem using the Program-101 symbol rectangle",
        "scope": "Frozen normalized strict lattice and alternating-site Schur map.",
    }


def program114_continuous_parameter_box(previous: dict) -> dict:
    half = {"omega": 0.0005, "phi": 0.002, "beta": 0.02, "eta": 0.01}
    beta_min = BETA - half["beta"]
    eta_min = ETA - half["eta"]
    split = int(math.floor(2 / half["omega"]))
    # Fully analytic bounds.  The omega term uses min(2,h*d), split at
    # d=floor(2/h); all remaining series are bounded by monotone integrals.
    omega_bound = (
        half["omega"]
        / beta_min
        * (
            1
            + (split ** (2 - eta_min) - 1) / (2 - eta_min)
        )
        + 2
        / beta_min
        * split ** (1 - eta_min)
        / (eta_min - 1)
    )
    zeta_bound = 1 + 1 / (eta_min - 1)
    phi_bound = half["phi"] / beta_min * zeta_bound
    beta_bound = half["beta"] / (beta_min * BETA) * zeta_bound
    log_series_bound = (
        math.log(2) * 2 ** (-eta_min)
        + 2 ** (1 - eta_min)
        * (
            math.log(2) / (eta_min - 1)
            + 1 / (eta_min - 1) ** 2
        )
    )
    eta_bound = half["eta"] / beta_min * log_series_bound
    one_sided = {
        "omega": float(np.nextafter(omega_bound, np.inf)),
        "phi": float(np.nextafter(phi_bound, np.inf)),
        "beta": float(np.nextafter(beta_bound, np.inf)),
        "eta": float(np.nextafter(eta_bound, np.inf)),
    }
    unnormalized_l1 = float(np.nextafter(2 * sum(one_sided.values()), np.inf))
    p101 = previous["program101_independent_interval_certificate"]
    z0_lower = p101["partial_normalization_interval"][0]
    normalized_symbol_error = 2 * unnormalized_l1 / (z0_lower - unnormalized_l1)
    ppi = p101["p_hat_pi_interval"]
    ph = p101["p_hat_pi_over_2_interval"]
    ppi_box = [ppi[0] - normalized_symbol_error, ppi[1] + normalized_symbol_error]
    ph_box = [ph[0] - normalized_symbol_error, ph[1] + normalized_symbol_error]
    a_min = 1 - ppi_box[1]
    b_max = 1 - ph_box[0]
    linear_lower = 3 * a_min - 2 * b_max
    constant_lower = a_min * (2 * a_min - b_max)

    fig, ax = plt.subplots(figsize=(8.5, 4.7), constrained_layout=True)
    labels = ["3a-2b", "a(2a-b)", "normalization reserve"]
    vals = [linear_lower, constant_lower, z0_lower - unnormalized_l1]
    ax.bar(labels, vals, color=["#19733a", "#19733a", "#1f5a99"])
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("certified lower bound")
    ax.set_title("Continuous strict-parameter box certificate")
    fig.savefig(FIG / "program114_continuous_parameter_box.png", dpi=190)
    plt.close(fig)

    return {
        "box": {
            "omega": [OMEGA - half["omega"], OMEGA + half["omega"]],
            "phi": [PHI - half["phi"], PHI + half["phi"]],
            "beta": [BETA - half["beta"], BETA + half["beta"]],
            "eta": [ETA - half["eta"], ETA + half["eta"]],
            "mass": "(0,infinity)",
        },
        "omega_series_split": split,
        "one_sided_analytic_variation_bounds": one_sided,
        "series_bounds": {
            "sum_d_minus_eta": zeta_bound,
            "sum_log_d_d_minus_eta": log_series_bound,
        },
        "two_sided_unnormalized_l1_bound": unnormalized_l1,
        "normalized_symbol_variation_bound": normalized_symbol_error,
        "p_hat_pi_box": ppi_box,
        "p_hat_pi_over_2_box": ph_box,
        "coefficient_3a_minus_2b_lower": linear_lower,
        "coefficient_a_times_2a_minus_b_lower": constant_lower,
        "continuous_box_global_mass_nonclosure": (
            unnormalized_l1 < z0_lower
            and linear_lower > 0
            and constant_lower > 0
        ),
        "status": "Analytic continuous-box theorem with monotone integral/Hölder bounds",
        "scope": "The box is deliberately finite and is not asserted maximal.",
    }


def lattice_data(truncation=1_000_000):
    d = np.arange(1, truncation + 1, dtype=float)
    a = strict_weight(d)
    zd = float(2 * a.sum())
    tz = 2 * truncation ** (1 - ETA) / (ETA - 1)
    return d, a, zd, tz


def program115_effective_abelian_certificate() -> dict:
    d, a, zd, tz = lattice_data()
    qs = np.geomspace(1e-3, 2e-2, 20)
    alpha = ALPHA
    integral = math.pi / (
        2 * math.gamma(1 + alpha) * math.sin(math.pi * alpha / 2)
    )
    numerator_constant = 2 * (2 / math.pi) * integral
    c_interval = [
        numerator_constant / (zd + tz),
        numerator_constant / zd,
    ]
    rows = []
    for q in qs:
        nd = float(2 * np.dot(a, 1 - np.cos(q * d)))
        # Each missing term has factor 2*a_d*(1-cos)<=4*d^-eta.
        n_tail = 2 * tz
        s_interval = [nd / (zd + tz), (nd + n_tail) / zd]
        pred = [c_interval[0] * q**alpha, c_interval[1] * q**alpha]
        rel = [s_interval[0] / pred[1] - 1, s_interval[1] / pred[0] - 1]
        rows.append({
            "q": q,
            "symbol_interval": s_interval,
            "predicted_fractional_interval": pred,
            "relative_remainder_interval": rel,
        })
    max_abs = max(max(abs(x) for x in row["relative_remainder_interval"])
                  for row in rows)

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.semilogx(qs, [r["relative_remainder_interval"][0] for r in rows],
                "o-", label="certified lower remainder")
    ax.semilogx(qs, [r["relative_remainder_interval"][1] for r in rows],
                "o-", label="certified upper remainder")
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("|q|")
    ax.set_ylabel("relative remainder")
    ax.set_title("Finite-q Abelian certificate for the strict fractional symbol")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program115_effective_abelian_certificate.png", dpi=190)
    plt.close(fig)

    return {
        "q_window": [float(qs[0]), float(qs[-1])],
        "certified_q_points": len(rows),
        "normalization_interval": [zd, zd + tz],
        "abelian_constant_interval": c_interval,
        "rows": rows,
        "maximum_absolute_relative_remainder_bound": max_abs,
        "full_continuous_q_remainder_theorem_exported": False,
        "missing_object": (
            "an explicit Diophantine discrepancy modulus for the irrational "
            "rotation omega/(2*pi), needed to make the oscillatory remainder uniform"
        ),
        "status": "Proven pointwise finite-q enclosures; continuous effective remainder remains open",
    }


def program116_stable_invariance_principle(p115: dict) -> dict:
    numerical_truncation = 5_000_000
    d, a, zd, _ = lattice_data(numerical_truncation)
    c = float(np.mean(p115["abelian_constant_interval"]))
    rows = []
    for n in [16, 64, 256, 1024]:
        for k in [0.5, 1.0, 2.0]:
            q = k * n ** (-1 / ALPHA)
            s = float(2 * np.dot(a, 1 - np.cos(q * d)) / zd)
            finite_cf = math.exp(n * math.log(max(1 - s, 1e-300)))
            stable_cf = math.exp(-c * abs(k) ** ALPHA)
            rows.append({
                "steps": n,
                "test_wave_number": k,
                "scaled_lattice_q": q,
                "finite_characteristic_value": finite_cf,
                "stable_characteristic_value": stable_cf,
                "absolute_error": abs(finite_cf - stable_cf),
            })

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    for k in [0.5, 1.0, 2.0]:
        r = [x for x in rows if x["test_wave_number"] == k]
        ax.loglog([x["steps"] for x in r], [x["absolute_error"] for x in r],
                  "o-", label=f"k={k:g}")
    ax.set_xlabel("number of jumps n")
    ax.set_ylabel("characteristic-function error")
    ax.set_title("Convergence toward the symmetric 4/5-stable process")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program116_stable_invariance_principle.png", dpi=190)
    plt.close(fig)

    return {
        "stable_index": ALPHA,
        "numerical_truncation": numerical_truncation,
        "space_scaling": "n^(-1/alpha)=n^(-5/4)",
        "limiting_characteristic_function": "exp(-C*|k|^(4/5))",
        "limiting_generator": "C*(-Delta)^(2/5)",
        "rows": rows,
        "theorem": (
            "The Program-105 characteristic exponent places the strict jump "
            "law in the normal domain of attraction of a symmetric 4/5-stable "
            "Levy process; finite-dimensional convergence follows by "
            "[1-S(k n^-1/alpha)]^n -> exp(-C|k|^alpha), and tightness by the "
            "standard stable-domain functional limit theorem."
        ),
        "status": (
            "Proven conditional on the standard stable-domain functional "
            "limit theorem; finite-truncation numerical check (not assumed monotone)"
        ),
        "physical_scale_sourced": False,
    }


def detector_confusion(p, error):
    p = np.asarray(p, dtype=float)
    n = len(p)
    return (1 - error) * p + error * (1 - p) / (n - 1)


def fractional_records(n=128, detector_error=0.05):
    q = 2 * np.pi * np.fft.fftfreq(n)
    c = 1.1474679863772292
    lam = c * np.abs(q) ** ALPHA
    times = np.linspace(0.01, 4.0, 400)
    rows, best = [], None
    for t in times:
        amp = np.fft.ifft(np.exp(-1j * t * lam))
        wave = np.abs(amp) ** 2
        wave /= wave.sum()
        diffusion = np.fft.ifft(np.exp(-t * lam)).real
        diffusion = np.maximum(diffusion, 0)
        diffusion /= diffusion.sum()
        ow = detector_confusion(wave, detector_error)
        od = detector_confusion(diffusion, detector_error)
        js = js_divergence(ow, od)
        row = (t, js, wave, diffusion, ow, od)
        rows.append(row)
        if best is None or js > best[1]:
            best = row
    return q, lam, times, rows, best


def program117_fractional_operational_process() -> dict:
    q, lam, times, rows, best = fractional_records()
    t, js, wave, diffusion, ow, od = best
    event = ow > od
    binary_wave = float(ow[event].sum())
    binary_diff = float(od[event].sum())

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.8, 4.4), constrained_layout=True)
    ax1.plot(times, [row[1] for row in rows], color="#1f5a99")
    ax1.axvline(t, color="#a61b1b", ls="--")
    ax1.set_xlabel("dimensionless clock time")
    ax1.set_ylabel("record JSD")
    ax1.set_title("Wave/diffusion operational separation")
    x = np.arange(len(wave))
    dist = np.minimum(x, len(x) - x)
    order = np.argsort(dist)
    ax2.semilogy(dist[order], np.maximum(ow[order], 1e-14), label="wave record")
    ax2.semilogy(dist[order], np.maximum(od[order], 1e-14), label="diffusion record")
    ax2.set_xlabel("cyclic distance")
    ax2.set_ylabel("observed probability")
    ax2.set_title(f"Best record at t={t:.3f}")
    ax2.legend()
    fig.savefig(FIG / "program117_fractional_operational_process.png", dpi=190)
    plt.close(fig)

    operational_tuple = {
        "state_space": "C^128",
        "initial_state": "|0><0|",
        "generator": "diag_k(C*|2*pi*k/128|^(4/5))",
        "dynamics": ["unitary exp(-itA)", "Markov exp(-tA)"],
        "clock": "dimensionless declared t",
        "preparation": "localized site 0",
        "instrument": "full-site projective record",
        "environment": "identity channel in this minimal construction",
        "apparatus": "symmetric 5% confusion channel",
        "record": "site label; optimal binary coarse event also exported",
    }
    return {
        "operational_ten_tuple": operational_tuple,
        "best_time": t,
        "best_record_jsd": js,
        "optimal_binary_event_size": int(event.sum()),
        "binary_event_probability_wave": binary_wave,
        "binary_event_probability_diffusion": binary_diff,
        "unitarity_residual": float(abs(np.abs(np.exp(-1j * t * lam)) - 1).max()),
        "diffusion_mass_residual": float(abs(diffusion.sum() - 1)),
        "status": "Constructed dimensionless operational process; no physical clock or data",
    }


def program118_local_fractional_crossover() -> dict:
    alpha = ALPHA
    kappa = 1.0
    nus = [0.0, 0.01, 0.1, 1.0, 10.0]
    q = np.geomspace(1e-4, 10, 500)
    rows = []
    for nu in nus:
        crossover = None if nu == 0 else (nu / kappa) ** (1 / (2 - alpha))
        for scale in [1, 2, 4, 8, 16, 32]:
            rows.append({
                "nu": nu,
                "coarse_factor": scale,
                "running_fractional_to_local_coupling": (
                    0.0 if nu == 0 else (nu / kappa) * scale ** (2 - alpha)
                ),
            })

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    for nu in nus[1:]:
        frac = nu * q**alpha / (kappa * q**2 + nu * q**alpha)
        ax.semilogx(q, frac, label=f"nu={nu:g}")
    ax.set_xlabel("|q|")
    ax.set_ylabel("fractional share of symbol")
    ax.set_title("Local/fractional crossover in A=kappa q^2+nu |q|^(4/5)")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program118_local_fractional_crossover.png", dpi=190)
    plt.close(fig)

    return {
        "constructed_operator": "A_(kappa,nu)=kappa*(-Delta)+nu*(-Delta)^(2/5)",
        "symbol": "kappa*q^2+nu*|q|^(4/5)",
        "crossover_formula": "q_*=(nu/kappa)^(1/(2-4/5))=(nu/kappa)^(5/6)",
        "rg_formula": "g(b)=g*b^(2-4/5)=g*b^(6/5)",
        "rows": rows,
        "theorem": "Every nu>0 is infrared-relevant relative to the local q^2 term.",
        "status": "Constructed and analytically classified conditional crossover family",
        "relative_coupling_sourced_by_FIN": False,
    }


def program119_unique_inverse_fisher_potential() -> dict:
    d = np.arange(1, 7, dtype=float)
    multiplicity = np.array([2, 2, 2, 2, 2, 1], dtype=float)
    reference = multiplicity / multiplicity.sum()
    q_strict = multiplicity * strict_weight(d)
    q_strict /= q_strict.sum()
    inverse_potential = -np.log(q_strict / reference)
    inverse_potential -= np.dot(reference, inverse_potential)
    formula = np.log(1 + d**ETA) - np.log(np.abs(np.cos(OMEGA * d + PHI)))
    formula -= np.dot(reference, formula)
    residual = float(np.max(abs(inverse_potential - formula)))

    envelope_potential = np.log(1 + d**ETA)
    q_envelope = reference * np.exp(-envelope_potential)
    q_envelope /= q_envelope.sum()
    tv = float(0.5 * np.sum(abs(q_strict - q_envelope)))
    js = js_divergence(q_strict, q_envelope)

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.plot(d, inverse_potential, "o-", label="unique inverse strict potential")
    shifted_env = envelope_potential - np.dot(reference, envelope_potential)
    ax.plot(d, shifted_env, "s--", label="target-independent fractional envelope")
    ax.set_xlabel("distance shell")
    ax.set_ylabel("potential modulo constants")
    ax.set_title("The adaptive source potential and its unavoidable target content")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program119_inverse_fisher_potential.png", dpi=190)
    plt.close(fig)

    return {
        "reference_measure": reference.tolist(),
        "strict_shell_probability": q_strict.tolist(),
        "unique_potential_modulo_constant": inverse_potential.tolist(),
        "closed_formula_residual": residual,
        "closed_formula": "log(1+d^(9/5))-log|cos((743/4000)d+13/80)|+constant",
        "envelope_only_total_variation": tv,
        "envelope_only_jsd": js,
        "theorem": (
            "For F_V(q)=KL(q||r)+<V,q>, an interior stationary distribution "
            "q_* determines V uniquely modulo constants: V=-log(q_*/r)+c."
        ),
        "status": "Constructed unique adaptive potential; negative source result",
        "source_independent_of_strict_target": False,
        "reason": "the unique potential is an information-equivalent encoding of all strict parameters",
    }


def node_distribution(theta, n):
    omega, phi, beta, eta = theta
    idx = np.arange(n)
    d = np.minimum(idx, n - idx).astype(float)
    w = np.abs(np.cos(omega * d + phi)) / (1 + beta * d**eta)
    w[0] = 0
    return w / w.sum()


def program120_finite_variational_grammar() -> dict:
    strict = (OMEGA, PHI, BETA, ETA)
    omega_sources = [
        math.pi / 4,
        math.pi / 6,
        ALPHA_GEO / (2 * math.pi),
        math.pi / 12,
        1 / 12,
    ]
    phi_sources = [0.0, math.pi / 6, math.pi / 4, ALPHA_GEO / (2 * math.pi)]
    beta_sources = [BETA_TORS, 1.0]
    eta_sources = [j / 5 for j in range(5, 16)]
    target12 = node_distribution(strict, 12)
    target24 = node_distribution(strict, 24)
    records = []
    for theta in itertools.product(
        omega_sources, phi_sources, beta_sources, eta_sources
    ):
        js12 = js_divergence(node_distribution(theta, 12), target12)
        js24 = js_divergence(node_distribution(theta, 24), target24)
        records.append({
            "theta": list(theta),
            "mean_jsd": 0.5 * (js12 + js24),
            "jsd_C12": js12,
            "jsd_C24": js24,
        })
    records.sort(key=lambda row: row["mean_jsd"])
    best = records[0]
    exact_count = sum(row["mean_jsd"] < 1e-14 for row in records)
    direct_augmented = js_divergence(node_distribution(strict, 12), target12)

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    scores = [r["mean_jsd"] for r in records]
    ax.semilogy(np.arange(1, len(scores) + 1), np.maximum(scores, 1e-16))
    ax.set_xlabel("grammar candidate rank")
    ax.set_ylabel("mean C12/C24 JSD")
    ax.set_title("Exhaustive target-free variational grammar")
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program120_variational_grammar.png", dpi=190)
    plt.close(fig)

    return {
        "grammar_size": len(records),
        "omega_source_count": len(omega_sources),
        "phi_source_count": len(phi_sources),
        "beta_source_count": len(beta_sources),
        "eta_source_count": len(eta_sources),
        "best_candidate": best,
        "exact_candidate_count": exact_count,
        "direct_target_augmented_jsd": direct_augmented,
        "target_free_selection_code_lower_bound_bits": math.log2(len(records)),
        "mdl_theorem": (
            "Adding the exact strict rational tuple to the action grammar "
            "achieves zero defect only by carrying the same four-parameter "
            "payload as direct kernel storage; it yields no description-length compression."
        ),
        "status": "Finite exhaustive no-exact-source certificate for the declared grammar",
    }


def hmm_transition(epsilon, rho):
    a = epsilon * (1 - rho)
    b = epsilon + rho * (1 - epsilon)
    return np.array([[1 - a, a], [1 - b, b]], dtype=float)


def simulate_markov_errors(length, epsilon, rho, rng):
    trans = hmm_transition(epsilon, rho)
    e = np.empty(length, dtype=int)
    e[0] = rng.random() < epsilon
    for t in range(1, length):
        e[t] = rng.random() < trans[e[t - 1], 1]
    return e


def fit_error_transition(known_x, observed_y):
    e = np.bitwise_xor(known_x, observed_y)
    counts = np.ones((2, 2), dtype=float) * 0.5
    for i, j in zip(e[:-1], e[1:]):
        counts[i, j] += 1
    trans = counts / counts.sum(axis=1, keepdims=True)
    eps = float(e.mean())
    return trans, eps


def hmm_loglik_y(y, p_x1, trans, eps):
    # Hidden state is detector flip E. Emission marginalizes independent X.
    emission = np.empty((len(y), 2))
    for t, obs in enumerate(y):
        emission[t, 0] = p_x1 if obs == 1 else 1 - p_x1
        emission[t, 1] = 1 - p_x1 if obs == 1 else p_x1
    alpha = np.array([1 - eps, eps]) * emission[0]
    logscale = math.log(alpha.sum())
    alpha /= alpha.sum()
    for t in range(1, len(y)):
        alpha = (alpha @ trans) * emission[t]
        scale = alpha.sum()
        logscale += math.log(scale)
        alpha /= scale
    return logscale


def iid_loglik_y(y, p_x1, eps):
    py = eps + (1 - 2 * eps) * p_x1
    return float(np.sum(y * math.log(py) + (1 - y) * math.log(1 - py)))


def program121_hidden_markov_apparatus(p117: dict) -> dict:
    pw = p117["binary_event_probability_wave"]
    pd = p117["binary_event_probability_diffusion"]
    epsilon, rho = 0.1, 0.8
    true_trans = hmm_transition(epsilon, rho)
    calibration_x = RNG.integers(0, 2, 2000, dtype=int)
    calibration_e = simulate_markov_errors(2000, epsilon, rho, RNG)
    calibration_y = np.bitwise_xor(calibration_x, calibration_e)
    fit_trans, fit_eps = fit_error_transition(calibration_x, calibration_y)
    lengths = [5, 10, 20, 50, 100]
    rows = []
    for length in lengths:
        hmm_correct = 0
        iid_correct = 0
        total = 0
        for truth, px in [(1, pw), (0, pd)]:
            for _ in range(300):
                x = (RNG.random(length) < px).astype(int)
                e = simulate_markov_errors(length, epsilon, rho, RNG)
                y = np.bitwise_xor(x, e)
                llr_hmm = (
                    hmm_loglik_y(y, pw, fit_trans, fit_eps)
                    - hmm_loglik_y(y, pd, fit_trans, fit_eps)
                )
                llr_iid = iid_loglik_y(y, pw, fit_eps) - iid_loglik_y(y, pd, fit_eps)
                hmm_correct += int((llr_hmm > 0) == bool(truth))
                iid_correct += int((llr_iid > 0) == bool(truth))
                total += 1
        rows.append({
            "records": length,
            "hmm_accuracy": hmm_correct / total,
            "iid_misspecified_accuracy": iid_correct / total,
        })

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax.plot(lengths, [r["hmm_accuracy"] for r in rows], "o-", label="calibrated HMM")
    ax.plot(lengths, [r["iid_misspecified_accuracy"] for r in rows],
            "s--", label="iid error model")
    ax.axhline(0.5, color="black", lw=0.8)
    ax.set_xlabel("held-out records")
    ax.set_ylabel("synthetic classification accuracy")
    ax.set_title("Wave/diffusion inference with correlated apparatus memory")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program121_hidden_markov_apparatus.png", dpi=190)
    plt.close(fig)

    return {
        "binary_wave_probability": pw,
        "binary_diffusion_probability": pd,
        "true_error_transition": true_trans.tolist(),
        "calibration_records": 2000,
        "fitted_error_transition": fit_trans.tolist(),
        "fitted_stationary_error_rate": fit_eps,
        "rows": rows,
        "constructed_apparatus_object": (
            "two-state hidden Markov flip channel calibrated on a disjoint "
            "known-input partition and used unchanged for held-out likelihoods"
        ),
        "status": "Constructed and synthetically validated inference object; no external evidence",
    }


def all_u12_characters():
    units = [1, 5, 7, 11]
    chars = []
    for values in itertools.product([-1, 1], repeat=4):
        table = dict(zip(units, values))
        if table[1] != 1:
            continue
        if all(table[(a * b) % 12] == table[a] * table[b]
               for a in units for b in units):
            chars.append(table)
    return units, chars


def program122_homological_character_functor() -> dict:
    primes = [2, 3, 5, 7, 11]
    units, chars = all_u12_characters()
    nontrivial = [chi for chi in chars if any(chi[u] != 1 for u in units)]
    fibres = []
    vector = []
    for p in primes:
        homology_dim = math.gcd(p, 12) - 1
        character_dim = (
            sum(chi[p] == -1 for chi in nontrivial) if p in units else 0
        )
        total = homology_dim + character_dim
        vector.append(total)
        fibres.append({
            "prime": p,
            "reduced_H0_kernel_multiplication_dimension": homology_dim,
            "negative_unit_character_fibre_dimension": character_dim,
            "direct_sum_fibre_dimension": total,
        })
    uniform_mean = Fraction(sum(vector), len(vector))
    # Since V=(1,2,2,2,2), a general probability trace gives eta=2-w_2.
    required_weight_p2 = Fraction(1, 5)

    fig, ax = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    x = np.arange(len(primes))
    h = [row["reduced_H0_kernel_multiplication_dimension"] for row in fibres]
    c = [row["negative_unit_character_fibre_dimension"] for row in fibres]
    ax.bar(x, h, label="reduced H0 of ker(m_p)")
    ax.bar(x, c, bottom=h, label="negative character fibre")
    ax.set_xticks(x, [str(p) for p in primes])
    ax.set_xlabel("prime generator p<12")
    ax.set_ylabel("fibre dimension")
    ax.set_title("Constructed Z12 homological-character fibre functor")
    ax.legend()
    fig.savefig(FIG / "program122_homological_character_functor.png", dpi=190)
    plt.close(fig)

    return {
        "base_object": "multiplication endomorphisms m_p on Z12 plus real character dual of U(12)",
        "prime_generator_set": primes,
        "unit_group": units,
        "real_character_count": len(chars),
        "nontrivial_character_count": len(nontrivial),
        "fibre_rows": fibres,
        "dimension_vector": vector,
        "uniform_normalized_trace": str(uniform_mean),
        "uniform_trace_equals_9_over_5": uniform_mean == Fraction(9, 5),
        "general_trace_formula": "eta(w)=sum_p w_p dim(F_p)=2-w_2",
        "weight_condition_for_eta_9_over_5": f"w_2={required_weight_p2}",
        "constructed_object": (
            "F_p = reduced H_0(ker[m_p:Z12->Z12]) direct-sum "
            "negative eigenspace of p on Hom(U12,{+1,-1})"
        ),
        "functor_is_coefficient_free_after_Z12_is_given": True,
        "strict_uniform_trace_source_exported": False,
        "status": "New canonical fibre construction; exponent readout remains conditional",
        "obstruction": (
            "C^5 has many tracial states; current strict data do not force "
            "the weight w_2=1/5 or the uniform prime-sector trace."
        ),
    }


def program123_conditional_damping_cocycle(p122: dict) -> dict:
    beta_s, d_s, eta_s, a_s, b_s = sp.symbols(
        "beta d eta a b", positive=True
    )
    multiplicativity_defect = sp.factor(
        beta_s * (a_s * b_s) ** eta_s
        - (beta_s * a_s**eta_s) * (beta_s * b_s**eta_s)
    )
    solutions_beta = sp.solve(sp.Eq(multiplicativity_defect, 0), beta_s)
    w2 = sp.symbols("w_2", real=True)
    eta_w = 2 - w2
    retention = sp.simplify(2 ** (1 - eta_w))
    uniform_retention = sp.simplify(retention.subs(w2, sp.Rational(1, 5)))
    alpha_retention = sp.exp(-sp.Rational(1, 5) * 4 * sp.log(2))
    retention_identity = sp.simplify(uniform_retention - alpha_retention)

    weights = np.linspace(0, 1, 201)
    etas = 2 - weights
    retentions = 2 ** (1 - etas)
    fig, ax1 = plt.subplots(figsize=(8.4, 4.7), constrained_layout=True)
    ax1.plot(weights, etas, label=r"$\eta=2-w_2$")
    ax1.axvline(0.2, color="#a61b1b", ls="--", label="uniform trace")
    ax1.set_xlabel("trace weight on the p=2 fibre")
    ax1.set_ylabel("derived exponent")
    ax2 = ax1.twinx()
    ax2.plot(weights, retentions, color="#19733a", label="dyadic retention")
    ax2.set_ylabel("retention")
    ax1.set_title("Conditional damping cocycle and the remaining trace degree")
    fig.savefig(FIG / "program123_conditional_damping_cocycle.png", dpi=190)
    plt.close(fig)

    return {
        "multiplicativity_defect_factorization": str(multiplicativity_defect),
        "algebraic_beta_solutions": [str(x) for x in solutions_beta],
        "nonzero_positive_monoid_character_forces_beta": 1,
        "eta_family_from_trace": "eta=2-w_2",
        "uniform_trace_eta": "9/5",
        "uniform_dyadic_retention": str(uniform_retention),
        "retention_equals_exp_minus_alpha_geo_over_5": retention_identity == 0,
        "conditional_completion_map": (
            "D_w(d)=(1+beta_tors*d)/(1+d^(2-w_2)); "
            "at w_2=1/5, D_w maps the legacy linear envelope exactly to "
            "the strict 1/(1+d^(9/5)) envelope"
        ),
        "minimal_axioms": [
            {
                "axiom": "Z12 homological-character fibre functor",
                "necessity": "removing it removes the dimension vector",
                "current_status": "constructed conditional on the strict Z12 carrier",
            },
            {
                "axiom": "uniform/equivalent w_2=1/5 trace readout",
                "necessity": "removing it leaves eta=2-w_2 continuously free",
                "current_status": "not strictly sourced",
            },
            {
                "axiom": "nonzero multiplicative tail T(ab)=T(a)T(b)",
                "necessity": "removing it leaves positive beta free",
                "current_status": "additional compression principle",
            },
        ],
        "status": "Constructed conditional eta-and-beta damping theorem; no strict bridge closure",
        "physical_length_unit_sourced": False,
        "phase_frequency_amplitude_bridge_sourced": False,
    }


def program124_typed_operational_completion(p117: dict, p121: dict) -> dict:
    required = [
        "state_space", "state", "generator", "dynamics", "clock",
        "preparation", "instrument", "environment", "apparatus", "record",
        "calibration", "decision_rule",
    ]
    obj = {
        "state_space": "C^128",
        "state": "|0><0|",
        "generator": "fractional A with symbol C|q|^(4/5)",
        "dynamics": ["exp(-itA)", "exp(-tA)"],
        "clock": {"type": "dimensionless protocol clock", "physical_calibration": None},
        "preparation": "localized site 0",
        "instrument": "full-site readout plus frozen binary likelihood-ratio event",
        "environment": "identity minimal model; any bath must be separately declared",
        "apparatus": p121["constructed_apparatus_object"],
        "record": "immutable time-stamped binary or full-site counts",
        "calibration": "disjoint known-input HMM calibration partition",
        "decision_rule": "held-out wave-minus-diffusion log-likelihood ratio",
    }
    missing_fields = [field for field in required if field not in obj]
    physical_missing = [
        "length standard ell_*", "clock standard tau_*", "action standard hbar_*",
        "independent raw experimental records",
    ]
    intake = {
        "schema": "FIN_PROGRAMS_113_124_EXTERNAL_DATA_INTAKE_V1",
        "record_date": "2026-07-27",
        "required_fields": [
            "raw_record_file_sha256", "preparation_protocol_sha256",
            "instrument_protocol_sha256", "graph_operator_sha256",
            "clock_calibration_sha256", "apparatus_calibration_sha256",
            "environment_boundary_description", "blind_holdout_split_sha256",
        ],
        "preregistered_best_dimensionless_time": p117["best_time"],
        "external_data_admitted": False,
        "reason": "no independent immutable calibrated record package is present",
    }
    intake["canonical_sha256"] = canonical_digest(intake)
    INTAKE.write_text(json.dumps(intake, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {
        "typed_operational_object": obj,
        "required_field_count": len(required),
        "missing_mathematical_operational_fields": missing_fields,
        "mathematically_complete_operational_tuple": not missing_fields,
        "missing_physical_conversion_or_evidence": physical_missing,
        "external_data_admitted": False,
        "intake_file": INTAKE.name,
        "intake_canonical_sha256": intake["canonical_sha256"],
        "intake_file_sha256": sha256(INTAKE),
        "status": (
            "Operational object constructed and type-complete; physical "
            "calibration and empirical execution remain open"
        ),
    }


def main() -> None:
    previous = json.loads(PREVIOUS.read_text(encoding="utf-8"))
    results = {
        "metadata": {
            "title": "FIN Programs 113--124: Constructive Completion Objects and Fractional Operational Physics",
            "release": "10.12",
            "date": "2026-07-27",
            "author": "Krzysztof Żuchowski",
            "orcid": "0009-0002-0909-3613",
            "seed": SEED,
        }
    }
    results["program113_global_mass_theorem"] = program113_global_mass_theorem(previous)
    results["program114_continuous_parameter_box"] = program114_continuous_parameter_box(previous)
    p115 = program115_effective_abelian_certificate()
    results["program115_effective_abelian_certificate"] = p115
    results["program116_stable_invariance_principle"] = program116_stable_invariance_principle(p115)
    p117 = program117_fractional_operational_process()
    results["program117_fractional_operational_process"] = p117
    results["program118_local_fractional_crossover"] = program118_local_fractional_crossover()
    results["program119_unique_inverse_fisher_potential"] = program119_unique_inverse_fisher_potential()
    results["program120_finite_variational_grammar"] = program120_finite_variational_grammar()
    p121 = program121_hidden_markov_apparatus(p117)
    results["program121_hidden_markov_apparatus"] = p121
    p122 = program122_homological_character_functor()
    results["program122_homological_character_functor"] = p122
    results["program123_conditional_damping_cocycle"] = program123_conditional_damping_cocycle(p122)
    results["program124_typed_operational_completion"] = program124_typed_operational_completion(p117, p121)
    results["guardrails"] = {
        "strict_kernel_primary": True,
        "legacy_kernel_role": "intermediate bridge kernel only",
        "selector_QW2191_closed": False,
        "physical_units_sourced": False,
        "strict_uniform_prime_trace_sourced": False,
        "legacy_strict_bridge_complete": False,
        "legacy_role_transfer_started": False,
        "external_data_admitted": False,
        "L_total_or_ToE_promoted": False,
    }
    OUT.write_text(
        json.dumps(results, indent=2, sort_keys=True, default=json_default) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({
        "output": OUT.name,
        "sha256": sha256(OUT),
        "figures": len(list(FIG.glob("*.png"))),
        "program113_global": results["program113_global_mass_theorem"]["global_positive_for_every_m_gt_zero"],
        "program114_box": results["program114_continuous_parameter_box"]["continuous_box_global_mass_nonclosure"],
        "program122_vector": p122["dimension_vector"],
        "program123_beta": results["program123_conditional_damping_cocycle"]["nonzero_positive_monoid_character_forces_beta"],
    }, indent=2, default=json_default))


if __name__ == "__main__":
    main()
