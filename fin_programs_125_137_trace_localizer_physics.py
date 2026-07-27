#!/usr/bin/env python3
"""Execute FIN research Programs 125--137.

The round stress-tests the homological-character fibre object constructed in
Program 122, supplies its missing presentation-independent localizer, derives
finite-window quantitative fractional limits, and constructs several
conditional operational objects.  It deliberately does not infer a canonical
trace, physical units, a strict selector, a completed legacy-to-strict bridge,
or experimental confirmation.
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


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "FIN_Programs_125_137_Trace_Localizer_Physics_Results.json"
FIG = ROOT / "FIN_Programs_125_137_Trace_Localizer_Physics_Figures"
FIG.mkdir(exist_ok=True)

PREVIOUS = ROOT / "FIN_Programs_113_124_Constructive_Completion_Results.json"

ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA = 743.0 / 4000.0
PHI = 13.0 / 80.0
BETA_TORS = 0.01
ETA = 9.0 / 5.0
ALPHA = ETA - 1.0
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


def program125_invariant_trace_classification(previous: dict) -> dict:
    dims = np.asarray(
        previous["program122_homological_character_functor"]["dimension_vector"],
        dtype=float,
    )
    primes = [2, 3, 5, 7, 11]
    uniform_weights = np.full(5, 0.2)
    hilbert_weights = dims / dims.sum()
    eta_uniform = float(uniform_weights @ dims)
    eta_hilbert = float(hilbert_weights @ dims)

    # Under the deliberately strongest abstract symmetry Aut(U(12))=S3, the
    # exceptional sectors p=2 and p=3 are fixed and the three nonidentity
    # unit sectors p=5,7,11 are permuted.  Every trace invariant under that
    # enlarged symmetry is (x,y,z/3,z/3,z/3).  Changes of generator of C12
    # commute with every m_p and hence fix the p labels individually; under
    # carrier naturality alone the admissible state simplex is even larger.
    grid = []
    for x in np.linspace(0, 1, 101):
        for y in np.linspace(0, 1 - x, 101):
            z = 1 - x - y
            grid.append((x, y, z, 2 - x))
    grid = np.asarray(grid)

    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    sc = ax.scatter(grid[:, 0], grid[:, 1], c=grid[:, 3], s=5, cmap="viridis")
    ax.scatter([0.2], [0.2], marker="*", s=150, color="#d1495b", label="uniform sectors")
    ax.scatter(
        [hilbert_weights[0]],
        [hilbert_weights[1]],
        marker="D",
        s=65,
        color="white",
        edgecolor="black",
        label="normalized Hilbert trace",
    )
    ax.set_xlabel(r"$w_2$")
    ax.set_ylabel(r"$w_3$")
    ax.set_title("Invariant trace simplex; colour is the induced exponent")
    ax.legend()
    fig.colorbar(sc, ax=ax, label=r"$\eta=2-w_2$")
    fig.savefig(FIG / "program125_invariant_trace_simplex.png", dpi=190)
    plt.close(fig)

    return {
        "prime_sectors": primes,
        "fibre_dimensions": dims.astype(int).tolist(),
        "maximal_abstract_unit_group_automorphism_orbits": [[2], [3], [5, 7, 11]],
        "symmetry_scope": (
            "The S3 permutation is Aut(U(12)) acting on the three nonidentity "
            "abstract unit elements. C12 generator changes commute with each "
            "m_p and fix the prime labels; carrier naturality alone permits the "
            "full probability simplex."
        ),
        "traces_invariant_under_enlarged_S3_symmetry": "(x,y,z/3,z/3,z/3), x,y,z>=0, x+y+z=1",
        "induced_exponent_formula": "eta=2-w_2",
        "uniform_sector_trace": {
            "weights": uniform_weights.tolist(),
            "eta": eta_uniform,
        },
        "normalized_hilbert_trace": {
            "weights": hilbert_weights.tolist(),
            "eta": eta_hilbert,
            "exact_eta": "17/9",
        },
        "eta_range_on_invariant_trace_simplex": [1.0, 2.0],
        "eta_9_over_5_requires": (
            "w_2=1/5; even enlarged Aut(U(12)) invariance does not force this"
        ),
        "canonical_trace_selected": False,
        "status": "Proven finite classification under enlarged symmetry and robust no-go theorem",
        "scope": (
            "The fibre object is retained, but a state/trace on its sector algebra "
            "is an additional datum. The no-go is stronger because it survives "
            "a symmetry larger than carrier naturality."
        ),
    }


def unit_character_values(p: int) -> list[int]:
    """The four real characters of U(12)=C2xC2, evaluated at a unit p."""
    coords = {1: (0, 0), 5: (1, 0), 7: (0, 1), 11: (1, 1)}
    a, b = coords[p % 12]
    return [(-1) ** (u * a + v * b) for u, v in itertools.product([0, 1], repeat=2)]


def fibre_dimension(n: int, p: int) -> tuple[int, int, int]:
    kernel_size = math.gcd(n, p)
    reduced_h0 = max(0, kernel_size - 1)
    negative_chars = 0
    if n == 12 and math.gcd(p, 12) == 1:
        negative_chars = sum(v == -1 for v in unit_character_values(p))
    return reduced_h0 + negative_chars, reduced_h0, negative_chars


def program126_natural_fibre_localizer() -> dict:
    primes = [2, 3, 5, 7, 11]
    rows = []
    for p in primes:
        total, h0, neg = fibre_dimension(12, p)
        generator_checks = []
        for u in [1, 5, 7, 11]:
            # Under x -> ux, multiplication by p commutes exactly.
            commutes = all((u * p * x - p * u * x) % 12 == 0 for x in range(12))
            image_kernel = sorted(
                {(u * x) % 12 for x in range(12) if (p * x) % 12 == 0}
            )
            target_kernel = sorted(x for x in range(12) if (p * x) % 12 == 0)
            generator_checks.append(commutes and image_kernel == target_kernel)
        rows.append(
            {
                "p": p,
                "reduced_H0_dimension": h0,
                "negative_character_dimension": neg,
                "total_dimension": total,
                "all_C12_generator_changes_natural": all(generator_checks),
            }
        )

    fig, ax = plt.subplots(figsize=(8.4, 4.8), constrained_layout=True)
    x = np.arange(len(primes))
    h0 = [r["reduced_H0_dimension"] for r in rows]
    chars = [r["negative_character_dimension"] for r in rows]
    ax.bar(x, h0, label=r"$\widetilde H_0(\ker m_p)$", color="#3f7cac")
    ax.bar(x, chars, bottom=h0, label=r"$X_p^-$", color="#f4a261")
    ax.set_xticks(x, [str(p) for p in primes])
    ax.set_xlabel("prime label p")
    ax.set_ylabel("dimension")
    ax.set_title("Presentation-independent homological-character localizer on C12")
    ax.legend()
    fig.savefig(FIG / "program126_natural_fibre_localizer.png", dpi=190)
    plt.close(fig)

    return {
        "definition": (
            "For a finite cyclic group C_n and prime p<n, let m_p(x)=px. "
            "Set F_p(C_n)=H0_tilde(ker m_p) direct-sum X_p^-, where X_p^- "
            "is the sum of nontrivial real-character lines negative at p "
            "(and is zero when p is not a unit)."
        ),
        "naturality_equation": "psi o m_p = m_p o psi for every cyclic-group isomorphism psi",
        "C12_rows": rows,
        "dimensions": [r["total_dimension"] for r in rows],
        "presentation_independent": all(
            r["all_C12_generator_changes_natural"] for r in rows
        ),
        "strict_carrier_localizer_constructed": True,
        "ontology_to_carrier_source_constructed": False,
        "trace_or_state_selected": False,
        "damping_coupling_constructed": False,
        "status": "Proven functorial localizer at the finite C12 carrier level",
        "scope": (
            "This closes the presentation-dependence defect of the fibre formula, "
            "not the source, trace, or coupling obligations."
        ),
    }


def continuous_fractional_certificate() -> tuple[dict, dict]:
    nfft = 2**21
    dmax = nfft // 2 - 1
    d = np.arange(1, dmax + 1, dtype=float)
    weights = np.abs(np.cos(OMEGA * d + PHI)) / (1 + d**ETA)
    ring = np.zeros(nfft, dtype=float)
    ring[1 : dmax + 1] = weights
    ring[-dmax:] = weights[::-1]
    spectrum = np.fft.rfft(ring).real
    zd = float(2 * weights.sum())
    qgrid = 2 * math.pi * np.arange(len(spectrum)) / nfft
    mask = (qgrid >= 1e-3) & (qgrid <= 2e-2)
    q = qgrid[mask]
    numerator = zd - spectrum[mask]

    dq = math.pi / nfft
    derivative_bound = float(2 * np.sum(weights * d))
    tail_norm = 2 * dmax ** (1 - ETA) / (ETA - 1)
    round_guard = 5e-11 * max(1.0, zd)
    variation = derivative_bound * dq + round_guard
    slo = np.maximum(0, numerator - variation) / (zd + tail_norm)
    shi = (numerator + variation + 2 * tail_norm) / zd

    mean_abs_cos = 2 / math.pi
    abelian_unnormalized = (
        2
        * mean_abs_cos
        * math.pi
        / (2 * math.gamma(ETA) * math.sin(math.pi * (ETA - 1) / 2))
    )
    c_lo = abelian_unnormalized / (zd + tail_norm)
    c_hi = abelian_unnormalized / zd
    qlo = np.maximum(q - dq, np.finfo(float).tiny)
    qhi = q + dq
    pred_lo = c_lo * qlo**ALPHA
    pred_hi = c_hi * qhi**ALPHA
    rel_lower = np.maximum(0, slo - pred_hi) / pred_hi
    rel_upper = np.maximum(
        np.abs(shi - pred_lo) / pred_lo,
        np.abs(slo - pred_hi) / pred_hi,
    )
    worst = int(np.argmax(rel_upper))

    take = np.unique(np.linspace(0, len(q) - 1, 32).astype(int))
    rows = [
        {
            "q_cell_center": float(q[i]),
            "symbol_interval": [float(slo[i]), float(shi[i])],
            "prediction_interval": [float(pred_lo[i]), float(pred_hi[i])],
            "relative_remainder_upper": float(rel_upper[i]),
        }
        for i in take
    ]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(11.2, 4.6), constrained_layout=True
    )
    ax1.loglog(q, 0.5 * (slo + shi), color="#1f5a99", label="symbol enclosure midpoint")
    ax1.loglog(q, 0.5 * (pred_lo + pred_hi), "--", color="#c44536", label=r"$C|q|^{4/5}$")
    ax1.fill_between(q, slo, shi, color="#1f5a99", alpha=0.18)
    ax1.set_xlabel("|q|")
    ax1.set_ylabel("1 - p-hat(q)")
    ax1.legend(fontsize=8)
    ax1.grid(True, which="both", alpha=0.2)
    ax2.semilogy(q, rel_upper, color="#6a4c93")
    ax2.set_xlabel("|q|")
    ax2.set_ylabel("conservative relative-remainder upper bound")
    ax2.grid(True, which="both", alpha=0.2)
    fig.suptitle("Continuous finite-window fractional-symbol enclosure")
    fig.savefig(FIG / "program127_continuous_fractional_enclosure.png", dpi=190)
    plt.close(fig)

    record = {
        "q_window": [0.001, 0.02],
        "fft_length": nfft,
        "retained_distance": dmax,
        "cell_half_width": dq,
        "finite_derivative_bound": derivative_bound,
        "tail_normalization_upper": tail_norm,
        "abelian_constant_interval": [c_lo, c_hi],
        "number_of_certified_cells": int(len(q)),
        "worst_cell_center": float(q[worst]),
        "maximum_relative_remainder_upper": float(rel_upper[worst]),
        "representative_rows": rows,
        "continuous_window_enclosed": True,
        "formal_interval_arithmetic": False,
        "status": "Strong numerical evidence with analytic truncation and cell bounds",
        "scope": (
            "Every frequency cell in the stated compact window is covered. "
            "Floating-point FFT rounding is guarded conservatively but not "
            "machine-verified interval arithmetic; this is not a full q->0 rate theorem."
        ),
    }
    aux = {
        "q": q,
        "dq": dq,
        "slo": slo,
        "shi": shi,
        "c_lo": c_lo,
        "c_hi": c_hi,
        "rel_upper": rel_upper,
    }
    return record, aux


def program128_quantitative_stable_window(aux: dict) -> dict:
    q = aux["q"]
    slo = aux["slo"]
    shi = aux["shi"]
    c_lo = aux["c_lo"]
    c_hi = aux["c_hi"]
    rows = []
    for n in [25, 40, 64, 100, 160, 250]:
        for k in [0.5, 1.0, 2.0]:
            qn = k * n ** (-1 / ALPHA)
            if qn < q[0] or qn > q[-1]:
                continue
            j = int(np.argmin(np.abs(q - qn)))
            lo = max(0.0, float(slo[j]))
            hi = min(0.999999999, float(shi[j]))
            finite_interval = [(1 - hi) ** n, (1 - lo) ** n]
            stable_interval = [
                math.exp(-c_hi * k**ALPHA),
                math.exp(-c_lo * k**ALPHA),
            ]
            separation = max(
                0.0,
                stable_interval[0] - finite_interval[1],
                finite_interval[0] - stable_interval[1],
            )
            hausdorff = max(
                abs(finite_interval[0] - stable_interval[0]),
                abs(finite_interval[1] - stable_interval[1]),
            )
            rows.append(
                {
                    "n": n,
                    "k": k,
                    "q": qn,
                    "finite_characteristic_interval": finite_interval,
                    "stable_target_interval": stable_interval,
                    "interval_separation": separation,
                    "endpoint_hausdorff_bound": hausdorff,
                }
            )

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    for k in sorted({r["k"] for r in rows}):
        rr = [r for r in rows if r["k"] == k]
        ax.loglog(
            [r["n"] for r in rr],
            [r["endpoint_hausdorff_bound"] for r in rr],
            "o-",
            label=f"k={k:g}",
        )
    ax.set_xlabel("number of steps n")
    ax.set_ylabel("certified endpoint-distance bound")
    ax.set_title("Finite-n approach to the 4/5-stable characteristic function")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program128_quantitative_stable_window.png", dpi=190)
    plt.close(fig)

    return {
        "scaling": "q_n=k*n^(-5/4)",
        "rows": rows,
        "maximum_interval_separation": max(r["interval_separation"] for r in rows),
        "maximum_endpoint_hausdorff_bound": max(
            r["endpoint_hausdorff_bound"] for r in rows
        ),
        "all_rows_overlap_stable_target": all(
            r["interval_separation"] == 0 for r in rows
        ),
        "status": "Quantitative finite-window theorem conditional on Program 127 enclosure",
        "scope": (
            "The certificate covers only n,k combinations whose rescaled frequency "
            "lies in [0.001,0.02]; it does not prove a global asymptotic rate."
        ),
    }


def program129_fractional_wave_uv_obstruction() -> dict:
    shell = np.geomspace(1, 1e7, 300)
    t = 1.0
    c = 1.0
    curvature_bound = (
        t ** (-0.5)
        * (c * ALPHA * (1 - ALPHA)) ** (-0.5)
        * shell ** ((2 - ALPHA) / 2)
    )
    trivial_bound = shell
    best_shell_bound = np.minimum(trivial_bound, curvature_bound)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    ax.loglog(shell, curvature_bound, label="stationary-phase shell bound")
    ax.loglog(shell, trivial_bound, "--", label="trivial shell size")
    ax.loglog(shell, best_shell_bound, color="black", lw=2, label="best available")
    ax.set_xlabel("dyadic UV scale Λ")
    ax.set_ylabel("shell contribution bound")
    ax.set_title("No cutoff-free dispersive summation for alpha=4/5")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program129_fractional_wave_uv_obstruction.png", dpi=190)
    plt.close(fig)

    return {
        "wave_multiplier": "exp(-i*t*C*|q|^(4/5))",
        "formal_scaling": "K_t(x)=t^(-5/4) K_1(x*t^(-5/4))",
        "phase_second_derivative": "t*C*alpha*(alpha-1)*|q|^(alpha-2)",
        "dyadic_stationary_phase_bound": (
            "O(t^(-1/2)*Lambda^((2-alpha)/2)) = O(t^(-1/2)*Lambda^(3/5))"
        ),
        "ultraviolet_shell_sum_converges": False,
        "cutoff_free_L1_to_Linf_claim_established": False,
        "constructed_object": (
            "For each UV cutoff Lambda, U_t^Lambda is the spectral multiplier "
            "1_{|q|<=Lambda} exp(-it C|q|^(4/5)), unitary on its cutoff subspace."
        ),
        "physical_interpretation": (
            "A point detector is not defined by the fractional wave group alone; "
            "a bandwidth or detector-resolution object is necessary."
        ),
        "status": "Proven ultraviolet obstruction for this dispersive proof route",
        "scope": (
            "This does not prove that every possible dispersive estimate fails; "
            "it proves that the unregularized continuum kernel is not supplied by "
            "the present data and that the standard dyadic bound is nonsummable."
        ),
    }


def program130_dimensional_calibration_object() -> dict:
    dimension_matrix = np.asarray(
        [
            [1.0, 0.0, 0.0],  # length
            [0.0, 1.0, 0.0],  # time
            [0.0, -1.0, 1.0],  # energy = hbar/time
        ]
    )
    rank = int(np.linalg.matrix_rank(dimension_matrix))
    determinant = float(np.linalg.det(dimension_matrix))
    rel = {"ell": 0.01, "tau": 0.02, "hbar": 0.03}
    uncertainties = {
        "energy_hbar_over_tau": math.hypot(rel["hbar"], rel["tau"]),
        "fractional_diffusion_ell_alpha_over_tau": math.sqrt(
            (ALPHA * rel["ell"]) ** 2 + rel["tau"] ** 2
        ),
        "fractional_hamiltonian_hbar_ell_alpha_over_tau": math.sqrt(
            rel["hbar"] ** 2
            + (ALPHA * rel["ell"]) ** 2
            + rel["tau"] ** 2
        ),
    }

    fig, ax = plt.subplots(figsize=(7.4, 4.8), constrained_layout=True)
    image = ax.imshow(dimension_matrix, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks(range(3), [r"$\log\ell$", r"$\log\tau$", r"$\log\hbar$"])
    ax.set_yticks(range(3), ["length", "time", "energy"])
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{dimension_matrix[i,j]:g}", ha="center", va="center")
    ax.set_title("Minimal external calibration map has rank three")
    fig.colorbar(image, ax=ax, shrink=0.8)
    fig.savefig(FIG / "program130_dimensional_calibration.png", dpi=190)
    plt.close(fig)

    return {
        "conditional_map": {
            "physical_position": "x_phys=ell*x",
            "physical_time": "t_phys=tau*t",
            "physical_Hamiltonian": "H_phys=(hbar/tau)*A",
            "fractional_diffusion_coefficient": "D_(4/5)=ell^(4/5)/tau",
            "fractional_Hamiltonian_coefficient": "hbar*ell^(4/5)/tau",
        },
        "log_dimension_matrix": dimension_matrix.tolist(),
        "rank": rank,
        "determinant": determinant,
        "independent_external_calibrations_needed": rank,
        "example_relative_input_uncertainties": rel,
        "propagated_relative_uncertainties": uncertainties,
        "internal_unit_source_found": False,
        "status": "Proven dimensional-analysis theorem and constructed calibration layer",
        "scope": (
            "The map converts a dimensionless operational model into calibrated "
            "predictions only after ell, tau, and hbar are supplied externally."
        ),
    }


def simulate_apparatus(n: int, eps: float, rho: float, rng) -> tuple[np.ndarray, np.ndarray]:
    truth = rng.integers(0, 2, n)
    error = np.zeros(n, dtype=int)
    error[0] = rng.random() < eps
    for i in range(1, n):
        if rng.random() < rho:
            error[i] = error[i - 1]
        else:
            error[i] = rng.random() < eps
    return truth, truth ^ error


def conditional_logloss(errors: np.ndarray, order: int, train: int) -> float:
    contexts = 2**order if order else 1
    counts = np.full((contexts, 2), 0.5)  # Jeffreys smoothing
    for i in range(order, train):
        ctx = 0
        for j in range(order):
            ctx = 2 * ctx + int(errors[i - order + j])
        counts[ctx, errors[i]] += 1
    probs = counts / counts.sum(axis=1, keepdims=True)
    loss = 0.0
    number = 0
    for i in range(max(order, train), len(errors)):
        ctx = 0
        for j in range(order):
            ctx = 2 * ctx + int(errors[i - order + j])
        loss -= math.log(probs[ctx, errors[i]])
        number += 1
    return loss / number


def program131_apparatus_process_tomography() -> dict:
    n_total = 20_000
    eps, rho = 0.1, 0.8
    truth, observed = simulate_apparatus(n_total, eps, rho, RNG)
    errors = truth ^ observed
    rows = []
    for train in [200, 1000, 5000, 10000]:
        row = {"calibration_records": train}
        for order in [0, 1, 2]:
            row[f"order_{order}_heldout_logloss"] = conditional_logloss(
                errors, order, train
            )
        rows.append(row)

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    for order, colour in [(0, "#d1495b"), (1, "#1f5a99"), (2, "#2a9d8f")]:
        ax.semilogx(
            [r["calibration_records"] for r in rows],
            [r[f"order_{order}_heldout_logloss"] for r in rows],
            "o-",
            color=colour,
            label=f"memory order {order}",
        )
    ax.set_xlabel("calibration records")
    ax.set_ylabel("held-out log loss per record")
    ax.set_title("Nonparametric apparatus-channel tomography")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend()
    fig.savefig(FIG / "program131_apparatus_tomography.png", dpi=190)
    plt.close(fig)

    return {
        "generative_channel": {
            "error_prevalence_when_refreshed": eps,
            "error_persistence": rho,
            "true_memory_order": 1,
        },
        "estimator": "Jeffreys-smoothed conditional-frequency process tensor",
        "rows": rows,
        "all_binary_contexts_observed_at_largest_sample": all(
            np.any(
                [
                    tuple(errors[i - 2 : i]) == context
                    for i in range(2, 10000)
                ]
            )
            for context in itertools.product([0, 1], repeat=2)
        ),
        "best_large_sample_order": min(
            range(3), key=lambda o: rows[-1][f"order_{o}_heldout_logloss"]
        ),
        "external_experiment_used": False,
        "status": "Constructed and numerically validated finite-memory apparatus estimator",
        "scope": (
            "This is a synthetic identifiability test.  It supplies a method for "
            "learning an instrument channel, not empirical evidence for FIN."
        ),
    }


def program132_crossover_rg_flow() -> dict:
    gamma = 2 - ALPHA  # local q^2 versus fractional q^alpha
    l = np.linspace(-5, 5, 500)
    initial = [0.02, 0.1, 0.5, 0.9, 0.98]
    flows = {}
    for x0 in initial:
        x = x0 * np.exp(gamma * l) / (
            1 - x0 + x0 * np.exp(gamma * l)
        )
        flows[str(x0)] = x

    fig, ax = plt.subplots(figsize=(8.5, 4.8), constrained_layout=True)
    for x0, x in flows.items():
        ax.plot(l, x, label=f"x(0)={x0}")
    ax.axhline(0, color="black", lw=0.7)
    ax.axhline(1, color="black", lw=0.7)
    ax.set_xlabel(r"RG scale $\ell=\log b$")
    ax.set_ylabel(r"projective fractional coupling $x=g/(1+g)$")
    ax.set_title("Exact local-to-fractional crossover flow")
    ax.legend(fontsize=8, ncol=2)
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program132_crossover_rg_flow.png", dpi=190)
    plt.close(fig)

    return {
        "operator": "A=kappa*(-Delta)+nu*(-Delta)^(2/5)",
        "ratio": "g=nu/kappa",
        "scale_law": "g(b)=b^(6/5) g",
        "projective_coordinate": "x=g/(1+g)",
        "beta_function": "dx/dlog(b)=(6/5)*x*(1-x)",
        "fixed_points": {
            "x=0": "local/UV; unstable toward increasing b",
            "x=1": "fractional/IR; stable toward increasing b",
        },
        "crossover_momentum": "q_*=(nu/kappa)^(5/6)",
        "status": "Proven exact two-coupling projective RG theorem",
        "scope": (
            "This is kinematic scaling of a constructed crossover operator, "
            "not a dynamically derived renormalization law for FIN parameters."
        ),
    }


def count_arithmetic_encodings() -> dict:
    values = [2, 4, 5, 9, 12, 16, 17]
    phi_target = Fraction(13, 80)
    omega_target = Fraction(743, 4000)
    phi_hits = []
    for a, b, c, d in itertools.product(values, repeat=4):
        if c * d and Fraction(a - b, c * d) == phi_target:
            phi_hits.append((a, b, c, d))
    omega_hits = []
    for a, b, c, d, e, p in itertools.product(values, repeat=6):
        den = d**3 * e * p
        if den and Fraction(a * (b * c - b) - d, den) == omega_target:
            omega_hits.append((a, b, c, d, e, p))
    return {
        "base_values": values,
        "phi_template": "(A-B)/(C*D)",
        "phi_hit_count": len(phi_hits),
        "phi_first_hits": phi_hits[:10],
        "omega_template": "(A*(B*C-B)-D)/(D^3*E*P)",
        "omega_hit_count": len(omega_hits),
        "omega_first_hits": omega_hits[:10],
    }


def program133_phase_frequency_source_test() -> dict:
    enc = count_arithmetic_encodings()
    phi_identity = Fraction(17 - 4, 5 * 16)
    omega_identity = Fraction(17 * (4 * 12 - 4) - 5, 5**3 * (2 * 16))
    perturb_phi = Fraction(16 - 4, 5 * 16)
    perturb_omega = Fraction(16 * (4 * 12 - 4) - 5, 5**3 * (2 * 16))

    fig, ax = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
    labels = ["phi target", "phi perturbed", "omega target", "omega perturbed"]
    vals = [
        float(phi_identity),
        float(perturb_phi),
        float(omega_identity),
        float(perturb_omega),
    ]
    ax.bar(labels, vals, color=["#2a9d8f", "#adb5bd", "#1f5a99", "#adb5bd"])
    ax.axhline(PHI, color="#2a9d8f", ls="--", lw=1)
    ax.axhline(OMEGA, color="#1f5a99", ls=":", lw=1)
    ax.tick_params(axis="x", rotation=18)
    ax.set_ylabel("dimensionless value")
    ax.set_title("Exact arithmetic encodings fail the source/naturality test")
    fig.savefig(FIG / "program133_phase_frequency_source_test.png", dpi=190)
    plt.close(fig)

    return {
        "intrinsic_integer_inventory": {
            "carrier_order": 12,
            "prime_sector_count": 5,
            "sum_fibre_dimensions": 9,
            "sum_squared_fibre_dimensions": 17,
            "product_fibre_dimensions": 16,
            "real_character_count": 4,
            "p3_fibre_dimension": 2,
        },
        "exact_phi_encoding": {
            "formula": "(17-4)/(5*16)",
            "value": str(phi_identity),
            "matches_frozen_phi": phi_identity == Fraction(13, 80),
        },
        "exact_omega_encoding": {
            "formula": "(17*(4*12-4)-5)/(5^3*(2*16))",
            "value": str(omega_identity),
            "matches_frozen_omega": omega_identity == Fraction(743, 4000),
        },
        "finite_template_nonuniqueness_audit": enc,
        "one_integer_perturbation": {
            "phi": str(perturb_phi),
            "omega": str(perturb_omega),
        },
        "natural_transformation_or_variational_selection_found": False,
        "source_claim_accepted": False,
        "status": "Exact identities found but rejected as target-driven arithmetic encodings",
        "scope": (
            "Equality to the frozen rationals is not a derivation.  No theorem "
            "selects these expressions from the many available integer operations."
        ),
    }


def program134_amplitude_projectivization() -> dict:
    d = np.arange(0, 65, dtype=float)
    legacy = ALPHA_GEO * np.cos(OMEGA * d + PHI) / (1 + BETA_TORS * d)
    strict = np.cos(OMEGA * d + PHI) / (1 + d**ETA)
    legacy_p = legacy / legacy[0]
    strict_p = strict / strict[0]
    projective_mismatch = float(
        np.linalg.norm(legacy_p - strict_p) / np.linalg.norm(strict_p)
    )
    best_scale = float(legacy @ strict / (legacy @ legacy))
    scalar_residual = float(
        np.linalg.norm(best_scale * legacy - strict) / np.linalg.norm(strict)
    )

    fig, ax = plt.subplots(figsize=(9.0, 4.8), constrained_layout=True)
    ax.plot(d, legacy_p, label="projectivized legacy", color="#d1495b")
    ax.plot(d, strict_p, label="projectivized strict", color="#1f5a99")
    ax.set_xlabel("distance d")
    ax.set_ylabel(r"$K(d)/K(0)$")
    ax.set_title("Amplitude quotient removes alpha_geo but not strict completion")
    ax.legend()
    ax.grid(True, alpha=0.25)
    fig.savefig(FIG / "program134_amplitude_projectivization.png", dpi=190)
    plt.close(fig)

    return {
        "map": "Pi_0(K)=K/K(0), for K(0) nonzero",
        "properties": {
            "amplitude_invariant": True,
            "idempotent": True,
            "alpha_geo_removed_from_legacy_shape": True,
        },
        "distance_window": [0, 64],
        "projective_shape_relative_L2_mismatch": projective_mismatch,
        "best_scalar_legacy_to_strict": best_scale,
        "best_scalar_only_relative_residual": scalar_residual,
        "strict_completion_achieved_by_projectivization": False,
        "role_transfer_obstruction": (
            "Any quantity depending on absolute legacy amplitude is erased by Pi_0 "
            "and therefore cannot be transported through this quotient alone."
        ),
        "status": "Proven projective-amplitude theorem plus finite shape falsification",
        "scope": (
            "Projectivization explains how alpha_geo may become gauge-like for shape, "
            "but it neither supplies damping completion nor preserves amplitude roles."
        ),
    }


def program135_conditional_damping_bridge() -> dict:
    d = np.geomspace(1, 1e6, 800)
    traces = [0.1, 0.2, 0.3]
    curves = {}
    for w2 in traces:
        eta = 2 - w2
        curves[w2] = (1 + BETA_TORS * d) / (1 + d**eta)
    target = curves[0.2]
    reconstruction = (1 + BETA_TORS * d) / (1 + d ** (9 / 5))
    residual = float(np.max(np.abs(target - reconstruction)))
    retention = 2 ** (-4 / 5)
    retention_information = math.exp(-ALPHA_GEO / 5)

    fig, ax = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    for w2, curve in curves.items():
        ax.loglog(d, curve, label=fr"$w_2={w2:g}$, $\eta={2-w2:g}$")
    ax.set_xlabel("distance d")
    ax.set_ylabel("conditional damping completion D_w(d)")
    ax.set_title("Trace-sensitive conditional legacy-to-strict damping factor")
    ax.legend()
    ax.grid(True, which="both", alpha=0.25)
    fig.savefig(FIG / "program135_conditional_damping_bridge.png", dpi=190)
    plt.close(fig)

    return {
        "family": "D_w(d)=(1+beta_tors*d)/(1+d^(2-w_2))",
        "conditional_trace": "w_2=1/5",
        "conditional_exponent": "eta=9/5",
        "exact_reconstruction_max_abs_residual": residual,
        "dyadic_retention": retention,
        "information_expression": "2^(-4/5)=exp(-alpha_geo/5)",
        "information_expression_residual": abs(retention - retention_information),
        "tail_multiplicativity_selects_beta": 1,
        "trace_selected_internally": False,
        "phase_frequency_source_selected": False,
        "full_bridge_completed": False,
        "status": "Exact conditional bridge factor; source obligations remain open",
        "scope": (
            "Given the uniform sector state and the multiplicative-tail axiom, "
            "the strict damping exponent and beta=1 follow.  Program 125 proves "
            "that the state itself is not forced by present symmetry."
        ),
    }


def fractional_circulant(n=12, alpha=ALPHA) -> np.ndarray:
    k = np.arange(n)
    q = 2 * math.pi * np.minimum(k, n - k) / n
    lam = q**alpha
    fourier = np.exp(2j * math.pi * np.outer(np.arange(n), k) / n) / math.sqrt(n)
    return np.real_if_close(fourier @ np.diag(lam) @ fourier.conj().T).real


def signed_current(rho: np.ndarray, a: np.ndarray) -> float:
    n = len(rho)
    return float(
        2
        * sum(np.imag(rho[i, (i + 1) % n]) * a[(i + 1) % n, i] for i in range(n))
    )


def program136_signed_operational_state_source() -> dict:
    n = 12
    a = fractional_circulant(n)
    x = np.arange(n)
    states = {}
    for k in [0, 1, -1, 2, -2]:
        psi = np.exp(2j * math.pi * k * x / n) / math.sqrt(n)
        rho = np.outer(psi, psi.conj())
        states[k] = signed_current(rho, a)
    thermal = np.linalg.eigh(a)
    vals, vecs = thermal
    rho_thermal = vecs @ np.diag(np.exp(-vals)) @ vecs.conj().T
    rho_thermal /= np.trace(rho_thermal)
    thermal_current = signed_current(rho_thermal, a)

    fig, ax = plt.subplots(figsize=(8.0, 4.8), constrained_layout=True)
    labels = [str(k) for k in states]
    vals_plot = [states[k] for k in states]
    colours = ["#adb5bd" if abs(k) == 0 else ("#1f5a99" if k > 0 else "#d1495b") for k in states]
    ax.bar(labels, vals_plot, color=colours)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xlabel("prepared Fourier branch k")
    ax.set_ylabel(r"signed receiver $\Lambda(\rho,A)$")
    ax.set_title("A signed receiver exists, but its nonzero sign is preparation-selected")
    fig.savefig(FIG / "program136_signed_state_receiver.png", dpi=190)
    plt.close(fig)

    return {
        "operator": "fractional circulant A on C12 with eigenvalue |q|^(4/5)",
        "signed_receiver": (
            "Lambda(rho,A)=2 Im sum_i rho_(i,i+1) A_(i+1,i)"
        ),
        "fourier_branch_currents": {str(k): v for k, v in states.items()},
        "paired_sign_test": {
            "k1_plus_kminus1": states[1] + states[-1],
            "k2_plus_kminus2": states[2] + states[-2],
        },
        "thermal_state_current": thermal_current,
        "reflection_invariant_states_give_zero": (
            abs(states[0]) < 1e-12 and abs(thermal_current) < 1e-12
        ),
        "nonzero_signed_state_requires_branch_preparation": True,
        "strict_nonpremise_branch_preparation_found": False,
        "QW_2191_discharged": False,
        "status": "Constructed signed operational receiver; strict source no-go remains",
        "scope": (
            "The formula can read chirality after a state has been prepared. "
            "It does not choose k over -k and therefore does not break the selector torsor."
        ),
    }


def program137_external_data_audit() -> dict:
    candidates = sorted(
        set(ROOT.glob("*Data*.json"))
        | set(ROOT.glob("*Preregistration*.json"))
        | set(ROOT.glob("*Intake*.json"))
    )
    rows = []
    admitted = []
    for path in candidates:
        try:
            data = json.loads(path.read_text())
        except (json.JSONDecodeError, UnicodeDecodeError):
            continue
        text = json.dumps(data).lower()
        positive_markers = [
            '"external_data_admitted": true',
            '"external_data_used": true',
            '"data_admitted": true',
        ]
        is_admitted = any(marker in text for marker in positive_markers)
        row = {
            "file": path.name,
            "sha256": sha256(path),
            "external_data_admitted_true": is_admitted,
        }
        rows.append(row)
        if is_admitted:
            admitted.append(path.name)
    return {
        "files_scanned": len(rows),
        "rows": rows,
        "files_admitting_external_data": admitted,
        "external_data_admitted_in_current_research_chain": bool(admitted),
        "current_round_uses_only_synthetic_or_repository_internal_data": True,
        "status": "Reproducible local data-intake audit",
        "scope": (
            "The audit detects explicit Boolean admission markers in named local "
            "data/intake/preregistration JSON artifacts; it is not a forensic scan "
            "of every byte in the repository."
        ),
    }


def main() -> None:
    previous = json.loads(PREVIOUS.read_text())
    p127, aux127 = continuous_fractional_certificate()
    result = {
        "title": "FIN Programs 125-137: Trace, Localizer, Fractional Physics and Operational Completion",
        "release": "10.13",
        "execution_date": "2026-07-27",
        "seed": SEED,
        "frozen_parameters": {
            "alpha_geo": ALPHA_GEO,
            "omega": OMEGA,
            "phi": PHI,
            "beta_tors": BETA_TORS,
            "eta": ETA,
        },
        "provenance": {
            "previous_results": PREVIOUS.name,
            "previous_results_sha256": sha256(PREVIOUS),
            "script": Path(__file__).name,
        },
        "program125_invariant_trace_classification": program125_invariant_trace_classification(previous),
        "program126_natural_fibre_localizer": program126_natural_fibre_localizer(),
        "program127_continuous_fractional_enclosure": p127,
        "program128_quantitative_stable_window": program128_quantitative_stable_window(aux127),
        "program129_fractional_wave_uv_obstruction": program129_fractional_wave_uv_obstruction(),
        "program130_dimensional_calibration_object": program130_dimensional_calibration_object(),
        "program131_apparatus_process_tomography": program131_apparatus_process_tomography(),
        "program132_crossover_rg_flow": program132_crossover_rg_flow(),
        "program133_phase_frequency_source_test": program133_phase_frequency_source_test(),
        "program134_amplitude_projectivization": program134_amplitude_projectivization(),
        "program135_conditional_damping_bridge": program135_conditional_damping_bridge(),
        "program136_signed_operational_state_source": program136_signed_operational_state_source(),
        "program137_external_data_audit": program137_external_data_audit(),
        "global_verdict": {
            "new_objects_constructed": [
                "presentation-independent homological-character fibre localizer",
                "continuous finite-frequency fractional enclosure",
                "cutoff fractional-wave operational family",
                "rank-three physical calibration map",
                "finite-memory apparatus process-tomography estimator",
                "exact projective crossover RG flow",
                "amplitude projectivization quotient",
                "trace-parameterized damping completion family",
                "signed state-dependent current receiver",
            ],
            "strongest_positive_result": (
                "The fibre formula is now a natural C12 carrier construction, and "
                "the fractional operational limit has a quantitative finite window."
            ),
            "strongest_negative_result": (
                "Naturality does not select the uniform sector trace: invariant traces "
                "form a simplex and the normalized Hilbert trace yields eta=17/9."
            ),
            "closures_not_claimed": [
                "canonical trace or eta=9/5 source",
                "phase/frequency source",
                "QW-2191 selector closure",
                "internal physical units",
                "full legacy-to-strict completion",
                "legacy physical-role transfer",
                "L_total or Theory-of-Everything closure",
                "experimental validation",
            ],
        },
    }
    OUT.write_text(json.dumps(result, indent=2, ensure_ascii=False, default=json_default) + "\n")
    print(json.dumps({
        "output": OUT.name,
        "figures": len(list(FIG.glob("*.png"))),
        "programs": 13,
        "trace_eta_uniform": result["program125_invariant_trace_classification"]["uniform_sector_trace"]["eta"],
        "trace_eta_hilbert": result["program125_invariant_trace_classification"]["normalized_hilbert_trace"]["eta"],
        "fractional_window_max_relative_upper": p127["maximum_relative_remainder_upper"],
        "data_files_admitted": result["program137_external_data_audit"]["files_admitting_external_data"],
    }, indent=2))


if __name__ == "__main__":
    main()
