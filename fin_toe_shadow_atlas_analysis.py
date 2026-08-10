#!/usr/bin/env python3
"""Counterfactual FIN-as-ToE atlas: exact transforms and null controls.

The ontological hypothesis H_TOE is never used as a numerical premise.  The
script audits which proposed shadows follow from the frozen strict operator,
which are generic to positive circulant Laplacians, and which require external
sector, scale, algebra or operational data.
"""

from __future__ import annotations

import json
import math
import platform
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy.linalg import expm
from scipy.optimize import minimize_scalar


ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "FIN_TOE_Shadow_Atlas_Results.json"
FIG_DIR = ROOT / "FIN_TOE_Shadow_Atlas_Figures"
SEED = 20260810
N = 12


def native(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [native(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    return value


def cyclic_distance(x: int, y: int, n: int = N) -> int:
    return min((x - y) % n, (y - x) % n)


def strict_profile(d: int) -> float:
    return math.cos(0.18575 * d + 0.16250) / (1.0 + d**1.8)


def radial_matrix(profile) -> np.ndarray:
    w = np.zeros((N, N))
    for x in range(N):
        for y in range(N):
            if x != y:
                w[x, y] = profile(cyclic_distance(x, y))
    return w


def laplacian(w: np.ndarray) -> np.ndarray:
    return np.diag(w.sum(axis=1)) - w


def symmetry_matrices() -> tuple[np.ndarray, np.ndarray]:
    shift = np.zeros((N, N))
    reflection = np.zeros((N, N))
    for x in range(N):
        shift[(x + 1) % N, x] = 1.0
        reflection[(-x) % N, x] = 1.0
    return shift, reflection


def exact_operator_shadow_audit() -> dict:
    w = radial_matrix(strict_profile)
    a = laplacian(w)
    eigenvalues, vectors = np.linalg.eigh(a)
    s = float(w.sum(axis=1)[0])
    t, z, beta = 0.55, 0.70, 0.80
    unitary = vectors @ np.diag(np.exp(-1j * t * eigenvalues)) @ vectors.T
    heat = vectors @ np.diag(np.exp(-t * eigenvalues)) @ vectors.T
    green = vectors @ np.diag(1.0 / (eigenvalues + z)) @ vectors.T
    gibbs_weights = np.exp(-beta * eigenvalues)
    gibbs_weights /= gibbs_weights.sum()
    gibbs = vectors @ np.diag(gibbs_weights) @ vectors.T
    wave = vectors @ np.diag(np.cos(t * np.sqrt(np.maximum(eigenvalues, 0)))) @ vectors.T

    unitary_spectrum = np.linalg.eigvals(unitary)
    unitary_reconstructed = np.sort(-np.angle(unitary_spectrum) / t)
    gibbs_spectrum = np.linalg.eigvalsh(gibbs)[::-1]
    gibbs_reconstructed_gaps = -np.log(gibbs_spectrum / gibbs_spectrum[0]) / beta

    reconstructed = {
        "from_heat": -np.log(np.linalg.eigvalsh(heat)[::-1]) / t,
        "from_green": 1.0 / np.linalg.eigvalsh(green)[::-1] - z,
        "from_unitary_principal_phase": unitary_reconstructed,
        "from_gibbs_gaps": gibbs_reconstructed_gaps,
    }
    # Sort methods whose matrix eigenvalue ordering is reversed by monotonicity.
    reconstructed["from_heat"] = np.sort(reconstructed["from_heat"])
    reconstructed["from_green"] = np.sort(reconstructed["from_green"])
    max_errors = {
        "heat": float(np.max(np.abs(reconstructed["from_heat"] - eigenvalues))),
        "green": float(np.max(np.abs(reconstructed["from_green"] - eigenvalues))),
        "unitary": float(np.max(np.abs(reconstructed["from_unitary_principal_phase"] - eigenvalues))),
        "gibbs_gaps": float(np.max(np.abs(reconstructed["from_gibbs_gaps"] - (eigenvalues - eigenvalues[0])))),
    }
    functions = {"A": a, "unitary": unitary, "heat": heat, "green": green, "wave": wave}
    commutators = {}
    names = list(functions)
    for i, left in enumerate(names):
        for right in names[i + 1 :]:
            commutators[f"{left}:{right}"] = float(np.linalg.norm(functions[left] @ functions[right] - functions[right] @ functions[left]))

    distinct = []
    for value in eigenvalues:
        if not distinct or abs(value - distinct[-1]) > 1e-10:
            distinct.append(float(value))
    shift, reflection = symmetry_matrices()
    symmetry_residuals = {
        name: {
            "translation": float(np.linalg.norm(matrix @ shift - shift @ matrix)),
            "reflection": float(np.linalg.norm(matrix @ reflection - reflection @ matrix)),
        }
        for name, matrix in functions.items()
    }

    mu = np.linalg.eigvalsh(w)
    resonance_energy_order_error = float(np.max(np.abs(np.sort(eigenvalues) - np.sort(s - mu))))
    markov = expm(-t * a)
    return {
        "strict_matrix_definition": "W_xx=0 and W_xy=K_strict(d(x,y)) for x!=y",
        "strict_row_sum_s": s,
        "strict_weights_by_distance": [strict_profile(d) for d in range(1, 7)],
        "eigenvalues": eigenvalues,
        "distinct_eigenvalue_count": len(distinct),
        "functional_algebra_dimension": len(distinct),
        "full_matrix_algebra_dimension": N * N,
        "commutative_dimension_deficit": N * N - len(distinct),
        "all_functional_commutator_norms": commutators,
        "symmetry_commutator_norms": symmetry_residuals,
        "common_spectrum_reconstruction_max_errors": max_errors,
        "unitary_reconstruction_premise": {
            "time": t,
            "t_lambda_max": float(t * eigenvalues[-1]),
            "principal_branch_bound_pi": math.pi,
            "known_zero_mode": True,
        },
        "highest_resonance_equals_lowest_A_order_error": resonance_energy_order_error,
        "markov_checks": {
            "minimum_entry": float(markov.min()),
            "maximum_row_sum_error": float(np.max(np.abs(markov.sum(axis=1) - 1))),
            "maximum_column_sum_error": float(np.max(np.abs(markov.sum(axis=0) - 1))),
            "detailed_balance_uniform_error": float(np.max(np.abs(markov - markov.T))),
        },
        "theorems": [
            "U_t, P_t, G_z, wave and Gibbs weights are functions of one spectral resolution.",
            "C*(A) is commutative and has dimension equal to the number of distinct eigenvalues; it cannot by itself supply a noncommutative observable algebra.",
            "Every f(A) inherits strict translation and reflection symmetry, so no orientation selector lies in functional calculus alone.",
            "For A=sI-W, maximizing an eigenvalue of W is exactly minimizing the corresponding eigenvalue of A.",
        ],
    }


def random_positive_circulant_laplacian(rng: np.random.Generator, target_s: float) -> np.ndarray:
    radial = rng.lognormal(mean=0.0, sigma=1.0, size=6)
    radial *= target_s / (2 * radial[:5].sum() + radial[5])
    return laplacian(radial_matrix(lambda d: radial[d - 1]))


def spectral_metrics(a: np.ndarray) -> dict:
    values = np.linalg.eigvalsh(a)
    positive = values[values > 1e-10]
    normalized = positive / positive.sum()
    return {
        "normalized_gap": float(positive[0] / positive.mean()),
        "condition_number": float(positive[-1] / positive[0]),
        "spectral_entropy": float(-np.sum(normalized * np.log(normalized))),
        "normalized_heat_trace_t1": float(np.sum(np.exp(-positive)) / len(positive)),
        "max_frequency_ratio": float(positive[-1] / positive.mean()),
    }


def null_ensemble_audit(samples: int = 10000) -> dict:
    rng = np.random.default_rng(SEED)
    strict_w = radial_matrix(strict_profile)
    strict_a = laplacian(strict_w)
    target_s = float(strict_w.sum(axis=1)[0])
    fin = spectral_metrics(strict_a)
    rows = []
    functional_shadow_passes = 0
    for _ in range(samples):
        a = random_positive_circulant_laplacian(rng, target_s)
        rows.append(spectral_metrics(a))
        t = 0.31
        p = expm(-t * a)
        if np.min(p) >= -1e-12 and np.max(np.abs(p.sum(axis=1) - 1)) < 1e-10:
            functional_shadow_passes += 1
    percentiles = {}
    distributions = {}
    for key in fin:
        values = np.asarray([row[key] for row in rows])
        distributions[key] = values
        percentiles[key] = float(100.0 * np.mean(values <= fin[key]))
    return {
        "sample_count": samples,
        "strict_metrics": fin,
        "strict_percentiles": percentiles,
        "analytic_unitary_heat_green_property_rate_in_declared_class": 1.0,
        "generic_markov_semigroup_pass_rate": functional_shadow_passes / samples,
        "interpretation": (
            "Functional-calculus transforms and Markov positivity occur throughout the sampled positive "
            "circulant-Laplacian null class (and follow analytically in that class). They establish structural "
            "coherence but are not physical derivations or FIN-specific signatures; percentiles depend on the "
            "declared lognormal measure."
        ),
        "_distributions": distributions,
    }


def legacy_strict_fit_audit() -> dict:
    d = np.arange(1, 7, dtype=float)
    target = np.asarray([strict_profile(int(x)) for x in d])

    def residual_for_log_beta(log_beta: float) -> float:
        beta = math.exp(log_beta)
        shape = np.cos(math.pi * d / 4 + math.pi / 6) / (1 + beta * d)
        amplitude = float(shape @ target / (shape @ shape))
        return float(np.linalg.norm(amplitude * shape - target) / np.linalg.norm(target))

    fit = minimize_scalar(residual_for_log_beta, bounds=(-12, 8), method="bounded", options={"xatol": 1e-14})
    beta = math.exp(float(fit.x))
    shape = np.cos(math.pi * d / 4 + math.pi / 6) / (1 + beta * d)
    amplitude = float(shape @ target / (shape @ shape))
    amplitude_only_shape = np.cos(math.pi * d / 4 + math.pi / 6) / (1 + 0.01 * d)
    amplitude_only = float(amplitude_only_shape @ target / (amplitude_only_shape @ amplitude_only_shape))
    amplitude_only_residual = float(np.linalg.norm(amplitude_only * amplitude_only_shape - target) / np.linalg.norm(target))
    return {
        "distance_domain": d,
        "strict_target": target,
        "best_legacy_star_amplitude": amplitude,
        "best_legacy_star_beta": beta,
        "best_relative_residual": float(fit.fun),
        "canonical_beta_0_01_best_amplitude": amplitude_only,
        "canonical_beta_0_01_relative_residual": amplitude_only_residual,
        "verdict": "No exact canonical-phase linear-damping bridge; a fit is not a typed completion or physical-role transfer theorem.",
    }


def theory_shadow_ledger() -> list[dict]:
    # Scores: mathematics, FIN-specificity, physical bridge readiness,
    # local falsifiability, and overinterpretation risk, each 0--5.
    return [
        {"theory": "spectral graph/operator theory", "relation": "E1 exact mathematics", "level": "L1 physical", "scores": [5, 4, 2, 5, 1]},
        {"theory": "information theory", "relation": "E1 for Shannon/data processing", "level": "L1", "scores": [5, 3, 1, 5, 3]},
        {"theory": "classical Markov diffusion", "relation": "E1 semigroup", "level": "L1", "scores": [5, 3, 2, 5, 2]},
        {"theory": "finite quantum kinematics", "relation": "E2 conditional", "level": "L1", "scores": [5, 2, 2, 5, 3]},
        {"theory": "Green/linear response", "relation": "E1 resolvent", "level": "L1", "scores": [5, 3, 2, 5, 2]},
        {"theory": "finite Gaussian graph model", "relation": "E2 conditional", "level": "L1", "scores": [4, 2, 1, 4, 3]},
        {"theory": "statistical mechanics", "relation": "E2 formal Gibbs", "level": "L1", "scores": [4, 2, 1, 4, 4]},
        {"theory": "open quantum systems", "relation": "E2 completion required", "level": "L1", "scores": [3, 2, 1, 4, 4]},
        {"theory": "DNLS/localized nonlinear waves", "relation": "E2 supplied nonlinear law", "level": "L1", "scores": [4, 3, 2, 4, 3]},
        {"theory": "Gaussian EFT/Schur reduction", "relation": "E1 elimination; E3 full RG", "level": "L1", "scores": [5, 3, 1, 5, 3]},
        {"theory": "classical Hamiltonian mechanics", "relation": "E3 analogy", "level": "L0", "scores": [2, 1, 0, 4, 3]},
        {"theory": "hydrodynamics/integrable solitons", "relation": "E3 analogy", "level": "L0", "scores": [2, 1, 0, 4, 4]},
        {"theory": "noncommutative geometry", "relation": "E3 incomplete spectral triple", "level": "L0-L1", "scores": [3, 2, 1, 4, 4]},
        {"theory": "tensor networks/holography", "relation": "E3 compression analogy", "level": "L0", "scores": [2, 1, 0, 4, 4]},
        {"theory": "causal sets", "relation": "E0 incompatible primitive", "level": "L0", "scores": [1, 0, 0, 5, 4]},
        {"theory": "Yang-Mills gauge theory", "relation": "E0 receiver only", "level": "L0", "scores": [1, 0, 0, 5, 5]},
        {"theory": "general relativity", "relation": "E3 spectral-geometry analogy", "level": "L0", "scores": [1, 0, 0, 4, 5]},
        {"theory": "Standard Model", "relation": "E0 no derived sector", "level": "L0", "scores": [0, 0, 0, 5, 5]},
    ]


def literature_ledger() -> list[dict]:
    return [
        {"area": "information", "work": "C. E. Shannon, A Mathematical Theory of Communication (1948)", "url": "https://doi.org/10.1002/j.1538-7305.1948.tb01338.x"},
        {"area": "information thermodynamics", "work": "R. Landauer, Irreversibility and Heat Generation in the Computing Process (1961)", "url": "https://doi.org/10.1147/rd.53.0183"},
        {"area": "statistical inference", "work": "E. T. Jaynes, Information Theory and Statistical Mechanics (1957)", "url": "https://doi.org/10.1103/PhysRev.106.620"},
        {"area": "open systems", "work": "G. Lindblad, On the Generators of Quantum Dynamical Semigroups (1976)", "url": "https://doi.org/10.1007/BF01608499"},
        {"area": "continuous-time quantum walk", "work": "E. Farhi and S. Gutmann, Quantum Computation and Decision Trees (1998)", "url": "https://arxiv.org/abs/quant-ph/9706062"},
        {"area": "solitons", "work": "N. Zabusky and M. Kruskal, Interaction of 'Solitons' in a Collisionless Plasma and the Recurrence of Initial States (1965)", "url": "https://doi.org/10.1103/PhysRevLett.15.240"},
        {"area": "renormalization", "work": "K. G. Wilson, The Renormalization Group: Critical Phenomena and the Kondo Problem (1975)", "url": "https://doi.org/10.1103/RevModPhys.47.773"},
        {"area": "gauge theory", "work": "C. N. Yang and R. Mills, Conservation of Isotopic Spin and Isotopic Gauge Invariance (1954)", "url": "https://doi.org/10.1103/PhysRev.96.191"},
        {"area": "electroweak theory", "work": "S. Weinberg, A Model of Leptons (1967)", "url": "https://doi.org/10.1103/PhysRevLett.19.1264"},
        {"area": "spectral action", "work": "A. H. Chamseddine and A. Connes, The Spectral Action Principle (1996/1997)", "url": "https://arxiv.org/abs/hep-th/9606001"},
        {"area": "causal sets", "work": "L. Bombelli et al., Space-Time as a Causal Set (1987)", "url": "https://doi.org/10.1103/PhysRevLett.59.521"},
        {"area": "tensor networks", "work": "G. Vidal, Entanglement Renormalization (2007)", "url": "https://doi.org/10.1103/PhysRevLett.99.220405"},
        {"area": "holography", "work": "B. Swingle, Entanglement Renormalization and Holography (2009)", "url": "https://arxiv.org/abs/0905.1317"},
        {"area": "emergent gravity", "work": "T. Jacobson, Thermodynamics of Spacetime (1995)", "url": "https://arxiv.org/abs/gr-qc/9504004"},
        {"area": "quantum mechanics", "work": "E. Schrodinger, Quantisierung als Eigenwertproblem (1926)", "url": "https://doi.org/10.1002/andp.19263840404"},
        {"area": "general relativity", "work": "A. Einstein, Die Grundlage der allgemeinen Relativitatstheorie (1916)", "url": "https://doi.org/10.1002/andp.19163540702"},
        {"area": "Euclidean QFT", "work": "K. Osterwalder and R. Schrader, Axioms for Euclidean Green's Functions (1973)", "url": "https://doi.org/10.1007/BF01645738"},
        {"area": "isospectral geometry", "work": "C. Gordon, D. Webb and S. Wolpert, One Cannot Hear the Shape of a Drum (1992)", "url": "https://doi.org/10.1090/S0273-0979-1992-00289-6"},
    ]


def make_figures(payload: dict) -> None:
    import warnings

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Unable to import Axes3D.*", category=UserWarning)
        import matplotlib.pyplot as plt

    FIG_DIR.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.22})
    eigen = np.asarray(payload["operator_audit"]["eigenvalues"])
    t, z, beta = 0.55, 0.70, 0.80
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.6))
    axes[0, 0].plot(eigen, np.cos(t * eigen), "o-", label="Re U")
    axes[0, 0].plot(eigen, -np.sin(t * eigen), "s-", label="Im U")
    axes[0, 0].set_title("unitary phase")
    axes[0, 0].legend()
    axes[0, 1].semilogy(eigen, np.exp(-t * eigen), "o-")
    axes[0, 1].set_title("heat/Markov decay")
    axes[1, 0].plot(eigen, 1 / (eigen + z), "o-")
    axes[1, 0].set_title("Green response")
    gibbs = np.exp(-beta * eigen); gibbs /= gibbs.sum()
    axes[1, 1].semilogy(eigen, gibbs, "o-")
    axes[1, 1].set_title("formal Gibbs weights")
    for ax in axes.flat:
        ax.set_xlabel("shared eigenvalue lambda")
    fig.tight_layout()
    fig.savefig(FIG_DIR / "common_spectrum_four_shadows.png", dpi=190)
    plt.close(fig)

    null = payload["null_ensemble"]
    distributions = null.pop("_distributions")
    metrics = list(null["strict_metrics"])
    fig, axes = plt.subplots(2, 3, figsize=(10.5, 6.5))
    for ax, key in zip(axes.flat, metrics):
        values = np.asarray(distributions[key])
        ax.hist(values, bins=55, color="#6baed6", alpha=0.85)
        ax.axvline(null["strict_metrics"][key], color="#cb181d", lw=2, label="FIN strict")
        ax.set_title(f"{key}\npercentile={null['strict_percentiles'][key]:.1f}")
    axes.flat[-1].axis("off")
    axes.flat[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "positive_circulant_null_ensemble.png", dpi=190)
    plt.close(fig)

    ledger = payload["theory_shadow_ledger"]
    labels = [row["theory"] for row in ledger]
    math_score = [row["scores"][0] for row in ledger]
    bridge = [row["scores"][2] for row in ledger]
    risk = [row["scores"][4] for row in ledger]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(9.0, 8.0))
    ax.barh(y + 0.22, math_score, height=0.22, label="mathematical support")
    ax.barh(y, bridge, height=0.22, label="physical bridge readiness")
    ax.barh(y - 0.22, risk, height=0.22, label="overinterpretation risk")
    ax.set_yticks(y, labels, fontsize=8)
    ax.set_xlim(0, 5.2)
    ax.invert_yaxis()
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "theory_shadow_maturity_map.png", dpi=190)
    plt.close(fig)


def main() -> None:
    null = null_ensemble_audit()
    payload = {
        "title": "FIN counterfactual ToE shadow atlas",
        "hypothesis": (
            "H_TOE is an organizing counterfactual, not an established conclusion. H_1A is the narrower, "
            "falsifiable ansatz that selected channels are functions of one frozen generator."
        ),
        "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "operator_audit": exact_operator_shadow_audit(),
        "null_ensemble": null,
        "legacy_strict_fit": legacy_strict_fit_audit(),
        "theory_shadow_ledger": theory_shadow_ledger(),
        "literature_ledger": literature_ledger(),
        "epistemic_boundary": (
            "Exact functional transforms are not physical shadows without sector, units, reduction, observable, "
            "instrument and record maps. H_TOE does not discharge QW-2191, source units, "
            "complete legacy-to-strict, transfer legacy roles, construct gauge/SM/GR sectors or provide laboratory data."
        ),
    }
    make_figures(payload)
    RESULTS.write_text(json.dumps(native(payload), indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"results": RESULTS.name, "null_samples": null["sample_count"], "figures": 3}, indent=2))


if __name__ == "__main__":
    main()
