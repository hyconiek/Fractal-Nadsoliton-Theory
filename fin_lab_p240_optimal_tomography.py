#!/usr/bin/env python3
"""FIN Laboratory Program 240: optimal heat-kernel tomography design.

This module is an executable specification, not evidence from a laboratory.
It freezes the strict twelve-node target, proves/implements the exact
dimensionless time and preparation-allocation design, evaluates a
distribution-free rectangular matrix-Bernstein bound, and performs a
deterministic finite-shot planning simulation.

Claim boundary
--------------
The code optimizes a conditional operational realization of
P_tau = exp(-tau A).  It does not generate an apparatus, a physical clock,
an independent record, a strict selector, or a physical law selecting the
heat rather than the unitary temporal ray.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import argparse
import hashlib
import json
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Lab_P240_242_Transfer_Figures"
RESULTS_PATH = ROOT / "FIN_Lab_P240_Optimal_Tomography_Results.json"
LOCK_PATH = ROOT / "FIN_Lab_P240_Design_Lock.json"

N = 12
OMEGA_STRICT = 0.18575
PHI_STRICT = 0.16250
BETA_STRICT = 1.0
ETA_STRICT = 1.8
ROW_SUM_REFERENCE = 1.660307278766099
SEED = 202607240


@dataclass(frozen=True)
class StrictTarget:
    n: int
    omega: float
    phi: float
    beta: float
    eta: float
    row_sum: float
    lambda_max: float
    spectral_gap: float
    tau_dimensionless: float


def canonical_digest(value: object) -> str:
    payload = json.dumps(
        value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_canonical(path: Path, core: dict) -> None:
    record = {"core": core, "canonical_core_sha256": canonical_digest(core)}
    path.write_text(
        json.dumps(record, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def strict_operator(n: int = N) -> tuple[np.ndarray, np.ndarray, float]:
    """Return (A, W, s) for the frozen strict cyclic-distance convention."""
    if n != 12:
        raise ValueError(
            "Program 240 freezes the certified C12 target only; "
            "naive fixed-formula refinement is not licensed."
        )
    w = np.zeros((n, n), dtype=float)
    for x in range(n):
        for y in range(n):
            if x == y:
                continue
            d = min((x - y) % n, (y - x) % n)
            w[x, y] = math.cos(OMEGA_STRICT * d + PHI_STRICT) / (
                1.0 + BETA_STRICT * d**ETA_STRICT
            )
    row_sums = w.sum(axis=1)
    if not np.allclose(row_sums, row_sums[0], atol=2e-15):
        raise ArithmeticError("strict matrix lost its constant-row-sum property")
    s = float(row_sums[0])
    a = s * np.eye(n) - w
    return a, w, s


def projective_signature_from_generator(a: np.ndarray) -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh((a + a.T) / 2.0)
    scale = float(eigenvalues[-1])
    if eigenvalues[0] < -1e-10 * max(1.0, scale) or scale <= 0:
        raise ValueError("generator is not positive semidefinite")
    positive = np.sort(eigenvalues[1:])
    return positive / positive[-1]


def symmetric_transition_estimate(counts: np.ndarray) -> np.ndarray:
    """Return the symmetric estimate used only for spectral reconstruction."""
    counts = np.asarray(counts, dtype=float)
    totals = counts.sum(axis=0)
    if np.any(totals <= 0):
        raise ValueError("every basis preparation needs at least one outcome")
    raw = counts / totals[None, :]
    return (raw + raw.T) / 2.0


def signature_from_counts(
    counts: np.ndarray, duration: float
) -> tuple[np.ndarray | None, dict]:
    estimate = symmetric_transition_estimate(counts)
    eigenvalues = np.linalg.eigvalsh(estimate)[::-1]
    diagnostics = {
        "minimum_transition_eigenvalue": float(eigenvalues[-1]),
        "maximum_transition_eigenvalue": float(eigenvalues[0]),
    }
    if eigenvalues[-1] <= 0 or eigenvalues[0] <= 0:
        return None, diagnostics
    relative = np.clip(eigenvalues / eigenvalues[0], 1e-15, 1.0)
    rates = -np.log(relative) / duration
    positive = np.sort(rates[1:])
    if positive[-1] <= 0:
        return None, diagnostics
    return positive / positive[-1], diagnostics


def sample_counts(
    transition: np.ndarray, shots_per_preparation: int, rng: np.random.Generator
) -> np.ndarray:
    counts = np.zeros_like(transition, dtype=np.int64)
    for column in range(transition.shape[1]):
        probabilities = np.clip(transition[:, column], 0.0, 1.0)
        probabilities /= probabilities.sum()
        counts[:, column] = rng.multinomial(shots_per_preparation, probabilities)
    return counts


def bernstein_operator_radius(
    n: int, shots_per_preparation: int, alpha: float
) -> float:
    """Distribution-free rectangular matrix-Bernstein operator radius.

    For independent one-hot outcomes in every column,

      P(||P_hat-P||_2 > eps) <= alpha,

    with variance proxy v <= n/(2m), increment bound R <= sqrt(2)/m,
    and dilation dimension 2n.
    """
    if shots_per_preparation <= 0:
        raise ValueError("shots_per_preparation must be positive")
    if not 0.0 < alpha < 1.0:
        raise ValueError("alpha must lie in (0,1)")
    m = float(shots_per_preparation)
    log_term = math.log(2.0 * n / alpha)
    variance_proxy = n / (2.0 * m)
    increment_bound = math.sqrt(2.0) / m
    return math.sqrt(2.0 * variance_proxy * log_term) + (
        2.0 * increment_bound * log_term / 3.0
    )


def log_generator_error_bound(
    transition_minimum_eigenvalue: float,
    transition_operator_radius: float,
    duration: float,
) -> float:
    """Bound max eigenrate error after normalization to the top eigenvalue.

    The factor two covers the perturbation of both mu_k and the empirical
    Perron eigenvalue used in log(mu_k/mu_0).
    """
    epsilon = transition_operator_radius
    mu_min = transition_minimum_eigenvalue
    if epsilon >= mu_min:
        return math.inf
    return 2.0 * epsilon / ((mu_min - epsilon) * duration)


def projective_signature_error_bound(
    absolute_rate_error: float, lambda_max: float
) -> float:
    if not math.isfinite(absolute_rate_error) or absolute_rate_error >= lambda_max:
        return math.inf
    return 2.0 * absolute_rate_error / (lambda_max - absolute_rate_error)


def shots_for_projective_bound(
    target_bound: float,
    n: int,
    alpha: float,
    transition_minimum_eigenvalue: float,
    duration: float,
    lambda_max: float,
) -> int:
    """Find the smallest equal column allocation meeting a conservative bound."""
    lower, upper = 1, 1
    while True:
        eps = bernstein_operator_radius(n, upper, alpha)
        rate = log_generator_error_bound(
            transition_minimum_eigenvalue, eps, duration
        )
        bound = projective_signature_error_bound(rate, lambda_max)
        if bound <= target_bound:
            break
        upper *= 2
        if upper > 2_000_000_000:
            raise RuntimeError("requested distribution-free bound is impractical")
    while lower < upper:
        middle = (lower + upper) // 2
        eps = bernstein_operator_radius(n, middle, alpha)
        rate = log_generator_error_bound(
            transition_minimum_eigenvalue, eps, duration
        )
        bound = projective_signature_error_bound(rate, lambda_max)
        if bound <= target_bound:
            upper = middle
        else:
            lower = middle + 1
    return lower


def scalar_inverse_log_amplification(x: float) -> float:
    """Worst-mode local inverse-log amplification after lambda_max=1 scaling."""
    if x <= 0:
        return math.inf
    return math.exp(x) / x


def analyze_design(
    shot_levels: tuple[int, ...] = (10_000, 25_000, 50_000),
    replicates: int = 300,
) -> dict:
    a, w, s = strict_operator()
    eigenvalues = np.linalg.eigvalsh(a)
    lambda_max = float(eigenvalues[-1])
    gap = float(eigenvalues[1])
    duration = 1.0 / lambda_max
    target_signature = eigenvalues[1:] / lambda_max
    transition = expm(-duration * a)
    transition_2 = transition @ transition
    mu_min = float(np.linalg.eigvalsh(transition)[0])
    mu_min_2 = float(np.linalg.eigvalsh(transition_2)[0])

    target = StrictTarget(
        n=N,
        omega=OMEGA_STRICT,
        phi=PHI_STRICT,
        beta=BETA_STRICT,
        eta=ETA_STRICT,
        row_sum=s,
        lambda_max=lambda_max,
        spectral_gap=gap,
        tau_dimensionless=duration,
    )

    # Exact analytic design: exp(x)/x is uniquely minimized at x=1.
    time_grid = np.linspace(0.15, 1.75, 321)
    amplification = np.asarray(
        [scalar_inverse_log_amplification(float(x)) for x in time_grid]
    )
    selected_index = int(np.argmin(amplification))
    selected_x = float(time_grid[selected_index])

    rng = np.random.default_rng(SEED)
    shot_rows = []
    distributions: dict[int, list[float]] = {}
    for shots in shot_levels:
        distances: list[float] = []
        spd_failures = 0
        for _ in range(replicates):
            counts = sample_counts(transition, shots, rng)
            signature, _ = signature_from_counts(counts, duration)
            if signature is None:
                spd_failures += 1
                continue
            distances.append(
                float(np.max(np.abs(signature - target_signature)))
            )
        distributions[shots] = distances
        epsilon = bernstein_operator_radius(N, shots, 0.05)
        rate_bound = log_generator_error_bound(mu_min, epsilon, duration)
        projective_bound = projective_signature_error_bound(
            rate_bound, lambda_max
        )
        shot_rows.append(
            {
                "shots_per_preparation": shots,
                "replicates": replicates,
                "spd_failures": spd_failures,
                "mean_projective_distance": float(np.mean(distances)),
                "p95_projective_distance": float(np.quantile(distances, 0.95)),
                "p214_0p02_pass_rate": float(
                    np.mean(np.asarray(distances) <= 0.02)
                ),
                "matrix_bernstein_operator_radius_alpha_0p05": epsilon,
                "distribution_free_rate_error_bound": rate_bound,
                "distribution_free_projective_bound": projective_bound,
            }
        )

    certified_shots = shots_for_projective_bound(
        target_bound=0.02,
        n=N,
        alpha=0.05,
        transition_minimum_eigenvalue=mu_min,
        duration=duration,
        lambda_max=lambda_max,
    )

    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.8), constrained_layout=True)
    axes[0].plot(time_grid, amplification, color="#1F5A99", lw=2.2)
    axes[0].axvline(1.0, color="#A61B1B", ls="--", label=r"$\tau\lambda_{\max}=1$")
    axes[0].scatter([1.0], [math.e], color="#A61B1B", zorder=4)
    axes[0].set_xlabel(r"dimensionless design time $x=\tau\lambda_{\max}$")
    axes[0].set_ylabel(r"worst-mode inverse-log factor $e^x/x$")
    axes[0].set_title("Exact time-design optimum")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.22)

    for shots in shot_levels:
        axes[1].hist(
            distributions[shots],
            bins=34,
            alpha=0.43,
            label=f"{shots:,} shots/input",
        )
    axes[1].axvline(0.02, color="black", ls="--", label="frozen fingerprint gate")
    axes[1].set_xlabel("maximum projective spectral error")
    axes[1].set_ylabel("simulation count")
    axes[1].set_title("Finite-shot planning (synthetic only)")
    axes[1].legend(fontsize=7)
    figure_path = FIGURE_DIR / "p240_optimal_tomography.png"
    fig.savefig(figure_path, dpi=210)
    plt.close(fig)

    lock_core = {
        "protocol_id": "FIN-LAB-P240-OPTIMAL-HEAT-TOMOGRAPHY-V1",
        "frozen_target": {
            **asdict(target),
            "strict_kernel": (
                "cos(0.18575*d+0.16250)/(1+d^1.8), "
                "zero diagonal, undirected cyclic distance on Z12"
            ),
            "row_sum_reference": ROW_SUM_REFERENCE,
        },
        "preparations": list(range(N)),
        "preparation_allocation": "uniform by C12 minimax symmetrization",
        "base_time_rule": "tau*lambda_max=1",
        "base_dimensionless_time": duration,
        "registered_process_times": {
            "calibration": "tau",
            "held_out": "2*tau",
            "optional_duality_diagnostic": ["tau/4", "tau/2", "tau"],
        },
        "recommended_shots_per_preparation_per_registered_time": 50_000,
        "registered_event_total": 2 * N * 50_000,
        "fingerprint_distance_threshold": 0.02,
        "primary_semigroup_test": (
            "maximum column total-variation distance between held-out "
            "P_2tau and P_tau@P_tau with a preregistered finite-count bound"
        ),
        "secondary_test": "projective strict spectral fingerprint",
        "physical_clock_required": True,
        "apparatus_required": True,
        "external_record_required": True,
        "program_241_gate_required": True,
        "synthetic_data_are_external_evidence": False,
        "strict_selector_or_unit_exported": False,
        "legacy_to_strict_bridge_completed": False,
        "external_validation_claimed": False,
    }
    write_canonical(LOCK_PATH, lock_core)

    result = {
        "program": 240,
        "constructed_object": "OptimalHeatKernelTomographyDesign",
        "strict_target": asdict(target),
        "strict_row_sum_reference_error": abs(s - ROW_SUM_REFERENCE),
        "strict_weight_minimum": float(
            np.min(w[np.triu_indices_from(w, k=1)])
        ),
        "generator_minimum_eigenvalue": float(eigenvalues[0]),
        "transition_minimum_eigenvalue_at_tau": mu_min,
        "transition_minimum_eigenvalue_at_2tau": mu_min_2,
        "exact_optimum_x": 1.0,
        "grid_optimum_x": selected_x,
        "uniform_allocation_theorem": (
            "For a convex cyclic-permutation-invariant design risk, averaging "
            "any allocation over C12 cannot increase risk and produces the "
            "uniform allocation. Thus a minimax design exists with equal "
            "shots for all twelve basis preparations."
        ),
        "time_optimum_theorem": (
            "After normalizing lambda_max=1, the local absolute-noise "
            "amplification of inverse-log recovery for the fastest mode is "
            "g(x)=exp(x)/x. Since g'(x)=exp(x)(x-1)/x^2, its unique positive "
            "minimum is x=1."
        ),
        "matrix_concentration_theorem": (
            "For m independent categorical shots in every one of n columns, "
            "the raw transition estimate obeys rectangular matrix Bernstein "
            "with v<=n/(2m), R<=sqrt(2)/m, and dilation dimension 2n. "
            "Symmetrization does not increase operator error. If epsilon is "
            "below lambda_min(P), Weyl plus log Lipschitz continuity yields "
            "the recorded nonasymptotic generator and projective bounds."
        ),
        "shot_rows": shot_rows,
        "shots_per_preparation_for_distribution_free_0p02_projective_bound": (
            certified_shots
        ),
        "design_lock": {
            "path": LOCK_PATH.name,
            "sha256": sha256_file(LOCK_PATH),
        },
        "figure": {
            "path": str(figure_path.relative_to(ROOT)),
            "sha256": sha256_file(figure_path),
        },
        "claim_status": {
            "analytic_time_optimum": "PROVEN_IN_DECLARED_LOCAL_NOISE_MODEL",
            "uniform_minimax_allocation": "PROVEN_BY_CYCLIC_SYMMETRIZATION",
            "finite_shot_rows": "SYNTHETIC_METHOD_EVIDENCE",
            "external_data_used": False,
            "external_validation_claimed": False,
        },
    }
    RESULTS_PATH.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--replicates",
        type=int,
        default=300,
        help="deterministic planning simulations per shot level",
    )
    args = parser.parse_args()
    result = analyze_design(replicates=args.replicates)
    print(json.dumps(result, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
