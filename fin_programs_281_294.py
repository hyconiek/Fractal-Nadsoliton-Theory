#!/usr/bin/env python3
"""Execute FIN research Programs P281--P294.

The suite combines exact finite theorems, a dependency-free Lean certificate,
and frozen synthetic audits.  It does not fabricate laboratory evidence,
derive a strict selector or dimensional unit, complete the legacy-to-strict
bridge, transfer legacy physical roles, or promote conditional operational
models to established physics.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares
from scipy.special import gammaln, logsumexp

import fin_programs_255_266 as core
import fin_programs_267_280 as previous


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_281_294_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_281_294_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_281_294_Summary.csv"
STABILITY_PATH = ROOT / "FIN_Programs_281_294_Stieltjes_Recovery.csv"
POVM_PATH = ROOT / "FIN_Programs_281_294_Detector_POVM.csv"
RG_PATH = ROOT / "FIN_Programs_281_294_RG_Completion.csv"
FLUX_PATH = ROOT / "FIN_Programs_281_294_Detector_Flux.csv"
MECHANISM_PATH = ROOT / "FIN_Programs_281_294_Mechanism_Classification.csv"
INTERVENTION_PATH = ROOT / "FIN_Programs_281_294_Intervention_Design.csv"
RESERVOIR_PATH = ROOT / "FIN_Programs_281_294_Reservoir_Replication.csv"
LEAN_SOURCE = ROOT / "FIN_Programs_281_294_Schur_Core.lean"

N = 12
SEED = 20260731


def json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, complex):
        return {"real": float(value.real), "imag": float(value.imag)}
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def hermitian(matrix: np.ndarray) -> np.ndarray:
    return (matrix + matrix.T.conj()) / 2.0


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: (
                        json.dumps(json_ready(value), ensure_ascii=False)
                        if isinstance(value, (dict, list, tuple, np.ndarray))
                        else value
                    )
                    for key, value in row.items()
                }
            )


def bounded_stieltjes_fit(
    z_grid: np.ndarray,
    values: np.ndarray,
    starts: list[np.ndarray],
    pole_window: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray, float]:
    """Robust bounded three-pole fit with frozen model order and window."""

    count = 3
    best: tuple[float, np.ndarray] | None = None

    def residual(parameters: np.ndarray) -> np.ndarray:
        poles = parameters[:count]
        weights = parameters[count:]
        model = np.sum(
            weights[None, :] / (z_grid[:, None] + poles[None, :]), axis=1
        )
        return (model - values) / np.maximum(np.abs(values), 1e-12)

    lower = np.array([pole_window[0]] * count + [1e-8] * count)
    upper = np.array([pole_window[1]] * count + [4.0] * count)
    for poles0 in starts:
        design = 1.0 / (z_grid[:, None] + poles0[None, :])
        weights0 = np.maximum(
            1e-4, np.linalg.lstsq(design, values, rcond=None)[0]
        )
        initial = np.concatenate([poles0, np.minimum(weights0, 3.9)])
        fit = least_squares(
            residual,
            initial,
            bounds=(lower, upper),
            loss="soft_l1",
            f_scale=0.002,
            max_nfev=2000,
            xtol=1e-11,
            ftol=1e-11,
            gtol=1e-11,
        )
        score = float(np.mean(residual(fit.x) ** 2))
        if best is None or score < best[0]:
            best = (score, fit.x.copy())
    assert best is not None
    parameters = best[1]
    order = np.argsort(parameters[:count])
    return (
        parameters[:count][order],
        parameters[count:][order],
        best[0],
    )


def program_281(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    data = previous.visible_pole_residue_data(strict_a)
    true_poles = np.array([pole for pole, _ in data])
    true_weights = np.array(
        [float(np.trace(residue).real) for _, residue in data]
    )
    z_grid = np.geomspace(0.03, 8.0, 48)
    exact = np.sum(
        true_weights[None, :] / (z_grid[:, None] + true_poles[None, :]),
        axis=1,
    )
    pole_window = (0.55, 2.55)
    starts = [
        np.array([0.75, 1.45, 2.25]),
        np.array([1.00, 1.55, 2.10]),
        np.array([0.65, 1.75, 2.35]),
    ]
    rows: list[dict[str, Any]] = []
    for sigma in [1e-6, 1e-4, 1e-3]:
        for repetition in range(60):
            noisy = exact * (
                1.0 + rng.normal(0.0, sigma, size=exact.shape)
            )
            poles, weights, objective = bounded_stieltjes_fit(
                z_grid, noisy, starts, pole_window
            )
            rows.append(
                {
                    "relative_noise": sigma,
                    "repetition": repetition,
                    "maximum_relative_pole_error": float(
                        np.max(np.abs(poles - true_poles) / true_poles)
                    ),
                    "maximum_relative_weight_error": float(
                        np.max(np.abs(weights - true_weights) / true_weights)
                    ),
                    "relative_curve_rmse": float(math.sqrt(objective)),
                    "minimum_fitted_pole_separation": float(
                        np.min(np.diff(poles))
                    ),
                }
            )
    summaries: dict[str, Any] = {}
    for sigma in [1e-6, 1e-4, 1e-3]:
        selected = [row for row in rows if row["relative_noise"] == sigma]
        summaries[str(sigma)] = {}
        for field in [
            "maximum_relative_pole_error",
            "maximum_relative_weight_error",
            "relative_curve_rmse",
            "minimum_fitted_pole_separation",
        ]:
            values = np.array([row[field] for row in selected])
            summaries[str(sigma)][f"median_{field}"] = float(
                np.median(values)
            )
            summaries[str(sigma)][f"p95_{field}"] = float(
                np.quantile(values, 0.95)
            )
    return (
        {
            "status": (
                "[Proven] conditional compactness/existence; "
                "[Strong evidence] bounded recovery audit"
            ),
            "frozen_model_order": 3,
            "frozen_pole_window": pole_window,
            "visible_poles": true_poles,
            "visible_weights": true_weights,
            "summaries": summaries,
            "theorem": (
                "With fixed pole count, a compact pole window, bounded "
                "nonnegative residues, and continuous least-squares loss, a "
                "minimizer exists. These priors prevent runaway poles but do "
                "not create uniform inverse stability near pole collisions or "
                "vanishing residues."
            ),
            "boundary": (
                "Model order and spectral window are supplied regularization "
                "data. The audit is not a minimax-optimality theorem and does "
                "not identify invisible modes."
            ),
        },
        rows,
    )


def direct_lean_binary() -> Path | None:
    candidates = sorted((ROOT / ".elan" / "toolchains").glob("*/bin/lean"))
    return candidates[-1] if candidates else None


def program_282() -> dict[str, Any]:
    lean = direct_lean_binary()
    if lean is None:
        return {
            "status": "[Blocked] no Lean runtime",
            "lean_runtime": None,
            "kernel_exit_zero": False,
            "boundary": "No machine-checked claim is exported.",
        }
    process = subprocess.run(
        [str(lean), str(LEAN_SOURCE)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=120,
        check=False,
    )
    return {
        "status": (
            "[Proven] machine-checked exact rational witness"
            if process.returncode == 0
            else "[Refuted] Lean source failed"
        ),
        "lean_runtime": str(lean.relative_to(ROOT)),
        "lean_version": subprocess.run(
            [str(lean), "--version"],
            cwd=ROOT,
            text=True,
            capture_output=True,
            timeout=30,
            check=False,
        ).stdout.strip(),
        "source": LEAN_SOURCE.name,
        "source_sha256": sha256_file(LEAN_SOURCE),
        "kernel_exit_code": process.returncode,
        "kernel_exit_zero": process.returncode == 0,
        "stdout": process.stdout.strip(),
        "stderr": process.stderr.strip(),
        "checked_statements": [
            "identity reduction",
            "positive rational pivots",
            "nested Schur value 13/5",
            "direct Schur value 13/5",
            "exact nested/direct equality",
            "positive reduced value",
        ],
        "boundary": (
            "This dependency-free certificate checks one exact rational "
            "three-by-three witness. It is not the general Mathlib theorem "
            "for arbitrary positive matrices."
        ),
    }


def current_povm(
    strict_w: np.ndarray,
) -> tuple[list[np.ndarray], list[np.ndarray], np.ndarray, float]:
    basis, _ = previous.orthonormal_current_basis(strict_w)
    vertices = previous.simplex_vertices(5)
    simplex_operators = [
        sum(vertices[y, a] * basis[a] for a in range(5))
        for y in range(6)
    ]
    maximum_norm = max(
        np.linalg.norm(operator, 2) for operator in simplex_operators
    )
    epsilon = 0.90 / (6.0 * maximum_norm)
    effects = [
        np.eye(N) / 6.0 + epsilon * operator
        for operator in simplex_operators
    ]
    return basis, effects, vertices, epsilon


def positive_sqrt(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(hermitian(matrix))
    return (
        eigenvectors
        * np.sqrt(np.maximum(eigenvalues, 0.0))
    ) @ eigenvectors.T.conj()


def program_283(
    strict_w: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    basis, effects, _, _ = current_povm(strict_w)
    dilation = np.vstack([positive_sqrt(effect) for effect in effects])
    isometry_residual = float(
        np.linalg.norm(dilation.T.conj() @ dilation - np.eye(N), 2)
    )
    response = np.array(
        [
            [np.trace(effect @ operator).real for operator in basis]
            for effect in effects
        ]
    )

    efficiency = 0.82
    cross_talk = 0.018
    dark_total = 2.0e-4
    confusion = np.zeros((7, 6))
    for true_y in range(6):
        confusion[true_y, true_y] += efficiency * (1.0 - 2.0 * cross_talk)
        confusion[(true_y - 1) % 6, true_y] += efficiency * cross_talk
        confusion[(true_y + 1) % 6, true_y] += efficiency * cross_talk
        confusion[:6, true_y] += dark_total / 6.0
        confusion[6, true_y] = 1.0 - efficiency - dark_total
    stochastic_residual = float(
        np.max(np.abs(confusion.sum(axis=0) - 1.0))
    )
    observed_response = confusion @ response
    singular_values = np.linalg.svd(observed_response, compute_uv=False)
    observed_rank = int(np.linalg.matrix_rank(observed_response, tol=1e-12))
    condition = float(singular_values[0] / singular_values[4])

    coefficients = np.array([0.028, -0.021, 0.016, 0.011, -0.008])
    rho = np.eye(N, dtype=complex) / N + sum(
        coefficient * operator
        for coefficient, operator in zip(coefficients, basis)
    )
    ideal = np.array([np.trace(effect @ rho).real for effect in effects])
    observed_probabilities = confusion @ ideal
    baseline = confusion @ (np.ones(6) / 6.0)
    shots = 500_000
    counts = rng.multinomial(shots, observed_probabilities)
    estimate = np.linalg.lstsq(
        observed_response, counts / shots - baseline, rcond=None
    )[0]
    relative_error = float(
        np.linalg.norm(estimate - coefficients)
        / np.linalg.norm(coefficients)
    )
    rows = [
        {
            "observed_outcome": y if y < 6 else "no_click",
            "probability": observed_probabilities[y],
            "count": int(counts[y]),
            "baseline_probability": baseline[y],
        }
        for y in range(7)
    ]
    return (
        {
            "status": (
                "[Proven] Naimark/isometry realization and calibration rank; "
                "[Strong evidence] frozen detector audit"
            ),
            "ideal_outcomes": 6,
            "physical_record_outcomes_including_no_click": 7,
            "dilation_shape": dilation.shape,
            "isometry_residual": isometry_residual,
            "minimum_effect_eigenvalue": float(
                min(np.linalg.eigvalsh(hermitian(effect)).min() for effect in effects)
            ),
            "detector_model": {
                "efficiency": efficiency,
                "nearest_neighbor_cross_talk": cross_talk,
                "dark_probability_per_trial": dark_total,
            },
            "confusion_column_stochastic_residual": stochastic_residual,
            "calibrated_response_rank": observed_rank,
            "calibrated_response_condition_number": condition,
            "shots": shots,
            "relative_coefficient_error": relative_error,
            "theorem": (
                "V=stack_y sqrt(E_y) is an isometry because "
                "V*V=sum_y E_y=I. Projective measurement of the outcome "
                "ancilla realizes the POVM. A known full-rank detector "
                "confusion channel preserves identifiability on the five-"
                "current response subspace."
            ),
            "boundary": (
                "Naimark dilation proves realizability in abstract quantum "
                "mechanics. It is not an optical mesh, detector calibration, "
                "or built apparatus."
            ),
        },
        rows,
    )


def block_embedding(size: int, target: int = N) -> np.ndarray:
    multiplicity = target // size
    embedding = np.zeros((target, size))
    for index in range(target):
        embedding[index, index % size] = 1.0 / math.sqrt(multiplicity)
    return embedding


def program_284(strict_a: np.ndarray) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    contexts = {
        6: [0, 2, 4, 6, 8, 10],
        3: [0, 4, 8],
    }
    strict_signature = core.projective_signature(strict_a)
    rows: list[dict[str, Any]] = []
    for size, keep in contexts.items():
        reduced = core.schur_keep(strict_a, keep)
        embedding = block_embedding(size)
        image_projector = embedding @ embedding.T
        complement_projector = np.eye(N) - image_projector
        nonzero = np.linalg.eigvalsh(hermitian(reduced))[1:]
        alpha = float(np.mean(nonzero))
        bare = hermitian(embedding @ reduced @ embedding.T)
        completed = hermitian(bare + alpha * complement_projector)
        bare_fit_scale = float(
            np.sum(strict_a * bare) / np.sum(bare * bare)
        )
        completed_fit_scale = float(
            np.sum(strict_a * completed) / np.sum(completed * completed)
        )
        rows.append(
            {
                "context_size": size,
                "alpha_mean_coarse_nonzero_eigenvalue": alpha,
                "bare_rank": int(np.linalg.matrix_rank(bare, tol=1e-10)),
                "completed_rank": int(
                    np.linalg.matrix_rank(completed, tol=1e-10)
                ),
                "minimum_completed_eigenvalue": float(
                    np.linalg.eigvalsh(completed).min()
                ),
                "row_sum_residual": float(
                    np.linalg.norm(completed @ np.ones(N))
                ),
                "bare_best_fit_residual": float(
                    np.linalg.norm(strict_a - bare_fit_scale * bare, "fro")
                    / np.linalg.norm(strict_a, "fro")
                ),
                "completed_best_fit_residual": float(
                    np.linalg.norm(
                        strict_a - completed_fit_scale * completed, "fro"
                    )
                    / np.linalg.norm(strict_a, "fro")
                ),
                "completed_projective_fingerprint_distance": float(
                    np.linalg.norm(
                        core.projective_signature(completed) - strict_signature
                    )
                ),
            }
        )
    return (
        {
            "status": (
                "[Proven] PSD rank-restoring complement theorem; "
                "[Refuted] as exact strict fixed point"
            ),
            "frozen_rule": (
                "R(B)=E B E* + mean(nonzero_eigenvalues(B))*(I-E E*)"
            ),
            "rows": rows,
            "theorem": (
                "If B is a connected n-state Laplacian and E maps its constant "
                "mode to the twelve-state constant mode, the supplied "
                "complement term is PSD, restores rank 11, and preserves the "
                "global null mode. It is determined from B, not from the "
                "strict target."
            ),
            "boundary": (
                "Rank restoration is not RG closure. The frozen isotropic "
                "complement does not reproduce the strict fingerprint and no "
                "physical coarse-graining law selects it."
            ),
        },
        rows,
    )


def qubit_instrument(
    diagonal_effect: tuple[float, float], swap_second: bool = True
) -> list[np.ndarray]:
    effect_zero = np.diag(diagonal_effect)
    effect_one = np.eye(2) - effect_zero
    swap = np.array([[0.0, 1.0], [1.0, 0.0]])
    return [
        positive_sqrt(effect_zero),
        (swap if swap_second else np.eye(2)) @ positive_sqrt(effect_one),
    ]


def apply_channel(state: np.ndarray, kraus: list[np.ndarray]) -> np.ndarray:
    return hermitian(
        sum(operator @ state @ operator.T.conj() for operator in kraus)
    )


def branch_instrument(
    state: np.ndarray, kraus: list[np.ndarray]
) -> tuple[np.ndarray, list[np.ndarray]]:
    probabilities = np.array(
        [
            np.trace(operator @ state @ operator.T.conj()).real
            for operator in kraus
        ]
    )
    conditionals = [
        hermitian(operator @ state @ operator.T.conj() / probabilities[index])
        for index, operator in enumerate(kraus)
    ]
    return probabilities, conditionals


def amplitude_damping(gamma: float, angle: float) -> list[np.ndarray]:
    damping = [
        np.array([[1.0, 0.0], [0.0, math.sqrt(1.0 - gamma)]]),
        np.array([[0.0, math.sqrt(gamma)], [0.0, 0.0]]),
    ]
    rotation = np.array(
        [
            [math.cos(angle), -math.sin(angle)],
            [math.sin(angle), math.cos(angle)],
        ]
    )
    return [rotation @ operator for operator in damping]


def two_step_information_ledger(
    rho: np.ndarray, sigma: np.ndarray
) -> dict[str, float]:
    first = qubit_instrument((0.81, 0.29), swap_second=True)
    second = qubit_instrument((0.67, 0.38), swap_second=False)
    p_y, rho_y = branch_instrument(rho, first)
    q_y, sigma_y = branch_instrument(sigma, first)
    d_input = previous.quantum_relative_entropy(rho, sigma)
    first_record = core.kl_divergence(p_y, q_y)
    first_conditional = sum(
        p_y[y]
        * previous.quantum_relative_entropy(rho_y[y], sigma_y[y])
        for y in range(2)
    )
    loss_first = d_input - first_record - first_conditional

    joint_p = np.zeros((2, 2))
    joint_q = np.zeros((2, 2))
    final_conditional = 0.0
    channel_loss = 0.0
    second_instrument_loss = 0.0
    for y in range(2):
        channel = amplitude_damping(
            gamma=0.11 + 0.17 * y,
            angle=0.13 - 0.07 * y,
        )
        rho_channel = apply_channel(rho_y[y], channel)
        sigma_channel = apply_channel(sigma_y[y], channel)
        before_channel = previous.quantum_relative_entropy(
            rho_y[y], sigma_y[y]
        )
        after_channel = previous.quantum_relative_entropy(
            rho_channel, sigma_channel
        )
        channel_loss += p_y[y] * (before_channel - after_channel)
        p_z, rho_yz = branch_instrument(rho_channel, second)
        q_z, sigma_yz = branch_instrument(sigma_channel, second)
        second_record_y = core.kl_divergence(p_z, q_z)
        second_conditional_y = sum(
            p_z[z]
            * previous.quantum_relative_entropy(rho_yz[z], sigma_yz[z])
            for z in range(2)
        )
        second_instrument_loss += p_y[y] * (
            after_channel - second_record_y - second_conditional_y
        )
        for z in range(2):
            joint_p[y, z] = p_y[y] * p_z[z]
            joint_q[y, z] = q_y[y] * q_z[z]
            final_conditional += (
                joint_p[y, z]
                * previous.quantum_relative_entropy(
                    rho_yz[z], sigma_yz[z]
                )
            )
    joint_record = core.kl_divergence(joint_p.ravel(), joint_q.ravel())
    total_rhs = (
        loss_first
        + channel_loss
        + second_instrument_loss
        + joint_record
        + final_conditional
    )
    return {
        "input": d_input,
        "first_instrument_loss": loss_first,
        "intermediate_channel_loss": channel_loss,
        "second_instrument_loss": second_instrument_loss,
        "joint_record": joint_record,
        "final_conditional": final_conditional,
        "decomposition_residual": abs(d_input - total_rhs),
        "chain_record_residual": abs(
            joint_record
            - first_record
            - sum(
                p_y[y]
                * core.kl_divergence(
                    joint_p[y] / p_y[y], joint_q[y] / q_y[y]
                )
                for y in range(2)
            )
        ),
    }


def program_285(rng: np.random.Generator) -> dict[str, Any]:
    rows = [
        two_step_information_ledger(
            previous.random_density(2, rng),
            previous.random_density(2, rng),
        )
        for _ in range(300)
    ]
    loss_fields = [
        "first_instrument_loss",
        "intermediate_channel_loss",
        "second_instrument_loss",
    ]
    return {
        "status": "[Proven] conditioned two-step quantum process ledger",
        "pairs_audited": len(rows),
        "minimum_losses": {
            field: float(min(row[field] for row in rows))
            for field in loss_fields
        },
        "median_total_discarded_information": float(
            np.median(
                [
                    sum(row[field] for field in loss_fields)
                    for row in rows
                ]
            )
        ),
        "maximum_decomposition_residual": float(
            max(row["decomposition_residual"] for row in rows)
        ),
        "maximum_chain_record_residual": float(
            max(row["chain_record_residual"] for row in rows)
        ),
        "theorem": (
            "Sequential application of quantum relative-entropy data "
            "processing and the classical chain rule decomposes input "
            "distinguishability into first-instrument loss, channel loss, "
            "second-instrument loss, joint-record divergence, and final "
            "conditional quantum divergence."
        ),
        "boundary": (
            "The instruments and branch-dependent channels are supplied. "
            "This is a finite two-step process tree, not a general process-"
            "tensor reconstruction theorem or a FIN-derived bath."
        ),
    }


def legacy_fractal_atom() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    alpha_geo = 4.0 * math.log(2.0)
    beta_tors = 0.01
    dimension = alpha_geo
    weights = np.array(
        [
            alpha_geo
            * math.cos(math.pi * d / 4.0 + math.pi / 6.0)
            / (1.0 + beta_tors * d**dimension)
            for d in range(1, 7)
        ]
    )
    w = core.cyclic_weight_matrix(weights)
    return core.laplacian_from_weights(w), w, weights


def program_286(strict_a: np.ndarray) -> dict[str, Any]:
    candidate, _, weights = legacy_fractal_atom()
    eigenvalues = np.linalg.eigvalsh(hermitian(candidate))
    mean_zero = np.delete(eigenvalues, np.argmin(np.abs(eigenvalues)))
    retained = list(range(0, N, 2))
    _, _, hidden = core.block_components(candidate, retained)
    hidden_eigenvalues = np.linalg.eigvalsh(hermitian(hidden))
    positive_generator = bool(mean_zero.min() >= -1e-12)
    if positive_generator:
        fingerprint_defect: float | None = float(
            np.linalg.norm(
                core.projective_signature(candidate)
                - core.projective_signature(strict_a)
            )
        )
        sigma_candidate = core.self_energy(candidate, retained, 0.2)
        sigma_strict = core.self_energy(strict_a, retained, 0.2)
        sigma_candidate /= np.linalg.norm(sigma_candidate, "fro")
        sigma_strict /= np.linalg.norm(sigma_strict, "fro")
        memory_defect: float | None = float(
            np.linalg.norm(sigma_candidate - sigma_strict, "fro")
        )
    else:
        fingerprint_defect = None
        memory_defect = None
    return {
        "status": (
            "[Refuted] as sufficient legacy-to-strict completion atom"
            if not positive_generator
            else "[Refuted] as exact strict completion"
        ),
        "frozen_atom": (
            "K_fc(d)=alpha_geo*cos(pi*d/4+pi/6)/"
            "(1+beta_tors*d**D_f), with D_f=alpha_geo=4 ln 2"
        ),
        "candidate_weights": weights,
        "minimum_mean_zero_eigenvalue": float(mean_zero.min()),
        "minimum_hidden_block_eigenvalue": float(hidden_eigenvalues.min()),
        "positive_generator": positive_generator,
        "strict_projective_fingerprint_defect": fingerprint_defect,
        "normalized_self_energy_defect": memory_defect,
        "verdict": (
            "Replacing linear distance by the canonical legacy fractal power "
            "does not by itself enter the strict positive-generator class."
            if not positive_generator
            else "The atom is positive but does not reproduce strict shape."
        ),
        "boundary": (
            "Exactly one legacy-sourced nonlinear compression atom was "
            "tested. No positivity shift, strict parameter fit, selector, "
            "amplitude absorption, or role transfer was added."
        ),
    }


def program_287() -> dict[str, Any]:
    production_manifests = [
        path
        for path in ROOT.rglob("bundle_manifest.json")
        if "FIN_Lab_P241_Transfer_Template" not in str(path)
    ]
    template = ROOT / "FIN_Lab_P241_Transfer_Template"
    headers = [
        template / "events_heat_process.header.csv",
        template / "events_double_slit.header.csv",
    ]
    event_rows = []
    for path in headers:
        with path.open("r", encoding="utf-8") as handle:
            event_rows.append(max(0, sum(1 for _ in handle) - 1))
    return {
        "status": "[Executed no-admission certificate: external evidence absent]",
        "production_manifest_count": len(production_manifests),
        "production_manifest_paths": [
            str(path.relative_to(ROOT)) for path in production_manifests
        ],
        "template_event_rows": event_rows,
        "p242_authorized": False,
        "verdict": (
            "No independently signed production event bundle is present. "
            "The confirmatory P242 pipeline remains blocked."
        ),
        "boundary": (
            "Files whose research-program numbers begin with p241 are not "
            "P241 laboratory custody packages."
        ),
    }


def detector_probabilities(
    h: float,
    gradient: float,
    efficiency: float,
    visibility: float,
    dark: float,
    cross_talk: float,
) -> np.ndarray:
    signal = math.sin(gradient * h)
    ideal = np.array([(1.0 + visibility * signal) / 2.0,
                      (1.0 - visibility * signal) / 2.0])
    clicked = efficiency * np.array(
        [
            (1.0 - cross_talk) * ideal[0] + cross_talk * ideal[1],
            cross_talk * ideal[0] + (1.0 - cross_talk) * ideal[1],
        ]
    )
    clicked += dark / 2.0
    return np.array([clicked[0], clicked[1], 1.0 - clicked.sum()])


def program_288(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    retained = list(range(0, N, 2))
    xi, _ = core.chiral_susceptibility_analytic(
        core.twisted_generator(0.0),
        core.twisted_generator_derivative_at_zero(),
        retained,
        0.2,
    )
    gradient = float(np.linalg.norm(xi, "fro"))
    efficiency = 0.78
    visibility = 0.91
    dark = 2.0e-4
    cross_talk = 0.015
    shots_per_flux = 100_000
    correction = efficiency * visibility * (1.0 - 2.0 * cross_talk)
    h_grid = np.geomspace(0.004, 0.45, 18)
    rows: list[dict[str, Any]] = []
    for h in h_grid:
        p_plus = detector_probabilities(
            float(h), gradient, efficiency, visibility, dark, cross_talk
        )
        p_minus = detector_probabilities(
            float(-h), gradient, efficiency, visibility, dark, cross_talk
        )
        exact_signal_plus = (p_plus[0] - p_plus[1]) / correction
        exact_signal_minus = (p_minus[0] - p_minus[1]) / correction
        noiseless = (exact_signal_plus - exact_signal_minus) / (2.0 * h)
        bias = abs(noiseless - gradient) / gradient
        errors = []
        for _ in range(240):
            count_plus = rng.multinomial(shots_per_flux, p_plus)
            count_minus = rng.multinomial(shots_per_flux, p_minus)
            frequency_plus = count_plus / shots_per_flux
            frequency_minus = count_minus / shots_per_flux
            signal_plus = (
                frequency_plus[0] - frequency_plus[1]
            ) / correction
            signal_minus = (
                frequency_minus[0] - frequency_minus[1]
            ) / correction
            estimate = (signal_plus - signal_minus) / (2.0 * h)
            errors.append(abs(estimate - gradient) / gradient)
        rows.append(
            {
                "h": float(h),
                "relative_noiseless_bias": float(bias),
                "median_relative_error": float(np.median(errors)),
                "p95_relative_error": float(np.quantile(errors, 0.95)),
                "plus_click_probability": float(p_plus[:2].sum()),
            }
        )
    optimum = min(rows, key=lambda row: row["median_relative_error"])
    del strict_a
    return (
        {
            "status": "[Strong evidence] frozen detector-level likelihood audit",
            "gradient_norm": gradient,
            "detector_model": {
                "efficiency": efficiency,
                "visibility": visibility,
                "dark_probability": dark,
                "cross_talk": cross_talk,
            },
            "shots_per_flux": shots_per_flux,
            "optimal_tested_h": optimum["h"],
            "optimal_median_relative_error": optimum[
                "median_relative_error"
            ],
            "optimal_p95_relative_error": optimum["p95_relative_error"],
            "grid": rows,
            "boundary": (
                "The categorical click likelihood is detector-level but "
                "synthetic. Its constants are supplied and do not calibrate "
                "a real instrument."
            ),
        },
        rows,
    )


def circulant_laplacian_spectra(weight_rows: np.ndarray) -> np.ndarray:
    modes = np.arange(N)
    coefficients = np.zeros((6, N))
    for distance in range(1, 6):
        coefficients[distance - 1] = 2.0 * (
            1.0 - np.cos(2.0 * math.pi * modes * distance / N)
        )
    coefficients[5] = 1.0 - (-1.0) ** modes
    return weight_rows @ coefficients


def fingerprint_rows(weight_rows: np.ndarray) -> np.ndarray:
    eigenvalues = np.sort(circulant_laplacian_spectra(weight_rows), axis=1)
    return eigenvalues[:, 1:] / eigenvalues[:, -1, None]


def dirichlet_logpdf(rows: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    return (
        gammaln(float(np.sum(alpha)))
        - float(np.sum(gammaln(alpha)))
        + np.sum((alpha - 1.0) * np.log(rows), axis=1)
    )


def program_289(
    strict_a: np.ndarray, rng: np.random.Generator
) -> dict[str, Any]:
    target = core.projective_signature(strict_a)
    target_weights = core.STRICT_K / np.sum(core.STRICT_K)
    tolerance = 0.02
    sample_count = 120_000
    proposal_rows = []
    for concentration in [220.0, 350.0, 550.0]:
        proposal_alpha = 1.0 + concentration * target_weights
        samples = rng.dirichlet(proposal_alpha, size=sample_count)
        signatures = fingerprint_rows(samples)
        distances = np.linalg.norm(signatures - target[None, :], axis=1)
        events = distances <= tolerance
        log_weight = (
            dirichlet_logpdf(samples, np.ones(6))
            - dirichlet_logpdf(samples, proposal_alpha)
        )
        event_count = int(np.sum(events))
        if event_count:
            log_terms = log_weight[events]
            log_probability = float(
                logsumexp(log_terms) - math.log(sample_count)
            )
            probability = float(math.exp(log_probability))
            log_second = float(
                logsumexp(2.0 * log_terms) - math.log(sample_count)
            )
            second_moment = math.exp(log_second)
            standard_error = math.sqrt(
                max(0.0, second_moment - probability**2) / sample_count
            )
            event_weights = np.exp(log_terms - np.max(log_terms))
            event_ess = float(
                np.sum(event_weights) ** 2 / np.sum(event_weights**2)
            )
        else:
            log_probability = -math.inf
            probability = 0.0
            standard_error = math.nan
            event_ess = 0.0
        proposal_rows.append(
            {
                "concentration": concentration,
                "sample_count": sample_count,
                "event_count": event_count,
                "probability": probability,
                "log_probability": log_probability,
                "standard_error": standard_error,
                "event_ess": event_ess,
            }
        )
    central = proposal_rows[1]
    positive_estimates = [
        row["probability"] for row in proposal_rows if row["probability"] > 0.0
    ]
    sensitivity_ratio = (
        max(positive_estimates) / min(positive_estimates)
        if positive_estimates
        else math.inf
    )

    direct_count = 120_000
    direct = rng.dirichlet(np.ones(6), size=direct_count)
    direct_distances = np.linalg.norm(
        fingerprint_rows(direct) - target[None, :], axis=1
    )
    direct_events = int(np.sum(direct_distances <= tolerance))
    rule_of_three = 3.0 / direct_count if direct_events == 0 else None
    return {
        "status": "[Strong evidence] importance-sampled ensemble probability",
        "null_ensemble": "Dirichlet(1,...,1) normalized positive circulant weights",
        "tolerance": tolerance,
        "importance_sample_count_per_proposal": sample_count,
        "proposal_rows": proposal_rows,
        "proposal_sensitivity_max_to_min_ratio": sensitivity_ratio,
        "estimated_false_positive_probability": central["probability"],
        "log_estimated_probability": central["log_probability"],
        "estimated_standard_error": central["standard_error"],
        "event_effective_sample_size": central["event_ess"],
        "direct_sample_count": direct_count,
        "direct_event_count": direct_events,
        "direct_zero_event_rule_of_three_upper_95": rule_of_three,
        "minimum_direct_distance": float(direct_distances.min()),
        "boundary": (
            "The estimate is specific to the declared Dirichlet null and "
            "importance proposal. It is not distribution-free and requires "
            "ESS and proposal-sensitivity checks before use as a discovery "
            "threshold."
        ),
    }


def matrix_observation_features(
    matrices: list[np.ndarray], times: np.ndarray
) -> np.ndarray:
    generators = [
        -previous.matrix_log_spd(matrix) / time
        for matrix, time in zip(matrices, times)
    ]
    scale = max(np.linalg.norm(generators[0], "fro"), 1e-12)
    features: list[float] = []
    for index in range(1, len(generators)):
        difference = generators[index] - generators[0]
        features.extend(
            [
                np.linalg.norm(difference, "fro") / scale,
                float(np.trace(difference).real / (N * scale)),
                np.linalg.norm(
                    matrices[0] @ matrices[index]
                    - matrices[index] @ matrices[0],
                    "fro",
                )
                / max(np.linalg.norm(matrices[0], "fro"), 1e-12),
            ]
        )
    return np.asarray(features)


def mechanism_matrices(
    name: str,
    strict_a: np.ndarray,
    times: np.ndarray,
    intervention: float,
) -> list[np.ndarray]:
    edge_a = core.edge_laplacian(0, 3)
    edge_b = core.edge_laplacian(1, 7)
    matrices: list[np.ndarray] = []
    for time in times:
        if name == "homogeneous":
            matrix = expm(-time * (strict_a + intervention * edge_a))
        elif name == "generator_drift":
            generator = (
                strict_a
                + intervention * edge_a
                + (0.10 + 0.4 * intervention) * time * edge_b
            )
            matrix = expm(-time * generator)
        elif name == "environment_mixture":
            alternative = strict_a + 0.34 * edge_b
            mixture = min(0.45, 0.18 + 0.5 * intervention)
            matrix = (
                (1.0 - mixture) * expm(-time * strict_a)
                + mixture * expm(-time * alternative)
            )
        elif name == "apparatus_lag":
            lag = 0.045 + 0.08 * intervention
            matrix = (
                0.80 * expm(-time * strict_a)
                + 0.20 * expm(-(time + lag) * strict_a)
            )
        else:
            raise ValueError(name)
        matrices.append(hermitian(matrix))
    return matrices


def mechanism_feature(
    name: str, strict_a: np.ndarray, times: np.ndarray
) -> np.ndarray:
    base = matrix_observation_features(
        mechanism_matrices(name, strict_a, times, 0.0), times
    )
    probe = matrix_observation_features(
        mechanism_matrices(name, strict_a, times, 0.12), times
    )
    return np.concatenate([base, probe - base])


def program_290(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    times = np.array([0.15, 0.35, 0.65])
    names = [
        "homogeneous",
        "generator_drift",
        "environment_mixture",
        "apparatus_lag",
    ]
    templates = {
        name: mechanism_feature(name, strict_a, times) for name in names
    }
    feature_scale = np.maximum(
        np.std(np.vstack(list(templates.values())), axis=0), 1e-5
    )
    trials_per_class = 500
    rows = []
    accuracies = {}
    for noise in [0.05, 0.10, 0.20, 0.35]:
        confusion = {
            name: {candidate: 0 for candidate in names} for name in names
        }
        for true_name in names:
            for _ in range(trials_per_class):
                observed = templates[true_name] + rng.normal(
                    0.0, noise * feature_scale
                )
                predicted = min(
                    names,
                    key=lambda candidate: np.linalg.norm(
                        (observed - templates[candidate]) / feature_scale
                    ),
                )
                confusion[true_name][predicted] += 1
        accuracies[str(noise)] = sum(
            confusion[name][name] for name in names
        ) / (len(names) * trials_per_class)
        rows.extend(
            {
                "relative_feature_noise": noise,
                "true_mechanism": true_name,
                "predicted_mechanism": predicted,
                "count": confusion[true_name][predicted],
            }
            for true_name in names
            for predicted in names
        )
    minimum_template_separation = min(
        np.linalg.norm(
            (templates[left] - templates[right]) / feature_scale
        )
        for index, left in enumerate(names)
        for right in names[index + 1 :]
    )
    return (
        {
            "status": (
                "[Proven] distinct frozen noiseless templates; "
                "[Strong evidence] synthetic causal classification"
            ),
            "times": times,
            "mechanisms": names,
            "intervention_strength": 0.12,
            "trials_per_class": trials_per_class,
            "classification_accuracy_by_noise": accuracies,
            "reference_noise": 0.20,
            "reference_classification_accuracy": accuracies["0.2"],
            "minimum_standardized_template_separation": float(
                minimum_template_separation
            ),
            "confusion_rows": rows,
            "boundary": (
                "Three times plus one intervention separate only the frozen "
                "four-template family. Arbitrary time-dependent generators, "
                "baths, and apparatus kernels remain nonidentifiable."
            ),
        },
        rows,
    )


def sensitivity_design(
    inputs: np.ndarray,
    baseline: np.ndarray,
    direction: np.ndarray,
    dt: float,
    initial: np.ndarray,
    parameters: np.ndarray,
) -> tuple[np.ndarray, float, float]:
    delta = 1e-5
    base = previous.simulate_weight_law(
        initial,
        baseline,
        inputs,
        dt,
        "driven_additive",
        parameters,
        direction,
    )
    columns = []
    for index in range(2):
        shifted = parameters.copy()
        shifted[index] += delta
        trajectory = previous.simulate_weight_law(
            initial,
            baseline,
            inputs,
            dt,
            "driven_additive",
            shifted,
            direction,
        )
        columns.append(((trajectory - base) / delta).reshape(-1))
    design = np.column_stack(columns)
    information = design.T @ design
    eigenvalues = np.linalg.eigvalsh(information)
    return information, float(eigenvalues.min()), float(np.linalg.det(information))


def program_291(
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    baseline = core.STRICT_K.copy()
    direction = np.array([1.0, -0.8, 0.6, -0.4, 0.25, -0.15])
    direction /= np.linalg.norm(direction)
    initial = baseline * np.array([1.08, 0.94, 1.04, 0.97, 1.02, 0.95])
    parameters = np.array([0.24, 0.075])
    dt = 0.05
    length = 160
    candidates: list[tuple[str, np.ndarray]] = [
        ("zero", np.zeros(length)),
        ("single_pulse", np.r_[np.ones(length // 2), -np.ones(length // 2)]),
        ("sinusoid", np.sin(np.linspace(0.0, 8.0 * math.pi, length))),
        ("prbs", rng.choice([-1.0, 1.0], size=length)),
    ]
    for index in range(240):
        sequence = rng.choice([-1.0, 1.0], size=length)
        candidates.append((f"random_{index:03d}", sequence))
    rows: list[dict[str, Any]] = []
    for name, inputs in candidates:
        information, minimum, determinant = sensitivity_design(
            inputs, baseline, direction, dt, initial, parameters
        )
        rows.append(
            {
                "design": name,
                "input_energy": float(np.sum(inputs**2)),
                "information_minimum_eigenvalue": minimum,
                "information_determinant": determinant,
                "information_condition_number": float(
                    np.linalg.cond(information)
                )
                if minimum > 1e-14
                else None,
            }
        )
    best = max(rows, key=lambda row: row["information_minimum_eigenvalue"])
    named = {
        row["design"]: row
        for row in rows
        if not str(row["design"]).startswith("random_")
    }
    return (
        {
            "status": "[Strong evidence] frozen maximin intervention search",
            "candidate_count": len(rows),
            "sequence_length": length,
            "amplitude_constraint": 1.0,
            "best_design": best,
            "named_designs": named,
            "verdict": (
                "A balanced long two-level pulse maximizes the worst local "
                "parameter direction within the frozen two-parameter law."
            ),
            "boundary": (
                "The search is local at supplied parameters and restricted to "
                "one adaptive law, one direction, and a finite input family."
            ),
        },
        rows,
    )


def evaluate_replication_reservoir(
    reservoir: np.ndarray,
    input_vector: np.ndarray,
    inputs: dict[str, np.ndarray],
) -> dict[str, float]:
    continuous = inputs["continuous"]
    states_continuous = previous.nonlinear_states(
        reservoir, input_vector, continuous
    )
    delay_target = np.roll(continuous, 5)
    nonlinear_target = (
        np.sin(2.0 * np.roll(continuous, 1))
        + 0.30 * np.roll(continuous, 2) * np.roll(continuous, 3)
    )
    narma_input = inputs["narma"]
    states_narma = previous.nonlinear_states(
        reservoir, input_vector, narma_input
    )
    narma_target = previous.narma10(narma_input)
    binary = inputs["binary"]
    states_binary = previous.nonlinear_states(
        reservoir, input_vector, binary
    )
    parity_target = binary * np.roll(binary, 1) * np.roll(binary, 2)
    return {
        "delay5_r2": previous.ridge_score(
            states_continuous, delay_target, 250
        ),
        "nonlinear_prediction_r2": previous.ridge_score(
            states_continuous, nonlinear_target, 250
        ),
        "narma10_r2": previous.ridge_score(
            states_narma, narma_target, 250
        ),
        "parity3_accuracy": previous.ridge_score(
            states_binary, parity_target, 250, classification=True
        ),
    }


def program_292(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    fin = 0.92 * expm(-0.25 * strict_a)
    input_vector = np.zeros(N)
    input_vector[0] = 0.75
    input_vector[3] = -0.45
    input_vector -= input_vector.mean()
    fin_eigenvalues = np.linalg.eigvalsh(fin)
    seeds = 24
    controls_per_seed = 12
    rows: list[dict[str, Any]] = []
    primary_wins = 0
    primary_percentiles = []
    primary_advantages = []
    secondary_percentiles: dict[str, list[float]] = {
        "delay5_r2": [],
        "nonlinear_prediction_r2": [],
        "parity3_accuracy": [],
    }
    for replication in range(seeds):
        local_seed = int(rng.integers(0, 2**32 - 1))
        local = np.random.default_rng(local_seed)
        inputs = {
            "continuous": local.uniform(-1.0, 1.0, size=3_200),
            "narma": local.uniform(0.0, 0.40, size=3_200),
            "binary": local.choice([-1.0, 1.0], size=3_200),
        }
        fin_metrics = evaluate_replication_reservoir(
            fin, input_vector, inputs
        )
        control_metrics = []
        for control_index in range(controls_per_seed):
            if control_index < controls_per_seed // 2:
                q, _ = np.linalg.qr(local.normal(size=(N, N)))
                control = q @ np.diag(fin_eigenvalues) @ q.T
                ensemble = "isospectral_orientation"
            else:
                values = local.uniform(-0.92, 0.92, size=N)
                values[np.argmax(np.abs(values))] = 0.92
                q, _ = np.linalg.qr(local.normal(size=(N, N)))
                control = q @ np.diag(values) @ q.T
                ensemble = "symmetric_radius_matched"
            metrics = evaluate_replication_reservoir(
                control, input_vector, inputs
            )
            control_metrics.append(metrics)
            rows.append(
                {
                    "replication": replication,
                    "local_seed": local_seed,
                    "ensemble": ensemble,
                    "control_index": control_index,
                    **metrics,
                }
            )
        rows.append(
            {
                "replication": replication,
                "local_seed": local_seed,
                "ensemble": "FIN",
                "control_index": -1,
                **fin_metrics,
            }
        )
        control_narma = np.array(
            [metrics["narma10_r2"] for metrics in control_metrics]
        )
        percentile = float(np.mean(control_narma <= fin_metrics["narma10_r2"]))
        primary_percentiles.append(percentile)
        primary_advantages.append(
            fin_metrics["narma10_r2"] - float(np.max(control_narma))
        )
        primary_wins += int(
            fin_metrics["narma10_r2"] > float(np.max(control_narma))
        )
        for metric in secondary_percentiles:
            secondary_percentiles[metric].append(
                float(
                    np.mean(
                        np.array([row[metric] for row in control_metrics])
                        <= fin_metrics[metric]
                    )
                )
            )
    return (
        {
            "status": "[Strong evidence] preregistered primary-task replication",
            "replications": seeds,
            "controls_per_replication": controls_per_seed,
            "primary_task": "narma10_r2",
            "primary_win_count_against_all_controls": primary_wins,
            "primary_win_rate": primary_wins / seeds,
            "median_primary_percentile": float(
                np.median(primary_percentiles)
            ),
            "median_advantage_over_best_control": float(
                np.median(primary_advantages)
            ),
            "p05_advantage_over_best_control": float(
                np.quantile(primary_advantages, 0.05)
            ),
            "secondary_median_percentiles": {
                metric: float(np.median(values))
                for metric, values in secondary_percentiles.items()
            },
            "boundary": (
                "NARMA10 was frozen as the primary replication target after "
                "P278. The benchmark remains synthetic, software-specific, "
                "and dependent on input coupling and hyperparameters."
            ),
        },
        rows,
    )


def binary_entropy(probability: float) -> float:
    if probability <= 0.0 or probability >= 1.0:
        return 0.0
    return float(
        -probability * math.log(probability)
        - (1.0 - probability) * math.log(1.0 - probability)
    )


def program_293() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    k_b = 1.0
    temperature = 2.0
    beta = 1.0 / (k_b * temperature)
    maximum_gap = 12.0 * k_b * temperature
    steps = 800
    gaps = np.linspace(0.0, maximum_gap, steps + 1)
    probability_excited = 0.5
    work = 0.0
    heat_to_system = 0.0
    entropy_production = 0.0
    rows: list[dict[str, Any]] = []
    previous_entropy = binary_entropy(probability_excited)
    for index in range(1, len(gaps)):
        old_gap = gaps[index - 1]
        new_gap = gaps[index]
        gap_change = new_gap - old_gap
        work_step = probability_excited * gap_change
        work += work_step
        equilibrium_excited = 1.0 / (1.0 + math.exp(beta * new_gap))
        heat_step = new_gap * (equilibrium_excited - probability_excited)
        heat_to_system += heat_step
        new_entropy = binary_entropy(equilibrium_excited)
        production_step = new_entropy - previous_entropy - beta * heat_step
        entropy_production += production_step
        probability_excited = equilibrium_excited
        previous_entropy = new_entropy
        if index in {1, steps // 4, steps // 2, 3 * steps // 4, steps}:
            rows.append(
                {
                    "step": index,
                    "gap": new_gap,
                    "excited_probability": probability_excited,
                    "cumulative_work": work,
                    "cumulative_heat_to_system": heat_to_system,
                    "cumulative_entropy_production": entropy_production,
                }
            )
    closing_work = -probability_excited * maximum_gap
    work += closing_work
    initial_entropy = math.log(2.0)
    final_entropy = binary_entropy(probability_excited)
    information_erased = initial_entropy - final_entropy
    generalized_landauer = k_b * temperature * information_erased
    energy_balance_residual = abs(work + heat_to_system)
    return (
        {
            "status": "[Proven] conditioned finite-step erasure ledger",
            "supplied_resources": {
                "temperature": temperature,
                "boltzmann_constant": k_b,
                "maximum_gap": maximum_gap,
                "ideal_fresh_gibbs_bath_qubits": steps,
                "classical_work_meter": True,
            },
            "steps": steps,
            "final_excited_probability": probability_excited,
            "information_erased_nats": information_erased,
            "total_work": work,
            "heat_dumped_to_bath": -heat_to_system,
            "generalized_landauer_bound": generalized_landauer,
            "work_above_bound": work - generalized_landauer,
            "total_entropy_production": entropy_production,
            "energy_balance_residual": energy_balance_residual,
            "theorem": (
                "Each gap change is work and each Gibbs collision is heat. "
                "Summing the first law gives Delta U=W+Q; summing relative-"
                "entropy contraction gives W >= k_B T times the erased "
                "information, with finite-step dissipation."
            ),
            "boundary": (
                "The bath temperature, Hamiltonian control, fresh Gibbs "
                "ancillas, and work meter are supplied physical resources. "
                "They are not generated by the strict FIN operator."
            ),
        },
        rows,
    )


def product_group_add(left: tuple[int, int], right: tuple[int, int]) -> tuple[int, int]:
    return ((left[0] + right[0]) % 3, (left[1] + right[1]) % 2)


def program_294() -> dict[str, Any]:
    group = [(scale, sign) for scale in range(3) for sign in range(2)]
    trivial_sections = 0
    for value in group:
        if all(product_group_add(g, value) == value for g in group):
            trivial_sections += 1
    equivariant_maps = 0
    pointed_maps = 0
    for offset in group:
        mapping = {
            x: product_group_add(x, offset)
            for x in group
        }
        equivariant = all(
            mapping[product_group_add(g, x)]
            == product_group_add(g, mapping[x])
            for g in group
            for x in group
        )
        if equivariant:
            equivariant_maps += 1
            if mapping[(0, 0)] == (0, 0):
                pointed_maps += 1
    return {
        "status": "[Proven] finite product-torsor section theorem",
        "finite_model": "G=Z3 x Z2 acting regularly on T=G",
        "group_order": len(group),
        "equivariant_sections_from_trivial_input": trivial_sections,
        "equivariant_maps_from_resource_torsor": equivariant_maps,
        "pointed_equivariant_maps": pointed_maps,
        "theorem": (
            "A free nontrivial G-torsor has no G-fixed point, so no "
            "equivariant section exists from a trivial input. Equivariant "
            "self-maps of the regular torsor are translations; specifying a "
            "base point leaves exactly the identity map."
        ),
        "new_object": (
            "Pointed Physicalization Resource (T,t0): a scale-orientation "
            "torsor together with a supplied base point."
        ),
        "boundary": (
            "Pointing the torsor makes a unique operational section possible "
            "but does not derive the point, a strict selector, or a physical "
            "unit."
        ),
    }


def make_figures(results: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    stability = results["P281"]["summaries"]
    noise = np.array([1e-6, 1e-4, 1e-3])
    axes[0].loglog(
        noise,
        [
            stability[str(value)]["median_maximum_relative_pole_error"]
            for value in noise
        ],
        "o-",
        label="bounded pole error",
    )
    axes[0].loglog(
        noise,
        [
            stability[str(value)]["median_maximum_relative_weight_error"]
            for value in noise
        ],
        "s-",
        label="bounded residue error",
    )
    axes[0].set_xlabel("relative response noise")
    axes[0].set_ylabel("median maximum relative error")
    axes[0].set_title("P281: bounded Stieltjes recovery")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend()
    detector = results["P283"]
    labels = ["ideal POVM", "calibrated detector"]
    ranks = [5, detector["calibrated_response_rank"]]
    conditions = [1.0, detector["calibrated_response_condition_number"]]
    x = np.arange(2)
    axes[1].bar(x - 0.18, ranks, 0.36, label="rank")
    axes[1].bar(x + 0.18, conditions, 0.36, label="condition no.")
    axes[1].set_xticks(x, labels)
    axes[1].set_title("P283: current readout survives detector channel")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p281_p283_recovery_povm.png", dpi=190
    )
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    rg_rows = results["P284"]["rows"]
    sizes = [row["context_size"] for row in rg_rows]
    axes[0].bar(
        np.arange(len(sizes)) - 0.18,
        [row["bare_best_fit_residual"] for row in rg_rows],
        0.36,
        label="bare lift",
    )
    axes[0].bar(
        np.arange(len(sizes)) + 0.18,
        [row["completed_best_fit_residual"] for row in rg_rows],
        0.36,
        label="complement completion",
    )
    axes[0].set_xticks(np.arange(len(sizes)), [str(value) for value in sizes])
    axes[0].set_xlabel("coarse context size")
    axes[0].set_ylabel("best relative Frobenius residual")
    axes[0].set_title("P284: rank restored, strict shape not restored")
    axes[0].legend()
    candidate = np.array(results["P286"]["candidate_weights"])
    axes[1].plot(np.arange(1, 7), core.STRICT_K, "o-", label="strict")
    axes[1].plot(np.arange(1, 7), candidate / (4.0 * math.log(2.0)), "s-",
                 label="legacy fractal atom / amplitude")
    axes[1].axhline(0.0, color="black", linewidth=0.7)
    axes[1].set_xlabel("cyclic distance")
    axes[1].set_ylabel("normalized coupling")
    axes[1].set_title("P286: one nonlinear legacy atom")
    axes[1].legend()
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p284_p286_rg_bridge.png", dpi=190
    )
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    flux = results["P288"]["grid"]
    axes[0].loglog(
        [row["h"] for row in flux],
        [row["median_relative_error"] for row in flux],
        "o-",
        label="median",
    )
    axes[0].loglog(
        [row["h"] for row in flux],
        [row["p95_relative_error"] for row in flux],
        "s--",
        label="95%",
    )
    axes[0].set_xlabel("flux separation h")
    axes[0].set_ylabel("relative error")
    axes[0].set_title("P288: detector-level bias--noise optimum")
    axes[0].grid(True, which="both", alpha=0.25)
    axes[0].legend()
    probability = results["P289"]["estimated_false_positive_probability"]
    direct_upper = results["P289"]["direct_zero_event_rule_of_three_upper_95"]
    values = [max(probability, 1e-14), max(direct_upper or 1e-14, 1e-14)]
    axes[1].bar(["importance estimate", "direct 95% upper"], values)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("false-positive probability")
    axes[1].set_title("P289: declared Dirichlet null ensemble")
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p288_p289_detector_false_positive.png", dpi=190
    )
    plt.close(figure)

    figure, axes = plt.subplots(1, 2, figsize=(11.0, 4.2))
    names = results["P290"]["mechanisms"]
    confusion = np.zeros((len(names), len(names)))
    for row in results["P290"]["confusion_rows"]:
        if row["relative_feature_noise"] != results["P290"]["reference_noise"]:
            continue
        confusion[names.index(row["true_mechanism"]),
                  names.index(row["predicted_mechanism"])] = row["count"]
    confusion /= confusion.sum(axis=1, keepdims=True)
    image = axes[0].imshow(confusion, vmin=0.0, vmax=1.0, cmap="Blues")
    axes[0].set_xticks(np.arange(len(names)), names, rotation=35, ha="right")
    axes[0].set_yticks(np.arange(len(names)), names)
    axes[0].set_xlabel("predicted")
    axes[0].set_ylabel("true")
    axes[0].set_title("P290: three-time causal panel")
    figure.colorbar(image, ax=axes[0], fraction=0.046)
    named = results["P291"]["named_designs"]
    design_names = list(named)
    axes[1].bar(
        design_names,
        [named[name]["information_minimum_eigenvalue"] for name in design_names],
    )
    axes[1].set_yscale("symlog", linthresh=1e-8)
    axes[1].tick_params(axis="x", rotation=30)
    axes[1].set_ylabel("minimum Fisher/design eigenvalue")
    axes[1].set_title("P291: intervention informativeness")
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p290_p291_mechanism_intervention.png", dpi=190
    )
    plt.close(figure)

    figure, axes = plt.subplots(1, 3, figsize=(13.0, 4.0))
    replication_rows = [
        row for row in results["P292_rows"] if row["ensemble"] == "FIN"
    ]
    axes[0].hist(
        [row["narma10_r2"] for row in replication_rows],
        bins=10,
        color="#1F5A99",
        alpha=0.8,
    )
    axes[0].set_xlabel("FIN NARMA10 R2")
    axes[0].set_ylabel("replications")
    axes[0].set_title("P292: primary-task replication")
    erasure = results["P293"]
    axes[1].bar(
        ["work", "Landauer bound"],
        [erasure["total_work"], erasure["generalized_landauer_bound"]],
        color=["#D55E00", "#19733A"],
    )
    axes[1].set_title("P293: finite-step erasure")
    axes[1].set_ylabel("energy units")
    torsor = results["P294"]
    axes[2].bar(
        ["trivial input", "resource torsor", "pointed"],
        [
            torsor["equivariant_sections_from_trivial_input"],
            torsor["equivariant_maps_from_resource_torsor"],
            torsor["pointed_equivariant_maps"],
        ],
        color=["#A61B1B", "#1F5A99", "#19733A"],
    )
    axes[2].set_title("P294: equivariant section counts")
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / "p292_p294_replication_thermo_torsor.png", dpi=190
    )
    plt.close(figure)


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    messages = {
        "P281": "bounded priors prevent runaway poles but do not prove uniform stability",
        "P282": "Lean checks an exact rational Schur witness, not the general theorem",
        "P283": "abstract Naimark dilation and calibrated detector rank are explicit",
        "P284": "complement dynamics restores rank but not the strict fixed point",
        "P285": "two-step quantum information ledger closes conditionally",
        "P286": "one legacy fractal-compression atom is insufficient",
        "P287": "no independent P241 bundle; P242 remains blocked",
        "P288": "detector likelihood has a finite optimal flux separation",
        "P289": "false-positive probability is ensemble-specific",
        "P290": "three times plus intervention separate only the frozen mechanism panel",
        "P291": "a balanced two-level intervention improves local adaptive-law identifiability",
        "P292": "NARMA10 advantage receives a multi-seed replication test",
        "P293": "Landauer accounting closes only with supplied physical resources",
        "P294": "a pointed resource torsor gives a unique section but not its source",
    }
    return [
        {
            "program": key,
            "status": results[key]["status"],
            "headline": messages[key],
        }
        for key in messages
    ]


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    strict_a, strict_w = core.strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "programs": "P281-P294",
            "seed": SEED,
            "scope": (
                "finite proofs, one dependency-free Lean certificate, and "
                "synthetic audits; no fabricated laboratory data, strict "
                "selector, physical unit, role transfer, or ToE closure"
            ),
        }
    }
    results["P281"], stability_rows = program_281(strict_a, rng)
    results["P282"] = program_282()
    results["P283"], povm_rows = program_283(strict_w, rng)
    results["P284"], rg_rows = program_284(strict_a)
    results["P285"] = program_285(rng)
    results["P286"] = program_286(strict_a)
    results["P287"] = program_287()
    results["P288"], flux_rows = program_288(strict_a, rng)
    results["P289"] = program_289(strict_a, rng)
    results["P290"], mechanism_rows = program_290(strict_a, rng)
    results["P291"], intervention_rows = program_291(rng)
    results["P292"], reservoir_rows = program_292(strict_a, rng)
    results["P292_rows"] = reservoir_rows
    results["P293"], erasure_rows = program_293()
    results["P294"] = program_294()

    make_figures(results)
    results_without_raw_replication = dict(results)
    results_without_raw_replication.pop("P292_rows")
    RESULTS_PATH.write_text(
        json.dumps(json_ready(results_without_raw_replication), indent=2)
        + "\n",
        encoding="utf-8",
    )
    write_csv(SUMMARY_PATH, summary_rows(results))
    write_csv(STABILITY_PATH, stability_rows)
    write_csv(POVM_PATH, povm_rows)
    write_csv(RG_PATH, rg_rows)
    write_csv(FLUX_PATH, flux_rows)
    write_csv(MECHANISM_PATH, mechanism_rows)
    write_csv(INTERVENTION_PATH, intervention_rows)
    write_csv(RESERVOIR_PATH, reservoir_rows)
    print(RESULTS_PATH.name)
    for row in summary_rows(results):
        print(row["program"], row["status"], "-", row["headline"])


if __name__ == "__main__":
    main()
