#!/usr/bin/env python3
"""Execute FIN research Programs P309--P322.

The round tests only hypothesis classes that remain admissible after Release
10.27.  It keeps legacy, repaired legacy, and strict operators typed
separately; labels synthetic calculations; and leaves external evidence gates
closed when custody, hardware, or preregistration are absent.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from math import gcd
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm, null_space
from scipy.optimize import least_squares, linprog, minimize_scalar, nnls

import fin_programs_255_266 as core
import fin_programs_267_280 as p267
import fin_programs_281_294 as p281
import fin_programs_295_308 as p295


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_309_322_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_309_322_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_309_322_Summary.csv"
MINIMAX_PATH = ROOT / "FIN_Programs_309_322_Minimax_Stieltjes.csv"
LOSS_PATH = ROOT / "FIN_Programs_309_322_Lossy_Mesh.csv"
PARENT_PATH = ROOT / "FIN_Programs_309_322_Common_Parent.csv"
REGULAR_PATH = ROOT / "FIN_Programs_309_322_Regular_Variation.csv"
OPERATIONAL_PATH = ROOT / "FIN_Programs_309_322_Operational_Distance.csv"
SIGNED_PATH = ROOT / "FIN_Programs_309_322_Signed_Path.csv"
DRIFT_PATH = ROOT / "FIN_Programs_309_322_Detector_Drift.csv"
CHAMBER_PATH = ROOT / "FIN_Programs_309_322_Spectral_Chambers.csv"
CLOCK_PATH = ROOT / "FIN_Programs_309_322_Clock_Design.csv"
ROLE_PATH = ROOT / "FIN_Programs_309_322_Role_Naturality.csv"
FORMAL_SOURCE = ROOT / "FIN_Programs_309_322_Formal_Core.lean"

N = 12
SEED = 20260802
ALPHA_GEO = 4.0 * math.log(2.0)
WEINBERG_BENCHMARK = ALPHA_GEO / 12.0


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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def hermitian(matrix: np.ndarray) -> np.ndarray:
    return (matrix + matrix.T.conj()) / 2.0


def positive_sqrt(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(hermitian(matrix))
    return (
        eigenvectors * np.sqrt(np.maximum(eigenvalues, 0.0))
    ) @ eigenvectors.T.conj()


def lean_binary() -> Path:
    candidates = [
        ROOT / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean",
        Path.home() / ".elan/toolchains/leanprover--lean4---v4.28.0/bin/lean",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("pinned Lean 4.28 binary not found")


# ---------------------------------------------------------------------------
# P309: minimax resolution modulus and E-optimal shared probes
# ---------------------------------------------------------------------------


def profiled_pole_information(
    z_grid: np.ndarray, poles: np.ndarray, weights: np.ndarray
) -> np.ndarray:
    jacobian = p295.stieltjes_jacobian(z_grid, poles, weights)
    order = len(poles)
    pole_jacobian = jacobian[:, :order]
    residue_jacobian = jacobian[:, order:]
    nuisance_projector = residue_jacobian @ np.linalg.pinv(residue_jacobian)
    residualizer = np.eye(len(jacobian)) - nuisance_projector
    return hermitian(pole_jacobian.T @ residualizer @ pole_jacobian)


def program_309(
    strict_a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    atoms = p267.visible_pole_residue_data(strict_a)
    poles = np.array([pole for pole, _ in atoms])
    residues = [residue for _, residue in atoms]
    z_grid = np.geomspace(0.03, 8.0, 42)

    candidates = [np.eye(6)[:, index] for index in range(6)]
    candidates.extend(
        vector / np.linalg.norm(vector)
        for vector in rng.normal(size=(90, 6))
    )

    selected: list[np.ndarray] = []
    remaining = list(candidates)
    design_rows: list[dict[str, Any]] = []
    total_information = np.zeros((3, 3))
    for probe_count in range(1, 7):
        best_index = -1
        best_score = -math.inf
        best_information = None
        best_weights = None
        for index, probe in enumerate(remaining):
            weights = np.array(
                [
                    [
                        float(probe.T @ residue @ probe)
                        for residue in residues
                    ]
                ]
            )
            information = profiled_pole_information(
                z_grid, poles, weights
            )
            combined = total_information + information
            score = float(np.linalg.eigvalsh(combined)[0])
            if score > best_score:
                best_score = score
                best_index = index
                best_information = information
                best_weights = weights[0]
        assert best_information is not None and best_weights is not None
        selected.append(remaining.pop(best_index))
        total_information += best_information
        eigenvalues = np.linalg.eigvalsh(total_information)
        design_rows.append(
            {
                "row_type": "e_optimal_probe",
                "probe_count": probe_count,
                "minimum_profiled_pole_information": float(eigenvalues[0]),
                "maximum_profiled_pole_information": float(eigenvalues[-1]),
                "profiled_condition_number": float(
                    eigenvalues[-1] / eigenvalues[0]
                ),
                "new_probe_residues": best_weights,
            }
        )

    collision_rows: list[dict[str, Any]] = []
    sigma = 1.0e-3
    center = 1.5
    residue = 0.5
    deltas = np.geomspace(1.0e-4, 0.20, 36)
    base = 2.0 * residue / (z_grid + center)
    for probe_count in (1, 4, 16):
        for delta in deltas:
            split = (
                residue / (z_grid + center - delta)
                + residue / (z_grid + center + delta)
            )
            normalized_difference = (split - base) / (sigma * base)
            kl = float(
                0.5
                * probe_count
                * np.sum(normalized_difference**2)
            )
            pinsker_tv_bound = min(1.0, math.sqrt(max(kl, 0.0) / 2.0))
            bayes_error_lower = 0.5 * (1.0 - pinsker_tv_bound)
            collision_rows.append(
                {
                    "row_type": "collision_lower_bound",
                    "probe_count": probe_count,
                    "relative_noise": sigma,
                    "pole_half_separation": float(delta),
                    "gaussian_kl": kl,
                    "bayes_error_lower_bound": bayes_error_lower,
                }
            )
    small = [
        row
        for row in collision_rows
        if row["probe_count"] == 1
        and row["pole_half_separation"] <= 0.01
        and row["gaussian_kl"] > 0.0
    ]
    kl_exponent = float(
        np.polyfit(
            np.log([row["pole_half_separation"] for row in small]),
            np.log([row["gaussian_kl"] for row in small]),
            1,
        )[0]
    )
    resolutions = {}
    for probe_count in (1, 4, 16):
        chosen = [
            row
            for row in collision_rows
            if row["probe_count"] == probe_count
            and row["bayes_error_lower_bound"] >= 0.25
        ]
        resolutions[str(probe_count)] = (
            max(row["pole_half_separation"] for row in chosen)
            if chosen
            else 0.0
        )
    rows = design_rows + collision_rows
    return (
        {
            "status": (
                "[Proven] two-pole minimax obstruction; "
                "[Strong evidence] finite E-optimal probe design"
            ),
            "true_poles": poles,
            "selected_probe_count": len(selected),
            "e_optimal_rows": design_rows,
            "collision_kl_power_exponent": kl_exponent,
            "quarter_error_resolution_by_probe_count": resolutions,
            "theorem": (
                "For an equal-residue pair mu±delta, the response differs "
                "from the merged pole only at order delta^2; Gaussian KL is "
                "therefore O(P*delta^4/sigma^2). Le Cam/Pinsker yields a "
                "nonzero testing-error lower bound at "
                "delta=O(sqrt(sigma)*P^(-1/4)). Fixed finite probe count "
                "cannot give a uniform Lipschitz inverse."
            ),
            "boundary": (
                "The lower bound is a worst-case equal-residue subfamily. "
                "The E-optimal design is local to the frozen strict atoms and "
                "candidate probe pool."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P310: completion-of-square Schur theorem and exact block witnesses
# ---------------------------------------------------------------------------


def fraction_inverse(matrix: list[list[Fraction]]) -> list[list[Fraction]]:
    dimension = len(matrix)
    augmented = [
        row[:] + [Fraction(int(i == j)) for j in range(dimension)]
        for i, row in enumerate(matrix)
    ]
    for column in range(dimension):
        pivot = next(
            row
            for row in range(column, dimension)
            if augmented[row][column] != 0
        )
        augmented[column], augmented[pivot] = (
            augmented[pivot],
            augmented[column],
        )
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(dimension):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                left - factor * right
                for left, right in zip(
                    augmented[row], augmented[column]
                )
            ]
    return [row[dimension:] for row in augmented]


def fraction_matmul(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [
            sum(
                (left[i][k] * right[k][j] for k in range(len(right))),
                Fraction(0),
            )
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def fraction_transpose(
    matrix: list[list[Fraction]],
) -> list[list[Fraction]]:
    return [list(row) for row in zip(*matrix)]


def fraction_add(
    left: list[list[Fraction]], right: list[list[Fraction]]
) -> list[list[Fraction]]:
    return [
        [a + b for a, b in zip(row_left, row_right)]
        for row_left, row_right in zip(left, right)
    ]


def fraction_quadratic(
    matrix: list[list[Fraction]], vector: list[Fraction]
) -> Fraction:
    return sum(
        (
            vector[i] * matrix[i][j] * vector[j]
            for i in range(len(vector))
            for j in range(len(vector))
        ),
        Fraction(0),
    )


def random_spd_fraction(
    dimension: int, rng: np.random.Generator
) -> list[list[Fraction]]:
    lower = [
        [
            Fraction(int(rng.integers(-3, 4))) if j <= i else Fraction(0)
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]
    for i in range(dimension):
        lower[i][i] = Fraction(int(rng.integers(1, 5)))
    product = fraction_matmul(lower, fraction_transpose(lower))
    return [
        [
            product[i][j] + (Fraction(1) if i == j else Fraction(0))
            for j in range(dimension)
        ]
        for i in range(dimension)
    ]


def program_310(rng: np.random.Generator) -> dict[str, Any]:
    exact_checks = 0
    for _ in range(300):
        n = int(rng.integers(1, 5))
        m = int(rng.integers(1, 4))
        reduced = random_spd_fraction(n, rng)
        inner = random_spd_fraction(m, rng)
        inner_inverse = fraction_inverse(inner)
        coupling = [
            [Fraction(int(rng.integers(-3, 4))) for _ in range(m)]
            for _ in range(n)
        ]
        correction = fraction_matmul(
            fraction_matmul(coupling, inner_inverse),
            fraction_transpose(coupling),
        )
        upper = fraction_add(reduced, correction)
        x = [Fraction(int(rng.integers(-3, 4))) for _ in range(n)]
        y = [Fraction(int(rng.integers(-3, 4))) for _ in range(m)]
        coupling_t_x = [
            sum(
                (coupling[i][j] * x[i] for i in range(n)),
                Fraction(0),
            )
            for j in range(m)
        ]
        shift = [
            sum(
                (
                    inner_inverse[i][j] * coupling_t_x[j]
                    for j in range(m)
                ),
                Fraction(0),
            )
            for i in range(m)
        ]
        completed = [y[i] + shift[i] for i in range(m)]
        full_value = (
            fraction_quadratic(upper, x)
            + 2
            * sum(
                (
                    x[i] * coupling[i][j] * y[j]
                    for i in range(n)
                    for j in range(m)
                ),
                Fraction(0),
            )
            + fraction_quadratic(inner, y)
        )
        decomposed = fraction_quadratic(reduced, x) + fraction_quadratic(
            inner, completed
        )
        if full_value != decomposed:
            raise AssertionError("exact block completion identity failed")
        exact_checks += 1

    process = subprocess.run(
        [str(lean_binary()), str(FORMAL_SOURCE)],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=False,
        timeout=90,
    )
    mathlib_available = (ROOT / ".lake/packages/mathlib").exists()
    return {
        "status": (
            "[Proven] general completion-square implication and exact block "
            "identities; [Blocked] full Mathlib matrix formalization"
        ),
        "exact_rational_block_checks": exact_checks,
        "lean_exit_code": process.returncode,
        "lean_stdout": process.stdout.strip(),
        "lean_stderr": process.stderr.strip(),
        "lean_source": FORMAL_SOURCE.name,
        "lean_source_sha256": sha256_file(FORMAL_SOURCE),
        "mathlib_available_in_repository_environment": mathlib_available,
        "theorem": (
            "If a block quadratic form has the exact Schur "
            "completion-of-square decomposition and the remainder is "
            "nonnegative, the Schur critical point is a global inner "
            "minimizer. Lean proves this implication for an abstract ordered "
            "additive value type."
        ),
        "boundary": (
            "The matrix expansion was checked exactly for 300 rational "
            "blocks. The local environment does not contain a repository "
            "Mathlib package proving the matrix identity and positivity from "
            "primitive finite-dimensional assumptions."
        ),
    }


# ---------------------------------------------------------------------------
# P311: loss/drift-completed complex-unitary current measurement
# ---------------------------------------------------------------------------


def reconstruct_perturbed_unitary(
    rotations: list[tuple[int, int, float, float, float]],
    diagonal: np.ndarray,
    parameter_sigma: float,
    rng: np.random.Generator,
) -> np.ndarray:
    phases = np.angle(np.diag(diagonal))
    if parameter_sigma:
        phases = phases + rng.normal(0.0, parameter_sigma, len(phases))
    work = np.diag(np.exp(1j * phases))
    for first, second, theta, phase_a, phase_b in reversed(rotations):
        if parameter_sigma:
            theta += float(rng.normal(0.0, parameter_sigma))
            phase_a += float(rng.normal(0.0, parameter_sigma))
            phase_b += float(rng.normal(0.0, parameter_sigma))
        cosine = math.cos(theta)
        sine = math.sin(theta)
        givens = np.array(
            [
                [
                    np.exp(-1j * phase_a) * cosine,
                    np.exp(-1j * phase_b) * sine,
                ],
                [
                    -np.exp(1j * phase_b) * sine,
                    np.exp(1j * phase_a) * cosine,
                ],
            ],
            dtype=complex,
        )
        row_first, row_second = givens.conj().T @ np.vstack(
            [work[first], work[second]]
        )
        work[first] = row_first
        work[second] = row_second
    return work


def detector_confusion_matrix() -> np.ndarray:
    efficiency = 0.86
    cross_talk = 0.012
    dark = 3.0e-4
    confusion = np.zeros((7, 7))
    for true_y in range(6):
        confusion[true_y, true_y] += efficiency * (1.0 - 2.0 * cross_talk)
        confusion[(true_y - 1) % 6, true_y] += efficiency * cross_talk
        confusion[(true_y + 1) % 6, true_y] += efficiency * cross_talk
        confusion[6, true_y] = 1.0 - efficiency
    confusion[:6, 6] = dark / 6.0
    confusion[6, 6] = 1.0 - dark
    return confusion


def effects_from_lossy_unitary(
    unitary: np.ndarray, transmissions: np.ndarray
) -> list[np.ndarray]:
    isometry = unitary[:, :N]
    effects = []
    for outcome in range(6):
        block = isometry[outcome * N : (outcome + 1) * N]
        attenuation = np.diag(
            transmissions[outcome * N : (outcome + 1) * N]
        )
        effects.append(hermitian(block.conj().T @ attenuation @ block))
    no_click = hermitian(np.eye(N) - sum(effects))
    return effects + [no_click]


def apply_confusion(
    effects: list[np.ndarray], confusion: np.ndarray
) -> list[np.ndarray]:
    return [
        hermitian(
            sum(
                confusion[observed, true] * effects[true]
                for true in range(7)
            )
        )
        for observed in range(7)
    ]


def program_311(
    strict_w: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    basis, ideal_effects_six, _, _ = p281.current_povm(strict_w)
    isometry = np.vstack(
        [p281.positive_sqrt(effect) for effect in ideal_effects_six]
    )
    complement = null_space(isometry.conj().T)
    ideal_unitary = np.column_stack([isometry, complement])
    rotations, diagonal = p295.givens_decompose_unitary(ideal_unitary)
    confusion = detector_confusion_matrix()
    ideal_effects = ideal_effects_six + [np.zeros((N, N))]
    ideal_observed = apply_confusion(ideal_effects, confusion)

    rows: list[dict[str, Any]] = []
    for mean_loss in (0.0, 0.01, 0.03, 0.08):
        for phase_sigma in (0.0, 1.0e-4, 1.0e-3):
            for repetition in range(12):
                perturbed = reconstruct_perturbed_unitary(
                    rotations, diagonal, phase_sigma, rng
                )
                if mean_loss == 0.0:
                    transmissions = np.ones(72)
                else:
                    transmissions = np.clip(
                        1.0
                        - rng.lognormal(
                            mean=math.log(mean_loss),
                            sigma=0.35,
                            size=72,
                        ),
                        0.0,
                        1.0,
                    )
                true_effects = effects_from_lossy_unitary(
                    perturbed, transmissions
                )
                observed_effects = apply_confusion(true_effects, confusion)
                response = np.array(
                    [
                        [
                            np.trace(effect @ operator).real
                            for operator in basis
                        ]
                        for effect in observed_effects
                    ]
                )
                singular_values = np.linalg.svd(
                    response, compute_uv=False
                )
                completeness = sum(observed_effects)
                minimum_eigenvalue = min(
                    float(np.linalg.eigvalsh(effect)[0])
                    for effect in observed_effects
                )
                rows.append(
                    {
                        "mean_mode_loss": mean_loss,
                        "phase_parameter_sigma": phase_sigma,
                        "repetition": repetition,
                        "response_rank": int(
                            np.linalg.matrix_rank(response, tol=1.0e-10)
                        ),
                        "response_minimum_singular_value": float(
                            singular_values[-1]
                        ),
                        "completeness_residual": float(
                            np.linalg.norm(completeness - np.eye(N), 2)
                        ),
                        "minimum_effect_eigenvalue": minimum_eigenvalue,
                        "maximum_effect_error_from_ideal_detector_model": float(
                            max(
                                np.linalg.norm(actual - reference, 2)
                                for actual, reference in zip(
                                    observed_effects, ideal_observed
                                )
                            )
                        ),
                        "mean_realized_transmission": float(
                            np.mean(transmissions)
                        ),
                    }
                )
    grouped = []
    for mean_loss in (0.0, 0.01, 0.03, 0.08):
        for phase_sigma in (0.0, 1.0e-4, 1.0e-3):
            selected = [
                row
                for row in rows
                if row["mean_mode_loss"] == mean_loss
                and row["phase_parameter_sigma"] == phase_sigma
            ]
            grouped.append(
                {
                    "mean_mode_loss": mean_loss,
                    "phase_parameter_sigma": phase_sigma,
                    "minimum_rank": min(row["response_rank"] for row in selected),
                    "median_response_minimum_singular_value": float(
                        np.median(
                            [
                                row["response_minimum_singular_value"]
                                for row in selected
                            ]
                        )
                    ),
                    "p95_maximum_effect_error": float(
                        np.quantile(
                            [
                                row[
                                    "maximum_effect_error_from_ideal_detector_model"
                                ]
                                for row in selected
                            ],
                            0.95,
                        )
                    ),
                }
            )
    return (
        {
            "status": (
                "[Proven] loss-completed POVM algebra; "
                "[Strong evidence] component-drift tolerance audit"
            ),
            "unitary_dimension": 72,
            "rotation_count": len(rotations),
            "detector_confusion_column_residual": float(
                np.max(np.abs(confusion.sum(axis=0) - 1.0))
            ),
            "grouped_results": grouped,
            "all_tested_response_ranks": sorted(
                {row["response_rank"] for row in rows}
            ),
            "maximum_completeness_residual": max(
                row["completeness_residual"] for row in rows
            ),
            "minimum_effect_eigenvalue": min(
                row["minimum_effect_eigenvalue"] for row in rows
            ),
            "boundary": (
                "Loss, phase drift, efficiency, crosstalk, and dark counts "
                "are supplied synthetic distributions. The no-click effect "
                "closes the algebraic POVM but no device was fabricated."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P312: cross-supported joint-feature parent across three carriers
# ---------------------------------------------------------------------------


def parent_for_commuting_psd(
    legacy: np.ndarray, strict: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    legacy_sqrt = positive_sqrt(legacy)
    strict_sqrt = positive_sqrt(strict)
    cross = legacy_sqrt @ strict_sqrt
    parent = hermitian(
        np.block([[legacy, cross], [cross.conj().T, strict]])
    )
    return parent, cross


def positive_projective_spectrum(matrix: np.ndarray) -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh(hermitian(matrix))
    positive = eigenvalues[eigenvalues > 1.0e-9]
    return positive / positive[-1]


def program_312() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    training_coefficients = None
    for size in (12, 16, 24):
        legacy_signed = p295.laplacian_from_profile(
            p295.legacy_weights_for_size(size)
        )
        legacy_abs = p295.laplacian_from_profile(
            np.abs(p295.legacy_weights_for_size(size))
        )
        strict = p295.laplacian_from_profile(
            p295.strict_weights_for_size(size)
        )
        parent, cross = parent_for_commuting_psd(legacy_abs, strict)
        legacy_spectrum = positive_projective_spectrum(legacy_abs)
        strict_spectrum = positive_projective_spectrum(strict)
        feature_matrix = np.column_stack(
            [
                legacy_spectrum,
                legacy_spectrum**2,
                legacy_spectrum**3,
            ]
        )
        if size == 12:
            training_coefficients = np.linalg.lstsq(
                feature_matrix, strict_spectrum, rcond=None
            )[0]
        assert training_coefficients is not None
        predicted = feature_matrix @ training_coefficients
        rows.append(
            {
                "carrier_size": size,
                "signed_legacy_minimum_eigenvalue": float(
                    np.linalg.eigvalsh(legacy_signed)[0]
                ),
                "parent_minimum_eigenvalue": float(
                    np.linalg.eigvalsh(parent)[0]
                ),
                "parent_rank": int(
                    np.linalg.matrix_rank(parent, tol=1.0e-9)
                ),
                "direct_sum_rank": int(
                    np.linalg.matrix_rank(legacy_abs, tol=1.0e-9)
                    + np.linalg.matrix_rank(strict, tol=1.0e-9)
                ),
                "legacy_compression_residual": float(
                    np.linalg.norm(parent[:size, :size] - legacy_abs, 2)
                ),
                "strict_compression_residual": float(
                    np.linalg.norm(parent[size:, size:] - strict, 2)
                ),
                "normalized_cross_support": float(
                    np.linalg.norm(cross, "fro")
                    / math.sqrt(
                        np.linalg.norm(legacy_abs, "fro")
                        * np.linalg.norm(strict, "fro")
                    )
                ),
                "c12_trained_cubic_spectral_rmse": float(
                    math.sqrt(np.mean((predicted - strict_spectrum) ** 2))
                ),
                "c12_trained_cubic_spectral_max_error": float(
                    np.max(np.abs(predicted - strict_spectrum))
                ),
            }
        )
    return (
        {
            "status": (
                "[Proven] nontrivial joint-feature parent after explicit "
                "legacy positivity repair; [Refuted] source/completion theorem"
            ),
            "rows": rows,
            "c12_trained_cubic_coefficients": training_coefficients,
            "construction": (
                "P=[[A_L, A_L^(1/2)A_S^(1/2)], "
                "[A_S^(1/2)A_L^(1/2), A_S]]. For commuting PSD circulants "
                "this is the Gram matrix of a shared mode-feature map and "
                "has nonzero cross support."
            ),
            "boundary": (
                "The parent imports the strict square root and an absolute-"
                "value repair of legacy; existence is therefore not a "
                "legacy-to-strict source law. A cubic map trained on C12 "
                "deteriorates on C16 and C24."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P313: regular-variation and positive hyperbolic-mixture obstruction
# ---------------------------------------------------------------------------


def program_313() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    beta_grid = np.geomspace(1.0e-4, 1.0e4, 241)
    train_d = np.arange(1, 7, dtype=float)
    target_train = 1.0 / (1.0 + train_d**p295.STRICT_ETA)
    design = 1.0 / (1.0 + train_d[:, None] * beta_grid[None, :])
    coefficients, residual = nnls(design, target_train)
    fitted_train = design @ coefficients
    tail_d = np.geomspace(20.0, 1.0e5, 180)
    fitted_tail = (
        1.0 / (1.0 + tail_d[:, None] * beta_grid[None, :])
    ) @ coefficients
    strict_tail = 1.0 / (1.0 + tail_d**p295.STRICT_ETA)
    fitted_slope = float(
        np.polyfit(np.log(tail_d[-60:]), np.log(fitted_tail[-60:]), 1)[0]
    )
    strict_slope = float(
        np.polyfit(np.log(tail_d[-60:]), np.log(strict_tail[-60:]), 1)[0]
    )
    rows = [
        {
            "row_type": "path_counting",
            "path_count_exponent_p": 1.6,
            "path_weight_exponent_q": q,
            "resulting_regular_variation_index": 1.6 - q,
            "interpretation": label,
        }
        for q, label in (
            (0.6, "historical stated pair: growth d^1"),
            (2.6, "required for legacy d^-1"),
            (3.4, "required for strict envelope d^-1.8"),
        )
    ]
    rows.extend(
        {
            "row_type": "tail_curve",
            "distance": float(distance),
            "positive_hyperbolic_mixture": float(mixture),
            "strict_envelope": float(strict_value),
        }
        for distance, mixture, strict_value in zip(
            tail_d, fitted_tail, strict_tail
        )
    )
    return (
        {
            "status": (
                "[Proven] positive hyperbolic-mixture tail obstruction; "
                "[Strong evidence] finite-fit extrapolation audit"
            ),
            "positive_nnls_train_relative_residual": float(
                np.linalg.norm(fitted_train - target_train)
                / np.linalg.norm(target_train)
            ),
            "positive_nnls_active_scales": int(
                np.sum(coefficients > 1.0e-10)
            ),
            "fitted_tail_log_slope": fitted_slope,
            "strict_tail_log_slope": strict_slope,
            "theorem": (
                "For every nonzero positive measure mu on beta>0, "
                "h(d)=int(1+beta*d)^(-1)dmu(beta) is not o(d^-1): some "
                "compact beta interval has positive mass and supplies a "
                "c/d lower bound. The strict envelope is asymptotic to "
                "d^-1.8=o(d^-1), so no positive mixture of legacy "
                "hyperbolic scales can equal it globally."
            ),
            "boundary": (
                "Signed measures and carrier-dependent cancellations evade "
                "the positive-mixture theorem and are quantified in P316."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P314: minimal nontorsion phase extension
# ---------------------------------------------------------------------------


def program_314() -> dict[str, Any]:
    strict_omega = Fraction("0.18575")
    strict_phi = Fraction("0.16250")
    denominator = math.lcm(strict_omega.denominator, strict_phi.denominator)
    omega_units = strict_omega.numerator * (
        denominator // strict_omega.denominator
    )
    phi_units = strict_phi.numerator * (
        denominator // strict_phi.denominator
    )
    common_gcd = gcd(abs(omega_units), abs(phi_units))
    base_numerator = common_gcd
    base_denominator = denominator
    legacy_phase_group_order = math.lcm(8, 12)
    return {
        "status": (
            "[Proven] torsion-to-nontorsion phase obstruction and minimal "
            "one-generator extension for the frozen decimal tuple"
        ),
        "legacy_generated_phase_group_order": legacy_phase_group_order,
        "strict_omega_exact_decimal_fraction": str(strict_omega),
        "strict_phi_exact_decimal_fraction": str(strict_phi),
        "strict_common_angular_unit_radians": (
            f"{base_numerator}/{base_denominator}"
        ),
        "strict_omega_integer_units": omega_units // common_gcd,
        "strict_phi_integer_units": phi_units // common_gcd,
        "integer_unit_gcd": gcd(
            omega_units // common_gcd, phi_units // common_gcd
        ),
        "minimal_extension": (
            "adjoin r=exp(i/4000); then exp(i*omega_S)=r^743 and "
            "exp(i*phi_S)=r^650"
        ),
        "theorem": (
            "The legacy phase subgroup generated by exp(i*pi/4) and "
            "exp(i*pi/6) is finite of order 24. Because pi is "
            "transcendental, exp(i/4000) has infinite order. A "
            "torsion-preserving homomorphism cannot map the finite legacy "
            "phase group onto the strict generator."
        ),
        "boundary": (
            "The exact 1/4000-radian generator depends on treating the frozen "
            "decimal tuple as exact. It is a target-coded extension resource, "
            "not a strict internal source or selector."
        ),
    }


# ---------------------------------------------------------------------------
# P315: operational legacy/strict prediction pseudometric
# ---------------------------------------------------------------------------


def prediction_distance_witness(
    strict_a: np.ndarray,
    legacy_a: np.ndarray,
    times: np.ndarray,
    dynamics: str,
    scale: float,
) -> tuple[float, dict[str, Any], list[dict[str, Any]]]:
    maximum = -1.0
    witness: dict[str, Any] = {}
    rows = []
    for time in times:
        if dynamics == "heat":
            strict_transition = expm(-time * strict_a)
            legacy_transition = expm(-time * scale * legacy_a)
        else:
            strict_transition = np.abs(expm(-1j * time * strict_a)) ** 2
            legacy_transition = (
                np.abs(expm(-1j * time * scale * legacy_a)) ** 2
            )
        time_max = 0.0
        for vertex in range(N):
            difference = strict_transition[vertex] - legacy_transition[vertex]
            tv = 0.5 * float(np.sum(np.abs(difference)))
            time_max = max(time_max, tv)
            if tv > maximum:
                maximum = tv
                witness = {
                    "time": float(time),
                    "preparation_vertex": vertex,
                    "positive_event_vertices": np.flatnonzero(
                        difference > 0.0
                    ).tolist(),
                    "event_probability_strict": float(
                        np.sum(strict_transition[vertex][difference > 0.0])
                    ),
                    "event_probability_legacy": float(
                        np.sum(legacy_transition[vertex][difference > 0.0])
                    ),
                }
        rows.append(
            {
                "dynamics": dynamics,
                "time": float(time),
                "maximum_vertex_tv": time_max,
            }
        )
    return maximum, witness, rows


def program_315(
    strict_a: np.ndarray,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    legacy_abs = p295.laplacian_from_profile(
        np.abs(p295.legacy_weights())
    )
    legacy_signed = p295.laplacian_from_profile(p295.legacy_weights())
    times = np.geomspace(0.02, 3.0, 42)
    results = {}
    all_rows: list[dict[str, Any]] = []
    for dynamics, legacy in (
        ("heat", legacy_abs),
        ("wave", legacy_signed),
    ):
        def objective(log_scale: float) -> float:
            distance, _, _ = prediction_distance_witness(
                strict_a,
                legacy,
                times,
                dynamics,
                math.exp(log_scale),
            )
            return distance

        optimum = minimize_scalar(
            objective,
            bounds=(-6.0, 4.0),
            method="bounded",
            options={"xatol": 1.0e-10},
        )
        scale = float(math.exp(optimum.x))
        distance, witness, rows = prediction_distance_witness(
            strict_a, legacy, times, dynamics, scale
        )
        for row in rows:
            row["optimized_legacy_scale"] = scale
        all_rows.extend(rows)
        alpha = 0.05
        conservative_shots = math.ceil(
            2.0 * math.log(2.0 / alpha) / distance**2
        )
        results[dynamics] = {
            "optimized_legacy_scale": scale,
            "finite_menu_operational_distance": distance,
            "witness": witness,
            "conservative_95pct_shot_bound": conservative_shots,
        }
    return (
        {
            "status": (
                "[Proven] finite-menu operational lower bounds; "
                "[Refuted] operational equivalence after scale alignment"
            ),
            "menu": {
                "preparations": "all 12 vertex states",
                "measurements": "12-outcome vertex measurement",
                "times": times,
                "heat_legacy": "absolute-value positivity repair",
                "wave_legacy": "historical signed self-adjoint generator",
            },
            "results": results,
            "boundary": (
                "The reported supremum is exact only on the declared finite "
                "time/preparation/measurement menu and is a lower bound on a "
                "larger operational or diamond/comb distance."
            ),
        },
        all_rows,
    )


# ---------------------------------------------------------------------------
# P316: minimum signed path-scale negativity
# ---------------------------------------------------------------------------


def program_316() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    distances = np.arange(7, dtype=float)
    target = 1.0 / (1.0 + distances**p295.STRICT_ETA)
    fourth_difference = float(np.diff(target, n=4)[0])
    analytic_lower_bound = max(0.0, -fourth_difference)
    rows = []
    best_solution = None
    for grid_size in (61, 121, 241):
        moment_nodes = np.linspace(0.0, 1.0, grid_size)
        moment_matrix = moment_nodes[None, :] ** distances[:, None]
        objective = np.concatenate(
            [np.zeros(grid_size), np.ones(grid_size)]
        )
        equality = np.column_stack([moment_matrix, -moment_matrix])
        result = linprog(
            objective,
            A_eq=equality,
            b_eq=target,
            bounds=[(0.0, None)] * (2 * grid_size),
            method="highs",
        )
        if not result.success:
            raise AssertionError(result.message)
        positive = result.x[:grid_size]
        negative = result.x[grid_size:]
        row = {
            "moment_grid_size": grid_size,
            "minimum_grid_negative_mass": float(result.fun),
            "positive_mass": float(np.sum(positive)),
            "negative_mass": float(np.sum(negative)),
            "moment_residual": float(
                np.max(np.abs(equality @ result.x - target))
            ),
            "analytic_continuum_negative_mass_lower_bound": (
                analytic_lower_bound
            ),
            "active_positive_atoms": int(np.sum(positive > 1.0e-10)),
            "active_negative_atoms": int(np.sum(negative > 1.0e-10)),
        }
        rows.append(row)
        if grid_size == 241:
            best_solution = (
                moment_nodes,
                positive,
                negative,
            )
    assert best_solution is not None
    nodes, positive, negative = best_solution
    atoms = [
        {
            "moment_node_t": float(nodes[index]),
            "signed_weight": float(positive[index] - negative[index]),
        }
        for index in range(len(nodes))
        if positive[index] > 1.0e-9 or negative[index] > 1.0e-9
    ]
    return (
        {
            "status": (
                "[Proven] positive path cone exclusion and analytic "
                "negativity lower bound; [Strong evidence] discretized "
                "minimum signed resource"
            ),
            "strict_moments_d0_to_d6": target,
            "fourth_forward_difference_at_zero": fourth_difference,
            "analytic_negative_mass_lower_bound": analytic_lower_bound,
            "grid_results": rows,
            "finest_grid_signed_atoms": atoms,
            "theorem": (
                "Every positive exponential-distance mixture is a Hausdorff "
                "moment sequence and has Delta^4 h(0)>=0. Strict has "
                "Delta^4 h(0)=-0.0715275581. For a normalized signed measure, "
                "|(t-1)^4|<=1 gives negative mass at least 0.0715275581."
            ),
            "boundary": (
                "The approximately 0.4067 optimum is exact only for each "
                "declared finite node grid; 0.07153 is the rigorous continuum "
                "lower bound. The signed measure is a resource certificate, "
                "not a physical mechanism."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P317: independent laboratory bundle gate
# ---------------------------------------------------------------------------


def external_gate(label: str) -> dict[str, Any]:
    manifests = sorted(ROOT.rglob("bundle_manifest.json"))
    admitted = []
    for path in manifests:
        relative = str(path.relative_to(ROOT))
        if "template" in relative.lower() or "example" in relative.lower():
            continue
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        required = {
            "provider",
            "registrar",
            "analyst",
            "holdout_hash",
            "event_file_hash",
        }
        if isinstance(payload, dict) and required.issubset(payload):
            admitted.append(relative)
    return {
        "status": (
            f"[Blocked by external evidence] no admitted {label} package"
        ),
        "manifest_count": len(manifests),
        "admitted_manifest_paths": admitted,
        "admitted": bool(admitted),
        "one_shot_pipeline_authorized": False,
        "boundary": (
            "Repository code cannot create independent provider/registrar/"
            "analyst roles, custody, apparatus calibration, or a previously "
            "frozen external hold-out."
        ),
    }


def program_317() -> dict[str, Any]:
    result = external_gate("P241 laboratory event")
    result["historical_internal_simulations_are_external_data"] = False
    return result


# ---------------------------------------------------------------------------
# P318: drift-aware sequential detector filter
# ---------------------------------------------------------------------------


def detector_jacobian(
    theta: np.ndarray, design: tuple[str, float]
) -> np.ndarray:
    steps = np.array([1.0e-5, 1.0e-6, 1.0e-6, 1.0e-7, 1.0e-7])
    jacobian = np.zeros((2, 5))
    for index, step in enumerate(steps):
        plus = theta.copy()
        minus = theta.copy()
        plus[index] += step
        minus[index] -= step
        jacobian[:, index] = (
            p295.detector_probabilities(plus, design)[:2]
            - p295.detector_probabilities(minus, design)[:2]
        ) / (2.0 * step)
    return jacobian


def ekf_update(
    mean: np.ndarray,
    covariance: np.ndarray,
    design: tuple[str, float],
    shots: int,
    counts: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    probability = p295.detector_probabilities(mean, design)
    predicted = probability[:2]
    observation = counts[:2] / shots
    jacobian = detector_jacobian(mean, design)
    covariance_observation = (
        np.diag(predicted) - np.outer(predicted, predicted)
    ) / shots
    innovation = (
        jacobian @ covariance @ jacobian.T
        + covariance_observation
        + 1.0e-10 * np.eye(2)
    )
    gain = covariance @ jacobian.T @ np.linalg.inv(innovation)
    updated_mean = mean + gain @ (observation - predicted)
    updated_covariance = hermitian(
        (np.eye(5) - gain @ jacobian) @ covariance
    )
    lower = np.array([0.05, 0.10, 0.10, 1.0e-7, 0.0])
    upper = np.array([2.50, 0.999, 0.999, 0.05, 0.20])
    return np.clip(updated_mean, lower, upper), updated_covariance


def run_detector_filter(
    rng: np.random.Generator,
    calibration_interval: int,
    dynamic: bool,
) -> tuple[float, float, float, list[dict[str, Any]]]:
    blocks = 70
    true = np.array([0.73, 0.81, 0.88, 0.0025, 0.035])
    mean = np.array([0.70, 0.79, 0.85, 0.0030, 0.030])
    covariance = np.diag([0.04, 0.02, 0.02, 3.0e-5, 0.003]) ** 2
    process_sd = np.array([0.0018, 0.0008, 0.0006, 2.0e-5, 0.00015])
    process_covariance = np.diag(process_sd**2) if dynamic else np.zeros((5, 5))
    errors = []
    covered = []
    rows = []
    calibration_designs = [
        ("dark", 0.0),
        ("crosstalk", 0.0),
        ("contrast", 0.0),
    ]
    for block in range(blocks):
        true += rng.normal(0.0, process_sd)
        true = np.clip(
            true,
            np.array([0.05, 0.10, 0.10, 1.0e-7, 0.0]),
            np.array([2.50, 0.999, 0.999, 0.05, 0.20]),
        )
        covariance = covariance + process_covariance
        designs = [("physics", 1.375), ("physics", -1.375)]
        if block % calibration_interval == 0:
            designs = calibration_designs + designs
        for design in designs:
            shots = 3500 if design[0] == "physics" else 2500
            counts = rng.multinomial(
                shots, p295.detector_probabilities(true, design)
            )
            mean, covariance = ekf_update(
                mean, covariance, design, shots, counts
            )
        error = float(mean[0] - true[0])
        sd = math.sqrt(max(float(covariance[0, 0]), 0.0))
        errors.append(error)
        covered.append(abs(error) <= 1.96 * sd)
        rows.append(
            {
                "calibration_interval": calibration_interval,
                "dynamic_filter": dynamic,
                "block": block,
                "true_gradient": float(true[0]),
                "estimated_gradient": float(mean[0]),
                "posterior_gradient_sd": sd,
                "calibration_block": block % calibration_interval == 0,
            }
        )
    burn = 10
    errors_array = np.array(errors[burn:])
    return (
        float(math.sqrt(np.mean(errors_array**2))),
        float(np.mean(np.abs(errors_array))),
        float(np.mean(covered[burn:])),
        rows,
    )


def program_318(
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    summaries = []
    for interval in (1, 5, 10, 20):
        for dynamic in (False, True):
            replicate_metrics = []
            for repetition in range(16):
                rmse, mae, coverage, trace = run_detector_filter(
                    rng, interval, dynamic
                )
                replicate_metrics.append((rmse, mae, coverage))
                if repetition == 0:
                    rows.extend(trace)
            summaries.append(
                {
                    "calibration_interval": interval,
                    "dynamic_filter": dynamic,
                    "median_gradient_rmse": float(
                        np.median([metric[0] for metric in replicate_metrics])
                    ),
                    "median_gradient_mae": float(
                        np.median([metric[1] for metric in replicate_metrics])
                    ),
                    "median_95pct_coverage": float(
                        np.median([metric[2] for metric in replicate_metrics])
                    ),
                    "replications": len(replicate_metrics),
                }
            )
    best_dynamic = min(
        (
            summary
            for summary in summaries
            if summary["dynamic_filter"]
        ),
        key=lambda summary: summary["median_gradient_rmse"],
    )
    matched_static = next(
        summary
        for summary in summaries
        if not summary["dynamic_filter"]
        and summary["calibration_interval"]
        == best_dynamic["calibration_interval"]
    )
    return (
        {
            "status": (
                "[Strong evidence] synthetic Bayesian/Laplace sequential "
                "drift audit"
            ),
            "summaries": summaries,
            "best_dynamic_design": best_dynamic,
            "matched_static_design": matched_static,
            "rmse_improvement_factor": float(
                matched_static["median_gradient_rmse"]
                / best_dynamic["median_gradient_rmse"]
            ),
            "boundary": (
                "The random-walk law, detector likelihood, process noise, "
                "and calibration cadence are supplied. Coverage and RMSE do "
                "not certify a real drifting apparatus."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P319: global spectral-order chamber complex and look-elsewhere rule
# ---------------------------------------------------------------------------


def mode_coefficient_matrix() -> np.ndarray:
    modes = np.arange(N)
    coefficients = np.zeros((6, N))
    for distance in range(1, 6):
        coefficients[distance - 1] = 2.0 * (
            1.0 - np.cos(2.0 * math.pi * modes * distance / N)
        )
    coefficients[5] = 1.0 - (-1.0) ** modes
    return coefficients[:, 1:7]


def fingerprint_jacobian_at(weights: np.ndarray) -> np.ndarray:
    coordinates = weights[:5]
    step = 1.0e-6
    jacobian = np.zeros((11, 5))
    for column in range(5):
        plus = coordinates.copy()
        minus = coordinates.copy()
        plus[column] += step
        minus[column] -= step
        plus_weights = np.r_[plus, 1.0 - np.sum(plus)]
        minus_weights = np.r_[minus, 1.0 - np.sum(minus)]
        jacobian[:, column] = (
            p281.fingerprint_rows(plus_weights[None, :])[0]
            - p281.fingerprint_rows(minus_weights[None, :])[0]
        ) / (2.0 * step)
    return jacobian


def program_319() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    mode_matrix = mode_coefficient_matrix()
    rows = []
    feasible = 0
    full_rank = 0
    determinants = []
    strict_weights = core.STRICT_K / np.sum(core.STRICT_K)
    strict_modes = strict_weights @ mode_matrix
    strict_order = tuple(np.argsort(strict_modes).tolist())
    strict_chamber_margin = None
    for permutation in itertools.permutations(range(6)):
        inequalities = []
        for left, right in zip(permutation[:-1], permutation[1:]):
            inequalities.append(
                np.r_[
                    mode_matrix[:, left] - mode_matrix[:, right],
                    1.0,
                ]
            )
        for index in range(6):
            row = np.zeros(7)
            row[index] = -1.0
            row[-1] = 1.0
            inequalities.append(row)
        result = linprog(
            np.r_[np.zeros(6), -1.0],
            A_ub=np.array(inequalities),
            b_ub=np.zeros(len(inequalities)),
            A_eq=np.r_[np.ones(6), 0.0][None, :],
            b_eq=[1.0],
            bounds=[(0.0, None)] * 7,
            method="highs",
        )
        if not result.success or result.x[-1] <= 1.0e-9:
            continue
        feasible += 1
        weights = result.x[:6]
        jacobian = fingerprint_jacobian_at(weights)
        rank = int(np.linalg.matrix_rank(jacobian, tol=1.0e-8))
        gram_determinant = float(np.linalg.det(jacobian.T @ jacobian))
        full_rank += int(rank == 5)
        determinants.append(gram_determinant)
        if permutation == strict_order:
            strict_chamber_margin = float(result.x[-1])
        rows.append(
            {
                "permutation": permutation,
                "maximum_simplex_and_order_margin": float(result.x[-1]),
                "center_jacobian_rank": rank,
                "center_gram_determinant": gram_determinant,
                "is_strict_order_chamber": permutation == strict_order,
            }
        )
    single_target_probability = 2.9e-8
    look_elsewhere = {
        str(target_count): min(
            1.0, target_count * single_target_probability
        )
        for target_count in (1, 10, 100, feasible)
    }
    return (
        {
            "status": (
                "[Proven] finite chamber upper bound and Bonferroni rule; "
                "[Strong evidence] exhaustive LP chamber feasibility/rank"
            ),
            "possible_mode_orders": math.factorial(6),
            "full_dimensional_feasible_chambers": feasible,
            "full_rank_chamber_centers": full_rank,
            "strict_mode_order": strict_order,
            "strict_order_chamber_maximum_margin": strict_chamber_margin,
            "minimum_center_gram_determinant": min(determinants),
            "maximum_center_gram_determinant": max(determinants),
            "p289_single_target_reference_probability": (
                single_target_probability
            ),
            "bonferroni_upper_bounds": look_elsewhere,
            "theorem": (
                "Six nonzero circulant Fourier mode values have at most 6!="
                "720 strict order chambers. For M preregistered targets, the "
                "family-wise false-positive probability is at most the sum "
                "of their individual probabilities, irrespective of target "
                "dependence."
            ),
            "boundary": (
                "LP feasibility and center ranks are floating-point "
                "certificates. The P289 probability belongs only to its "
                "declared Dirichlet null and cannot be assigned to arbitrary "
                "targets without a separate estimate."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P320: restricted-clock identifiability and E-optimal time design
# ---------------------------------------------------------------------------


def clock_transition(
    strict_a: np.ndarray,
    time: float,
    theta: np.ndarray,
) -> np.ndarray:
    log_scale, quadratic, cubic = theta
    tau = time + quadratic * time**2 + cubic * time**3
    return expm(-math.exp(log_scale) * tau * strict_a)


def clock_fisher(
    strict_a: np.ndarray,
    times: tuple[float, ...],
    theta: np.ndarray,
    shots_per_vertex_time: int,
) -> np.ndarray:
    steps = np.array([1.0e-5, 1.0e-6, 1.0e-6])
    fisher = np.zeros((3, 3))
    for time in times:
        probability = clock_transition(strict_a, time, theta)
        derivatives = []
        for index, step in enumerate(steps):
            plus = theta.copy()
            minus = theta.copy()
            plus[index] += step
            minus[index] -= step
            derivatives.append(
                (
                    clock_transition(strict_a, time, plus)
                    - clock_transition(strict_a, time, minus)
                )
                / (2.0 * step)
            )
        for vertex in range(N):
            p = np.maximum(probability[vertex], 1.0e-12)
            jacobian = np.column_stack(
                [derivative[vertex] for derivative in derivatives]
            )
            fisher += (
                shots_per_vertex_time
                * jacobian.T
                @ np.diag(1.0 / p)
                @ jacobian
            )
    return hermitian(fisher)


def program_320(
    strict_a: np.ndarray,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    theta = np.array([0.0, 0.08, -0.03])
    time_grid = np.linspace(0.05, 1.20, 12)
    rows = []
    best = None
    for times in itertools.combinations(time_grid, 4):
        fisher = clock_fisher(strict_a, times, theta, 500)
        eigenvalues = np.linalg.eigvalsh(fisher)
        score = float(eigenvalues[0])
        condition = float(eigenvalues[-1] / eigenvalues[0])
        row = {
            "times": times,
            "minimum_fisher_eigenvalue": score,
            "condition_number": condition,
            "cramer_rao_sd_log_scale": float(
                math.sqrt(np.linalg.inv(fisher)[0, 0])
            ),
            "cramer_rao_sd_quadratic_clock": float(
                math.sqrt(np.linalg.inv(fisher)[1, 1])
            ),
            "cramer_rao_sd_cubic_clock": float(
                math.sqrt(np.linalg.inv(fisher)[2, 2])
            ),
        }
        rows.append(row)
        if best is None or score > best[0]:
            best = (score, row, fisher)
    assert best is not None
    reference_times = (0.15, 0.35, 0.65, 1.0)
    reference_fisher = clock_fisher(
        strict_a, reference_times, theta, 500
    )
    reference_eigenvalues = np.linalg.eigvalsh(reference_fisher)

    # Without the calibration tau'(0)=1, predictions depend on only the
    # three products kappa*a, kappa*b, kappa*c. The four-parameter
    # coefficient sensitivity has an exact one-dimensional null direction.
    uncalibrated_sensitivity = np.array(
        [
            [
                time + theta[1] * time**2 + theta[2] * time**3,
                time,
                time**2,
                time**3,
            ]
            for time in time_grid
        ]
    )
    uncalibrated_rank = int(
        np.linalg.matrix_rank(uncalibrated_sensitivity, tol=1.0e-10)
    )
    return (
        {
            "status": (
                "[Proven] scale/clock rank obstruction; "
                "[Strong evidence] restricted-clock E-optimal time design"
            ),
            "clock_family": (
                "tau(t)=t+b*t^2+c*t^3 with externally calibrated tau'(0)=1"
            ),
            "best_design": best[1],
            "reference_times": reference_times,
            "reference_minimum_fisher_eigenvalue": float(
                reference_eigenvalues[0]
            ),
            "e_optimal_improvement_factor": float(
                best[0] / reference_eigenvalues[0]
            ),
            "uncalibrated_four_parameter_sensitivity_rank": (
                uncalibrated_rank
            ),
            "uncalibrated_parameter_count": 4,
            "boundary": (
                "Identifiability holds only inside the preregistered cubic "
                "clock family after an external slope calibration. Higher-"
                "order warps can reproduce a finite design as in P304."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P321: physical reservoir/hardware admission gate
# ---------------------------------------------------------------------------


def program_321() -> dict[str, Any]:
    result = external_gate("physical reservoir hardware")
    result.update(
        {
            "calibrated_clock_record": False,
            "detector_transfer_record": False,
            "paired_control_runs": False,
            "frozen_primary_endpoint": False,
        }
    )
    return result


# ---------------------------------------------------------------------------
# P322: single electroweak-role naturality certificate
# ---------------------------------------------------------------------------


def spectral_role_candidates(size: int) -> dict[str, float]:
    strict = p295.laplacian_from_profile(
        p295.strict_weights_for_size(size)
    )
    eigenvalues = np.linalg.eigvalsh(strict)
    positive = eigenvalues[eigenvalues > 1.0e-9]
    normalized = positive / np.sum(positive)
    return {
        "normalized_gap": float(positive[0] / positive[-1]),
        "trace_square_ratio": float(
            np.sum(positive**2) / np.sum(positive) ** 2
        ),
        "spectral_shannon_entropy_normalized": float(
            -np.sum(normalized * np.log(normalized))
            / math.log(len(normalized))
        ),
        "zero_mode_fraction": 1.0 / size,
        "alpha_geo_over_carrier_size": ALPHA_GEO / size,
    }


def program_322() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    rows = []
    by_candidate: dict[str, list[float]] = {}
    for size in (12, 16, 24):
        candidates = spectral_role_candidates(size)
        for name, value in candidates.items():
            rows.append(
                {
                    "carrier_size": size,
                    "candidate": name,
                    "value": value,
                    "absolute_error_from_legacy_weinberg_benchmark": abs(
                        value - WEINBERG_BENCHMARK
                    ),
                }
            )
            by_candidate.setdefault(name, []).append(value)
    stability = {
        name: {
            "range_across_carriers": max(values) - min(values),
            "maximum_benchmark_error": max(
                abs(value - WEINBERG_BENCHMARK) for value in values
            ),
        }
        for name, values in by_candidate.items()
    }
    certificate = {
        "typed_completion_naturality": False,
        "ordered_U1_SU2_sector_functor": False,
        "strict_coupling_observables_gprime_g": False,
        "renormalization_scale_semantics": False,
        "carrier_and_representation_invariance": False,
        "independent_experimental_holdout": False,
        "certificate_pass": False,
    }
    return (
        {
            "status": (
                "[Proven] two-sector pointing obstruction; "
                "[Refuted] current electroweak role-transfer certificate"
            ),
            "legacy_benchmark": WEINBERG_BENCHMARK,
            "candidate_stability": stability,
            "role_transfer_certificate": certificate,
            "theorem": (
                "An unlabelled spectral operator supplies no ordered choice "
                "between two gauge-sector roles. Sector swap is a free "
                "two-point action with no invariant point. Since "
                "g'^2/(g^2+g'^2) changes to 1-x under swapping, a value other "
                "than 1/2 requires sector pointing and coupling semantics. "
                "The kernel alone does not provide them."
            ),
            "factorization_test": (
                "Any physical role R that factors as Rhat∘K must be constant "
                "on kernel fibers. The historical alpha_geo/12 map has no "
                "exported strict Rhat and changes its carrier-denominator "
                "version on C12/C16/C24."
            ),
            "boundary": (
                "This obstructs an unpointed, spectrum-only electroweak role. "
                "It does not prove that no axiom-augmented sector theory or "
                "future strict sector source can define such a role."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# Figures, summary, and main
# ---------------------------------------------------------------------------


def make_figures(results: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    plt.style.use("seaborn-v0_8-whitegrid")

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    design = results["P309"]["e_optimal_rows"]
    axes[0].semilogy(
        [row["probe_count"] for row in design],
        [row["minimum_profiled_pole_information"] for row in design],
        "o-",
    )
    collision = [
        row
        for row in results["_P309_rows"]
        if row["row_type"] == "collision_lower_bound"
    ]
    for count in (1, 4, 16):
        selected = [
            row for row in collision if row["probe_count"] == count
        ]
        axes[1].semilogx(
            [row["pole_half_separation"] for row in selected],
            [row["bayes_error_lower_bound"] for row in selected],
            label=f"{count} probes",
        )
    axes[0].set(
        xlabel="greedy probe count",
        ylabel="minimum profiled pole information",
    )
    axes[1].set(
        xlabel="pole half-separation",
        ylabel="Bayes-error lower bound",
        ylim=(-0.02, 0.52),
    )
    axes[1].legend()
    fig.suptitle("P309 — local probe design versus minimax pole collision")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p309_minimax_stieltjes.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    grouped = results["P311"]["grouped_results"]
    for phase in (0.0, 1.0e-4, 1.0e-3):
        selected = [
            row
            for row in grouped
            if row["phase_parameter_sigma"] == phase
        ]
        axes[0].plot(
            [row["mean_mode_loss"] for row in selected],
            [row["p95_maximum_effect_error"] for row in selected],
            "o-",
            label=f"phase σ={phase:g}",
        )
        axes[1].plot(
            [row["mean_mode_loss"] for row in selected],
            [
                row["median_response_minimum_singular_value"]
                for row in selected
            ],
            "o-",
            label=f"phase σ={phase:g}",
        )
    axes[0].set(xlabel="mean mode loss", ylabel="95% max effect error")
    axes[1].set(
        xlabel="mean mode loss",
        ylabel="median minimum response singular value",
    )
    axes[0].legend()
    fig.suptitle("P311 — loss and component-drift completed current readout")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p311_lossy_mesh.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    parent_rows = results["P312"]["rows"]
    axes[0].plot(
        [row["carrier_size"] for row in parent_rows],
        [row["c12_trained_cubic_spectral_rmse"] for row in parent_rows],
        "o-",
        color="#9f4f42",
    )
    tail_rows = [
        row
        for row in results["_P313_rows"]
        if row["row_type"] == "tail_curve"
    ]
    axes[1].loglog(
        [row["distance"] for row in tail_rows],
        [row["positive_hyperbolic_mixture"] for row in tail_rows],
        label="positive legacy-scale mixture",
    )
    axes[1].loglog(
        [row["distance"] for row in tail_rows],
        [row["strict_envelope"] for row in tail_rows],
        label="strict envelope",
    )
    axes[0].set(
        xlabel="carrier size",
        ylabel="C12-trained cubic spectral RMSE",
    )
    axes[1].set(xlabel="distance", ylabel="envelope")
    axes[1].legend()
    fig.suptitle("P312/P313 — common parent existence versus tail no-go")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p312_p313_parent_tail.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    for dynamics in ("heat", "wave"):
        selected = [
            row
            for row in results["_P315_rows"]
            if row["dynamics"] == dynamics
        ]
        axes[0].semilogx(
            [row["time"] for row in selected],
            [row["maximum_vertex_tv"] for row in selected],
            label=dynamics,
        )
    signed_rows = results["P316"]["grid_results"]
    axes[1].plot(
        [row["moment_grid_size"] for row in signed_rows],
        [row["minimum_grid_negative_mass"] for row in signed_rows],
        "o-",
        label="grid optimum",
    )
    axes[1].axhline(
        results["P316"]["analytic_negative_mass_lower_bound"],
        ls="--",
        color="black",
        label="analytic lower bound",
    )
    axes[0].set(xlabel="time", ylabel="maximum vertex TV")
    axes[1].set(xlabel="moment grid size", ylabel="negative mass")
    axes[0].legend()
    axes[1].legend()
    fig.suptitle("P315/P316 — operational separation and signed resource")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p315_p316_distance_negativity.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    drift_rows = [
        row
        for row in results["_P318_rows"]
        if row["calibration_interval"] == 5 and row["dynamic_filter"]
    ]
    axes[0].plot(
        [row["block"] for row in drift_rows],
        [row["true_gradient"] for row in drift_rows],
        label="true",
    )
    axes[0].plot(
        [row["block"] for row in drift_rows],
        [row["estimated_gradient"] for row in drift_rows],
        label="filter",
    )
    chamber_rows = results["_P319_rows"]
    axes[1].hist(
        [row["maximum_simplex_and_order_margin"] for row in chamber_rows],
        bins=24,
        color="#557a95",
    )
    axes[0].set(xlabel="block", ylabel="gradient")
    axes[1].set(
        xlabel="maximum chamber margin",
        ylabel="number of mode-order chambers",
    )
    axes[0].legend()
    fig.suptitle("P318/P319 — detector drift and spectral chamber complex")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p318_p319_drift_chambers.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2))
    clock_rows = results["_P320_rows"]
    axes[0].hist(
        [row["minimum_fisher_eigenvalue"] for row in clock_rows],
        bins=28,
        color="#688f4e",
    )
    role_rows = results["_P322_rows"]
    candidates = sorted({row["candidate"] for row in role_rows})
    for candidate in candidates:
        selected = [row for row in role_rows if row["candidate"] == candidate]
        axes[1].plot(
            [row["carrier_size"] for row in selected],
            [row["value"] for row in selected],
            "o-",
            label=candidate.replace("_", " "),
        )
    axes[1].axhline(
        WEINBERG_BENCHMARK,
        color="black",
        ls="--",
        label="legacy benchmark",
    )
    axes[0].set(
        xlabel="minimum Fisher eigenvalue",
        ylabel="four-time designs",
    )
    axes[1].set(xlabel="carrier size", ylabel="candidate scalar")
    axes[1].legend(fontsize=7)
    fig.suptitle("P320/P322 — clock design and role naturality failure")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p320_p322_clock_role.png", dpi=180)
    plt.close(fig)


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    headlines = {
        "P309": "finite probes cannot remove the pole-collision minimax modulus",
        "P310": "Schur minimization follows from a certified square completion",
        "P311": "loss and no-click effects complete the synthetic current readout",
        "P312": "a cross-supported parent exists but imports both target operators",
        "P313": "positive legacy-scale mixtures cannot have the strict tail",
        "P314": "strict phase needs one new infinite-order generator",
        "P315": "legacy and strict remain operationally separable after scale fit",
        "P316": "strict attenuation requires a nonzero signed path resource",
        "P317": "the independent laboratory bundle gate remains closed",
        "P318": "a tuned dynamic filter improves the supplied drift model",
        "P319": "the positive simplex contains all 720 spectral order chambers",
        "P320": "clock identifiability requires a restricted family and calibration",
        "P321": "the physical reservoir hardware gate remains closed",
        "P322": "an unpointed spectrum cannot transfer the electroweak role",
    }
    return [
        {
            "program": program,
            "status": results[program]["status"],
            "headline": headline,
        }
        for program, headline in headlines.items()
    ]


def main() -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    strict_a, strict_w = core.strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "programs": "P309-P322",
            "release": "10.28",
            "seed": SEED,
            "scope": (
                "exact finite theorems, dependency-free Lean, synthetic "
                "tolerance/design audits, and external-evidence gates; no "
                "silent legacy/strict identity, selector closure, dimensional "
                "source, role transfer, laboratory claim, L_total, or ToE"
            ),
            "new_theoretical_objects": {
                "O62_minimax_spectral_resolution_modulus": (
                    "the pole-separation scale below which finite noisy probes "
                    "cannot uniformly distinguish a split spectral atom"
                ),
                "O63_loss_completed_measurement_package": (
                    "unitary mesh + mode transmission + no-click effect + "
                    "detector confusion + calibration record"
                ),
                "O64_joint_feature_parent": (
                    "a cross-supported Gram parent whose typed compressions "
                    "are repaired legacy and strict"
                ),
                "O65_regular_variation_bridge_class": (
                    "path-count and attenuation laws classified by their "
                    "Mellin/regular-variation index"
                ),
                "O66_minimal_nontorsion_phase_extension": (
                    "one infinite-order generator r=exp(i/4000) for the frozen "
                    "strict omega/phi tuple"
                ),
                "O67_operational_kernel_pseudometric": (
                    "supremal outcome-law TV over a declared finite menu"
                ),
                "O68_signed_path_resource": (
                    "Jordan negative mass required outside the positive "
                    "path-scale cone"
                ),
                "O69_drift_aware_sequential_instrument": (
                    "detector likelihood + process prior + calibration cadence "
                    "+ posterior stopping/coverage rule"
                ),
                "O70_spectral_chamber_complex": (
                    "the finite stratification of positive radial weights by "
                    "Fourier-mode order"
                ),
                "O71_clock_identifiability_certificate": (
                    "restricted clock family + external slope calibration + "
                    "full-rank operational Fisher information"
                ),
                "O72_two_sector_pointing_obstruction": (
                    "a free sector-swap torsor preventing an ordered physical "
                    "coupling role from an unlabelled spectrum"
                ),
                "O73_minimal_bridge_resource_budget": (
                    "(signed path negativity, nontorsion phase generator, "
                    "independently sourced cross-parent law, operational "
                    "pointing, dimensional scale); every component remains "
                    "separately typed"
                ),
            },
        }
    }
    results["P309"], p309_rows = program_309(strict_a, rng)
    results["P310"] = program_310(rng)
    results["P311"], p311_rows = program_311(strict_w, rng)
    results["P312"], p312_rows = program_312()
    results["P313"], p313_rows = program_313()
    results["P314"] = program_314()
    results["P315"], p315_rows = program_315(strict_a)
    results["P316"], p316_rows = program_316()
    results["P317"] = program_317()
    results["P318"], p318_rows = program_318(rng)
    results["P319"], p319_rows = program_319()
    results["P320"], p320_rows = program_320(strict_a)
    results["P321"] = program_321()
    results["P322"], p322_rows = program_322()

    results["_P309_rows"] = p309_rows
    results["_P313_rows"] = p313_rows
    results["_P315_rows"] = p315_rows
    results["_P318_rows"] = p318_rows
    results["_P319_rows"] = p319_rows
    results["_P320_rows"] = p320_rows
    results["_P322_rows"] = p322_rows
    make_figures(results)

    public_results = {
        key: value for key, value in results.items() if not key.startswith("_")
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(public_results), indent=2) + "\n",
        encoding="utf-8",
    )
    write_csv(SUMMARY_PATH, summary_rows(results))
    write_csv(MINIMAX_PATH, p309_rows)
    write_csv(LOSS_PATH, p311_rows)
    write_csv(PARENT_PATH, p312_rows)
    write_csv(REGULAR_PATH, p313_rows)
    write_csv(OPERATIONAL_PATH, p315_rows)
    write_csv(SIGNED_PATH, p316_rows)
    write_csv(DRIFT_PATH, p318_rows)
    write_csv(CHAMBER_PATH, p319_rows)
    write_csv(CLOCK_PATH, p320_rows)
    write_csv(ROLE_PATH, p322_rows)
    for row in summary_rows(results):
        print(row["program"], row["status"], "-", row["headline"])
    print(RESULTS_PATH.name)


if __name__ == "__main__":
    main()
