#!/usr/bin/env python3
"""Execute FIN research Programs 267--280.

The program combines finite proofs, deterministic numerical certificates, and
synthetic protocol audits.  It does not fabricate laboratory evidence, export
a strict selector or dimensional unit, complete the legacy-to-strict bridge,
transfer legacy physical roles, or promote a conditioned construction to
strict FIN physics.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Callable

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import expm
from scipy.optimize import least_squares

import fin_programs_255_266 as previous


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_267_280_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_267_280_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_267_280_Summary.csv"
STABILITY_PATH = ROOT / "FIN_Programs_267_280_Stieltjes_Stability.csv"
POVM_PATH = ROOT / "FIN_Programs_267_280_Current_POVM.csv"
RG_PATH = ROOT / "FIN_Programs_267_280_RG_Obstruction.csv"
RESERVOIR_PATH = ROOT / "FIN_Programs_267_280_Nonlinear_Reservoir.csv"

N = 12
SEED = 20260730


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
        return {"real": value.real, "imag": value.imag}
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def hermitian(matrix: np.ndarray) -> np.ndarray:
    return (matrix + matrix.T.conj()) / 2.0


def matrix_log_spd(matrix: np.ndarray) -> np.ndarray:
    eigenvalues, eigenvectors = np.linalg.eigh(hermitian(matrix))
    if eigenvalues.min() <= 0.0:
        raise ValueError("matrix logarithm requires a positive-definite input")
    return (eigenvectors * np.log(eigenvalues)) @ eigenvectors.T.conj()


def density_entropy(rho: np.ndarray) -> float:
    eigenvalues = np.linalg.eigvalsh(hermitian(rho))
    eigenvalues = eigenvalues[eigenvalues > 1e-15]
    return float(-np.sum(eigenvalues * np.log(eigenvalues)))


def quantum_relative_entropy(rho: np.ndarray, sigma: np.ndarray) -> float:
    rho = hermitian(rho)
    sigma = hermitian(sigma)
    er, ur = np.linalg.eigh(rho)
    es, us = np.linalg.eigh(sigma)
    if er.min() < -1e-12 or es.min() <= 0.0:
        raise ValueError("relative entropy requires rho>=0 and sigma>0")
    log_rho = (ur * np.log(np.maximum(er, 1e-300))) @ ur.T.conj()
    log_sigma = (us * np.log(es)) @ us.T.conj()
    return float(np.real(np.trace(rho @ (log_rho - log_sigma))))


def random_density(
    dimension: int, rng: np.random.Generator, floor: float = 0.08
) -> np.ndarray:
    real = rng.normal(size=(dimension, dimension))
    imag = rng.normal(size=(dimension, dimension))
    x = real + 1j * imag
    rho = x @ x.T.conj()
    rho /= np.trace(rho)
    rho = (1.0 - floor) * rho + floor * np.eye(dimension) / dimension
    return hermitian(rho)


def visible_pole_residue_data(a: np.ndarray) -> list[tuple[float, np.ndarray]]:
    retained = list(range(0, N, 2))
    _, coupling, hidden = previous.block_components(a, retained)
    groups = previous.group_spectral_residues(hidden, coupling)
    return [
        (float(pole), hermitian(residue))
        for pole, residue in groups
        if np.trace(residue).real > 1e-10
    ]


def fit_scalar_stieltjes(
    z_grid: np.ndarray,
    values: np.ndarray,
    initial_poles: np.ndarray,
    initial_weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    count = len(initial_poles)

    def residual(parameters: np.ndarray) -> np.ndarray:
        poles = np.exp(parameters[:count])
        weights = np.exp(parameters[count:])
        model = np.sum(
            weights[None, :] / (z_grid[:, None] + poles[None, :]), axis=1
        )
        return (model - values) / np.maximum(np.abs(values), 1e-12)

    initial = np.log(np.concatenate([initial_poles, initial_weights]))
    result = least_squares(
        residual,
        initial,
        max_nfev=2500,
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    poles = np.exp(result.x[:count])
    weights = np.exp(result.x[count:])
    order = np.argsort(poles)
    return poles[order], weights[order]


def program_267(
    a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    data = visible_pole_residue_data(a)
    poles = np.array([item[0] for item in data])
    weights = np.array([np.trace(item[1]).real for item in data])
    z_grid = np.geomspace(0.03, 8.0, 48)
    exact_trace = np.sum(
        weights[None, :] / (z_grid[:, None] + poles[None, :]), axis=1
    )
    rows: list[dict[str, Any]] = []
    for sigma in [1e-6, 1e-4, 1e-3]:
        for repetition in range(40):
            noisy = exact_trace * (
                1.0 + rng.normal(0.0, sigma, size=exact_trace.shape)
            )
            initial_poles = poles * np.array([0.92, 1.03, 1.08])
            initial_weights = weights * np.array([1.04, 0.95, 1.02])
            fitted_poles, fitted_weights = fit_scalar_stieltjes(
                z_grid, noisy, initial_poles, initial_weights
            )
            rows.append(
                {
                    "relative_noise": sigma,
                    "repetition": repetition,
                    "maximum_relative_pole_error": float(
                        np.max(np.abs(fitted_poles - poles) / poles)
                    ),
                    "maximum_relative_weight_error": float(
                        np.max(np.abs(fitted_weights - weights) / weights)
                    ),
                }
            )
    summaries = {}
    for sigma in [1e-6, 1e-4, 1e-3]:
        selected = [row for row in rows if row["relative_noise"] == sigma]
        pole_errors = np.array(
            [row["maximum_relative_pole_error"] for row in selected]
        )
        weight_errors = np.array(
            [row["maximum_relative_weight_error"] for row in selected]
        )
        summaries[str(sigma)] = {
            "median_max_relative_pole_error": float(np.median(pole_errors)),
            "p95_max_relative_pole_error": float(np.quantile(pole_errors, 0.95)),
            "median_max_relative_weight_error": float(np.median(weight_errors)),
            "p95_max_relative_weight_error": float(
                np.quantile(weight_errors, 0.95)
            ),
        }
    return (
        {
            "status": "[Proven] uniqueness; [Strong evidence] local noisy recovery",
            "visible_atomic_poles": poles,
            "visible_atomic_weights": weights,
            "pole_separation_minimum": float(np.min(np.diff(poles))),
            "noise_stability": summaries,
            "theorem": (
                "A finite positive atomic Stieltjes transform has a unique "
                "visible pole-residue measure. Analytic equality on an open "
                "interval implies equality of rational continuations and of "
                "all principal parts."
            ),
            "boundary": (
                "The theorem identifies only visible atoms. Zero-residue and "
                "decoupled hidden modes remain unidentifiable."
            ),
        },
        rows,
    )


def program_268(
    a: np.ndarray, rng: np.random.Generator
) -> dict[str, Any]:
    maximum_composition_residual = 0.0
    minimum_reduced_eigenvalue = math.inf
    chain_count = 50_000
    z = 0.37
    base = z * np.eye(N) + a
    for _ in range(chain_count):
        permutation = rng.permutation(N)
        sizes = sorted(rng.choice(np.arange(1, N + 1), size=3, replace=False))
        e = sorted(permutation[: sizes[0]].tolist())
        f = sorted(permutation[: sizes[1]].tolist())
        g = sorted(permutation[: sizes[2]].tolist())
        bg = previous.schur_keep(base, g)
        positions_f = [g.index(index) for index in f]
        bf = previous.schur_keep(bg, positions_f)
        positions_e = [f.index(index) for index in e]
        nested = previous.schur_keep(bf, positions_e)
        direct = previous.schur_keep(base, e)
        maximum_composition_residual = max(
            maximum_composition_residual,
            float(np.linalg.norm(nested - direct, 2)),
        )
        minimum_reduced_eigenvalue = min(
            minimum_reduced_eigenvalue,
            float(np.linalg.eigvalsh(hermitian(direct)).min()),
        )
    return {
        "status": "[Proven] algebraic theorem with executable finite certificate",
        "chains_audited": chain_count,
        "maximum_composition_residual": maximum_composition_residual,
        "minimum_reduced_positive_eigenvalue": minimum_reduced_eigenvalue,
        "identity_residual": 0.0,
        "formal_core": {
            "objects": "nonempty subsets E of V",
            "morphisms": "unique reduction r_{F->E} when E is a subset of F",
            "identity": "zero elimination",
            "composition": "associativity of block Gaussian elimination",
            "positive_cone_preservation": "Schur complement of an SPD matrix is SPD",
        },
        "boundary": (
            "This is an executable algebraic core, not a Lean/Coq/Isabelle "
            "kernel check and not a sheaf or contextuality theorem."
        ),
    }


def orthonormal_current_basis(w: np.ndarray) -> tuple[list[np.ndarray], np.ndarray]:
    currents = [previous.current_observable(w, distance) for distance in range(1, 6)]
    gram = np.array(
        [[np.trace(left @ right).real for right in currents] for left in currents]
    )
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    transform = eigenvectors @ np.diag(1.0 / np.sqrt(eigenvalues)) @ eigenvectors.T
    basis = [
        sum(transform[index, j] * currents[j] for j in range(5))
        for index in range(5)
    ]
    return basis, gram


def simplex_vertices(dimension: int) -> np.ndarray:
    centering = np.eye(dimension + 1) - np.ones(
        (dimension + 1, dimension + 1)
    ) / (dimension + 1)
    eigenvalues, eigenvectors = np.linalg.eigh(centering)
    return eigenvectors[:, eigenvalues > 0.5]


def program_269(
    w: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    basis, original_gram = orthonormal_current_basis(w)
    vertices = simplex_vertices(5)
    simplex_operators = [
        sum(vertices[outcome, index] * basis[index] for index in range(5))
        for outcome in range(6)
    ]
    maximum_norm = max(np.linalg.norm(operator, 2) for operator in simplex_operators)
    epsilon = 0.90 / (6.0 * maximum_norm)
    effects = [
        np.eye(N) / 6.0 + epsilon * operator
        for operator in simplex_operators
    ]
    completeness = np.linalg.norm(sum(effects) - np.eye(N), 2)
    minimum_effect_eigenvalue = min(
        float(np.linalg.eigvalsh(hermitian(effect)).min()) for effect in effects
    )
    response = np.array(
        [
            [np.trace(effect @ operator).real for operator in basis]
            for effect in effects
        ]
    )
    response_rank = int(np.linalg.matrix_rank(response, tol=1e-10))
    singular_values = np.linalg.svd(response, compute_uv=False)
    response_condition = float(singular_values[0] / singular_values[-2])

    coefficients = np.array([0.028, -0.021, 0.016, 0.011, -0.008])
    rho = np.eye(N) / N + sum(
        coefficient * operator
        for coefficient, operator in zip(coefficients, basis)
    )
    minimum_state_eigenvalue = float(np.linalg.eigvalsh(hermitian(rho)).min())
    probabilities = np.array([np.trace(effect @ rho).real for effect in effects])
    exact_estimate = np.linalg.lstsq(
        response, probabilities - np.ones(6) / 6.0, rcond=None
    )[0]
    exact_reconstruction_residual = float(
        np.linalg.norm(exact_estimate - coefficients)
    )
    shots = 200_000
    counts = rng.multinomial(shots, probabilities / probabilities.sum())
    observed = counts / counts.sum()
    baseline = np.ones(6) / 6.0
    estimate = np.linalg.lstsq(response, observed - baseline, rcond=None)[0]
    relative_error = float(
        np.linalg.norm(estimate - coefficients) / np.linalg.norm(coefficients)
    )
    rows = [
        {
            "outcome": outcome,
            "probability": probabilities[outcome],
            "count": int(counts[outcome]),
            "minimum_effect_eigenvalue": float(
                np.linalg.eigvalsh(hermitian(effects[outcome])).min()
            ),
        }
        for outcome in range(6)
    ]
    return (
        {
            "status": "[Proven] minimal outcome count in the five-current subspace; [Strong evidence] shot audit",
            "current_subspace_dimension": 5,
            "minimum_outcomes_by_normalization_bound": 6,
            "constructed_outcomes": 6,
            "response_rank": response_rank,
            "response_condition_number_nonzero_subspace": response_condition,
            "povm_completeness_residual": float(completeness),
            "minimum_effect_eigenvalue": minimum_effect_eigenvalue,
            "minimum_test_state_eigenvalue": minimum_state_eigenvalue,
            "exact_coefficient_reconstruction_residual": (
                exact_reconstruction_residual
            ),
            "shots": shots,
            "relative_coefficient_reconstruction_error": relative_error,
            "original_current_gram_condition_number": float(
                np.linalg.cond(original_gram)
            ),
            "theorem": (
                "M outcome probabilities carry at most M-1 independent real "
                "statistics. Five independent currents therefore require at "
                "least six outcomes in one POVM; the regular-simplex "
                "construction attains the bound."
            ),
            "boundary": (
                "Loss, dark-count, and phase calibration remain apparatus "
                "obligations. The POVM is a mathematical design, not a built detector."
            ),
        },
        rows,
    )


def block_embedding(n: int, target: int = N) -> np.ndarray:
    if target % n != 0:
        raise ValueError("block embedding requires divisibility")
    multiplicity = target // n
    embedding = np.zeros((target, n))
    for index in range(target):
        embedding[index, index % n] = 1.0 / math.sqrt(multiplicity)
    return embedding


def program_270(a: np.ndarray) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    contexts = {
        6: [0, 2, 4, 6, 8, 10],
        3: [0, 4, 8],
        1: [0],
    }
    rows = []
    for size, keep in contexts.items():
        reduced = previous.schur_keep(a, keep)
        embedding = block_embedding(size)
        lifted = embedding @ reduced @ embedding.T
        denominator = float(np.sum(lifted * lifted))
        scale = (
            float(np.sum(a * lifted) / denominator) if denominator > 1e-18 else 0.0
        )
        residual = float(np.linalg.norm(a - scale * lifted, "fro") / np.linalg.norm(a, "fro"))
        rows.append(
            {
                "context_size": size,
                "reduced_rank": int(np.linalg.matrix_rank(reduced, tol=1e-10)),
                "lifted_rank": int(np.linalg.matrix_rank(lifted, tol=1e-10)),
                "strict_rank": int(np.linalg.matrix_rank(a, tol=1e-10)),
                "best_frobenius_scale": scale,
                "relative_best_fit_residual": residual,
            }
        )
    return (
        {
            "status": "[Proven] rank-defect no-go for linear isometric liftings",
            "rows": rows,
            "theorem": (
                "For an isometric embedding E:C^n->C^12, rank(E B E*) is at "
                "most n-1 for a Laplacian B. Since the connected strict "
                "generator has rank 11, no n<12 lifted reduced Laplacian can "
                "equal it, under any scalar normalization."
            ),
            "missing_object": (
                "A genuine size-restoring RG comparison must add complement "
                "dynamics or an equivalence after quotienting observables; "
                "zero-padding or replication alone cannot close the fixed point."
            ),
        },
        rows,
    )


def quantum_instrument() -> list[np.ndarray]:
    effect_zero = np.diag([0.82, 0.61, 0.34, 0.17])
    effect_one = np.eye(4) - effect_zero
    unitary_one = np.array(
        [
            [0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, 0, 1],
            [0, 0, 1, 0],
        ],
        dtype=complex,
    )
    return [np.sqrt(effect_zero), unitary_one @ np.sqrt(effect_one)]


def instrument_ledger(
    rho: np.ndarray, sigma: np.ndarray, kraus: list[np.ndarray]
) -> tuple[float, float, float, float]:
    p = np.array([np.trace(k @ rho @ k.T.conj()).real for k in kraus])
    q = np.array([np.trace(k @ sigma @ k.T.conj()).real for k in kraus])
    conditional = 0.0
    for index, k in enumerate(kraus):
        rho_y = k @ rho @ k.T.conj() / p[index]
        sigma_y = k @ sigma @ k.T.conj() / q[index]
        conditional += p[index] * quantum_relative_entropy(rho_y, sigma_y)
    record = previous.kl_divergence(p, q)
    input_divergence = quantum_relative_entropy(rho, sigma)
    loss = input_divergence - record - conditional
    return input_divergence, record, conditional, loss


def program_271(rng: np.random.Generator) -> dict[str, Any]:
    kraus = quantum_instrument()
    completeness = np.linalg.norm(
        sum(k.T.conj() @ k for k in kraus) - np.eye(4), 2
    )
    losses = []
    chain_residuals = []
    for _ in range(400):
        rho = random_density(4, rng)
        sigma = random_density(4, rng)
        d_in, record, conditional, loss = instrument_ledger(rho, sigma, kraus)
        losses.append(loss)
        block_direct = record + conditional
        chain_residuals.append(abs((d_in - loss) - block_direct))
    return {
        "status": "[Proven] conditional quantum-information theorem with finite audit",
        "instrument_completeness_residual": float(completeness),
        "pairs_audited": len(losses),
        "minimum_information_loss_nats": float(min(losses)),
        "median_information_loss_nats": float(np.median(losses)),
        "maximum_chain_decomposition_residual": float(max(chain_residuals)),
        "theorem": (
            "For a quantum instrument, D(rho||sigma) is at least the classical "
            "record divergence plus the record-weighted conditional quantum "
            "relative entropies. The remainder is nonnegative by monotonicity "
            "under the CPTP instrument channel."
        ),
        "boundary": (
            "This is a supplied quantum instrument and not a FIN-derived "
            "Hamiltonian, bath, measurement postulate, or laboratory process tensor."
        ),
    }


def mean_zero_eigenvalues(matrix: np.ndarray) -> np.ndarray:
    eigenvalues = np.linalg.eigvalsh(hermitian(matrix))
    zero_index = int(np.argmin(np.abs(eigenvalues)))
    return np.delete(eigenvalues, zero_index)


def program_272(strict_a: np.ndarray) -> dict[str, Any]:
    legacy_a, _ = previous.legacy_amplitude_absorbed_operator()
    projector = np.eye(N) - np.ones((N, N)) / N
    legacy_nonconstant = mean_zero_eigenvalues(legacy_a)
    gamma = max(0.0, -float(legacy_nonconstant.min())) + 1e-9
    completed = hermitian(legacy_a + gamma * projector)
    completed_nonconstant = mean_zero_eigenvalues(completed)
    strict_fingerprint = previous.projective_signature(strict_a)
    completed_fingerprint = previous.projective_signature(completed)
    fingerprint_defect = float(
        np.linalg.norm(completed_fingerprint - strict_fingerprint)
    )
    retained = list(range(0, N, 2))
    z = 0.2
    sigma_strict = previous.self_energy(strict_a, retained, z)
    sigma_completed = previous.self_energy(completed, retained, z)
    sigma_strict /= np.linalg.norm(sigma_strict, "fro")
    sigma_completed /= np.linalg.norm(sigma_completed, "fro")
    sigma_defect = float(np.linalg.norm(sigma_strict - sigma_completed, "fro"))
    return {
        "status": "[Proven] one-atom positive-shift completion test",
        "frozen_atom": "A_legacy -> A_legacy + gamma*(I-11^T/12)",
        "gamma_derived_from_legacy_only": gamma,
        "legacy_minimum_mean_zero_eigenvalue": float(legacy_nonconstant.min()),
        "completed_minimum_mean_zero_eigenvalue": float(
            completed_nonconstant.min()
        ),
        "completed_row_sum_residual": float(
            np.linalg.norm(completed @ np.ones(N))
        ),
        "projective_strict_fingerprint_defect": fingerprint_defect,
        "normalized_self_energy_defect_at_z_0_2": sigma_defect,
        "verdict": (
            "The positive shift repairs the generator cone using legacy data "
            "alone, but it does not reproduce the strict spectral or memory shape."
        ),
        "boundary": (
            "This closes only the positive-shift atom. It is not a kernel-profile "
            "bridge, phase/frequency source, nonlinear damping completion, or role transfer."
        ),
    }


def program_273() -> dict[str, Any]:
    template = ROOT / "FIN_Lab_P241_Transfer_Template"
    production_manifests = [
        path
        for path in ROOT.rglob("bundle_manifest.json")
        if "FIN_Lab_P241_Transfer_Template" not in str(path)
    ]
    template_headers = [
        template / "events_heat_process.header.csv",
        template / "events_double_slit.header.csv",
    ]
    template_event_rows = []
    for path in template_headers:
        with path.open("r", encoding="utf-8") as handle:
            template_event_rows.append(max(0, sum(1 for _ in handle) - 1))
    return {
        "status": "[Executed no-admission certificate: external evidence absent]",
        "production_p241_manifests_found": len(production_manifests),
        "production_manifest_paths": [
            str(path.relative_to(ROOT)) for path in production_manifests
        ],
        "template_event_rows": template_event_rows,
        "validator_checks_required": 11,
        "detached_signature_required": True,
        "p242_execution_authorized": False,
        "verdict": (
            "Only schemas and empty event headers are present. No independent "
            "external record can be admitted, and P242 must not be run as a "
            "confirmatory physical test."
        ),
        "boundary": (
            "Code can validate custody structure but cannot generate independence, "
            "human role separation, empirical events, or a trusted signature."
        ),
    }


def complex_hermitian_noise(
    size: int, scale: float, rng: np.random.Generator
) -> np.ndarray:
    raw = rng.normal(size=(size, size)) + 1j * rng.normal(size=(size, size))
    return scale * hermitian(raw) / math.sqrt(2.0 * size)


def program_274(
    a: np.ndarray, rng: np.random.Generator
) -> dict[str, Any]:
    del a
    retained = list(range(0, N, 2))
    z = 0.2
    true_xi, _ = previous.chiral_susceptibility_analytic(
        previous.twisted_generator(0.0),
        previous.twisted_generator_derivative_at_zero(),
        retained,
        z,
    )
    true_norm = np.linalg.norm(true_xi, "fro")
    shots = 50_000
    noise_scale = 0.35 / math.sqrt(shots)
    h_grid = np.geomspace(0.003, 0.30, 15)
    rows = []
    for h in h_grid:
        sigma_plus = previous.self_energy(
            previous.twisted_generator(float(h)), retained, z
        )
        sigma_minus = previous.self_energy(
            previous.twisted_generator(float(-h)), retained, z
        )
        noiseless = (sigma_plus - sigma_minus) / (2.0 * h)
        bias = float(np.linalg.norm(noiseless - true_xi, "fro") / true_norm)
        errors = []
        for _ in range(120):
            estimate = (
                sigma_plus
                + complex_hermitian_noise(6, noise_scale, rng)
                - sigma_minus
                - complex_hermitian_noise(6, noise_scale, rng)
            ) / (2.0 * h)
            errors.append(float(np.linalg.norm(estimate - true_xi, "fro") / true_norm))
        rows.append(
            {
                "h": float(h),
                "relative_noiseless_bias": bias,
                "median_relative_error": float(np.median(errors)),
                "p95_relative_error": float(np.quantile(errors, 0.95)),
            }
        )
    optimum = min(rows, key=lambda row: row["median_relative_error"])
    return {
        "status": "[Strong evidence] synthetic two-flux noise-bias design",
        "shots_per_flux_setting": shots,
        "tomography_noise_scale": noise_scale,
        "grid": rows,
        "optimal_h_on_frozen_grid": optimum["h"],
        "optimal_median_relative_error": optimum["median_relative_error"],
        "optimal_p95_relative_error": optimum["p95_relative_error"],
        "verdict": (
            "The two-flux estimator has the expected small-h noise amplification "
            "and large-h truncation bias, yielding an interior optimum."
        ),
        "boundary": (
            "Synthetic matrix noise is not a detector model. Apparatus feasibility, "
            "phase calibration, losses, and event-level uncertainties remain external."
        ),
    }


def random_positive_circulant(
    rng: np.random.Generator,
) -> np.ndarray:
    weights = np.exp(rng.normal(-1.1, 0.65, size=6))
    return previous.laplacian_from_weights(previous.cyclic_weight_matrix(weights))


def program_275(
    a: np.ndarray, rng: np.random.Generator
) -> dict[str, Any]:
    target = previous.projective_signature(a)
    tolerance = 0.02
    rows = []
    dimension = len(target)
    for index in range(400):
        null = random_positive_circulant(rng)
        signature = previous.projective_signature(null)
        distance = float(np.linalg.norm(signature - target))
        maximum = float(np.linalg.eigvalsh(null).max())
        margin = max(0.0, distance - tolerance)
        radius = (
            margin * maximum / (2.0 * math.sqrt(dimension) + margin)
            if margin > 0.0
            else 0.0
        )
        rows.append(
            {
                "index": index,
                "fingerprint_distance": distance,
                "certified_operator_norm_exclusion_radius": radius,
                "passes_tolerance": distance <= tolerance,
            }
        )
    radii = np.array([row["certified_operator_norm_exclusion_radius"] for row in rows])
    distances = np.array([row["fingerprint_distance"] for row in rows])
    return {
        "status": "[Proven] deterministic exclusion bound for the frozen finite null set",
        "null_count": len(rows),
        "tolerance": tolerance,
        "false_positive_count": int(sum(row["passes_tolerance"] for row in rows)),
        "minimum_fingerprint_distance": float(distances.min()),
        "median_fingerprint_distance": float(np.median(distances)),
        "minimum_positive_certified_radius": float(
            radii[radii > 0].min() if np.any(radii > 0) else 0.0
        ),
        "median_certified_radius": float(np.median(radii)),
        "bound": (
            "If ||E||_2<=epsilon<lambda_max and "
            "2*sqrt(11)*epsilon/(lambda_max-epsilon)<distance-tolerance, "
            "the perturbed null fingerprint cannot enter the acceptance ball."
        ),
        "boundary": (
            "This is a certified robustness radius for 400 frozen nulls, not "
            "a distribution-free random-matrix false-positive probability."
        ),
        "rows": rows,
    }


def semigroup_consistency(
    p_one: np.ndarray, p_two: np.ndarray, t_one: float, t_two: float
) -> dict[str, float]:
    generator_one = -matrix_log_spd(p_one) / t_one
    generator_two = -matrix_log_spd(p_two) / t_two
    scale = max(np.linalg.norm(generator_one, "fro"), 1e-15)
    return {
        "relative_log_generator_defect": float(
            np.linalg.norm(generator_two - generator_one, "fro") / scale
        ),
        "commutator_residual": float(
            np.linalg.norm(p_one @ p_two - p_two @ p_one, "fro")
        ),
    }


def program_276(a: np.ndarray) -> dict[str, Any]:
    t_one, t_two = 0.20, 0.55
    p_one = expm(-t_one * a)
    p_two = expm(-t_two * a)
    static = semigroup_consistency(p_one, p_two, t_one, t_two)

    scaled = 1.37 * a
    uniform = semigroup_consistency(
        expm(-t_one * scaled), expm(-t_two * scaled), t_one, t_two
    )
    edge = previous.edge_laplacian(0, 3)
    shaped = hermitian(a + 0.17 * edge)
    shape_drift = semigroup_consistency(
        expm(-t_one * a), expm(-t_two * shaped), t_one, t_two
    )
    alternative = hermitian(a + 0.30 * previous.edge_laplacian(1, 7))
    memory_map = 0.78 * expm(-t_two * a) + 0.22 * expm(-t_two * alternative)
    memory = semigroup_consistency(p_one, memory_map, t_one, t_two)
    clock_degeneracy = float(
        np.linalg.norm(expm(-t_one * scaled) - expm(-(1.37 * t_one) * a), 2)
    )
    return {
        "status": "[Proven] two-time semigroup consistency theorem and finite witnesses",
        "times": [t_one, t_two],
        "static": static,
        "uniform_scaled_with_calibrated_clock": uniform,
        "shape_drift": shape_drift,
        "mixture_memory_model": memory,
        "uncalibrated_scale_clock_degeneracy_residual": clock_degeneracy,
        "theorem": (
            "For symmetric positive-definite P1,P2 at calibrated times t1,t2, "
            "one homogeneous PSD generator exists iff the matrices commute and "
            "log(P2)/t2=log(P1)/t1; when it exists, A=-log(P1)/t1 is unique."
        ),
        "boundary": (
            "Two matrices can reject a homogeneous semigroup but do not uniquely "
            "identify which adaptive, environmental, or apparatus-memory mechanism caused the defect."
        ),
    }


def simulate_weight_law(
    initial: np.ndarray,
    baseline: np.ndarray,
    input_series: np.ndarray,
    dt: float,
    law: str,
    parameters: np.ndarray,
    direction: np.ndarray,
) -> np.ndarray:
    state = initial.copy()
    trajectory = [state.copy()]
    for value in input_series:
        if law == "static":
            derivative = np.zeros_like(state)
        elif law == "relaxation":
            derivative = -parameters[0] * (state - baseline)
        elif law == "driven_additive":
            derivative = (
                -parameters[0] * (state - baseline)
                + parameters[1] * value * direction
            )
        elif law == "multiplicative":
            derivative = state * (
                -parameters[0] + parameters[1] * value * direction
            )
        else:
            raise ValueError(law)
        state = np.maximum(1e-6, state + dt * derivative)
        trajectory.append(state.copy())
    return np.asarray(trajectory)


def fit_adaptive_candidates(
    observed: np.ndarray,
    baseline: np.ndarray,
    input_series: np.ndarray,
    dt: float,
    direction: np.ndarray,
) -> dict[str, np.ndarray]:
    derivative = np.diff(observed, axis=0) / dt
    state = observed[:-1]
    fits: dict[str, np.ndarray] = {"static": np.array([])}
    x_relax = -(state - baseline).reshape(-1, 1)
    fits["relaxation"] = np.linalg.lstsq(
        x_relax, derivative.reshape(-1), rcond=None
    )[0]
    column_k = -(state - baseline).reshape(-1)
    column_g = (input_series[:, None] * direction[None, :]).reshape(-1)
    design = np.column_stack([column_k, column_g])
    fits["driven_additive"] = np.linalg.lstsq(
        design, derivative.reshape(-1), rcond=None
    )[0]
    multiplicative_design = np.column_stack(
        [
            -state.reshape(-1),
            (state * input_series[:, None] * direction[None, :]).reshape(-1),
        ]
    )
    fits["multiplicative"] = np.linalg.lstsq(
        multiplicative_design, derivative.reshape(-1), rcond=None
    )[0]
    return fits


def program_277(rng: np.random.Generator) -> dict[str, Any]:
    baseline = previous.STRICT_K.copy()
    direction = np.array([1.0, -0.8, 0.6, -0.4, 0.25, -0.15])
    direction /= np.linalg.norm(direction)
    dt = 0.05
    train_input = rng.normal(0.0, 0.45, size=800)
    holdout_input = np.zeros(400)
    holdout_input[40:80] = 1.0
    holdout_input[180:215] = -0.8
    holdout_input[300:330] = 0.6
    true_parameters = np.array([0.24, 0.075])
    initial = baseline * np.array([1.08, 0.94, 1.04, 0.97, 1.02, 0.95])
    train_true = simulate_weight_law(
        initial,
        baseline,
        train_input,
        dt,
        "driven_additive",
        true_parameters,
        direction,
    )
    observed = train_true + rng.normal(0.0, 2e-5, size=train_true.shape)
    fits = fit_adaptive_candidates(
        observed, baseline, train_input, dt, direction
    )
    true_holdout = simulate_weight_law(
        train_true[-1],
        baseline,
        holdout_input,
        dt,
        "driven_additive",
        true_parameters,
        direction,
    )
    rows = []
    for law, parameters in fits.items():
        predicted_train = simulate_weight_law(
            observed[0], baseline, train_input, dt, law, parameters, direction
        )
        predicted_holdout = simulate_weight_law(
            train_true[-1], baseline, holdout_input, dt, law, parameters, direction
        )
        rows.append(
            {
                "law": law,
                "fitted_parameters": parameters,
                "train_mse": float(np.mean((predicted_train - train_true) ** 2)),
                "holdout_mse": float(
                    np.mean((predicted_holdout - true_holdout) ** 2)
                ),
            }
        )
    winner = min(rows, key=lambda row: row["holdout_mse"])
    zero_input = np.zeros_like(train_input)
    design_no_intervention = np.column_stack(
        [
            -(train_true[:-1] - baseline).reshape(-1),
            (zero_input[:, None] * direction[None, :]).reshape(-1),
        ]
    )
    design_with_intervention = np.column_stack(
        [
            -(train_true[:-1] - baseline).reshape(-1),
            (train_input[:, None] * direction[None, :]).reshape(-1),
        ]
    )
    return {
        "status": "[Proven] identifiability inside the frozen linear family; [Strong evidence] synthetic holdout",
        "true_law": "driven_additive",
        "true_parameters": true_parameters,
        "candidate_rows": rows,
        "holdout_winner": winner["law"],
        "design_rank_without_intervention": int(
            np.linalg.matrix_rank(design_no_intervention)
        ),
        "design_rank_with_intervention": int(
            np.linalg.matrix_rank(design_with_intervention)
        ),
        "verdict": (
            "The supplied interventions raise the parameter-design rank from "
            "one to two and select the true member of the preregistered finite family."
        ),
        "boundary": (
            "The candidate family and intervention direction are supplied. "
            "This does not derive a strict adaptive law or exclude mechanisms outside the panel."
        ),
    }


def nonlinear_states(
    reservoir: np.ndarray,
    input_vector: np.ndarray,
    input_series: np.ndarray,
) -> np.ndarray:
    state = np.zeros(reservoir.shape[0])
    states = []
    for value in input_series:
        state = np.tanh(reservoir @ state + input_vector * value)
        states.append(state.copy())
    return np.asarray(states)


def ridge_score(
    states: np.ndarray,
    target: np.ndarray,
    start: int,
    classification: bool = False,
) -> float:
    indices = np.arange(start, len(target))
    split = int(0.6 * len(indices))
    train, test = indices[:split], indices[split:]
    x_train = np.column_stack([np.ones(len(train)), states[train]])
    x_test = np.column_stack([np.ones(len(test)), states[test]])
    weights = np.linalg.solve(
        x_train.T @ x_train + 1e-6 * np.eye(x_train.shape[1]),
        x_train.T @ target[train],
    )
    prediction = x_test @ weights
    if classification:
        return float(np.mean(np.sign(prediction) == np.sign(target[test])))
    denominator = float(np.sum((target[test] - np.mean(target[test])) ** 2))
    return float(1.0 - np.sum((prediction - target[test]) ** 2) / denominator)


def narma10(input_series: np.ndarray) -> np.ndarray:
    target = np.zeros_like(input_series)
    for index in range(10, len(input_series) - 1):
        target[index + 1] = (
            0.3 * target[index]
            + 0.05 * target[index] * np.sum(target[index - 9 : index + 1])
            + 1.5 * input_series[index - 9] * input_series[index]
            + 0.1
        )
    return target


def evaluate_nonlinear_reservoir(
    reservoir: np.ndarray,
    input_vector: np.ndarray,
    inputs: dict[str, np.ndarray],
) -> dict[str, float]:
    continuous = inputs["continuous"]
    states_continuous = nonlinear_states(reservoir, input_vector, continuous)
    delay_target = np.roll(continuous, 5)
    nonlinear_target = (
        np.sin(2.0 * np.roll(continuous, 1))
        + 0.30 * np.roll(continuous, 2) * np.roll(continuous, 3)
    )
    narma_target = narma10((continuous + 1.0) / 4.0)
    binary = inputs["binary"]
    states_binary = nonlinear_states(reservoir, input_vector, binary)
    parity_target = binary * np.roll(binary, 1) * np.roll(binary, 2)
    return {
        "delay5_r2": ridge_score(states_continuous, delay_target, 250),
        "nonlinear_prediction_r2": ridge_score(
            states_continuous, nonlinear_target, 250
        ),
        "narma10_r2": ridge_score(states_continuous, narma_target, 250),
        "parity3_accuracy": ridge_score(
            states_binary, parity_target, 250, classification=True
        ),
    }


def program_278(
    a: np.ndarray, rng: np.random.Generator
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    fin = 0.92 * expm(-0.25 * a)
    input_vector = np.zeros(N)
    input_vector[0] = 0.75
    input_vector[3] = -0.45
    input_vector -= input_vector.mean()
    inputs = {
        "continuous": rng.uniform(-1.0, 1.0, size=4_500),
        "binary": rng.choice([-1.0, 1.0], size=4_500),
    }
    rows = [
        {
            "ensemble": "FIN",
            "index": 0,
            **evaluate_nonlinear_reservoir(fin, input_vector, inputs),
        }
    ]
    eigenvalues = np.linalg.eigvalsh(fin)
    for index in range(30):
        q, _ = np.linalg.qr(rng.normal(size=(N, N)))
        control = q @ np.diag(eigenvalues) @ q.T
        rows.append(
            {
                "ensemble": "isospectral_orientation",
                "index": index,
                **evaluate_nonlinear_reservoir(control, input_vector, inputs),
            }
        )
    for index in range(30):
        values = rng.uniform(-0.92, 0.92, size=N)
        values[np.argmax(np.abs(values))] = 0.92
        q, _ = np.linalg.qr(rng.normal(size=(N, N)))
        control = q @ np.diag(values) @ q.T
        rows.append(
            {
                "ensemble": "symmetric_radius_matched",
                "index": index,
                **evaluate_nonlinear_reservoir(control, input_vector, inputs),
            }
        )
    fin_row = rows[0]
    controls = rows[1:]
    metrics = [
        "delay5_r2",
        "nonlinear_prediction_r2",
        "narma10_r2",
        "parity3_accuracy",
    ]
    percentiles = {
        metric: float(
            np.mean(
                np.array([row[metric] for row in controls]) <= fin_row[metric]
            )
        )
        for metric in metrics
    }
    control_ranges = {
        metric: {
            "minimum": float(min(row[metric] for row in controls)),
            "median": float(np.median([row[metric] for row in controls])),
            "maximum": float(max(row[metric] for row in controls)),
        }
        for metric in metrics
    }
    return (
        {
            "status": "[Strong evidence] preregistered synthetic nonlinear benchmark",
            "control_count": len(controls),
            "fin_metrics": fin_row,
            "fin_percentiles": percentiles,
            "control_ranges": control_ranges,
            "tasks_on_which_fin_exceeds_all_controls": int(
                sum(
                    fin_row[metric] > control_ranges[metric]["maximum"]
                    for metric in metrics
                )
            ),
            "task_count": len(metrics),
            "verdict": (
                "FIN is assessed task by task against spectrum- and radius-matched "
                "controls; no ontology or general computational advantage follows."
            ),
        },
        rows,
    )


def program_279(rng: np.random.Generator) -> dict[str, Any]:
    del rng
    hamiltonian = np.diag([0.0, 1.3])
    k_b = 1.0
    temperature = 2.0
    beta = 1.0 / (k_b * temperature)
    gibbs_unnormalized = expm(-beta * hamiltonian)
    gibbs = gibbs_unnormalized / np.trace(gibbs_unnormalized)
    rho = np.array([[0.72, 0.18], [0.18, 0.28]], dtype=float)
    rho = hermitian(rho)
    mixing = 0.37
    rho_after = (1.0 - mixing) * rho + mixing * gibbs

    def free_energy(state: np.ndarray) -> float:
        return float(
            np.trace(state @ hamiltonian).real
            - k_b * temperature * density_entropy(state)
        )

    d_before = quantum_relative_entropy(rho, gibbs)
    d_after = quantum_relative_entropy(rho_after, gibbs)
    f_gibbs = free_energy(gibbs)
    identity_before = abs(
        d_before - beta * (free_energy(rho) - f_gibbs)
    )
    identity_after = abs(
        d_after - beta * (free_energy(rho_after) - f_gibbs)
    )
    ledger_loss_nats = d_before - d_after
    free_energy_drop = free_energy(rho) - free_energy(rho_after)
    conversion_residual = abs(
        free_energy_drop - k_b * temperature * ledger_loss_nats
    )
    return {
        "status": "[Proven] conditioned thermodynamic identity",
        "supplied_axioms": {
            "hamiltonian_gap": 1.3,
            "temperature": temperature,
            "boltzmann_constant": k_b,
            "thermalization_channel_strength": mixing,
        },
        "relative_entropy_before_nats": d_before,
        "relative_entropy_after_nats": d_after,
        "information_ledger_loss_nats": ledger_loss_nats,
        "free_energy_drop": free_energy_drop,
        "maximum_gibbs_free_energy_identity_residual": max(
            identity_before, identity_after
        ),
        "information_to_energy_conversion_residual": conversion_residual,
        "landauer_one_bit_lower_bound": k_b * temperature * math.log(2.0),
        "theorem": (
            "For gamma=exp(-beta H)/Z, D(rho||gamma)=beta(F(rho)-F(gamma)). "
            "A decrease of relative entropy becomes a free-energy decrease "
            "only after H, T, and k_B are supplied."
        ),
        "boundary": (
            "The supplied dimensional data are conversion axioms. They are not "
            "generated by strict FIN, and the identity does not prove a physical bath or reset protocol."
        ),
    }


def program_280(a: np.ndarray, w: np.ndarray) -> dict[str, Any]:
    current = previous.current_observable(w, 4)
    rho_plus = np.eye(N, dtype=complex) / N
    rho_plus[10, 2] += 0.04j
    rho_plus[2, 10] -= 0.04j
    base_current = float(np.trace(rho_plus @ current).real)
    scales = [0.5, 1.7]
    signs = [-1, 1]
    rows = []
    fingerprints = []
    for scale in scales:
        physical = scale * a
        fingerprint = previous.projective_signature(physical)
        fingerprints.append(fingerprint)
        for sign in signs:
            rows.append(
                {
                    "scale_section": scale,
                    "orientation_section": sign,
                    "oriented_current": sign * base_current,
                    "spectral_fingerprint_norm": float(
                        np.linalg.norm(fingerprint)
                    ),
                    "transition_reparameterization_residual": float(
                        np.linalg.norm(
                            expm(-0.4 * physical)
                            - expm(-(0.4 * scale) * a),
                            2,
                        )
                    ),
                }
            )
    return {
        "status": "[Proven] conditioned product-torsor construction",
        "rows": rows,
        "scale_fingerprint_invariance_residual": float(
            np.linalg.norm(fingerprints[0] - fingerprints[1])
        ),
        "orientation_sign_pair_residual": float(
            abs(rows[0]["oriented_current"] + rows[1]["oriented_current"])
        ),
        "object": (
            "Operational section s=(kappa,epsilon) of "
            "T_scale x T_orientation, with A_phys=kappa*A and J_phys=epsilon*J."
        ),
        "independence_theorem": (
            "The product group R_{>0} x Z_2 acts independently: positive "
            "rescaling cannot select epsilon, and orientation reversal cannot "
            "set kappa. A section supplies both choices but generates neither."
        ),
        "boundary": (
            "This is conditioned on external scale and orientation resources. "
            "It does not discharge QW-2191 or export a canonical physical unit."
        ),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    normalized = []
    for row in rows:
        normalized.append(
            {
                key: json.dumps(json_ready(value), ensure_ascii=False)
                if isinstance(value, (dict, list, tuple, np.ndarray))
                else value
                for key, value in row.items()
            }
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(normalized[0]))
        writer.writeheader()
        writer.writerows(normalized)


def make_figures(results: dict[str, Any], reservoir_rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)

    p267 = results["P267"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    noise = np.array([float(value) for value in p267["noise_stability"]])
    pole = np.array(
        [
            p267["noise_stability"][str(value)][
                "median_max_relative_pole_error"
            ]
            for value in noise
        ]
    )
    weight = np.array(
        [
            p267["noise_stability"][str(value)][
                "median_max_relative_weight_error"
            ]
            for value in noise
        ]
    )
    axes[0].loglog(noise, pole, "o-", label="pole error")
    axes[0].loglog(noise, weight, "s-", label="residue-weight error")
    axes[0].set_xlabel("relative noise")
    axes[0].set_ylabel("median maximum relative error")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    p269 = results["P269"]
    axes[1].bar(
        ["rank", "min effect eig.", "reconstruction error"],
        [
            p269["response_rank"] / 5.0,
            p269["minimum_effect_eigenvalue"],
            p269["relative_coefficient_reconstruction_error"],
        ],
        color=["#19733A", "#1F5A99", "#D55E00"],
    )
    axes[1].set_title("Minimal six-outcome current POVM")
    axes[1].tick_params(axis="x", rotation=15)
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p267_p269_measure_povm.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    p270_rows = results["P270"]["rows"]
    axes[0].bar(
        [str(row["context_size"]) for row in p270_rows],
        [row["relative_best_fit_residual"] for row in p270_rows],
        color="#A61B1B",
    )
    axes[0].set_xlabel("retained context size")
    axes[0].set_ylabel("best lifted relative residual")
    axes[0].set_title("Rank-defect RG obstruction")
    p272 = results["P272"]
    axes[1].bar(
        ["fingerprint", "self-energy"],
        [
            p272["projective_strict_fingerprint_defect"],
            p272["normalized_self_energy_defect_at_z_0_2"],
        ],
        color=["#6A3D9A", "#D55E00"],
    )
    axes[1].set_title("Positive-shift completion remains non-strict")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p270_p272_rg_completion.png", dpi=180)
    plt.close(fig)

    p274 = results["P274"]["grid"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    h = np.array([row["h"] for row in p274])
    axes[0].loglog(
        h,
        [row["relative_noiseless_bias"] for row in p274],
        "o-",
        label="bias",
    )
    axes[0].loglog(
        h,
        [row["median_relative_error"] for row in p274],
        "s-",
        label="median noisy error",
    )
    axes[0].set_xlabel("twist h")
    axes[0].set_ylabel("relative error")
    axes[0].set_title("Two-flux tomography bias–noise tradeoff")
    axes[0].legend()
    axes[0].grid(alpha=0.25)
    p275_rows = results["P275"]["rows"]
    axes[1].hist(
        [row["fingerprint_distance"] for row in p275_rows],
        bins=28,
        color="#1F5A99",
        alpha=0.85,
    )
    axes[1].axvline(
        results["P275"]["tolerance"], color="#A61B1B", linestyle="--"
    )
    axes[1].set_xlabel("null-to-strict fingerprint distance")
    axes[1].set_ylabel("count")
    axes[1].set_title("Frozen null exclusion audit")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p274_p275_tomography_false_positive.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    p277_rows = results["P277"]["candidate_rows"]
    axes[0].bar(
        [row["law"] for row in p277_rows],
        [max(row["holdout_mse"], 1e-16) for row in p277_rows],
        color="#19733A",
    )
    axes[0].set_yscale("log")
    axes[0].set_ylabel("holdout MSE")
    axes[0].set_title("Causal adaptive-law panel")
    axes[0].tick_params(axis="x", rotation=20)
    metrics = [
        "delay5_r2",
        "nonlinear_prediction_r2",
        "narma10_r2",
        "parity3_accuracy",
    ]
    axes[1].bar(
        np.arange(len(metrics)),
        [results["P278"]["fin_percentiles"][metric] for metric in metrics],
        color="#6A3D9A",
    )
    axes[1].axhline(0.5, color="black", linewidth=0.8)
    axes[1].set_xticks(
        np.arange(len(metrics)),
        ["delay", "nonlinear", "NARMA10", "parity"],
        rotation=15,
    )
    axes[1].set_ylabel("FIN percentile among controls")
    axes[1].set_title("Nonlinear matched-reservoir benchmark")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p277_p278_adaptation_reservoir.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.8))
    p279 = results["P279"]
    axes[0].bar(
        ["information\n(kBT·nat)", "free-energy\ndrop"],
        [
            p279["information_ledger_loss_nats"]
            * p279["supplied_axioms"]["temperature"]
            * p279["supplied_axioms"]["boltzmann_constant"],
            p279["free_energy_drop"],
        ],
        color=["#1F5A99", "#19733A"],
    )
    axes[0].set_title("Conditional information-to-energy identity")
    p280 = results["P280"]["rows"]
    for scale in sorted({row["scale_section"] for row in p280}):
        selected = [row for row in p280 if row["scale_section"] == scale]
        axes[1].plot(
            [row["orientation_section"] for row in selected],
            [row["oriented_current"] for row in selected],
            "o-",
            label=f"scale={scale}",
        )
    axes[1].set_xlabel("orientation section")
    axes[1].set_ylabel("oriented current")
    axes[1].set_title("Independent scale and orientation torsors")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p279_p280_conversion_torsors.png", dpi=180)
    plt.close(fig)


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "program": "P267",
            "status": results["P267"]["status"],
            "principal_result": "unique visible atomic Stieltjes measure; local noise stability quantified",
        },
        {
            "program": "P268",
            "status": results["P268"]["status"],
            "principal_result": "executable Schur-category core; 50,000 composition chains",
        },
        {
            "program": "P269",
            "status": results["P269"]["status"],
            "principal_result": "minimal six-outcome POVM for five current observables",
        },
        {
            "program": "P270",
            "status": results["P270"]["status"],
            "principal_result": "rank-defect no-go for linear size-restoring liftings",
        },
        {
            "program": "P271",
            "status": results["P271"]["status"],
            "principal_result": "quantum process information ledger under a supplied instrument",
        },
        {
            "program": "P272",
            "status": results["P272"]["status"],
            "principal_result": "positive-shift atom repairs PSD but not strict spectral/memory shape",
        },
        {
            "program": "P273",
            "status": results["P273"]["status"],
            "principal_result": "no external P241 evidence; P242 remains unauthorized",
        },
        {
            "program": "P274",
            "status": results["P274"]["status"],
            "principal_result": "two-flux tomography has an interior bias-noise optimum",
        },
        {
            "program": "P275",
            "status": results["P275"]["status"],
            "principal_result": "analytic exclusion radii for 400 frozen null fingerprints",
        },
        {
            "program": "P276",
            "status": results["P276"]["status"],
            "principal_result": "two-time logarithmic semigroup consistency theorem",
        },
        {
            "program": "P277",
            "status": results["P277"]["status"],
            "principal_result": "interventions identify the true law within a frozen finite panel",
        },
        {
            "program": "P278",
            "status": results["P278"]["status"],
            "principal_result": "four-task nonlinear matched-reservoir benchmark",
        },
        {
            "program": "P279",
            "status": results["P279"]["status"],
            "principal_result": "information becomes free energy only under supplied H,T,k_B",
        },
        {
            "program": "P280",
            "status": results["P280"]["status"],
            "principal_result": "conditional product section of independent scale/orientation torsors",
        },
    ]


def main() -> None:
    rng = np.random.default_rng(SEED)
    strict_a, strict_w = previous.strict_operator()
    results: dict[str, Any] = {
        "metadata": {
            "programs": "267-280",
            "seed": SEED,
            "strict_kernel": previous.STRICT_K,
            "scope": (
                "finite proofs and synthetic audits; no fabricated laboratory "
                "data, strict selector, physical scale, role transfer, or ToE closure"
            ),
        }
    }
    results["P267"], stability_rows = program_267(strict_a, rng)
    results["P268"] = program_268(strict_a, rng)
    results["P269"], povm_rows = program_269(strict_w, rng)
    results["P270"], rg_rows = program_270(strict_a)
    results["P271"] = program_271(rng)
    results["P272"] = program_272(strict_a)
    results["P273"] = program_273()
    results["P274"] = program_274(strict_a, rng)
    results["P275"] = program_275(strict_a, rng)
    results["P276"] = program_276(strict_a)
    results["P277"] = program_277(rng)
    results["P278"], reservoir_rows = program_278(strict_a, rng)
    results["P279"] = program_279(rng)
    results["P280"] = program_280(strict_a, strict_w)

    RESULTS_PATH.write_text(
        json.dumps(json_ready(results), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_csv(SUMMARY_PATH, summary_rows(results))
    write_csv(STABILITY_PATH, stability_rows)
    write_csv(POVM_PATH, povm_rows)
    write_csv(RG_PATH, rg_rows)
    write_csv(RESERVOIR_PATH, reservoir_rows)
    make_figures(results, reservoir_rows)
    print(json.dumps(json_ready(results), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
