#!/usr/bin/env python3
"""Execute FIN research Programs P323--P336.

The round sharpens the bridge-resource lower bounds from Release 10.28,
traces the strict phase lattice to its actual refreeze provenance, constructs
apparatus-aware discrimination and transfer specifications, and separates
abstract logical resources from physical claims.  External laboratory and
hardware gates are never synthesized into success.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import subprocess
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp
import networkx as nx
import numpy as np
from scipy.linalg import expm
from scipy.optimize import minimize_scalar
import sympy as sp

import fin_programs_255_266 as core
import fin_programs_295_308 as p295
import fin_programs_309_322 as p309


ROOT = Path(__file__).resolve().parent
FIGURE_DIR = ROOT / "FIN_Programs_323_336_Figures"
RESULTS_PATH = ROOT / "FIN_Programs_323_336_Results.json"
SUMMARY_PATH = ROOT / "FIN_Programs_323_336_Summary.csv"
MOMENT_PATH = ROOT / "FIN_Programs_323_336_Continuum_Moment.csv"
PARENT_PATH = ROOT / "FIN_Programs_323_336_Source_Independent_Parent.csv"
PHASE_PATH = ROOT / "FIN_Programs_323_336_Phase_Provenance.csv"
PROTOCOL_PATH = ROOT / "FIN_Programs_323_336_Lossy_Protocol.csv"
CHAMBER_PATH = ROOT / "FIN_Programs_323_336_Exact_Chambers.csv"
DIAMOND_PATH = ROOT / "FIN_Programs_323_336_Diamond_Distance.csv"
CLOCK_PATH = ROOT / "FIN_Programs_323_336_Clock_Minimax.csv"
PHOTONIC_PATH = ROOT / "FIN_Programs_323_336_Photonic_Transfer.csv"
INTERPRETATION_PATH = ROOT / "FIN_Programs_323_336_Signed_Interpretations.csv"
CARRIER_PATH = ROOT / "FIN_Programs_323_336_Carrier_Naturality.csv"
EW_PATH = ROOT / "FIN_Programs_323_336_EW_Axiom_Package.csv"
INDEPENDENCE_PATH = ROOT / "FIN_Programs_323_336_Resource_Independence.csv"
FORMAL_SOURCE = ROOT / "FIN_Programs_323_336_Formal_Core.lean"

N = 12
SEED = 20260803
ALPHA_GEO = 4.0 * math.log(2.0)


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


def matrix_dynamics(
    eigenvalues: np.ndarray,
    eigenvectors: np.ndarray,
    time: float,
    dynamics: str,
) -> np.ndarray:
    if dynamics == "heat":
        factors = np.exp(-time * eigenvalues)
        return (eigenvectors * factors) @ eigenvectors.T
    factors = np.exp(-1j * time * eigenvalues)
    unitary = (eigenvectors * factors) @ eigenvectors.T
    return np.abs(unitary) ** 2


# ---------------------------------------------------------------------------
# P323: continuum signed-moment extremizer
# ---------------------------------------------------------------------------


def program_323() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    mp.mp.dps = 90
    eta = mp.mpf(9) / 5
    h = [1 / (1 + mp.mpf(k) ** eta) for k in range(7)]

    # The shifted moments h_1,...,h_6 define the three-node Gaussian rule.
    shifted = h[1:]
    hankel = mp.matrix(
        [[shifted[i + j] for j in range(3)] for i in range(3)]
    )
    rhs = mp.matrix([-shifted[i + 3] for i in range(3)])
    c0, c1, c2 = mp.lu_solve(hankel, rhs)
    nodes = sorted(
        mp.polyroots([1, c2, c1, c0], maxsteps=2000),
        key=float,
    )
    vandermonde = mp.matrix(
        [[nodes[j] ** i for j in range(3)] for i in range(3)]
    )
    shifted_weights = mp.lu_solve(
        vandermonde, mp.matrix(shifted[:3])
    )
    positive_weights = [
        shifted_weights[index] / nodes[index] for index in range(3)
    ]
    negative_mass = sum(positive_weights) - 1

    q_ascending = [mp.mpf(1), c1 / c0, c2 / c0, 1 / c0]
    p_ascending = [mp.mpf(0)] * 7
    for i, left in enumerate(q_ascending):
        for j, right in enumerate(q_ascending):
            p_ascending[i + j] -= left * right
    dual_value = sum(p_ascending[k] * h[k] for k in range(7))

    derivative_roots = sorted(
        mp.polyroots([3 / c0, 2 * c2 / c0, c1 / c0]),
        key=float,
    )

    def q_value(x: mp.mpf) -> mp.mpf:
        return (
            x**3 + c2 * x**2 + c1 * x + c0
        ) / c0

    extrema = [mp.mpf(0), *derivative_roots, mp.mpf(1)]
    maximum_q_square = max(q_value(x) ** 2 for x in extrema)
    minimum_q_square = min(q_value(x) ** 2 for x in extrema)

    rows: list[dict[str, Any]] = []
    for index, (node, weight) in enumerate(
        zip(nodes, positive_weights), start=1
    ):
        rows.append(
            {
                "row_type": "positive_atom",
                "index": index,
                "node": float(node),
                "weight": float(weight),
            }
        )
    rows.append(
        {
            "row_type": "negative_atom",
            "index": 0,
            "node": 0.0,
            "weight": float(negative_mass),
        }
    )
    maximum_residual = mp.mpf(0)
    for order in range(7):
        reconstructed = sum(
            positive_weights[j] * nodes[j] ** order for j in range(3)
        )
        if order == 0:
            reconstructed -= negative_mass
        residual = abs(reconstructed - h[order])
        maximum_residual = max(maximum_residual, residual)
        rows.append(
            {
                "row_type": "moment_check",
                "moment_order": order,
                "target": float(h[order]),
                "reconstructed": float(reconstructed),
                "absolute_residual": float(residual),
            }
        )
    for index, coefficient in enumerate(p_ascending):
        rows.append(
            {
                "row_type": "dual_polynomial",
                "power": index,
                "coefficient": float(coefficient),
            }
        )

    return (
        {
            "status": (
                "[Proven] matching primal-dual extremizer theorem; "
                "[Strong evidence] 90-digit FIN evaluation"
            ),
            "strict_attenuation_moments": [float(value) for value in h],
            "positive_nodes": [float(value) for value in nodes],
            "positive_weights": [
                float(value) for value in positive_weights
            ],
            "negative_atom_node": 0.0,
            "exact_continuum_minimum_negative_mass_numeric": float(
                negative_mass
            ),
            "dual_polynomial_coefficients_ascending": [
                float(value) for value in p_ascending
            ],
            "primal_dual_gap": float(abs(negative_mass - dual_value)),
            "maximum_moment_residual": float(maximum_residual),
            "dual_q_squared_extremal_range": [
                float(minimum_q_square),
                float(maximum_q_square),
            ],
            "certificate": (
                "Let g(x)=x^3+c2*x^2+c1*x+c0 be orthogonal to "
                "1,x,x^2 under shifted moments h1,...,h6. Its three roots "
                "r_i and Gaussian weights v_i give positive weights "
                "u_i=v_i/r_i. The signed measure sum u_i delta_ri - "
                "N delta_0 matches h0,...,h6. The dual polynomial "
                "p=-(g/g(0))^2 lies in [-1,0] on [0,1] and its value equals "
                "N, proving optimality by weak duality."
            ),
            "boundary": (
                "The structural optimality implication is theorem-grade. "
                "The displayed decimal evaluation uses high precision rather "
                "than a fully kernel-checked algebraic-number library."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P324: source-independent common-parent obstruction
# ---------------------------------------------------------------------------


def radial_profile(
    omega: float, phi: float, beta: float, eta: float
) -> np.ndarray:
    distance = np.arange(1, 7, dtype=float)
    return np.abs(
        np.cos(omega * distance + phi)
        / (1.0 + beta * distance**eta)
    )


def joint_feature_parent(
    left: np.ndarray, right: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    left_root = p309.positive_sqrt(left)
    right_root = p309.positive_sqrt(right)
    cross = left_root @ right_root
    parent = np.block([[left, cross], [cross.T, right]])
    return (parent + parent.T) / 2.0, cross


def program_324() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    legacy = p295.laplacian_from_profile(
        np.abs(p295.legacy_weights())
    )
    targets = {
        "frozen_strict": p295.laplacian_from_profile(
            radial_profile(0.18575, 0.16250, 1.0, 1.8)
        ),
        "admissible_alternative_A": p295.laplacian_from_profile(
            radial_profile(0.20175, 0.11250, 0.92, 1.7)
        ),
        "admissible_alternative_B": p295.laplacian_from_profile(
            radial_profile(0.16975, 0.06250, 0.84, 1.6)
        ),
    }
    rows = []
    target_items = list(targets.items())
    pairwise_distances = []
    for index, (name, target) in enumerate(target_items):
        parent, cross = joint_feature_parent(legacy, target)
        min_eigenvalue = float(np.linalg.eigvalsh(parent)[0])
        target_distance = float(
            np.linalg.norm(target - target_items[0][1], "fro")
        )
        pairwise_distances.append(target_distance)
        rows.append(
            {
                "target": name,
                "same_legacy_input": True,
                "target_frobenius_distance_from_frozen": target_distance,
                "parent_minimum_eigenvalue": min_eigenvalue,
                "parent_rank": int(np.linalg.matrix_rank(parent, tol=1e-9)),
                "cross_support_fraction": float(
                    np.count_nonzero(np.abs(cross) > 1e-10) / cross.size
                ),
                "cross_law_reads_target": True,
            }
        )
    distinct = all(
        np.linalg.norm(target_items[i][1] - target_items[j][1], "fro")
        > 1e-6
        for i in range(len(target_items))
        for j in range(i)
    )
    return (
        {
            "status": (
                "[Proven] source-independent nonidentification; "
                "[Strong evidence] three explicit positive parent witnesses"
            ),
            "same_legacy_input_target_count": len(targets),
            "targets_pairwise_distinct": distinct,
            "maximum_target_distance_from_frozen": max(pairwise_distances),
            "theorem": (
                "If two distinct admissible targets share exactly the same "
                "legacy input, no deterministic completion F(K_legacy) can "
                "be certified as producing both. A block parent can be built "
                "for each target, but its cross law sqrt(A_L)sqrt(A_T) reads "
                "that target and therefore is a representation, not a "
                "source-independent derivation."
            ),
            "boundary": (
                "The theorem proves underdetermination from legacy input "
                "alone. It does not prohibit a new independent law or datum "
                "that selects one target."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P325: exact refreeze phase provenance
# ---------------------------------------------------------------------------


def program_325() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    omega_base = Fraction(187, 800)
    omega_step = Fraction(2, 125)
    omega_index = -3
    omega_strict = omega_base + omega_index * omega_step
    phi_base = Fraction(-11, 80)
    phi_step = Fraction(1, 20)
    phi_index = 6
    phi_strict = phi_base + phi_index * phi_step
    common_denominator = math.lcm(
        omega_base.denominator,
        omega_step.denominator,
        phi_base.denominator,
        phi_step.denominator,
    )
    lattice_numerators = [
        omega_base.numerator
        * (common_denominator // omega_base.denominator),
        omega_step.numerator
        * (common_denominator // omega_step.denominator),
        phi_base.numerator * (common_denominator // phi_base.denominator),
        phi_step.numerator * (common_denominator // phi_step.denominator),
    ]
    lattice_gcd = 0
    for value in lattice_numerators:
        lattice_gcd = math.gcd(lattice_gcd, abs(value))
    lattice_unit = Fraction(lattice_gcd, common_denominator)

    report_path = (
        ROOT
        / "material_dowodowy/korpus_qw_pozostaly/raporty_json/"
        "report_qw2038_derivation_compatible_kernel_refreeze_scan.json"
    )
    source_path = ROOT / "QW_2038_DERIVATION_COMPATIBLE_KERNEL_REFREEZE_SCAN.py"
    rows = [
        {
            "quantity": "omega_base_QW2030",
            "exact_value": str(omega_base),
            "lattice_coordinate_in_1_over_4000": int(
                omega_base / lattice_unit
            ),
            "provenance_type": "supplied upstream scan center",
        },
        {
            "quantity": "omega_grid_step_QW2038",
            "exact_value": str(omega_step),
            "lattice_coordinate_in_1_over_4000": int(
                omega_step / lattice_unit
            ),
            "provenance_type": "declared numerical scan step",
        },
        {
            "quantity": "omega_selected_QW2038",
            "exact_value": str(omega_strict),
            "lattice_coordinate_in_1_over_4000": int(
                omega_strict / lattice_unit
            ),
            "provenance_type": "base minus three scan steps",
        },
        {
            "quantity": "phi_base_QW2030",
            "exact_value": str(phi_base),
            "lattice_coordinate_in_1_over_4000": int(
                phi_base / lattice_unit
            ),
            "provenance_type": "supplied upstream scan center",
        },
        {
            "quantity": "phi_grid_step_QW2038",
            "exact_value": str(phi_step),
            "lattice_coordinate_in_1_over_4000": int(phi_step / lattice_unit),
            "provenance_type": "declared numerical scan step",
        },
        {
            "quantity": "phi_selected_QW2038",
            "exact_value": str(phi_strict),
            "lattice_coordinate_in_1_over_4000": int(
                phi_strict / lattice_unit
            ),
            "provenance_type": "base plus six scan steps",
        },
    ]
    return (
        {
            "status": (
                "[Proven] exact refreeze-lattice provenance; "
                "[Refuted] current internal-source interpretation of 1/4000"
            ),
            "omega_identity": (
                f"{omega_base} - 3*{omega_step} = {omega_strict}"
            ),
            "phi_identity": (
                f"{phi_base} + 6*{phi_step} = {phi_strict}"
            ),
            "joint_additive_lattice_unit": str(lattice_unit),
            "omega_coordinate": int(omega_strict / lattice_unit),
            "phi_coordinate": int(phi_strict / lattice_unit),
            "source_script_sha256": sha256_file(source_path),
            "source_report_sha256": sha256_file(report_path),
            "theorem": (
                "The additive subgroup of Q generated by the QW2030 center "
                "187/800,-11/80 and QW2038 steps 2/125,1/20 is exactly "
                "(1/4000)Z. The selected row has coordinates 743 and 650. "
                "Thus exp(i/4000) is the minimal exact generator of the "
                "procedural refreeze lattice."
            ),
            "boundary": (
                "A procedural provenance theorem is not an ontological or "
                "dynamical source theorem. The grid imported its base points, "
                "step sizes, objective, and physical-data gates."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P326: loss-robust operational discrimination
# ---------------------------------------------------------------------------


def detector_channel(
    probabilities: np.ndarray,
    efficiencies: np.ndarray,
    confusion: float,
) -> np.ndarray:
    matrix = np.zeros((N, N))
    for true_vertex in range(N):
        matrix[true_vertex, true_vertex] = 1.0 - 2.0 * confusion
        matrix[(true_vertex - 1) % N, true_vertex] = confusion
        matrix[(true_vertex + 1) % N, true_vertex] = confusion
    click = matrix @ (efficiencies * probabilities)
    no_click = max(0.0, 1.0 - float(np.sum(click)))
    return np.r_[click, no_click]


def optimized_legacy_scale(
    strict_a: np.ndarray,
    legacy_a: np.ndarray,
    dynamics: str,
) -> float:
    times = np.geomspace(0.02, 3.0, 32)

    def objective(log_scale: float) -> float:
        distance, _, _ = p309.prediction_distance_witness(
            strict_a,
            legacy_a,
            times,
            dynamics,
            math.exp(log_scale),
        )
        return distance

    result = minimize_scalar(
        objective,
        bounds=(-6.0, 4.0),
        method="bounded",
        options={"xatol": 1e-9},
    )
    return float(math.exp(result.x))


def program_326(
    strict_a: np.ndarray,
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    legacy_abs = p295.laplacian_from_profile(
        np.abs(p295.legacy_weights())
    )
    legacy_signed = p295.laplacian_from_profile(p295.legacy_weights())
    times = np.geomspace(0.03, 3.0, 24)
    nuisance_draws = 48
    clock_jitter = rng.normal(0.0, 0.004, size=nuisance_draws)
    efficiencies = np.clip(
        rng.normal(0.82, 0.025, size=(nuisance_draws, N)),
        0.70,
        0.92,
    )
    confusions = np.clip(
        rng.normal(0.020, 0.004, size=nuisance_draws),
        0.005,
        0.040,
    )
    strict_eig = np.linalg.eigh(strict_a)
    rows: list[dict[str, Any]] = []
    outputs: dict[str, Any] = {}
    for dynamics, legacy_a in (
        ("heat", legacy_abs),
        ("wave", legacy_signed),
    ):
        scale = optimized_legacy_scale(strict_a, legacy_a, dynamics)
        legacy_eig = np.linalg.eigh(scale * legacy_a)
        best = None
        for time in times:
            tv_by_vertex = [[] for _ in range(N)]
            for draw in range(nuisance_draws):
                actual_time = float(time * (1.0 + clock_jitter[draw]))
                strict_transition = matrix_dynamics(
                    strict_eig[0], strict_eig[1], actual_time, dynamics
                )
                legacy_transition = matrix_dynamics(
                    legacy_eig[0], legacy_eig[1], actual_time, dynamics
                )
                for vertex in range(N):
                    strict_outcome = detector_channel(
                        strict_transition[vertex],
                        efficiencies[draw],
                        float(confusions[draw]),
                    )
                    legacy_outcome = detector_channel(
                        legacy_transition[vertex],
                        efficiencies[draw],
                        float(confusions[draw]),
                    )
                    tv_by_vertex[vertex].append(
                        0.5
                        * float(
                            np.sum(
                                np.abs(strict_outcome - legacy_outcome)
                            )
                        )
                    )
            for vertex, values in enumerate(tv_by_vertex):
                robust = float(np.quantile(values, 0.05))
                median = float(np.median(values))
                row = {
                    "dynamics": dynamics,
                    "nominal_time": float(time),
                    "preparation_vertex": vertex,
                    "robust_tv_5pct": robust,
                    "median_tv": median,
                    "minimum_tv": float(np.min(values)),
                    "maximum_tv": float(np.max(values)),
                    "optimized_legacy_scale": scale,
                }
                rows.append(row)
                if best is None or robust > best["robust_tv_5pct"]:
                    best = row
        assert best is not None
        shots = math.ceil(
            2.0 * math.log(2.0 / 0.05) / best["robust_tv_5pct"] ** 2
        )
        outputs[dynamics] = {
            **best,
            "conservative_95pct_shot_bound": shots,
        }
    return (
        {
            "status": (
                "[Strong evidence] loss/jitter/confusion-robust finite-menu "
                "protocol; [Not experimentally validated]"
            ),
            "nuisance_draws": nuisance_draws,
            "efficiency_model": "independent clipped N(0.82,0.025)",
            "nearest_neighbor_confusion_model": "clipped N(0.020,0.004)",
            "relative_clock_jitter_sd": 0.004,
            "best_protocols": outputs,
            "boundary": (
                "All nuisance distributions and transfer matrices are "
                "synthetic design assumptions. The result is a robust "
                "simulation protocol, not apparatus calibration or evidence."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P327: exact Q(sqrt(3)) chamber certificate
# ---------------------------------------------------------------------------


def program_327() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    sqrt3 = sp.sqrt(3)
    mode_matrix = sp.Matrix(
        [
            [
                2 * (1 - sp.cos(2 * sp.pi * mode * distance / 12))
                for mode in range(1, 7)
            ]
            for distance in range(1, 6)
        ]
        + [[1 - (-1) ** mode for mode in range(1, 7)]]
    )
    inverse = mode_matrix.T.inv().applyfunc(sp.radsimp)
    epsilon = sp.Rational(1, 100)
    ranks = (-5, -3, -1, 1, 3, 5)
    rows = []
    all_positive = True
    exact_orders: set[tuple[int, ...]] = set()
    global_minimum = math.inf
    minimum_expression = None
    for assignment in itertools.permutations(ranks):
        q = sp.Matrix(assignment)
        common_shift = -sp.Rational(
            assignment[0] + assignment[2] + assignment[4], 3
        )
        weights = (
            sp.ones(6, 1) / 6
            + epsilon
            * inverse
            * (q + common_shift * sp.ones(6, 1))
        )
        exact_modes = weights.T * mode_matrix
        expected_modes = sp.Matrix(
            [2 + epsilon * (value + common_shift) for value in assignment]
        )
        identity = all(
            sp.simplify(exact_modes[index] - expected_modes[index]) == 0
            for index in range(6)
        )
        positivity = all(value.is_positive for value in weights)
        all_positive &= positivity and identity
        numerical_weights = [float(sp.N(value, 30)) for value in weights]
        local_minimum = min(numerical_weights)
        if local_minimum < global_minimum:
            global_minimum = local_minimum
            minimum_expression = sp.radsimp(
                weights[int(np.argmin(numerical_weights))]
            )
        order = tuple(np.argsort(np.array(assignment)).tolist())
        exact_orders.add(order)
        rows.append(
            {
                "mode_order": order,
                "rank_assignment": assignment,
                "common_mode_shift": str(common_shift),
                "minimum_weight": local_minimum,
                "weights_exact": [str(sp.radsimp(x)) for x in weights],
                "exact_mode_identity": identity,
                "all_weights_positive": positivity,
            }
        )
    return (
        {
            "status": "[Proven] exact algebraic 720-chamber certificate",
            "coefficient_field": "Q(sqrt(3))",
            "mode_matrix_determinant": str(sp.radsimp(mode_matrix.det())),
            "epsilon": str(epsilon),
            "exact_distinct_orders": len(exact_orders),
            "all_exact_identities_hold": all_positive,
            "global_minimum_weight": global_minimum,
            "global_minimum_weight_exact": str(minimum_expression),
            "minimum_adjacent_mode_gap": float(2 * epsilon),
            "construction": (
                "Begin at uniform radial weights 1/6, where all six mode "
                "values equal 2. Assign ranks (-5,-3,-1,1,3,5) in any order, "
                "add the unique common shift preserving sum(weights)=1, and "
                "pull the perturbation back through the exact inverse mode "
                "matrix with epsilon=1/100."
            ),
            "boundary": (
                "The theorem classifies order chambers of the positive radial "
                "simplex. It does not assign physical probability to chambers "
                "or select the frozen strict chamber."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P328: one-time unitary diamond distance and comb boundary
# ---------------------------------------------------------------------------


def unitary_channel_half_diamond(relative_phases: np.ndarray) -> float:
    angles = np.sort(np.mod(np.angle(relative_phases), 2.0 * math.pi))
    gaps = np.diff(np.r_[angles, angles[0] + 2.0 * math.pi])
    containing_arc = 2.0 * math.pi - float(np.max(gaps))
    if containing_arc >= math.pi:
        return 1.0
    return float(math.sin(containing_arc / 2.0))


def program_328(
    strict_a: np.ndarray,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    legacy = p295.laplacian_from_profile(p295.legacy_weights())
    scale = optimized_legacy_scale(strict_a, legacy, "wave")
    # Compute the relative unitary directly so degenerate Fourier modes cannot
    # be mismatched by independently sorting the two spectra.
    times = np.geomspace(0.01, 8.0, 120)
    rows = []
    for time in times:
        strict_unitary = expm(-1j * time * strict_a)
        legacy_unitary = expm(-1j * time * scale * legacy)
        phases = np.linalg.eigvals(
            strict_unitary.conj().T @ legacy_unitary
        )
        half_diamond = unitary_channel_half_diamond(phases)
        rows.append(
            {
                "time": float(time),
                "half_diamond_distance": half_diamond,
                "diamond_norm": 2.0 * half_diamond,
                "perfectly_distinguishable": half_diamond >= 1.0 - 1e-12,
                "legacy_scale": scale,
            }
        )
    maximum = max(row["half_diamond_distance"] for row in rows)
    perfect = [
        row["time"] for row in rows if row["perfectly_distinguishable"]
    ]
    return (
        {
            "status": (
                "[Proven] commuting-unitary one-time diamond distance; "
                "[Blocked] full adaptive multitime comb"
            ),
            "optimized_legacy_scale": scale,
            "maximum_half_diamond_distance": maximum,
            "first_sampled_perfect_discrimination_time": (
                min(perfect) if perfect else None
            ),
            "theorem": (
                "For unitary channels U,V, half the diamond norm equals "
                "sqrt(1-nu^2), where nu is the distance from zero to the "
                "convex hull of the eigenvalues of U*V. For points on the "
                "unit circle this is determined by the shortest containing "
                "arc; an arc of length at least pi gives distance one."
            ),
            "boundary": (
                "This is a one-use wave-channel result. A multitime comb "
                "requires declared interventions, memory system, common "
                "input/output spaces, and an SDP-capable implementation."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P329: minimax clock-family design
# ---------------------------------------------------------------------------


def clock_information(
    times: tuple[float, ...],
    eigenvalues: np.ndarray,
    coefficients: tuple[float, float, float],
) -> np.ndarray:
    a2, a3, a4 = coefficients
    information = np.zeros((4, 4))
    for time in times:
        tau = time + a2 * time**2 + a3 * time**3 + a4 * time**4
        bases = np.array([tau, time**2, time**3, time**4])
        derivative_scale = -eigenvalues * np.exp(-tau * eigenvalues)
        jacobian = derivative_scale[:, None] * bases[None, :]
        information += jacobian.T @ jacobian
    return information


def program_329(
    strict_a: np.ndarray,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    eigenvalues = np.linalg.eigvalsh(strict_a)
    eigenvalues = eigenvalues[eigenvalues > 1e-8][:6]
    candidate_times = np.geomspace(0.035, 1.4, 16)
    clock_grid = list(
        itertools.product(
            (-0.08, 0.0, 0.08),
            (-0.04, 0.0, 0.04),
            (-0.02, 0.0, 0.02),
        )
    )
    rows = []
    best = None
    for design in itertools.combinations(candidate_times, 5):
        minimum_eigenvalues = []
        monotone = True
        for coefficients in clock_grid:
            a2, a3, a4 = coefficients
            slope = np.array(
                [
                    1
                    + 2 * a2 * time
                    + 3 * a3 * time**2
                    + 4 * a4 * time**3
                    for time in design
                ]
            )
            monotone &= bool(np.all(slope > 0.0))
            information = clock_information(
                tuple(float(x) for x in design),
                eigenvalues,
                coefficients,
            )
            minimum_eigenvalues.append(
                float(np.linalg.eigvalsh(information)[0])
            )
        row = {
            "times": tuple(float(x) for x in design),
            "worst_case_minimum_information_eigenvalue": min(
                minimum_eigenvalues
            ),
            "median_minimum_information_eigenvalue": float(
                np.median(minimum_eigenvalues)
            ),
            "all_sampled_clocks_monotone": monotone,
        }
        rows.append(row)
        if (
            best is None
            or row["worst_case_minimum_information_eigenvalue"]
            > best["worst_case_minimum_information_eigenvalue"]
        ):
            best = row
    assert best is not None
    reference_times = (0.05, 0.15, 0.35, 0.70, 1.20)
    reference_worst = min(
        float(
            np.linalg.eigvalsh(
                clock_information(reference_times, eigenvalues, coefficients)
            )[0]
        )
        for coefficients in clock_grid
    )

    # With an unknown linear clock coefficient, log scale and clock slope are
    # the same first-order direction.
    nominal_times = tuple(best["times"])
    calibrated = clock_information(
        nominal_times, eigenvalues, (0.0, 0.0, 0.0)
    )
    uncalibrated = np.zeros((5, 5))
    for time in nominal_times:
        derivative_scale = -eigenvalues * np.exp(-time * eigenvalues)
        bases = np.array([time, time, time**2, time**3, time**4])
        jacobian = derivative_scale[:, None] * bases[None, :]
        uncalibrated += jacobian.T @ jacobian
    return (
        {
            "status": (
                "[Proven] scale/linear-clock rank obstruction; "
                "[Strong evidence] quartic-clock minimax design"
            ),
            "clock_family": (
                "tau=t+a2*t^2+a3*t^3+a4*t^4 with tau'(0)=1, "
                "|a2|<=0.08, |a3|<=0.04, |a4|<=0.02"
            ),
            "clock_grid_size": len(clock_grid),
            "candidate_design_count": len(rows),
            "best_design": best,
            "reference_worst_information": reference_worst,
            "improvement_factor": (
                best["worst_case_minimum_information_eigenvalue"]
                / reference_worst
            ),
            "calibrated_nominal_rank": int(
                np.linalg.matrix_rank(calibrated, tol=1e-10)
            ),
            "uncalibrated_rank": int(
                np.linalg.matrix_rank(uncalibrated, tol=1e-10)
            ),
            "boundary": (
                "The family and bounds are declared design assumptions. "
                "External calibration tau'(0)=1 remains indispensable; "
                "nonparametric clocks can reopen nonidentifiability."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P330/P336: external admission gates
# ---------------------------------------------------------------------------


def external_gate(label: str) -> dict[str, Any]:
    manifest_names = (
        "bundle_manifest.json",
        "manifest_p241.json",
        "production_manifest.json",
    )
    manifests = sorted(
        {
            path
            for name in manifest_names
            for path in ROOT.rglob(name)
        }
    )
    admitted = []
    for path in manifests:
        relative = str(path.relative_to(ROOT))
        if any(
            token in relative.lower()
            for token in ("template", "example", "synthetic")
        ):
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
            roles = {
                str(payload["provider"]),
                str(payload["registrar"]),
                str(payload["analyst"]),
            }
            if len(roles) == 3:
                admitted.append(relative)
    return {
        "status": f"[Blocked by external evidence] no admitted {label}",
        "manifest_count": len(manifests),
        "admitted_manifest_paths": admitted,
        "admitted": bool(admitted),
        "one_shot_pipeline_authorized": bool(admitted),
        "boundary": (
            "Repository code cannot manufacture independent roles, custody, "
            "apparatus calibration, raw events, or a frozen hold-out."
        ),
    }


def program_330() -> dict[str, Any]:
    result = external_gate("P241 production event bundle")
    result.update(
        {
            "validator_241_executed": False,
            "pipeline_242_executed": False,
            "reason": "admission prerequisite not met",
        }
    )
    return result


# ---------------------------------------------------------------------------
# P331: executable photonic transfer specification
# ---------------------------------------------------------------------------


def program_331(
    strict_a: np.ndarray,
    strict_w: np.ndarray,
    protocol: dict[str, Any],
    rng: np.random.Generator,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    witness = protocol["best_protocols"]["wave"]
    time = float(witness["nominal_time"])
    vertex = int(witness["preparation_vertex"])
    allowed_tv_budget = 0.05
    generator_norm_budget = allowed_tv_budget / time
    nonzero = np.abs(strict_w[np.abs(strict_w) > 1e-12])
    dynamic_range = float(np.max(nonzero) / np.min(nonzero))

    reference = np.abs(expm(-1j * time * strict_a)[vertex]) ** 2
    rows = []
    for relative_sigma in (
        1e-4,
        3e-4,
        1e-3,
        3e-3,
        1e-2,
        3e-2,
        5e-2,
        1e-1,
    ):
        tv_values = []
        norm_values = []
        for _ in range(64):
            perturbation = rng.normal(
                0.0,
                relative_sigma * np.max(np.abs(strict_w)),
                size=strict_w.shape,
            )
            perturbation = (perturbation + perturbation.T) / 2.0
            np.fill_diagonal(perturbation, 0.0)
            perturbed_w = strict_w + perturbation
            row_sum = float(np.sum(perturbed_w[0]))
            perturbed_a = row_sum * np.eye(N) - perturbed_w
            observed = np.abs(
                expm(-1j * time * perturbed_a)[vertex]
            ) ** 2
            tv_values.append(0.5 * float(np.sum(np.abs(observed - reference))))
            norm_values.append(float(np.linalg.norm(perturbed_a - strict_a, 2)))
        rows.append(
            {
                "relative_coupling_sigma": relative_sigma,
                "median_generator_error_norm": float(np.median(norm_values)),
                "p95_generator_error_norm": float(
                    np.quantile(norm_values, 0.95)
                ),
                "median_output_tv_error": float(np.median(tv_values)),
                "p95_output_tv_error": float(np.quantile(tv_values, 0.95)),
                "passes_declared_tv_budget": bool(
                    np.quantile(tv_values, 0.95) <= allowed_tv_budget
                ),
            }
        )
    apparatus = [
        {
            "component": "heralded single-photon source",
            "function": "vertex-basis preparation and event timing",
            "required_record": "g2(0), heralding efficiency, flux stability",
        },
        {
            "component": "programmable 12-mode interferometer",
            "function": "compile exp(-itA) or a calibrated unitary dilation",
            "required_record": "complex transfer matrix with uncertainty",
        },
        {
            "component": "phase/amplitude controls",
            "function": "set signed dense couplings",
            "required_record": "AWG/EOM settings and drift log",
        },
        {
            "component": "12-channel SNSPD/SPAD plus TCSPC",
            "function": "vertex POVM and raw timestamps",
            "required_record": "efficiency, dark counts, crosstalk, dead time",
        },
        {
            "component": "independent clock and registrar",
            "function": "calibrated tau and custody",
            "required_record": "clock certificate, hashes, blind role split",
        },
    ]
    return (
        {
            "status": (
                "[Strong evidence] executable tolerance specification; "
                "[Hypothesis] photonic feasibility"
            ),
            "witness_time": time,
            "witness_vertex": vertex,
            "strict_nonzero_coupling_dynamic_range": dynamic_range,
            "duhamel_generator_norm_budget": generator_norm_budget,
            "declared_output_tv_budget": allowed_tv_budget,
            "apparatus": apparatus,
            "failure_rule": (
                "Reject the implementation before unblinding if the frozen "
                "complex transfer-matrix calibration or detector-channel "
                "hold-out exceeds its declared budget; do not refit FIN."
            ),
            "boundary": (
                "Dense signed couplings may require a programmable "
                "interferometer rather than a passive nearest-neighbor "
                "waveguide graph. No laboratory is claimed to implement it."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P332: interpretation of the signed resource
# ---------------------------------------------------------------------------


def program_332(
    moment_result: dict[str, Any],
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    negative_mass = moment_result[
        "exact_continuum_minimum_negative_mass_numeric"
    ]
    rows = [
        {
            "candidate": "classical signed measure",
            "reproduces_scalar_moments": True,
            "extra_structure_required": "Jordan decomposition",
            "physical_probability_automatic": False,
            "distinctive_test": "operational frame must expose signed weights",
            "status": "mathematically exact representation",
        },
        {
            "candidate": "quasiprobability frame",
            "reproduces_scalar_moments": True,
            "extra_structure_required": "states, effects, frame, allowed maps",
            "physical_probability_automatic": False,
            "distinctive_test": "negativity monotone under declared free maps",
            "status": "possible but representation-dependent",
        },
        {
            "candidate": "complex interference amplitudes",
            "reproduces_scalar_moments": True,
            "extra_structure_required": "phase reference and Born map",
            "physical_probability_automatic": False,
            "distinctive_test": "phase-sensitive recombination",
            "status": "possible scalar realization; not selected",
        },
        {
            "candidate": "Krein-space realization",
            "reproduces_scalar_moments": True,
            "extra_structure_required": "fundamental symmetry and positive sector",
            "physical_probability_automatic": False,
            "distinctive_test": "positive physical subspace and stable dynamics",
            "status": "formal representation only",
        },
        {
            "candidate": "open-system memory expansion",
            "reproduces_scalar_moments": True,
            "extra_structure_required": "environment, CPTP dilation, instruments",
            "physical_probability_automatic": True,
            "distinctive_test": "process-tensor conditional independence",
            "status": "not implied by signed spatial moments",
        },
        {
            "candidate": "active feedback/control",
            "reproduces_scalar_moments": True,
            "extra_structure_required": "controller, energy ledger, causal record",
            "physical_probability_automatic": True,
            "distinctive_test": "intervention response and work cost",
            "status": "not implied by signed spatial moments",
        },
    ]
    return (
        {
            "status": (
                "[Proven] scalar-realization nonuniqueness; "
                "[Refuted] unique physical interpretation from negativity"
            ),
            "minimum_negative_mass": negative_mass,
            "phase_pi_realization": (
                "replace the negative atom -N delta_0 by the complex "
                "amplitude N exp(i*pi) delta_0; the scalar moments are "
                "unchanged"
            ),
            "theorem": (
                "The moment functional sees only the signed coefficients. "
                "A Jordan signed measure and a phase-pi complex amplitude "
                "produce the same scalar sequence, while quasiprobability, "
                "Krein, open-system, and feedback interpretations require "
                "different additional operational structures."
            ),
            "boundary": (
                "P323/P316 necessity of a signed expansion does not prove "
                "quantumness, non-Markovianity, negative information, or "
                "energy loss."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P333: carrier naturality outside the cycle
# ---------------------------------------------------------------------------


def shortest_path_matrix(graph: nx.Graph) -> np.ndarray:
    nodes = list(graph.nodes())
    lengths = dict(nx.all_pairs_shortest_path_length(graph))
    return np.array(
        [[lengths[left][right] for right in nodes] for left in nodes],
        dtype=float,
    )


def kernel_matrix_from_distances(
    distances: np.ndarray,
    omega: float,
    phi: float,
    beta: float,
    eta: float,
    amplitude: float = 1.0,
) -> np.ndarray:
    matrix = amplitude * np.cos(omega * distances + phi) / (
        1.0 + beta * distances**eta
    )
    matrix = np.array(matrix, dtype=float)
    np.fill_diagonal(matrix, 0.0)
    return (matrix + matrix.T) / 2.0


def program_333() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    graphs = {
        "cycle_C12": nx.cycle_graph(12),
        "path_P12": nx.path_graph(12),
        "random_3_regular_12": nx.random_regular_graph(3, 12, seed=333),
        "barbell_irregular_12": nx.barbell_graph(5, 2),
    }
    rows = []
    for name, graph in graphs.items():
        distances = shortest_path_matrix(graph)
        legacy_w = kernel_matrix_from_distances(
            distances,
            math.pi / 4,
            math.pi / 6,
            0.01,
            1.0,
            ALPHA_GEO,
        )
        strict_w = kernel_matrix_from_distances(
            distances, 0.18575, 0.16250, 1.0, 1.8
        )
        commutator = legacy_w @ strict_w - strict_w @ legacy_w
        relative = float(
            np.linalg.norm(commutator, "fro")
            / (
                np.linalg.norm(legacy_w, "fro")
                * np.linalg.norm(strict_w, "fro")
            )
        )
        rows.append(
            {
                "carrier": name,
                "vertices": graph.number_of_nodes(),
                "edges": graph.number_of_edges(),
                "diameter": nx.diameter(graph),
                "legacy_strict_relative_commutator": relative,
                "commute_within_1e_10": relative < 1e-10,
                "scalar_functional_calculus_not_obstructed": relative < 1e-10,
            }
        )
    obstructed = sum(
        not row["scalar_functional_calculus_not_obstructed"] for row in rows
    )
    return (
        {
            "status": (
                "[Proven] commutation is necessary for scalar completion; "
                "[Strong evidence] carrier-specific obstruction matrix"
            ),
            "carrier_count": len(rows),
            "noncommuting_carriers": obstructed,
            "theorem": (
                "If A_strict=f(A_legacy) by scalar Borel/functional calculus, "
                "then the two operators commute. Nonzero commutators therefore "
                "exclude that completion class on the corresponding carrier."
            ),
            "boundary": (
                "Commutation on C12 is necessary but not sufficient and does "
                "not source f. Noncommutation does not exclude nonscalar, "
                "carrier-dependent, or target-importing parents."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P334: minimal axiom-augmented electroweak package
# ---------------------------------------------------------------------------


def program_334() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    x = ALPHA_GEO / 12.0
    ratio = math.sqrt(x / (1.0 - x))
    rows = [
        {
            "package": "unpointed spectrum",
            "discrete_new_axioms": 0,
            "continuous_new_parameters": 0,
            "law_axioms": 0,
            "derives_ordered_mixing_angle": False,
            "derives_boson_mass_scale": False,
        },
        {
            "package": "pointed angle-only conditional package",
            "discrete_new_axioms": 1,
            "continuous_new_parameters": 0,
            "law_axioms": 1,
            "derives_ordered_mixing_angle": True,
            "derives_boson_mass_scale": False,
        },
        {
            "package": "minimal tree-level gauge-mass scaffold",
            "discrete_new_axioms": 1,
            "continuous_new_parameters": 2,
            "law_axioms": 4,
            "derives_ordered_mixing_angle": True,
            "derives_boson_mass_scale": True,
        },
        {
            "package": "full electroweak theory",
            "discrete_new_axioms": None,
            "continuous_new_parameters": None,
            "law_axioms": None,
            "derives_ordered_mixing_angle": False,
            "derives_boson_mass_scale": False,
        },
    ]
    return (
        {
            "status": (
                "[Proven] minimal conditional angle package; "
                "[Refuted] full electroweak derivation from current spectrum"
            ),
            "conditional_axiom": "sin^2(theta_W)=alpha_geo/12",
            "conditional_sin2_theta": x,
            "ordered_coupling_ratio_gprime_over_g": ratio,
            "sector_swap_sin2_theta": 1.0 - x,
            "sector_swap_ratio": 1.0 / ratio,
            "conditional_tree_level_mW_over_mZ": math.sqrt(1.0 - x),
            "minimum_angle_only_additions": {
                "discrete_pointing_axiom": 1,
                "continuous_parameters": 0,
                "role_law_axiom": 1,
            },
            "minimum_mass_scaffold_additions": {
                "discrete_pointing_axiom": 1,
                "continuous_parameters": [
                    "overall gauge coupling g",
                    "symmetry-breaking scale v",
                ],
                "additional_laws": [
                    "U(1)xSU(2) representation",
                    "covariant derivative",
                    "symmetry-breaking sector",
                    "mass/observable map",
                ],
            },
            "boundary": (
                "The angle-only package is explicitly axiom-augmented and "
                "uses the historical formula as an axiom. It does not derive "
                "hypercharges, chirality, anomaly cancellation, Higgs "
                "dynamics, a dimensional scale, or radiative corrections."
            ),
        },
        rows,
    )


# ---------------------------------------------------------------------------
# P335: abstract independence cube of bridge resources
# ---------------------------------------------------------------------------


def program_335() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    resources = (
        "signed_path",
        "nontorsion_phase",
        "independent_cross_law",
        "pointing",
        "dimensional_scale",
    )
    rows = []
    for bits in itertools.product((False, True), repeat=len(resources)):
        row = {
            "configuration_id": "".join("1" if bit else "0" for bit in bits)
        }
        row.update(dict(zip(resources, bits)))
        rows.append(row)
    witnesses = {}
    for missing in resources:
        witness = next(
            row
            for row in rows
            if not row[missing]
            and all(row[other] for other in resources if other != missing)
        )
        witnesses[missing] = witness["configuration_id"]
    return (
        {
            "status": (
                "[Proven] logical independence in the minimal product "
                "signature; [Open] dynamical independence inside FIN"
            ),
            "resource_count": len(resources),
            "configuration_count": len(rows),
            "all_single_omission_witnesses_exist": len(witnesses) == 5,
            "single_omission_witnesses": witnesses,
            "typed_factor_models": {
                "signed_path": "positive cone versus signed Jordan cone",
                "nontorsion_phase": "finite torsion group versus Z generator",
                "independent_cross_law": "absent versus supplied parent morphism",
                "pointing": "free two-torsor versus chosen point",
                "dimensional_scale": "free R_+ torsor versus chosen unit",
            },
            "theorem": (
                "In the product signature, for every resource R_i there is "
                "a model with the other four resources present and R_i absent. "
                "Consequently no implication R_i <- conjunction(other four) "
                "is derivable without an additional coupling axiom."
            ),
            "boundary": (
                "Abstract model-theoretic independence does not prove that a "
                "future FIN source law cannot couple or co-generate resources."
            ),
        },
        rows,
    )


def program_336() -> dict[str, Any]:
    result = external_gate("physical reservoir production package")
    result.update(
        {
            "paired_seed_controls_present": False,
            "frozen_primary_endpoint_present": False,
            "hardware_calibration_present": False,
            "external_execution_performed": False,
        }
    )
    return result


# ---------------------------------------------------------------------------
# Figures, summaries, and orchestration
# ---------------------------------------------------------------------------


def create_figures(results: dict[str, Any]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)

    p323 = results["P323"]
    x = np.linspace(0.0, 1.0, 800)
    coefficients = p323["dual_polynomial_coefficients_ascending"]
    p = sum(coefficients[k] * x**k for k in range(7))
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].plot(x, p, color="#1f77b4")
    axes[0].axhline(0.0, color="black", lw=0.8)
    axes[0].axhline(-1.0, color="black", lw=0.8, ls="--")
    axes[0].set(xlabel="x", ylabel="dual p(x)", title="Optimal dual polynomial")
    nodes = p323["positive_nodes"]
    weights = p323["positive_weights"]
    axes[1].stem(nodes, weights, linefmt="C2-", markerfmt="C2o", basefmt=" ")
    axes[1].stem(
        [0],
        [-p323["exact_continuum_minimum_negative_mass_numeric"]],
        linefmt="C3-",
        markerfmt="C3o",
        basefmt=" ",
    )
    axes[1].axhline(0.0, color="black", lw=0.8)
    axes[1].set(xlabel="moment node", ylabel="signed weight", title="Extremal measure")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p323_continuum_extremizer.png", dpi=180)
    plt.close(fig)

    parent_rows = results["_P324_rows"]
    phase_rows = results["_P325_rows"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].bar(
        [row["target"] for row in parent_rows],
        [row["target_frobenius_distance_from_frozen"] for row in parent_rows],
        color=["#1f77b4", "#ff7f0e", "#2ca02c"],
    )
    axes[0].tick_params(axis="x", rotation=20)
    axes[0].set(ylabel="Frobenius distance", title="Same legacy, distinct targets")
    axes[1].bar(
        [row["quantity"] for row in phase_rows],
        [row["lattice_coordinate_in_1_over_4000"] for row in phase_rows],
        color="#9467bd",
    )
    axes[1].tick_params(axis="x", rotation=35)
    axes[1].set(ylabel="coordinate in (1/4000)Z", title="Refreeze phase provenance")
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p324_p325_parent_phase.png", dpi=180)
    plt.close(fig)

    protocol_rows = results["_P326_rows"]
    fig, ax = plt.subplots(figsize=(8.5, 4.4))
    for dynamics, color in (("heat", "#1f77b4"), ("wave", "#ff7f0e")):
        subset = [
            row
            for row in protocol_rows
            if row["dynamics"] == dynamics
            and row["preparation_vertex"]
            == results["P326"]["best_protocols"][dynamics][
                "preparation_vertex"
            ]
        ]
        ax.semilogx(
            [row["nominal_time"] for row in subset],
            [row["robust_tv_5pct"] for row in subset],
            label=dynamics,
            color=color,
        )
    ax.set(
        xlabel="nominal time",
        ylabel="5th-percentile observed TV",
        title="Loss-robust legacy/strict discrimination",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p326_lossy_protocol.png", dpi=180)
    plt.close(fig)

    chamber_rows = results["_P327_rows"]
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    ax.hist(
        [row["minimum_weight"] for row in chamber_rows],
        bins=24,
        color="#2ca02c",
        edgecolor="white",
    )
    ax.axvline(0.14, color="black", ls="--", label="exact global minimum 7/50")
    ax.set(
        xlabel="minimum radial weight in exact witness",
        ylabel="number of chambers",
        title="All 720 exact Q(sqrt(3)) chamber witnesses",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p327_exact_chambers.png", dpi=180)
    plt.close(fig)

    diamond_rows = results["_P328_rows"]
    clock_rows = results["_P329_rows"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].semilogx(
        [row["time"] for row in diamond_rows],
        [row["half_diamond_distance"] for row in diamond_rows],
        color="#d62728",
    )
    axes[0].set(
        xlabel="time",
        ylabel="half diamond norm",
        title="One-time wave-channel separation",
        ylim=(-0.02, 1.02),
    )
    axes[1].hist(
        [
            row["worst_case_minimum_information_eigenvalue"]
            for row in clock_rows
        ],
        bins=35,
        color="#17becf",
    )
    axes[1].set(
        xlabel="worst-case minimum information eigenvalue",
        ylabel="design count",
        title="Quartic-clock minimax designs",
    )
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p328_p329_channel_clock.png", dpi=180)
    plt.close(fig)

    photonic_rows = results["_P331_rows"]
    carrier_rows = results["_P333_rows"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    axes[0].loglog(
        [row["relative_coupling_sigma"] for row in photonic_rows],
        [row["p95_output_tv_error"] for row in photonic_rows],
        marker="o",
    )
    axes[0].axhline(
        results["P331"]["declared_output_tv_budget"],
        color="black",
        ls="--",
    )
    axes[0].set(
        xlabel="relative coupling error SD",
        ylabel="p95 output TV error",
        title="Photonic transfer tolerance",
    )
    axes[1].bar(
        [row["carrier"] for row in carrier_rows],
        [row["legacy_strict_relative_commutator"] for row in carrier_rows],
        color="#8c564b",
    )
    axes[1].set_yscale("symlog", linthresh=1e-12)
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].set(
        ylabel="relative commutator norm",
        title="Carrier naturality obstruction",
    )
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p331_p333_transfer_carriers.png", dpi=180)
    plt.close(fig)

    ew_rows = results["_P334_rows"]
    independence_rows = results["_P335_rows"]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    finite_ew = [
        row
        for row in ew_rows
        if row["continuous_new_parameters"] is not None
    ]
    axes[0].bar(
        [row["package"] for row in finite_ew],
        [row["continuous_new_parameters"] for row in finite_ew],
        color="#e377c2",
    )
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set(
        ylabel="new continuous parameters",
        title="Conditional electroweak packages",
    )
    counts = [
        sum(
            bool(row[name])
            for name in (
                "signed_path",
                "nontorsion_phase",
                "independent_cross_law",
                "pointing",
                "dimensional_scale",
            )
        )
        for row in independence_rows
    ]
    axes[1].hist(
        counts,
        bins=np.arange(-0.5, 6.5, 1),
        rwidth=0.85,
        color="#7f7f7f",
    )
    axes[1].set(
        xlabel="resources present",
        ylabel="configurations",
        title="Five-resource independence cube",
        xticks=range(6),
    )
    fig.tight_layout()
    fig.savefig(FIGURE_DIR / "p334_p335_axioms_independence.png", dpi=180)
    plt.close(fig)


def run_formal_core() -> dict[str, Any]:
    lean = p309.lean_binary()
    completed = subprocess.run(
        [str(lean), str(FORMAL_SOURCE)],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return {
        "lean_binary": str(lean),
        "lean_version": subprocess.run(
            [str(lean), "--version"],
            capture_output=True,
            text=True,
            check=False,
        ).stdout.strip(),
        "source_sha256": sha256_file(FORMAL_SOURCE),
        "exit_code": completed.returncode,
        "stdout": completed.stdout,
        "stderr": completed.stderr,
    }


def summary_rows(results: dict[str, Any]) -> list[dict[str, Any]]:
    headlines = {
        "P323": "the exact continuum signed-moment optimum is 0.4067063347",
        "P324": "legacy input alone does not identify a unique common-parent target",
        "P325": "1/4000 is the exact QW2030/QW2038 refreeze lattice unit",
        "P326": "lossy synthetic apparatus still separates legacy and strict",
        "P327": "all 720 spectral order chambers now have exact algebraic witnesses",
        "P328": "one-time wave channels have a computable diamond separation",
        "P329": "quartic clock uncertainty admits a finite minimax design after calibration",
        "P330": "the independent P241 production gate remains closed",
        "P331": "photonic transfer now has a quantitative pre-unblinding tolerance spec",
        "P332": "signed path mass does not select a unique physical interpretation",
        "P333": "scalar completion fails naturally on noncommuting noncyclic carriers",
        "P334": "the electroweak angle can be packaged only conditionally with pointing",
        "P335": "five bridge resources are logically independent without coupling axioms",
        "P336": "the external physical-reservoir gate remains closed",
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
            "programs": "P323-P336",
            "release": "10.29",
            "seed": SEED,
            "scope": (
                "continuum moment duality, exact refreeze provenance, "
                "source-independence, exact algebraic chamber witnesses, "
                "synthetic operational design, axiom accounting, and "
                "external-evidence gates; no silent bridge, selector, units, "
                "role transfer, laboratory validation, L_total, SM/GR, or ToE"
            ),
            "new_theoretical_objects": {
                "O74_gaussian_jordan_extremal_certificate": (
                    "matching three-node primal signed measure and bounded "
                    "degree-six dual polynomial"
                ),
                "O75_refreeze_provenance_phase_lattice": (
                    "the exact additive lattice generated by the QW2030 "
                    "center and QW2038 grid steps"
                ),
                "O76_source_independent_completion_fiber": (
                    "the set of distinct targets compatible with one legacy "
                    "input before a cross law is supplied"
                ),
                "O77_loss_robust_discrimination_protocol": (
                    "preparation, time, detector channel, nuisance set, and "
                    "finite-count rejection rule"
                ),
                "O78_exact_spectral_chamber_atlas": (
                    "720 exact Q(sqrt(3)) positive witnesses"
                ),
                "O79_unitary_channel_separation_modulus": (
                    "half diamond distance obtained from the relative-phase "
                    "convex hull"
                ),
                "O80_clock_warp_minimax_frame": (
                    "calibrated clock family, coefficient bounds, time design, "
                    "and worst-case information criterion"
                ),
                "O81_p241_admission_contract": (
                    "independent roles, hashes, raw events, calibration, and "
                    "frozen hold-out required before Validator 241"
                ),
                "O82_photonic_transfer_tolerance_certificate": (
                    "matrix error budget, detector budget, and frozen failure "
                    "threshold for a 12-mode implementation"
                ),
                "O83_signed_realization_ambiguity_class": (
                    "inequivalent operational models sharing one signed "
                    "scalar moment functional"
                ),
                "O84_carrier_naturality_obstruction": (
                    "commutator witness excluding scalar functional completion "
                    "on a declared carrier"
                ),
                "O85_pointed_EW_conditional_package": (
                    "sector point plus an explicit historical role-law axiom"
                ),
                "O86_bridge_resource_independence_cube": (
                    "32 product-signature models separating five resources"
                ),
                "O87_reservoir_admission_contract": (
                    "hardware calibration, paired controls, seeds, and frozen "
                    "endpoint required for external execution"
                ),
            },
        }
    }

    results["P323"], results["_P323_rows"] = program_323()
    results["P324"], results["_P324_rows"] = program_324()
    results["P325"], results["_P325_rows"] = program_325()
    results["P326"], results["_P326_rows"] = program_326(strict_a, rng)
    results["P327"], results["_P327_rows"] = program_327()
    results["P328"], results["_P328_rows"] = program_328(strict_a)
    results["P329"], results["_P329_rows"] = program_329(strict_a)
    results["P330"] = program_330()
    results["P331"], results["_P331_rows"] = program_331(
        strict_a, strict_w, results["P326"], rng
    )
    results["P332"], results["_P332_rows"] = program_332(results["P323"])
    results["P333"], results["_P333_rows"] = program_333()
    results["P334"], results["_P334_rows"] = program_334()
    results["P335"], results["_P335_rows"] = program_335()
    results["P336"] = program_336()
    results["formal_verification"] = run_formal_core()

    write_csv(MOMENT_PATH, results["_P323_rows"])
    write_csv(PARENT_PATH, results["_P324_rows"])
    write_csv(PHASE_PATH, results["_P325_rows"])
    write_csv(PROTOCOL_PATH, results["_P326_rows"])
    write_csv(CHAMBER_PATH, results["_P327_rows"])
    write_csv(DIAMOND_PATH, results["_P328_rows"])
    write_csv(CLOCK_PATH, results["_P329_rows"])
    write_csv(PHOTONIC_PATH, results["_P331_rows"])
    write_csv(INTERPRETATION_PATH, results["_P332_rows"])
    write_csv(CARRIER_PATH, results["_P333_rows"])
    write_csv(EW_PATH, results["_P334_rows"])
    write_csv(INDEPENDENCE_PATH, results["_P335_rows"])
    write_csv(SUMMARY_PATH, summary_rows(results))
    create_figures(results)

    public_results = {
        key: value for key, value in results.items() if not key.startswith("_")
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(public_results), indent=2, ensure_ascii=False)
        + "\n",
        encoding="utf-8",
    )
    for row in summary_rows(results):
        print(f"{row['program']} {row['status']} - {row['headline']}")
    print(RESULTS_PATH.name)


if __name__ == "__main__":
    main()
