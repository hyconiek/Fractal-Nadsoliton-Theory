#!/usr/bin/env python3
"""Local FIN checkpoint batch P435, P436, and P440.

P435 gives an exact one-slot process-tester SDP certificate, constructs the
two-slot comb instance, and records the local multi-slot solver obstruction.
P436 proves a positive heralded-erasure code gain in one declared cell using
exact rational interval arithmetic and exact characteristic-polynomial root
isolation.
P440 derives a detector-envelope minimax Jordan sampler and a finite-sample
confidence ledger.

No routine reads network resources, laboratory records, or external audits.
"""

from __future__ import annotations

import csv
from fractions import Fraction
from math import comb, factorial, isqrt
import importlib.util
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import sympy as sp


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_435_436_440"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P435_PATH = ROOT / f"{PREFIX}_P435_Comb_Certificate.csv"
P436_PATH = ROOT / f"{PREFIX}_P436_Erasure_Intervals.csv"
P440_PATH = ROOT / f"{PREFIX}_P440_Detector_Sampler.csv"
COMB_INSTANCE_PATH = ROOT / f"{PREFIX}_P435_Two_Slot_Comb_Instance.npz"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p435_p436_p440_certificates.png"
P429_CSV = ROOT / "FIN_Programs_428_430_P429_Krawczyk.csv"
P417_CSV = ROOT / "FIN_Programs_411_427_Noisy_Comb_Gap.csv"


def json_ready(value: Any) -> Any:
    if isinstance(value, Fraction):
        # Exact characteristic-polynomial isolators can have denominators with
        # many thousands of digits.  The executable certificate retains those
        # rationals in memory; result tables use compact decimal endpoints once
        # the exact representation would cease to be a useful interchange
        # format.
        if max(value.numerator.bit_length(), value.denominator.bit_length()) <= 512:
            return str(value)
        return format(float(value), ".17g")
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    return value


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: (
                json.dumps(json_ready(row.get(key)), ensure_ascii=False)
                if isinstance(row.get(key), (list, tuple, dict))
                else json_ready(row.get(key, ""))
            ) for key in keys})


# ---------------------------------------------------------------------------
# P435: process-tester SDP admission certificate and two-slot instance
# ---------------------------------------------------------------------------


def partial_trace_square(matrix: np.ndarray, dims: tuple[int, ...], axis: int) -> np.ndarray:
    """Trace one subsystem from a square operator with the declared ordering."""

    shaped = matrix.reshape(*dims, *dims)
    traced = np.trace(shaped, axis1=axis, axis2=axis + len(dims))
    remaining = int(np.prod([dim for index, dim in enumerate(dims) if index != axis]))
    return traced.reshape(remaining, remaining)


def qubit_phase_dephasing_choi(coherence: float, theta: float, sign: int) -> np.ndarray:
    """Unnormalised Choi matrix in output-input ordering."""

    choi = np.zeros((4, 4), dtype=complex)
    choi[0, 0] = 1.0
    choi[3, 3] = 1.0
    choi[0, 3] = coherence * np.exp(-1j * sign * theta)
    choi[3, 0] = np.conjugate(choi[0, 3])
    return choi


def program_435() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    coherence = 4.0 / 5.0
    theta = math.pi / 8.0
    j_plus = qubit_phase_dephasing_choi(coherence, theta, +1)
    j_minus = qubit_phase_dephasing_choi(coherence, theta, -1)

    # Exact one-slot tester: |+><+| input and the Helstrom output POVM.
    rho = np.full((2, 2), 0.5, dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    m_plus = (np.eye(2) + sigma_y) / 2
    m_minus = (np.eye(2) - sigma_y) / 2
    t_plus = np.kron(m_plus, rho.T)
    t_minus = np.kron(m_minus, rho.T)
    primal = 0.5 * np.trace(t_plus @ j_plus).real + 0.5 * np.trace(t_minus @ j_minus).real

    # Exact common-upper-bound dual: Y=(J+ + J- + |J+ - J-|)/4.
    difference = j_plus - j_minus
    evals, evecs = np.linalg.eigh(difference)
    absolute_difference = (evecs * np.abs(evals)) @ evecs.conj().T
    dual_y = (j_plus + j_minus + absolute_difference) / 4
    reduced_y = partial_trace_square(dual_y, (2, 2), 0)
    dual = float(np.max(np.linalg.eigvalsh(reduced_y)))
    exact_formula = "(1 + (4/5)*sin(pi/8))/2"
    exact_value = (1 + coherence * math.sin(theta)) / 2

    one_slot_checks = {
        "choi_plus_causality_defect": float(np.linalg.norm(partial_trace_square(j_plus, (2, 2), 0) - np.eye(2))),
        "choi_minus_causality_defect": float(np.linalg.norm(partial_trace_square(j_minus, (2, 2), 0) - np.eye(2))),
        "tester_normalization_defect": float(np.linalg.norm(t_plus + t_minus - np.kron(np.eye(2), rho.T))),
        "minimum_tester_eigenvalue": float(min(np.linalg.eigvalsh(t_plus).min(), np.linalg.eigvalsh(t_minus).min())),
        "minimum_dual_slack_eigenvalue": float(min(np.linalg.eigvalsh(dual_y - j_plus / 2).min(), np.linalg.eigvalsh(dual_y - j_minus / 2).min())),
        "dual_reduction_scalar_defect": float(np.linalg.norm(reduced_y - exact_value * np.eye(2))),
        "primal_dual_gap": float(dual - primal),
    }

    # The intended n=2 memoryless process instance in order Y2,X2,Y1,X1.
    r_plus = np.kron(j_plus, j_plus)
    r_minus = np.kron(j_minus, j_minus)
    causality_plus = partial_trace_square(r_plus, (2, 2, 2, 2), 0)
    causality_minus = partial_trace_square(r_minus, (2, 2, 2, 2), 0)
    target_plus = np.kron(np.eye(2), j_plus)
    target_minus = np.kron(np.eye(2), j_minus)
    two_slot_checks = {
        "comb_plus_causality_defect": float(np.linalg.norm(causality_plus - target_plus)),
        "comb_minus_causality_defect": float(np.linalg.norm(causality_minus - target_minus)),
        "process_dimension": 16,
        "primal_real_variables_before_equalities": 580,
        "dual_real_variables_before_inequality_reduction": 273,
    }
    np.savez_compressed(
        COMB_INSTANCE_PATH,
        J_plus=j_plus,
        J_minus=j_minus,
        R_plus=r_plus,
        R_minus=r_minus,
        coherence=np.array([coherence]),
        theta=np.array([theta]),
    )

    installed = {
        module: importlib.util.find_spec(module) is not None
        for module in ("cvxpy", "cvxopt", "clarabel", "scs")
    }
    with P417_CSV.open(encoding="utf-8") as handle:
        selected = next(
            row for row in csv.DictReader(handle)
            if row["uses"] == "2" and row["coherence"] == "0.8"
            and row["threshold_fraction"] == "0.5"
        )
    two_slot_lower = float(selected["optimized_parallel_lower"])
    two_slot_upper = float(selected["adaptive_hybrid_upper"])

    rows = [
        {
            "scope": "one_slot_exact",
            "uses": 1,
            "coherence": coherence,
            "theta": theta,
            "primal_success": primal,
            "dual_success": dual,
            "gap": dual - primal,
            "certificate": exact_formula,
            "status": "Proven",
        },
        {
            "scope": "two_slot_comb_instance",
            "uses": 2,
            "coherence": coherence,
            "theta": theta,
            "primal_success": (1 + two_slot_lower) / 2,
            "dual_success": (1 + two_slot_upper) / 2,
            "gap": (two_slot_upper - two_slot_lower) / 2,
            "certificate": "finite comb instance exported; no installed local SDP solver",
            "status": "Open / locally blocked",
        },
    ]
    return ({
        "status": (
            "[Proven] exact one-slot process-tester SDP primal/dual equality at a nonideal cell; "
            "[Constructed] explicit two-slot memoryless comb and primal/dual SDP; "
            "[Blocked locally] no installed SDP solver or checker closes the two-slot gap"
        ),
        "one_slot_exact_formula": exact_formula,
        "one_slot_optimal_success": exact_value,
        "one_slot_checks": one_slot_checks,
        "two_slot_checks": two_slot_checks,
        "two_slot_parallel_half_distance_lower": two_slot_lower,
        "two_slot_hybrid_half_distance_upper": two_slot_upper,
        "two_slot_success_gap": (two_slot_upper - two_slot_lower) / 2,
        "solver_inventory": installed,
        "two_slot_primal_constraints": [
            "T_plus,T_minus >= 0",
            "T_plus+T_minus = I_Y2 tensor B2",
            "Tr_X2(B2) = I_Y1 tensor rho",
            "rho >= 0 and Tr(rho)=1",
        ],
        "two_slot_dual_constraints": [
            "Y >= R_plus/2 and Y >= R_minus/2",
            "Tr_Y2(Y) <= I_X2 tensor Z",
            "Tr_Y1(Z) <= lambda I_X1",
        ],
        "boundary": (
            "The exact certificate proves only the one-slot reduced qubit channel. "
            "The exported two-slot instance is a genuine comb specification, but no matching "
            "multi-slot dual certificate is claimed. It is not the full twelve-mode FIN channel."
        ),
        "new_object": "O158 Process-Tester SDP Admission Certificate",
    }, rows)


# ---------------------------------------------------------------------------
# P436: rational interval certificate for one O149 gain cell
# ---------------------------------------------------------------------------


Interval = tuple[Fraction, Fraction]


def iv_add(left: Interval, right: Interval) -> Interval:
    return left[0] + right[0], left[1] + right[1]


def iv_mul(left: Interval, right: Interval) -> Interval:
    values = [left[i] * right[j] for i in (0, 1) for j in (0, 1)]
    return min(values), max(values)


def iv_scale(value: Interval, scalar: Fraction | int) -> Interval:
    return iv_mul(value, (Fraction(scalar), Fraction(scalar)))


def rational_sqrt_interval(value: Fraction, decimal_places: int = 24) -> Interval:
    if value < 0:
        raise ValueError("square root of a negative rational")
    if value == 0:
        return Fraction(0), Fraction(0)
    scale = 10**decimal_places
    root = isqrt((value.numerator * scale * scale) // value.denominator)
    while Fraction((root + 1) ** 2, scale * scale) <= value:
        root += 1
    while Fraction(root**2, scale * scale) > value:
        root -= 1
    return Fraction(root, scale), Fraction(root + 1, scale)


def atan_inverse_interval(denominator: int, odd_index: int) -> Interval:
    """Alternating-series enclosure of atan(1/denominator)."""

    if odd_index % 2 != 1:
        raise ValueError("the lower partial-sum index must be odd")
    x = Fraction(1, denominator)

    def partial(last_index: int) -> Fraction:
        return sum(
            (-1) ** index * x ** (2 * index + 1) / Fraction(2 * index + 1)
            for index in range(last_index + 1)
        )

    return partial(odd_index), partial(odd_index + 1)


def pi_interval() -> Interval:
    """Machin enclosure pi=16 atan(1/5)-4 atan(1/239)."""

    a5 = atan_inverse_interval(5, 59)
    a239 = atan_inverse_interval(239, 19)
    return 16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0]


PI_INTERVAL = pi_interval()


def sine_point_interval(value: Fraction, odd_index: int = 19) -> Interval:
    """Alternating sine enclosure for 0 <= value <= 2*pi/5."""

    def partial(last_index: int) -> Fraction:
        return sum(
            (-1) ** index * value ** (2 * index + 1) / Fraction(factorial(2 * index + 1))
            for index in range(last_index + 1)
        )

    return partial(odd_index), partial(odd_index + 1)


def p436_sine_interval(weight_difference: int) -> Interval:
    """Enclose sin((2*pi/15)*weight_difference), |difference|<=3."""

    if weight_difference == 0:
        return Fraction(0), Fraction(0)
    if weight_difference < 0:
        positive = p436_sine_interval(-weight_difference)
        return -positive[1], -positive[0]
    low_x = PI_INTERVAL[0] * Fraction(2 * weight_difference, 15)
    high_x = PI_INTERVAL[1] * Fraction(2 * weight_difference, 15)
    return sine_point_interval(low_x)[0], sine_point_interval(high_x)[1]


def skew_matrix_intervals(probabilities: list[Fraction], survivors: int) -> list[list[Interval]]:
    uses = 3
    lost = uses - survivors
    dimension = 2**survivors
    coherence = Fraction(4, 5)
    matrix = [[(Fraction(0), Fraction(0)) for _ in range(dimension)] for _ in range(dimension)]
    for left_index in range(dimension):
        left_bits = [(left_index >> bit) & 1 for bit in range(survivors)]
        for right_index in range(dimension):
            right_bits = [(right_index >> bit) & 1 for bit in range(survivors)]
            hamming = sum(left != right for left, right in zip(left_bits, right_bits))
            difference = sum(left_bits) - sum(right_bits)
            amplitude_sum: Interval = (Fraction(0), Fraction(0))
            for lost_value in range(2**lost):
                lost_weight = lost_value.bit_count()
                left_weight = sum(left_bits) + lost_weight
                right_weight = sum(right_bits) + lost_weight
                radicand = (
                    probabilities[left_weight]
                    * probabilities[right_weight]
                    / Fraction(comb(uses, left_weight) * comb(uses, right_weight))
                )
                amplitude_sum = iv_add(amplitude_sum, rational_sqrt_interval(radicand))
            matrix[left_index][right_index] = iv_scale(
                iv_mul(amplitude_sum, p436_sine_interval(difference)),
                2 * coherence**hamming,
            )
    return matrix


def exact_midpoint_trace_distance_bounds(
    intervals: list[list[Interval]], decimal_places: int = 16
) -> tuple[Fraction, Fraction, dict[str, Any]]:
    """Bound one half nuclear norm of an interval real skew matrix.

    A rounded rational midpoint is handled by exact characteristic-polynomial
    root isolation for K^T K.  The interval remainder is added with the
    conservative bound ||E||_*/2 <= d^(3/2) max|E_ij|/2 <= d^2 max|E_ij|/2.
    """

    dimension = len(intervals)
    scale = 10**decimal_places
    midpoint: list[list[Fraction]] = []
    maximum_radius = Fraction(0)
    for row in intervals:
        midpoint_row = []
        for low, high in row:
            center = (low + high) / 2
            rounded = Fraction(round(float(center) * scale), scale)
            maximum_radius = max(maximum_radius, abs(rounded - low), abs(high - rounded))
            midpoint_row.append(rounded)
        midpoint.append(midpoint_row)

    rational_matrix = sp.Matrix([
        [sp.Rational(value.numerator, value.denominator) for value in row]
        for row in midpoint
    ])
    gram = rational_matrix.T * rational_matrix
    polynomial = gram.charpoly().as_poly()
    isolated = sp.intervals(polynomial, eps=sp.Rational(1, 10**20))
    singular_intervals: list[Interval] = []
    for (low_high, multiplicity) in isolated:
        low_root, high_root = low_high
        low_root = max(low_root, sp.Rational(0))
        high_root = max(high_root, sp.Rational(0))
        low_fraction = Fraction(int(low_root.p), int(low_root.q))
        high_fraction = Fraction(int(high_root.p), int(high_root.q))
        square_root = (
            rational_sqrt_interval(low_fraction, 20)[0],
            rational_sqrt_interval(high_fraction, 20)[1],
        )
        singular_intervals.extend([square_root] * multiplicity)
    if len(singular_intervals) != dimension:
        raise RuntimeError("root isolation did not return the full multiplicity")

    midpoint_lower = sum(value[0] for value in singular_intervals) / 2
    midpoint_upper = sum(value[1] for value in singular_intervals) / 2
    perturbation = Fraction(dimension * dimension, 2) * maximum_radius
    return midpoint_lower - perturbation, midpoint_upper + perturbation, {
        "dimension": dimension,
        "isolating_intervals": len(isolated),
        "root_multiplicity": len(singular_intervals),
        "maximum_entry_radius": maximum_radius,
        "nuclear_half_norm_perturbation_bound": perturbation,
    }


def program_436() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    laws = {
        "rational_O149_candidate": [
            Fraction(227692, 10**6), Fraction(272308, 10**6),
            Fraction(272308, 10**6), Fraction(227692, 10**6),
        ],
        "product": [Fraction(1, 8), Fraction(3, 8), Fraction(3, 8), Fraction(1, 8)],
        "GHZ": [Fraction(1, 2), Fraction(0), Fraction(0), Fraction(1, 2)],
    }
    totals: dict[str, Interval] = {}
    rows: list[dict[str, Any]] = []
    for name, probabilities in laws.items():
        total: Interval = (Fraction(0), Fraction(0))
        for survivors in range(1, 4):
            trace_bounds = exact_midpoint_trace_distance_bounds(
                skew_matrix_intervals(probabilities, survivors)
            )
            low, high, detail = trace_bounds
            mask_probability = (
                Fraction(comb(3, survivors))
                * Fraction(4, 5) ** survivors
                * Fraction(1, 5) ** (3 - survivors)
            )
            total = iv_add(total, iv_scale((low, high), mask_probability))
            rows.append({
                "law": name,
                "survivors": survivors,
                "matrix_dimension": 2**survivors,
                "mask_probability_total": mask_probability,
                "trace_distance_lower": low,
                "trace_distance_upper": high,
                "interval_width": high - low,
                **detail,
            })
        totals[name] = total
        rows.append({
            "law": name,
            "survivors": "weighted_total",
            "matrix_dimension": "orthogonal herald blocks",
            "mask_probability_total": Fraction(1),
            "trace_distance_lower": total[0],
            "trace_distance_upper": total[1],
            "interval_width": total[1] - total[0],
        })

    candidate = totals["rational_O149_candidate"]
    baseline_upper = max(totals["product"][1], totals["GHZ"][1])
    certified_gain = candidate[0] - baseline_upper
    if certified_gain <= 0:
        raise RuntimeError("P436 failed to certify a positive gain")
    return ({
        "status": (
            "[Computer-assisted proof] the declared rational symmetric code strictly "
            "beats both product and GHZ baselines in the n=3, q=eta=4/5, theta=2*pi/15 cell"
        ),
        "cell": {"uses": 3, "coherence": "4/5", "survival": "4/5", "theta": "2*pi/15"},
        "candidate_probabilities": laws["rational_O149_candidate"],
        "candidate_total_interval": candidate,
        "product_total_interval": totals["product"],
        "ghz_total_interval": totals["GHZ"],
        "certified_gain_lower": certified_gain,
        "pi_interval_width": PI_INTERVAL[1] - PI_INTERVAL[0],
        "proof_method": (
            "Machin and alternating Taylor rational intervals; exact partial traces; "
            "exact rational characteristic-polynomial root isolation for K^T K; "
            "nuclear-norm perturbation enclosure"
        ),
        "boundary": (
            "This proves feasibility and strict improvement over two named baseline families "
            "in one reduced three-use cell. It does not prove global simplex optimality, "
            "adaptive-comb optimality, apparatus implementability, or a twelve-mode result."
        ),
        "new_object": "O159 Rational Heralded-Code Gain Certificate",
    }, rows)


# ---------------------------------------------------------------------------
# P440: detector-envelope minimax Jordan sampler
# ---------------------------------------------------------------------------


def parse_fraction(text: str) -> Fraction:
    return Fraction(text)


def p429_atom_midpoints() -> tuple[np.ndarray, np.ndarray]:
    variables: dict[str, Fraction] = {}
    with P429_CSV.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["row_type"] != "certified_variable":
                continue
            name = row["variable"]
            if name.startswith("x") or name.startswith("w"):
                variables[name] = (parse_fraction(row["box_lower"]) + parse_fraction(row["box_upper"])) / 2
    nodes = [variables[f"x{index}"] for index in range(6)] + [Fraction(1)]
    weights = [variables[f"w{index}"] for index in range(7)]
    return np.array([float(value) for value in nodes]), np.array([float(value) for value in weights])


def program_440() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    nodes, weights = p429_atom_midpoints()
    features = np.array([[node**order for order in range(12)] for node in nodes])
    radii = np.linalg.norm(features, axis=1)
    moments = weights @ features

    # A declared local nuisance envelope, not a claim about any real detector.
    efficiency_low = 0.65 + 0.25 * nodes
    efficiency_high = np.minimum(0.99, efficiency_low + 0.05)
    dark_upper = 0.01 + 0.02 * (1 - nodes)
    worst_second_moment_factor = (
        1 / efficiency_low
        + dark_upper * (1 - dark_upper) / efficiency_low**2
    )

    baseline = np.abs(weights)
    baseline /= baseline.sum()
    p422 = np.abs(weights) * radii
    p422 /= p422.sum()
    minimax = np.abs(weights) * radii * np.sqrt(worst_second_moment_factor)
    minimax /= minimax.sum()

    def worst_mse_coefficient(probabilities: np.ndarray) -> float:
        second = np.sum(
            weights**2 * radii**2 * worst_second_moment_factor / probabilities
        )
        return float(second - np.sum(moments**2))

    coefficients = {
        "absolute_weight_baseline": worst_mse_coefficient(baseline),
        "P422_spectral_radius": worst_mse_coefficient(p422),
        "P440_detector_minimax": worst_mse_coefficient(minimax),
    }

    # Uniform finite-sample coordinatewise Hoeffding ledger.  The estimator
    # uses supplied calibrated epsilon_i,d_i; the design is robust over the
    # declared envelope.  B,D in {0,1} and X=B+D-d_i.
    alpha = 0.05
    coordinate_tolerance = 0.02
    coordinate_widths = []
    for order in range(12):
        outcomes = []
        for index in range(7):
            factor = (
                weights[index] * nodes[index] ** order
                / (minimax[index] * efficiency_low[index])
            )
            # The affine dependence on d means extrema over [0,d_upper]
            # occur at an endpoint.  Include both endpoints; using only
            # d_upper would underestimate the largest positive outcome.
            for dark_probability in (0.0, dark_upper[index]):
                outcomes.extend([
                    factor * (-dark_probability),
                    factor * (1 - dark_probability),
                    factor * (2 - dark_probability),
                ])
        coordinate_widths.append(max(outcomes) - min(outcomes))
    maximum_width = max(coordinate_widths)
    sufficient_samples = math.ceil(
        maximum_width**2 * math.log(24 / alpha) / (2 * coordinate_tolerance**2)
    )

    rows = []
    for index in range(7):
        rows.append({
            "atom": index,
            "node": nodes[index],
            "weight": weights[index],
            "feature_radius": radii[index],
            "efficiency_lower": efficiency_low[index],
            "efficiency_upper": efficiency_high[index],
            "dark_probability_upper": dark_upper[index],
            "worst_second_moment_factor": worst_second_moment_factor[index],
            "absolute_weight_probability": baseline[index],
            "P422_probability": p422[index],
            "P440_minimax_probability": minimax[index],
        })
    return ({
        "status": (
            "[Proven conditional theorem] detector-envelope minimax sampling law for a supplied "
            "calibrated Bernoulli-efficiency plus subtracted-dark-count model; "
            "[Computed design] P429 atom allocation and a uniform finite-sample ledger"
        ),
        "detector_model": (
            "given atom i, X_i=B_i+D_i-d_i with B_i~Bernoulli(epsilon_i), "
            "D_i~Bernoulli(d_i), independent; estimate w_i f_i X_i/(q_i epsilon_i)"
        ),
        "nuisance_envelope": (
            "epsilon_i in [0.65+0.25*x_i, min(0.99,0.70+0.25*x_i)], "
            "d_i in [0,0.01+0.02*(1-x_i)]"
        ),
        "minimax_law": (
            "q_i proportional to |w_i| ||f_i|| sqrt(1/epsilon_low_i + "
            "d_upper_i*(1-d_upper_i)/epsilon_low_i^2)"
        ),
        "proof": (
            "The worst second-moment factor is decreasing in epsilon and increasing in d on "
            "d<=1/2; Cauchy-Schwarz minimizes sum_i a_i/q_i on the simplex."
        ),
        "worst_mse_coefficients": coefficients,
        "reduction_vs_absolute_weight": 1 - coefficients["P440_detector_minimax"] / coefficients["absolute_weight_baseline"],
        "reduction_vs_P422": 1 - coefficients["P440_detector_minimax"] / coefficients["P422_spectral_radius"],
        "confidence_ledger": {
            "familywise_alpha": alpha,
            "coordinate_tolerance": coordinate_tolerance,
            "maximum_single_draw_coordinate_range": maximum_width,
            "hoeffding_sufficient_attempts": sufficient_samples,
            "coordinates": 12,
        },
        "boundary": (
            "The nuisance envelope is a transparent synthetic design assumption. The theorem "
            "does not supply detector calibration, dark-count stationarity, event independence, "
            "apparatus data, custody, or a physical FIN realization."
        ),
        "new_object": "O160 Detector-Envelope Minimax Jordan Sampler",
    }, rows)


def make_figure(p435_rows: list[dict[str, Any]], p436: dict[str, Any], p440: dict[str, Any], p440_rows: list[dict[str, Any]]) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.25))

    one = p435_rows[0]
    axes[0].bar(["primal", "dual"], [one["primal_success"], one["dual_success"]], color=["#3b82f6", "#f59e0b"])
    axes[0].set_ylim(0.64, 0.67)
    axes[0].set_title("P435 exact one-slot SDP")
    axes[0].set_ylabel("equal-prior success")
    axes[0].text(0.5, 0.642, "two-slot comb exported;\nmatching dual still open", ha="center", va="bottom", fontsize=9)

    labels = ["candidate", "product", "GHZ"]
    intervals = [
        p436["candidate_total_interval"],
        p436["product_total_interval"],
        p436["ghz_total_interval"],
    ]
    centers = [sum(float(Fraction(v)) for v in interval) / 2 for interval in intervals]
    errors = [(float(Fraction(interval[1])) - float(Fraction(interval[0]))) / 2 for interval in intervals]
    axes[1].bar(labels, centers, color=["#10b981", "#64748b", "#a855f7"])
    axes[1].errorbar(range(3), centers, yerr=errors, fmt="none", ecolor="black", capsize=3)
    axes[1].set_title("P436 certified heralded gain")
    axes[1].set_ylabel("half trace distance")

    atoms = [row["atom"] for row in p440_rows]
    axes[2].plot(atoms, [row["absolute_weight_probability"] for row in p440_rows], "o--", label="|weight|")
    axes[2].plot(atoms, [row["P422_probability"] for row in p440_rows], "s--", label="P422")
    axes[2].plot(atoms, [row["P440_minimax_probability"] for row in p440_rows], "^-", label="P440")
    axes[2].set_title("P440 detector-robust allocation")
    axes[2].set_xlabel("certified atom index")
    axes[2].set_ylabel("sampling probability")
    axes[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def main() -> None:
    p435, rows435 = program_435()
    p436, rows436 = program_436()
    p440, rows440 = program_440()
    payload = {
        "metadata": {
            "programs": ["P435", "P436", "P440"],
            "checkpoint": "P435-P440 bounded batch",
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
        },
        "P435": p435,
        "P436": p436,
        "P440": p440,
    }
    RESULTS_PATH.write_text(json.dumps(json_ready(payload), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    write_csv(P435_PATH, rows435)
    write_csv(P436_PATH, rows436)
    write_csv(P440_PATH, rows440)
    write_csv(SUMMARY_PATH, [
        {"program": "P435", "status": p435["status"], "new_object": p435["new_object"], "boundary": p435["boundary"]},
        {"program": "P436", "status": p436["status"], "new_object": p436["new_object"], "boundary": p436["boundary"]},
        {"program": "P440", "status": p440["status"], "new_object": p440["new_object"], "boundary": p440["boundary"]},
    ])
    make_figure(rows435, p436, p440, rows440)
    print(json.dumps({
        "p435_one_slot_gap": p435["one_slot_checks"]["primal_dual_gap"],
        "p435_two_slot_success_gap": p435["two_slot_success_gap"],
        "p436_certified_gain_lower": float(Fraction(p436["certified_gain_lower"])),
        "p440_reduction_vs_P422": p440["reduction_vs_P422"],
        "p440_hoeffding_attempts": p440["confidence_ledger"]["hoeffding_sufficient_attempts"],
    }, indent=2))


if __name__ == "__main__":
    main()
