#!/usr/bin/env python3
"""FIN local research Programs P448--P450.

P448 constructs a concave fine-grained-erasure majorant and certifies a
global upper bound on the complete P446 three-use probability simplex.
P449 derives the exact recursive support cone of causal n-slot tester
normalizers and audits the first three-slot cone.
P450 tests representation independence of O163 with an exact finite-
difference null cycle.

All inputs and computations are local.  No physical record is consumed.
"""

from __future__ import annotations

import csv
from fractions import Fraction
import heapq
import itertools
import json
import math
from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

import fin_programs_435_436_440 as p435
import fin_programs_445_447 as p445


ROOT = Path(__file__).resolve().parent
PREFIX = "FIN_Programs_448_450"
RESULTS_PATH = ROOT / f"{PREFIX}_Results.json"
SUMMARY_PATH = ROOT / f"{PREFIX}_Summary.csv"
P448_PATH = ROOT / f"{PREFIX}_P448_Full_Simplex_Majorant.csv"
P449_PATH = ROOT / f"{PREFIX}_P449_Three_Slot_Recursion.csv"
P449_WITNESS_PATH = ROOT / f"{PREFIX}_P449_Three_Slot_Witness.npz"
P450_PATH = ROOT / f"{PREFIX}_P450_Null_Cycle.csv"
FIGURE_DIR = ROOT / f"{PREFIX}_Figures"
FIGURE_PATH = FIGURE_DIR / "p448_p450_global_majorant_and_gauge_obstruction.png"

Interval = tuple[Fraction, Fraction]


def json_ready(value: Any) -> Any:
    if isinstance(value, Fraction):
        if max(value.numerator.bit_length(), value.denominator.bit_length()) <= 512:
            return str(value)
        return format(float(value), ".17g")
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
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
            writer.writerow({
                key: (
                    json.dumps(json_ready(row.get(key)), ensure_ascii=False)
                    if isinstance(row.get(key), (list, tuple, dict, np.ndarray))
                    else json_ready(row.get(key, ""))
                )
                for key in keys
            })


# ---------------------------------------------------------------------------
# P448: concave fine-grained-erasure majorant
# ---------------------------------------------------------------------------


def fine_interval_skew(
    probabilities: list[Interval], survivors: int, lost_weight: int
) -> list[list[Interval]]:
    """One lost-computational-value contribution, before coarse graining."""

    dimension = 2**survivors
    result = [
        [(Fraction(0), Fraction(0)) for _ in range(dimension)]
        for _ in range(dimension)
    ]
    coherence = Fraction(4, 5)
    for left in range(dimension):
        for right in range(dimension):
            wl0 = left.bit_count()
            wr0 = right.bit_count()
            wl = wl0 + lost_weight
            wr = wr0 + lost_weight
            hamming = (left ^ right).bit_count()
            radicand = p445.iv_scale(
                p445.iv_mul(probabilities[wl], probabilities[wr]),
                Fraction(1, math.comb(3, wl) * math.comb(3, wr)),
            )
            amplitude = p445.iv_sqrt(radicand, 25)
            result[left][right] = p445.iv_scale(
                p445.iv_mul(amplitude, p435.p436_sine_interval(wl0 - wr0)),
                2 * coherence**hamming,
            )
    return result


def fine_component_interval(
    probabilities: list[Interval], survivors: int, lost_weight: int
) -> Interval:
    matrix = fine_interval_skew(probabilities, survivors, lost_weight)
    if survivors == 1:
        return p445.iv_abs(matrix[0][1])
    if survivors == 2:
        return p445.skew_three_distance(p445.interval_project(matrix, p445.U2))
    symmetric = p445.skew_four_distance(p445.interval_project(matrix, p445.U3))
    standard = p445.interval_project(matrix, p445.V3)
    return p445.iv_add(
        symmetric,
        p445.iv_scale(p445.iv_abs(standard[0][1]), 2),
    )


def p448_majorant_interval(probabilities: list[Interval]) -> Interval:
    terms: list[Interval] = []
    for survivors, mask_probability in (
        (1, Fraction(12, 125)),
        (2, Fraction(48, 125)),
        (3, Fraction(64, 125)),
    ):
        lost = 3 - survivors
        for lost_weight in range(lost + 1):
            multiplicity = math.comb(lost, lost_weight)
            terms.append(p445.iv_scale(
                fine_component_interval(probabilities, survivors, lost_weight),
                mask_probability * multiplicity,
            ))
    return p445.iv_sum(terms)


def fine_float_skew(
    probabilities: list[p445.FInterval], survivors: int, lost_weight: int
) -> list[list[p445.FInterval]]:
    """Outward-rounded binary64 version used by the finite cover."""

    dimension = 2**survivors
    result = [[(0.0, 0.0) for _ in range(dimension)] for _ in range(dimension)]
    for left in range(dimension):
        for right in range(dimension):
            wl0 = left.bit_count()
            wr0 = right.bit_count()
            wl = wl0 + lost_weight
            wr = wr0 + lost_weight
            hamming = (left ^ right).bit_count()
            denominator = math.comb(3, wl) * math.comb(3, wr)
            radicand = p445.fscale(
                p445.fmul(probabilities[wl], probabilities[wr]),
                (1 / denominator, 1 / denominator),
            )
            amplitude = p445.fsqrt(radicand)
            q_power = (
                p445.Q_FI[0] ** hamming,
                p445.Q_FI[1] ** hamming,
            )
            result[left][right] = p445.fscale(
                p445.fmul(
                    p445.fmul(amplitude, p445.SINES_FI[wl0 - wr0]),
                    q_power,
                ),
                (2.0, 2.0),
            )
    return result


def fine_float_component(
    probabilities: list[p445.FInterval], survivors: int, lost_weight: int
) -> p445.FInterval:
    matrix = fine_float_skew(probabilities, survivors, lost_weight)
    if survivors == 1:
        return p445.fabsiv(matrix[0][1])
    if survivors == 2:
        block = p445.fproject(matrix, p445.U2_FI)
        return p445.fsqrt(p445.fsum([
            p445.fsquare(block[i][j])
            for i in range(3) for j in range(i + 1, 3)
        ]))
    block = p445.fproject(matrix, p445.U3_FI)
    squared = p445.fsum([
        p445.fsquare(block[i][j])
        for i in range(4) for j in range(i + 1, 4)
    ])
    pfaffian = p445.fadd(
        p445.fsub(
            p445.fmul(block[0][1], block[2][3]),
            p445.fmul(block[0][2], block[1][3]),
        ),
        p445.fmul(block[0][3], block[1][2]),
    )
    symmetric = p445.fsqrt(p445.fadd(
        squared, p445.fscale(p445.fabsiv(pfaffian), (2.0, 2.0))
    ))
    standard = p445.fproject(matrix, p445.V3_FI)
    return p445.fadd(
        symmetric,
        p445.fscale(p445.fabsiv(standard[0][1]), (2.0, 2.0)),
    )


def p448_majorant_float_interval(
    probabilities: list[p445.FInterval],
) -> p445.FInterval:
    terms: list[p445.FInterval] = []
    for survivors, mask_probability in (
        (1, 12 / 125), (2, 48 / 125), (3, 64 / 125)
    ):
        lost = 3 - survivors
        for lost_weight in range(lost + 1):
            coefficient = mask_probability * math.comb(lost, lost_weight)
            terms.append(p445.fscale(
                fine_float_component(probabilities, survivors, lost_weight),
                (coefficient, coefficient),
            ))
    return p445.fsum(terms)


def fine_component_numeric(
    probabilities: np.ndarray, survivors: int, lost_weight: int
) -> float:
    dimension = 2**survivors
    matrix = np.zeros((dimension, dimension), dtype=complex)
    q = 0.8
    theta = 2 * math.pi / 15
    for left in range(dimension):
        for right in range(dimension):
            wl0 = left.bit_count()
            wr0 = right.bit_count()
            wl = wl0 + lost_weight
            wr = wr0 + lost_weight
            amplitude = math.sqrt(
                max(0.0, probabilities[wl] * probabilities[wr])
                / (math.comb(3, wl) * math.comb(3, wr))
            )
            matrix[left, right] = (
                2j * q ** ((left ^ right).bit_count())
                * math.sin(theta * (wl0 - wr0)) * amplitude
            )
    return float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(matrix))))


def p448_majorant_numeric(probabilities: np.ndarray) -> float:
    total = 0.0
    for survivors, mask_probability in ((1, 12 / 125), (2, 48 / 125), (3, 64 / 125)):
        lost = 3 - survivors
        for lost_weight in range(lost + 1):
            total += (
                mask_probability * math.comb(lost, lost_weight)
                * fine_component_numeric(probabilities, survivors, lost_weight)
            )
    return total


def palindromic_probabilities(interval_a: Interval) -> list[Interval]:
    interval_b = (
        Fraction(1, 2) - interval_a[1],
        Fraction(1, 2) - interval_a[0],
    )
    return [interval_a, interval_b, interval_b, interval_a]


def p448_branch_certificate(tolerance: float = 1e-3) -> dict[str, Any]:
    """Outward-rounded interval cover of the palindromic majorant."""

    def probabilities(interval_a: p445.FInterval) -> list[p445.FInterval]:
        interval_b = (
            p445.fdown(0.5 - interval_a[1]),
            p445.fup(0.5 - interval_a[0]),
        )
        return [interval_a, interval_b, interval_b, interval_a]

    initial = (0.0, 0.5)
    initial_value = p448_majorant_float_interval(probabilities(initial))
    heap: list[tuple[float, int, p445.FInterval, p445.FInterval]] = []
    counter = 0
    heapq.heappush(heap, (-initial_value[1], counter, initial, initial_value))

    seed = 0.22
    seed_value = p448_majorant_float_interval(probabilities((seed, seed)))
    incumbent_lower = seed_value[0]
    incumbent_a = seed
    pruned_upper = incumbent_lower
    evaluated = 0
    terminal_boxes = 0

    while heap and -heap[0][0] > incumbent_lower + tolerance:
        _, _, interval_a, value_interval = heapq.heappop(heap)
        midpoint = (interval_a[0] + interval_a[1]) / 2
        midpoint_value = p448_majorant_float_interval(probabilities((midpoint, midpoint)))
        evaluated += 1
        if midpoint_value[0] > incumbent_lower:
            incumbent_lower = midpoint_value[0]
            incumbent_a = midpoint

        for child in (
            (interval_a[0], midpoint),
            (midpoint, interval_a[1]),
        ):
            child_value = p448_majorant_float_interval(probabilities(child))
            counter += 1
            if child_value[1] <= incumbent_lower:
                pruned_upper = max(pruned_upper, child_value[1])
                terminal_boxes += 1
            else:
                heapq.heappush(heap, (-child_value[1], counter, child, child_value))
        if evaluated > 50000:
            raise RuntimeError("P448 majorant certificate exceeded finite budget")

    live_upper = -heap[0][0] if heap else incumbent_lower
    global_upper = max(pruned_upper, live_upper)
    live_boxes = [(box, value) for _, _, box, value in heap]
    relevant = [
        (box, value) for box, value in live_boxes
        if value[1] >= incumbent_lower
    ]
    if relevant:
        hull = (
            min(box[0] for box, _ in relevant),
            max(box[1] for box, _ in relevant),
        )
    else:
        hull = (incumbent_a, incumbent_a)
    return {
        "incumbent_a": incumbent_a,
        "incumbent_lower": incumbent_lower,
        "global_upper": global_upper,
        "majorant_gap": global_upper - incumbent_lower,
        "maximizer_hull": hull,
        "evaluated_boxes": evaluated,
        "live_boxes": len(heap),
        "terminal_boxes": terminal_boxes,
    }


def program_448() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    certificate = p448_branch_certificate()
    original_lower = 0.46327828319203235
    full_simplex_gap = certificate["global_upper"] - original_lower

    rng = np.random.default_rng(448)
    maximum_monotonicity_defect = -math.inf
    maximum_symmetrization_defect = -math.inf
    for probabilities in rng.dirichlet(np.full(4, 0.35), size=256):
        original = p445.p446_numeric_objective(probabilities)
        majorant = p448_majorant_numeric(probabilities)
        symmetric = (probabilities + probabilities[::-1]) / 2
        symmetric_majorant = p448_majorant_numeric(symmetric)
        maximum_monotonicity_defect = max(
            maximum_monotonicity_defect, original - majorant
        )
        maximum_symmetrization_defect = max(
            maximum_symmetrization_defect, majorant - symmetric_majorant
        )

    rows = [{
        "incumbent_a": certificate["incumbent_a"],
        "majorant_lower": certificate["incumbent_lower"],
        "majorant_global_upper": certificate["global_upper"],
        "majorant_certificate_gap": certificate["majorant_gap"],
        "original_certified_lower": original_lower,
        "original_full_simplex_upper_gap": full_simplex_gap,
        "maximizer_hull_lower": certificate["maximizer_hull"][0],
        "maximizer_hull_upper": certificate["maximizer_hull"][1],
        "evaluated_boxes": certificate["evaluated_boxes"],
        "live_boxes": certificate["live_boxes"],
        "terminal_boxes": certificate["terminal_boxes"],
    }]

    return ({
        "status": (
            "[Computer-assisted proof] a global full-simplex upper bound via a "
            "concave fine-grained-erasure majorant; [Open] exact full-simplex globality"
        ),
        "majorization_theorem": (
            "Coarse graining over the unobserved lost computational values is CPTP, so "
            "the original heralded trace distance is no larger than the fine-grained "
            "distance. Each fine block has form Tr|P^(1/2) H P^(1/2)|/2 and is concave "
            "in diagonal P by the nuclear-norm variational formula. Reversal invariance "
            "then places the majorant maximum on p=(a,1/2-a,1/2-a,a)."
        ),
        "certificate": certificate,
        "original_certified_lower": original_lower,
        "original_full_simplex_global_upper": certificate["global_upper"],
        "original_full_simplex_gap": full_simplex_gap,
        "maximum_random_monotonicity_defect": maximum_monotonicity_defect,
        "maximum_random_symmetrization_defect": maximum_symmetrization_defect,
        "proof_boundary": (
            "The majorant reveals the lost computational value and can be strictly larger "
            "than the operational coarse-grained objective. Therefore the certificate does "
            "not prove that the P436 witness is the exact full-simplex maximizer."
        ),
        "new_object": "O164 Concave Fine-Grained-Erasure Majorant",
    }, rows)


# ---------------------------------------------------------------------------
# P449: exact recursive support cone for causal testers
# ---------------------------------------------------------------------------


def support_affine_dimension(slots: int) -> int:
    if slots < 1:
        raise ValueError("at least one slot required")
    dimension = 1
    for level in range(2, slots + 1):
        dimension += (2 ** (level - 1)) ** 2
    return dimension


def central_support_normalizer(slots: int) -> np.ndarray:
    return np.eye(2**slots, dtype=complex) / 2**slots


def recursive_three_slot_sample(seed: int, radius: float = 0.20) -> np.ndarray:
    """A feasible N3=B0 direct-sum B1 with B0+B1=I4/4."""

    rng = np.random.default_rng(seed)
    raw = rng.normal(size=(4, 4)) + 1j * rng.normal(size=(4, 4))
    hermitian = (raw + raw.conj().T) / 2
    hermitian /= np.linalg.norm(hermitian, 2)
    difference = radius * hermitian
    inherited = central_support_normalizer(2)
    block_zero = (inherited + difference) / 2
    block_one = (inherited - difference) / 2
    result = np.zeros((8, 8), dtype=complex)
    result[:4, :4] = block_zero
    result[4:, 4:] = block_one
    return result


def compressed_process_difference(slots: int, q: float, theta: float) -> np.ndarray:
    plus = p435.qubit_phase_dephasing_choi(q, theta, +1)
    minus = p435.qubit_phase_dephasing_choi(q, theta, -1)
    process_plus = plus
    process_minus = minus
    for _ in range(slots - 1):
        process_plus = np.kron(process_plus, plus)
        process_minus = np.kron(process_minus, minus)
    support = []
    for history in itertools.product((0, 3), repeat=slots):
        index = 0
        for digit in history:
            index = 4 * index + digit
        support.append(index)
    delta = process_plus - process_minus
    return delta[np.ix_(support, support)]


def tester_distance(normalizer: np.ndarray, delta: np.ndarray) -> float:
    eigenvalues, eigenvectors = np.linalg.eigh(normalizer)
    root = (
        eigenvectors * np.sqrt(np.maximum(eigenvalues, 0))
    ) @ eigenvectors.conj().T
    return float(0.5 * np.sum(np.abs(np.linalg.eigvalsh(root @ delta @ root))))


def p449_sine_interval(weight_difference: int) -> p435.Interval:
    """Enclose sin(pi*weight_difference/8), |weight_difference|<=3."""

    if weight_difference == 0:
        return Fraction(0), Fraction(0)
    if weight_difference < 0:
        positive = p449_sine_interval(-weight_difference)
        return -positive[1], -positive[0]
    low = p435.PI_INTERVAL[0] * Fraction(weight_difference, 8)
    high = p435.PI_INTERVAL[1] * Fraction(weight_difference, 8)
    return p435.sine_point_interval(low)[0], p435.sine_point_interval(high)[1]


def p449_echo_witness_interval() -> tuple[Interval, Interval]:
    """Certified half-distance intervals for the echo witness and GHZ."""

    probabilities = [
        Fraction(value, 40) for value in (9, 1, 9, 1, 1, 9, 1, 9)
    ]
    q = Fraction(4, 5)
    matrix: list[list[Interval]] = []
    for left in range(8):
        row: list[Interval] = []
        for right in range(8):
            amplitude = p435.rational_sqrt_interval(
                probabilities[left] * probabilities[right], 28
            )
            hamming = (left ^ right).bit_count()
            difference = left.bit_count() - right.bit_count()
            row.append(p435.iv_scale(
                p435.iv_mul(amplitude, p449_sine_interval(difference)),
                2 * q**hamming,
            ))
        matrix.append(row)
    witness_low, witness_high, _ = p435.exact_midpoint_trace_distance_bounds(
        matrix, 18
    )
    ghz = p435.iv_scale(p449_sine_interval(3), q**3)
    return (witness_low, witness_high), ghz


def p449_full_echo_extension() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Explicit full Xi3 and tester sum extending the support witness."""

    inherited = np.eye(4) / 4
    difference = np.diag([1, -1, 1, -1]) / 5
    block_zero = (inherited + difference) / 2
    block_one = (inherited - difference) / 2
    support_normalizer = np.zeros((8, 8), dtype=complex)
    support_normalizer[:4, :4] = block_zero
    support_normalizer[4:, 4:] = block_one

    # Xi3 acts on X3,Y2,X2,Y1,X1.  Its X3 blocks must sum to
    # I_Y2 tensor Xi2=I_16/4.  Embed B0 on the equality support
    # Y2=X2 and Y1=X1; the complement is filled by I/8.
    equality_support = [
        np.ravel_multi_index((b, b, c, c), (2, 2, 2, 2))
        for b, c in itertools.product((0, 1), repeat=2)
    ]
    full_zero = np.eye(16, dtype=complex) / 8
    for index, support_index in enumerate(equality_support):
        full_zero[support_index, support_index] = block_zero[index, index]
    full_one = np.eye(16, dtype=complex) / 4 - full_zero
    xi3 = np.zeros((32, 32), dtype=complex)
    xi3[:16, :16] = full_zero
    xi3[16:, 16:] = full_one
    tester_sum = np.kron(np.eye(2), xi3)
    return support_normalizer, xi3, tester_sum


def program_449() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    q = 0.8
    theta = math.pi / 8
    delta = compressed_process_difference(3, q, theta)
    ghz = np.zeros((8, 8), dtype=complex)
    ghz[0, 0] = ghz[-1, -1] = 0.5
    ghz_distance = tester_distance(ghz, delta)

    echo, xi3, tester_sum = p449_full_echo_extension()
    echo_distance = tester_distance(echo, delta)
    echo_interval, ghz_interval = p449_echo_witness_interval()
    certified_advantage = echo_interval[0] - ghz_interval[1]

    # Independent full-space causality and compression checks.
    xi3_trace = p435.partial_trace_square(xi3, (2, 2, 2, 2, 2), 0)
    xi2 = np.eye(8) / 4
    xi2_trace = p435.partial_trace_square(xi2, (2, 2, 2), 0)
    full_support = []
    for a, b, c in itertools.product((0, 1), repeat=3):
        full_support.append(np.ravel_multi_index(
            (a, a, b, b, c, c), (2, 2, 2, 2, 2, 2)
        ))
    compressed_tester_sum = tester_sum[np.ix_(full_support, full_support)]
    full_causality_residual = float(np.linalg.norm(xi3_trace - np.eye(16) / 4))
    inherited_causality_residual = float(np.linalg.norm(xi2_trace - np.eye(4) / 2))
    compression_residual = float(np.linalg.norm(compressed_tester_sum - echo))
    full_minimum_eigenvalue = float(np.linalg.eigvalsh(xi3)[0])
    np.savez_compressed(
        P449_WITNESS_PATH,
        support_normalizer=echo,
        Xi3=xi3,
        tester_sum=tester_sum,
        process_difference_support=delta,
        history_probabilities=np.array([9, 1, 9, 1, 1, 9, 1, 9]) / 40,
        q=np.array([q]),
        theta=np.array([theta]),
    )

    rows: list[dict[str, Any]] = []
    maximum_random_distance = 0.0
    maximum_recursion_residual = 0.0
    minimum_eigenvalue = math.inf
    maximum_offdiagonal_norm = 0.0
    for seed in range(32):
        normalizer = recursive_three_slot_sample(449000 + seed)
        block_zero = normalizer[:4, :4]
        block_one = normalizer[4:, 4:]
        inherited = block_zero + block_one
        recursion_residual = float(np.linalg.norm(
            inherited - central_support_normalizer(2)
        ))
        eigenvalue = float(np.linalg.eigvalsh(normalizer)[0])
        distance = tester_distance(normalizer, delta)
        offdiagonal = float(np.linalg.norm(
            normalizer - np.diag(np.diag(normalizer))
        ))
        maximum_random_distance = max(maximum_random_distance, distance)
        maximum_recursion_residual = max(maximum_recursion_residual, recursion_residual)
        minimum_eigenvalue = min(minimum_eigenvalue, eigenvalue)
        maximum_offdiagonal_norm = max(maximum_offdiagonal_norm, offdiagonal)
        rows.append({
            "seed": seed,
            "minimum_eigenvalue": eigenvalue,
            "recursion_residual": recursion_residual,
            "offdiagonal_frobenius_norm": offdiagonal,
            "three_slot_half_distance": distance,
            "beats_ghz": distance > ghz_distance + 1e-12,
        })

    dimensions = {
        str(slots): support_affine_dimension(slots) for slots in range(1, 6)
    }
    return ({
        "status": (
            "[Proven] exact recursive characterization of the compressed causal tester "
            "normalizer cone; [Computer-assisted proof] an explicit rational echo-history "
            "tester beats three-slot GHZ; [Open] the global three-slot optimum"
        ),
        "recursion": (
            "N_n=B_0 direct-sum B_1, B_0>=0, B_1>=0, and "
            "B_0+B_1 belongs to C_(n-1); C_1={diag(r,1-r):0<=r<=1}."
        ),
        "affine_dimension_recurrence": "d_1=1; d_n=d_(n-1)+4^(n-1)",
        "affine_dimensions_n1_to_n5": dimensions,
        "three_slot_affine_dimension": support_affine_dimension(3),
        "three_slot_diagonal_simplex_dimension": 7,
        "three_slot_coherence_excess_dimension": support_affine_dimension(3) - 7,
        "ghz_half_distance": ghz_distance,
        "ghz_half_distance_interval": ghz_interval,
        "echo_history_law": ["9/40", "1/40", "9/40", "1/40", "1/40", "9/40", "1/40", "9/40"],
        "echo_rule": "weight 9/40 when the first and third history bits agree, 1/40 otherwise; the middle bit is free",
        "echo_half_distance": echo_distance,
        "echo_half_distance_interval": echo_interval,
        "certified_echo_advantage_over_ghz": certified_advantage,
        "echo_success_probability_lower": (Fraction(1) + echo_interval[0]) / 2,
        "ghz_success_probability_upper": (Fraction(1) + ghz_interval[1]) / 2,
        "full_xi3_causality_residual": full_causality_residual,
        "inherited_xi2_causality_residual": inherited_causality_residual,
        "full_support_compression_residual": compression_residual,
        "full_xi3_minimum_eigenvalue": full_minimum_eigenvalue,
        "maximum_random_feasible_distance": maximum_random_distance,
        "maximum_recursion_residual": maximum_recursion_residual,
        "minimum_random_normalizer_eigenvalue": minimum_eigenvalue,
        "maximum_random_offdiagonal_norm": maximum_offdiagonal_norm,
        "optimization_boundary": (
            "The exact cone has 21 affine dimensions, fourteen more than the diagonal "
            "history simplex. A diagonal correlated echo-history law already refutes GHZ "
            "optimality. The P445 one-coherence square completion does not survive as a "
            "scalar identity; no matching three-slot dual or global optimum is claimed."
        ),
        "new_object": "O165 Recursive Causal-Support Cone",
    }, rows)


# ---------------------------------------------------------------------------
# P450: representation dependence and an exact null cycle
# ---------------------------------------------------------------------------


def p429_midpoint_fractions() -> tuple[list[Fraction], list[Fraction]]:
    values: dict[str, Fraction] = {}
    with p445.P429_CSV.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["row_type"] != "certified_variable":
                continue
            name = row["variable"]
            if name.startswith("x") or name.startswith("w"):
                values[name] = (
                    Fraction(row["box_lower"]) + Fraction(row["box_upper"])
                ) / 2
    nodes = [values[f"x{index}"] for index in range(6)] + [Fraction(1)]
    weights = [values[f"w{index}"] for index in range(7)]
    return nodes, weights


def null_cycle() -> tuple[list[Fraction], list[Fraction]]:
    nodes = [Fraction(index, 12) for index in range(13)]
    weights = [Fraction((-1) ** index * math.comb(12, index)) for index in range(13)]
    return nodes, weights


def aggregate_measure(
    nodes: list[Fraction], weights: list[Fraction]
) -> tuple[list[Fraction], list[Fraction]]:
    aggregate: dict[Fraction, Fraction] = {}
    for node, weight in zip(nodes, weights):
        aggregate[node] = aggregate.get(node, Fraction(0)) + weight
    pairs = sorted((node, weight) for node, weight in aggregate.items() if weight)
    return [node for node, _ in pairs], [weight for _, weight in pairs]


def detector_score(node: float, weight: float) -> float:
    radius = math.sqrt(sum(node ** (2 * order) for order in range(12)))
    efficiency = 0.65 + 0.25 * node
    dark = 0.03 - 0.02 * node
    hbar = 1 / efficiency + dark * (1 - dark) / efficiency**2
    return abs(weight) * radius * math.sqrt(hbar)


def measure_risk(nodes: list[Fraction], weights: list[Fraction]) -> tuple[float, float]:
    scores = [detector_score(float(x), float(w)) for x, w in zip(nodes, weights)]
    moments = [
        sum(float(w) * float(x) ** order for x, w in zip(nodes, weights))
        for order in range(12)
    ]
    score_sum = sum(scores)
    return score_sum**2 - sum(moment * moment for moment in moments), score_sum


def program_450() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    base_nodes, base_weights = p429_midpoint_fractions()
    cycle_nodes, cycle_weights = null_cycle()
    exact_residuals = [
        sum(weight * node**order for node, weight in zip(cycle_nodes, cycle_weights))
        for order in range(12)
    ]
    if any(exact_residuals):
        raise AssertionError("finite-difference null cycle failed exact moment test")

    base_risk, base_score = measure_risk(base_nodes, base_weights)
    rows: list[dict[str, Any]] = []
    risks = []
    for epsilon in (
        Fraction(0), Fraction(1, 10**8), Fraction(1, 10**6),
        Fraction(1, 10**4), Fraction(1, 10**3), Fraction(1, 100),
    ):
        nodes, weights = aggregate_measure(
            base_nodes + cycle_nodes,
            base_weights + [epsilon * weight for weight in cycle_weights],
        )
        risk, score = measure_risk(nodes, weights)
        moment_defects = [
            sum(weight * node**order for node, weight in zip(nodes, weights))
            - sum(weight * node**order for node, weight in zip(base_nodes, base_weights))
            for order in range(12)
        ]
        maximum_exact_defect = max(abs(value) for value in moment_defects)
        risks.append(risk)
        rows.append({
            "epsilon": epsilon,
            "atom_count_after_aggregation": len(nodes),
            "maximum_exact_moment_defect": maximum_exact_defect,
            "allocation_score_sum": score,
            "minimax_risk_coefficient": risk,
            "risk_ratio_to_reduced_representation": risk / base_risk,
        })

    positive_split_invariance = True
    split_index = 3
    split_fraction = Fraction(2, 5)
    split_nodes = list(base_nodes)
    split_weights = list(base_weights)
    node = split_nodes.pop(split_index)
    weight = split_weights.pop(split_index)
    split_nodes.extend([node, node])
    split_weights.extend([split_fraction * weight, (1 - split_fraction) * weight])
    split_risk, split_score = measure_risk(split_nodes, split_weights)
    positive_split_invariance &= abs(split_risk - base_risk) < 1e-12
    positive_split_invariance &= abs(split_score - base_score) < 1e-12

    return ({
        "status": (
            "[Proven] O163 is invariant under same-sign splitting at an identical node but "
            "is not invariant under exact moment-null-cycle augmentation"
        ),
        "null_cycle": (
            "nu_12=sum_(j=0)^12 (-1)^j binom(12,j) delta_(j/12); "
            "its moments of orders 0 through 11 vanish exactly."
        ),
        "maximum_exact_null_moment_residual": max(abs(value) for value in exact_residuals),
        "same_sign_identical_node_split_invariant": positive_split_invariance,
        "base_risk": base_risk,
        "largest_tested_risk": max(risks),
        "largest_tested_risk_ratio": max(risks) / base_risk,
        "unboundedness_theorem": (
            "For mu_epsilon=mu+epsilon*nu_12 the twelve target moments are fixed, while "
            "the O163 score sum grows at least linearly in |epsilon| because nu_12 has "
            "nonzero support. Hence the declared minimax risk is unbounded along one exact "
            "representation gauge orbit."
        ),
        "necessity": (
            "A preparation-level representation rule is necessary. Minimal support, reduced "
            "Jordan form, minimum total variation, or an explicit instrument must be supplied; "
            "the moment functional alone cannot determine the sampler."
        ),
        "physical_boundary": (
            "The null cycle is a mathematical counterexample, not an apparatus process. "
            "Selecting a representation remains an operational premise."
        ),
        "new_object": "O166 Moment-Null Representation Gauge",
    }, rows)


def make_figure(
    p448: dict[str, Any], p449_rows: list[dict[str, Any]], p450_rows: list[dict[str, Any]]
) -> None:
    FIGURE_DIR.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.25))

    grid = np.linspace(0, 0.5, 180)
    original = [
        p445.p446_numeric_objective(np.array([a, 0.5 - a, 0.5 - a, a]))
        for a in grid
    ]
    majorant = [
        p448_majorant_numeric(np.array([a, 0.5 - a, 0.5 - a, a]))
        for a in grid
    ]
    axes[0].plot(grid, original, label="coarse objective", color="#2563eb")
    axes[0].plot(grid, majorant, label="fine majorant", color="#dc2626")
    axes[0].axhline(
        float(Fraction(p448["original_full_simplex_global_upper"])),
        color="#dc2626", linestyle="--", linewidth=0.9,
    )
    axes[0].set_title("P448 full-simplex majorant")
    axes[0].set_xlabel("palindromic a")
    axes[0].set_ylabel("half distance")
    axes[0].legend(fontsize=8)

    axes[1].hist(
        [row["three_slot_half_distance"] for row in p449_rows],
        bins=10, color="#0f766e", alpha=0.8,
    )
    axes[1].axvline(p449_rows[0].get("ghz_half_distance", 0), alpha=0)
    axes[1].set_title("P449 feasible 21D cone audit")
    axes[1].set_xlabel("three-slot half distance")
    axes[1].set_ylabel("count")

    epsilons = [float(Fraction(row["epsilon"])) for row in p450_rows]
    ratios = [row["risk_ratio_to_reduced_representation"] for row in p450_rows]
    axes[2].semilogx(epsilons[1:], ratios[1:], "o-", color="#7c3aed")
    axes[2].axhline(1, color="#64748b", linestyle="--", linewidth=0.9)
    axes[2].set_title("P450 exact null-cycle gauge")
    axes[2].set_xlabel("epsilon")
    axes[2].set_ylabel("risk ratio")

    fig.tight_layout()
    fig.savefig(FIGURE_PATH, dpi=220)
    plt.close(fig)


def main() -> None:
    p448_result, p448_rows = program_448()
    p449_result, p449_rows = program_449()
    for row in p449_rows:
        row["ghz_half_distance"] = p449_result["ghz_half_distance"]
    p450_result, p450_rows = program_450()

    payload = {
        "metadata": {
            "programs": ["P448", "P449", "P450"],
            "date": "2026-08-01",
            "local_only": True,
            "external_physical_evidence": False,
            "selector_discharged": False,
            "dimensional_source_exported": False,
            "legacy_strict_bridge_complete": False,
        },
        "P448": p448_result,
        "P449": p449_result,
        "P450": p450_result,
    }
    RESULTS_PATH.write_text(
        json.dumps(json_ready(payload), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    write_csv(P448_PATH, p448_rows)
    write_csv(P449_PATH, p449_rows)
    write_csv(P450_PATH, p450_rows)
    write_csv(SUMMARY_PATH, [
        {"program": "P448", "status": p448_result["status"], "new_object": p448_result["new_object"], "boundary": p448_result["proof_boundary"]},
        {"program": "P449", "status": p449_result["status"], "new_object": p449_result["new_object"], "boundary": p449_result["optimization_boundary"]},
        {"program": "P450", "status": p450_result["status"], "new_object": p450_result["new_object"], "boundary": p450_result["physical_boundary"]},
    ])
    make_figure(p448_result, p449_rows, p450_rows)
    print(json.dumps({
        "p448_full_simplex_upper": float(Fraction(p448_result["original_full_simplex_global_upper"])),
        "p448_full_simplex_gap": float(Fraction(p448_result["original_full_simplex_gap"])),
        "p449_dimension": p449_result["three_slot_affine_dimension"],
        "p449_ghz": p449_result["ghz_half_distance"],
        "p450_risk_ratio": p450_result["largest_tested_risk_ratio"],
    }, indent=2))


if __name__ == "__main__":
    main()
