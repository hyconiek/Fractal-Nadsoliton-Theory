#!/usr/bin/env python3
"""Scratch probe: exact fifth-mode symbolic histogram selector proof for d5.

The distance-5 pair-count probe isolated d5 from the single h_5 channel.  This
probe reconnects that core to the earlier "highest resonance" language by
computing the fifth Fourier-mode coherent power P_5 exactly over the 35 feasible
five-node pair-distance histograms on Z_12.

Finite result:
- For a five-node support with histogram h, P_5(h) is exactly
  5 + 2*sum_d h_d*cos(5*pi*d/6), an element of Q(sqrt(3)).
- The d5 histogram h=(0,3,2,1,4,0) is the unique maximizer with
  P_5 = 7 + 4*sqrt(3).
- The closest non-d5 histogram h=(1,3,2,1,3,0) has
  P_5 = 7 + 2*sqrt(3), so the exact symbolic gap is 2*sqrt(3).
- This is a proof-style replay of fifth-mode selection without floating-point
  spectral numerics.

No false pass: the fifth-mode channel is supplied as a selector/readout premise;
this does not derive why strict geometry must choose k=5, and it does not
discharge QW-2191.
"""
from __future__ import annotations

import json
from collections import Counter
from functools import cmp_to_key
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_fifth_mode_symbolic_histogram_selector_proof_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_fifth_mode_symbolic_histogram_selector_proof_report.md"

N = 12
ACTIVE_COUNT = 5
TARGET_MODE = 5
TARGET_STEP = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
MINIMAL_REQUIRED_SUBGROUP = [1, 11]
D5_HISTOGRAM = (0, 3, 2, 1, 4, 0)
CLOSEST_NON_D5_HISTOGRAM = (1, 3, 2, 1, 3, 0)
# cos(5*pi*d/6) for folded distances d=1..6 as rational + sqrt(3)*coefficient.
FIFTH_MODE_COS_COEFFICIENTS_D1_TO_D6 = (
    (Fraction(0), Fraction(-1, 2)),
    (Fraction(1, 2), Fraction(0)),
    (Fraction(0), Fraction(0)),
    (Fraction(-1, 2), Fraction(0)),
    (Fraction(0), Fraction(1, 2)),
    (Fraction(-1), Fraction(0)),
)

Support = tuple[int, ...]
Histogram = tuple[int, int, int, int, int, int]
Quadratic = tuple[Fraction, Fraction]  # a + b*sqrt(3)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def quadratic_text(value: Quadratic) -> str:
    rational, sqrt3_coeff = value
    if sqrt3_coeff == 0:
        return fraction_text(rational)
    if rational == 0:
        if sqrt3_coeff == 1:
            return "sqrt3"
        if sqrt3_coeff == -1:
            return "-sqrt3"
        return f"{fraction_text(sqrt3_coeff)}*sqrt3"
    sign = "+" if sqrt3_coeff > 0 else "-"
    coeff_abs = abs(sqrt3_coeff)
    coeff_text = "sqrt3" if coeff_abs == 1 else f"{fraction_text(coeff_abs)}*sqrt3"
    return f"{fraction_text(rational)} {sign} {coeff_text}"


def quadratic_json(value: Quadratic) -> dict[str, str]:
    return {
        "rational_part": fraction_text(value[0]),
        "sqrt3_coefficient": fraction_text(value[1]),
        "expression": quadratic_text(value),
    }


def sign_quadratic(value: Quadratic) -> int:
    """Return sign of a + b*sqrt(3) exactly."""
    rational, sqrt3_coeff = value
    if sqrt3_coeff == 0:
        return (rational > 0) - (rational < 0)
    if rational == 0:
        return (sqrt3_coeff > 0) - (sqrt3_coeff < 0)
    if rational > 0 and sqrt3_coeff > 0:
        return 1
    if rational < 0 and sqrt3_coeff < 0:
        return -1
    # Mixed signs: compare |rational| with |sqrt3_coeff|*sqrt(3).
    square_difference = 3 * sqrt3_coeff * sqrt3_coeff - rational * rational
    if sqrt3_coeff > 0 and rational < 0:
        return (square_difference > 0) - (square_difference < 0)
    # sqrt3_coeff < 0 and rational > 0
    return (square_difference < 0) - (square_difference > 0)


def compare_quadratic(left: Quadratic, right: Quadratic) -> int:
    return sign_quadratic((left[0] - right[0], left[1] - right[1]))


def canonical_support(nodes: set[int] | tuple[int, ...] | list[int]) -> Support:
    return tuple(sorted(node % N for node in nodes))


def all_supports() -> list[Support]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def teacher_support(step: int, translate: int) -> Support:
    return canonical_support({translate + index * step for index in range(ACTIVE_COUNT)})


def teacher_orbit(step: int) -> list[Support]:
    return sorted({teacher_support(step, translate) for translate in range(N)})


def folded_distance(left: int, right: int) -> int:
    distance = (right - left) % N
    return min(distance, N - distance)


def distance_histogram(support: Support) -> Histogram:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded_distance(left, right) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def fifth_mode_power_from_histogram(histogram: Histogram) -> Quadratic:
    rational = Fraction(ACTIVE_COUNT)
    sqrt3_coeff = Fraction(0)
    for count, (cos_rational, cos_sqrt3_coeff) in zip(histogram, FIFTH_MODE_COS_COEFFICIENTS_D1_TO_D6):
        rational += 2 * count * cos_rational
        sqrt3_coeff += 2 * count * cos_sqrt3_coeff
    return rational, sqrt3_coeff


def folded_mode(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def fifth_mode_channel_stabilizer() -> list[int]:
    return [unit for unit in AUT_UNITS if folded_mode(unit * TARGET_MODE) == TARGET_MODE]


def histogram_power_rows(histogram_counter: Counter[Histogram]) -> list[dict[str, Any]]:
    rows = []
    for histogram, support_count in histogram_counter.items():
        power = fifth_mode_power_from_histogram(histogram)
        rows.append(
            {
                "distance_histogram_d1_to_d6": list(histogram),
                "support_count": support_count,
                "p5_power_exact": quadratic_json(power),
            }
        )
    return sorted(
        rows,
        key=cmp_to_key(
            lambda left, right: -compare_quadratic(
                (Fraction(left["p5_power_exact"]["rational_part"]), Fraction(left["p5_power_exact"]["sqrt3_coefficient"])),
                (Fraction(right["p5_power_exact"]["rational_part"]), Fraction(right["p5_power_exact"]["sqrt3_coefficient"])),
            )
        ),
    )


def main() -> None:
    supports = all_supports()
    teacher = teacher_orbit(TARGET_STEP)
    teacher_set = set(teacher)
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    rows = histogram_power_rows(histogram_counter)
    top = rows[0]
    second = rows[1]
    top_power = fifth_mode_power_from_histogram(D5_HISTOGRAM)
    second_power = fifth_mode_power_from_histogram(CLOSEST_NON_D5_HISTOGRAM)
    gap = (top_power[0] - second_power[0], top_power[1] - second_power[1])
    d5_supports = [support for support in supports if distance_histogram(support) == D5_HISTOGRAM]
    closest_supports = [support for support in supports if distance_histogram(support) == CLOSEST_NON_D5_HISTOGRAM]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_FIFTH_MODE_SYMBOLIC_HISTOGRAM_SELECTOR_PROOF_PROBE__NOT_A_THEOREM",
        "status": "fifth-mode-symbolically-selects-d5-conditionally-not-origin-theorem",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "target_mode": TARGET_MODE,
            "teacher_step": TARGET_STEP,
            "teacher_orbit_size": len(teacher),
            "automorphism_units": AUT_UNITS,
            "minimal_required_subgroup_from_previous_probe": MINIMAL_REQUIRED_SUBGROUP,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "symbolic_power_formula": {
            "formula": "P_5(h)=5+2*sum_d h_d*cos(5*pi*d/6)",
            "coherent_power_field": "Q(sqrt(3))",
            "cos_coefficients_d1_to_d6_as_rational_plus_sqrt3_coeff": [
                {"rational_part": fraction_text(r), "sqrt3_coefficient": fraction_text(s)}
                for r, s in FIFTH_MODE_COS_COEFFICIENTS_D1_TO_D6
            ],
            "comparison_method": "Exact sign comparison for a+b*sqrt(3), no floating-point ranking.",
            "fifth_mode_channel_stabilizer": fifth_mode_channel_stabilizer(),
        },
        "selector_certificate": {
            "top_histogram": top,
            "second_histogram": second,
            "top_histogram_is_d5": top["distance_histogram_d1_to_d6"] == list(D5_HISTOGRAM),
            "second_histogram_is_closest_non_d5": second["distance_histogram_d1_to_d6"] == list(CLOSEST_NON_D5_HISTOGRAM),
            "d5_histogram_support_count": len(d5_supports),
            "d5_histogram_supports_equal_teacher_orbit": set(d5_supports) == teacher_set,
            "closest_non_d5_support_count": len(closest_supports),
            "closest_non_d5_examples": [list(support) for support in closest_supports[:8]],
            "symbolic_gap_top_minus_second": quadratic_json(gap),
            "positive_symbolic_gap": sign_quadratic(gap) > 0,
        },
        "relation_to_previous_core": {
            "distance5_pair_count_core": "h_5 uniquely selects d5 as a supplied objective.",
            "fifth_mode_readout": "P_5 also selects d5 exactly, but includes additional h_1,h_2,h_4,h_6 terms with Q(sqrt(3)) coefficients.",
            "finite_gain": "This bridges the h_5 core selector back to a symbolic fifth-resonance readout without numerical Fourier approximations.",
        },
        "candidate_source_interpretation": {
            "finite_gain": "The fifth-mode coherent-power selector can be proven exactly over histogram classes in Q(sqrt(3)).",
            "conditional_positive_result": "If a fifth-mode resonance readout is supplied, the unique maximizer is the d5 teacher orbit.",
            "honest_limit": "This probe does not derive why strict geometry must choose k=5 rather than another channel.",
            "relation_to_previous_probe": "The h_5 core probe supplied a minimal pair-count selector; this probe supplies a symbolic fifth-mode resonance replay.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to carry a finite fifth-mode self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose k=5.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives fifth-mode channel selection as a strict source.",
            "The symbolic P_5 selector is conditional on supplying k=5 as the relevant resonance readout.",
            "This is a conditional symbolic histogram selector proof, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the k=5 readout internally; otherwise keep fifth-mode selection as a conditional symbolic resonance certificate.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian fifth-mode symbolic histogram selector proof probe\n\n"
        "Status: fifth-mode power selects d5 symbolically; not an origin theorem.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- P5 formula: `P_5(h)=5+2*sum_d h_d*cos(5*pi*d/6)` in `Q(sqrt(3))`.\n"
        f"- Top histogram: `{top}`.\n"
        f"- Second histogram: `{second}`.\n"
        f"- Symbolic gap: `{quadratic_text(gap)}`; positive `{sign_quadratic(gap) > 0}`.\n"
        f"- Fifth-mode channel stabilizer: `{fifth_mode_channel_stabilizer()}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: exact fifth-mode readout selects d5, but k=5 remains a supplied selector channel.\n"
        "- No false pass: no strict fifth-mode source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
