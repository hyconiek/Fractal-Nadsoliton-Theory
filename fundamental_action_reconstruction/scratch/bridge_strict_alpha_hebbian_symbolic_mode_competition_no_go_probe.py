#!/usr/bin/env python3
"""Scratch probe: exact symbolic mode-competition no-go for d5.

The fifth-mode symbolic probe proved that a supplied k=5 readout selects the d5
histogram exactly in Q(sqrt(3)).  This probe checks the stronger and more risky
claim: does unconstrained highest coherent Fourier power over non-DC modes select
d5 without specifying k=5?

Finite result:
- Each non-DC folded mode k=1..6 has an exact top histogram certificate.
- k=5 selects d5 with P_5=7+4*sqrt(3) and gap 2*sqrt(3).
- k=1 selects the contiguous/unit-related mirror class with the same exact power.
- Unconstrained mode competition is won by the Nyquist/parity channel k=6 with
  P_6=25, strictly above the fifth-mode d5 power.
- Therefore "highest resonance" alone is not enough; the fifth-mode channel lock
  remains an explicit selector premise.

No false pass: this is an exact symbolic obstruction to an unconstrained
highest-resonance theorem, not a strict-core d5 source theorem and not a
QW-2191 discharge.
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
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_symbolic_mode_competition_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_symbolic_mode_competition_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
D5_HISTOGRAM = (0, 3, 2, 1, 4, 0)
CONTIGUOUS_HISTOGRAM = (4, 3, 2, 1, 0, 0)
NYQUIST_HISTOGRAM = (0, 4, 0, 4, 0, 2)
# cos(k*pi*d/6) for folded distance d=1..6 represented as rational + sqrt(3)*coefficient.
COS_COEFFICIENTS_BY_MODE: dict[int, tuple[tuple[Fraction, Fraction], ...]] = {
    1: ((0, Fraction(1, 2)), (Fraction(1, 2), 0), (0, 0), (Fraction(-1, 2), 0), (0, Fraction(-1, 2)), (-1, 0)),
    2: ((Fraction(1, 2), 0), (Fraction(-1, 2), 0), (-1, 0), (Fraction(-1, 2), 0), (Fraction(1, 2), 0), (1, 0)),
    3: ((0, 0), (-1, 0), (0, 0), (1, 0), (0, 0), (-1, 0)),
    4: ((Fraction(-1, 2), 0), (Fraction(-1, 2), 0), (1, 0), (Fraction(-1, 2), 0), (Fraction(-1, 2), 0), (1, 0)),
    5: ((0, Fraction(-1, 2)), (Fraction(1, 2), 0), (0, 0), (Fraction(-1, 2), 0), (0, Fraction(1, 2)), (-1, 0)),
    6: ((-1, 0), (1, 0), (-1, 0), (1, 0), (-1, 0), (1, 0)),
}

Support = tuple[int, ...]
Histogram = tuple[int, int, int, int, int, int]
Quadratic = tuple[Fraction, Fraction]


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
    rational, sqrt3_coeff = value
    if sqrt3_coeff == 0:
        return (rational > 0) - (rational < 0)
    if rational == 0:
        return (sqrt3_coeff > 0) - (sqrt3_coeff < 0)
    if rational > 0 and sqrt3_coeff > 0:
        return 1
    if rational < 0 and sqrt3_coeff < 0:
        return -1
    square_difference = 3 * sqrt3_coeff * sqrt3_coeff - rational * rational
    if sqrt3_coeff > 0 and rational < 0:
        return (square_difference > 0) - (square_difference < 0)
    return (square_difference < 0) - (square_difference > 0)


def compare_quadratic(left: Quadratic, right: Quadratic) -> int:
    return sign_quadratic((left[0] - right[0], left[1] - right[1]))


def canonical_support(nodes: set[int] | tuple[int, ...] | list[int]) -> Support:
    return tuple(sorted(node % N for node in nodes))


def all_supports() -> list[Support]:
    return [tuple(support) for support in combinations(range(N), ACTIVE_COUNT)]


def folded_distance(left: int, right: int) -> int:
    distance = (right - left) % N
    return min(distance, N - distance)


def distance_histogram(support: Support) -> Histogram:
    counts = [0] * (N // 2)
    for left, right in combinations(support, 2):
        counts[folded_distance(left, right) - 1] += 1
    return tuple(counts)  # type: ignore[return-value]


def folded_mode(value: int) -> int:
    residue = value % N
    return min(residue, N - residue)


def mode_channel_stabilizer(mode: int) -> list[int]:
    return [unit for unit in AUT_UNITS if folded_mode(unit * mode) == mode]


def mode_power_from_histogram(mode: int, histogram: Histogram) -> Quadratic:
    rational = Fraction(ACTIVE_COUNT)
    sqrt3_coeff = Fraction(0)
    for count, (cos_rational, cos_sqrt3_coeff) in zip(histogram, COS_COEFFICIENTS_BY_MODE[mode]):
        rational += 2 * count * Fraction(cos_rational)
        sqrt3_coeff += 2 * count * Fraction(cos_sqrt3_coeff)
    return rational, sqrt3_coeff


def histogram_rows_for_mode(mode: int, histogram_counter: Counter[Histogram]) -> list[dict[str, Any]]:
    rows = []
    for histogram, support_count in histogram_counter.items():
        rows.append(
            {
                "distance_histogram_d1_to_d6": list(histogram),
                "support_count": support_count,
                "power_exact": quadratic_json(mode_power_from_histogram(mode, histogram)),
            }
        )
    return sorted(
        rows,
        key=cmp_to_key(
            lambda left, right: -compare_quadratic(
                (Fraction(left["power_exact"]["rational_part"]), Fraction(left["power_exact"]["sqrt3_coefficient"])),
                (Fraction(right["power_exact"]["rational_part"]), Fraction(right["power_exact"]["sqrt3_coefficient"])),
            )
        ),
    )


def mode_certificate(mode: int, histogram_counter: Counter[Histogram]) -> dict[str, Any]:
    rows = histogram_rows_for_mode(mode, histogram_counter)
    top = rows[0]
    second = rows[1]
    top_value = mode_power_from_histogram(mode, tuple(top["distance_histogram_d1_to_d6"]))
    second_value = mode_power_from_histogram(mode, tuple(second["distance_histogram_d1_to_d6"]))
    gap = (top_value[0] - second_value[0], top_value[1] - second_value[1])
    return {
        "mode": mode,
        "channel_stabilizer": mode_channel_stabilizer(mode),
        "top_histogram": top,
        "second_histogram": second,
        "top_minus_second_gap_exact": quadratic_json(gap),
        "positive_gap": sign_quadratic(gap) > 0,
        "top_is_d5": top["distance_histogram_d1_to_d6"] == list(D5_HISTOGRAM),
        "top_is_contiguous": top["distance_histogram_d1_to_d6"] == list(CONTIGUOUS_HISTOGRAM),
        "top_is_nyquist": top["distance_histogram_d1_to_d6"] == list(NYQUIST_HISTOGRAM),
    }


def main() -> None:
    supports = all_supports()
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    mode_certificates = {str(mode): mode_certificate(mode, histogram_counter) for mode in range(1, 7)}
    top_competitors = [
        (
            mode_power_from_histogram(mode, tuple(certificate["top_histogram"]["distance_histogram_d1_to_d6"])),
            mode,
            certificate,
        )
        for mode, certificate in ((int(key), value) for key, value in mode_certificates.items())
    ]
    top_competitors_sorted = sorted(top_competitors, key=cmp_to_key(lambda left, right: -compare_quadratic(left[0], right[0])))
    global_winner_value, global_winner_mode, global_winner_certificate = top_competitors_sorted[0]
    fifth_value = mode_power_from_histogram(5, D5_HISTOGRAM)
    global_minus_fifth = (global_winner_value[0] - fifth_value[0], global_winner_value[1] - fifth_value[1])

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_SYMBOLIC_MODE_COMPETITION_NO_GO_PROBE__NOT_A_THEOREM",
        "status": "unconstrained-highest-resonance-selects-nyquist-not-d5",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "modes_scanned": [1, 2, 3, 4, 5, 6],
            "automorphism_units": AUT_UNITS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "mode_certificates": mode_certificates,
        "mode_competition_certificate": {
            "global_winner_mode": global_winner_mode,
            "global_winner_power_exact": quadratic_json(global_winner_value),
            "global_winner_histogram": global_winner_certificate["top_histogram"],
            "global_winner_is_nyquist": global_winner_mode == 6 and global_winner_certificate["top_is_nyquist"],
            "fifth_mode_d5_power_exact": quadratic_json(fifth_value),
            "global_winner_minus_fifth_mode_d5_exact": quadratic_json(global_minus_fifth),
            "global_winner_strictly_exceeds_fifth_mode_d5": sign_quadratic(global_minus_fifth) > 0,
            "fifth_mode_selects_d5_conditionally": mode_certificates["5"]["top_is_d5"],
            "unconstrained_highest_resonance_selects_d5": global_winner_mode == 5 and global_winner_certificate["top_is_d5"],
        },
        "candidate_source_interpretation": {
            "finite_gain": "Exact symbolic mode competition separates the conditional fifth-mode d5 selector from the unconstrained highest-resonance claim.",
            "conditional_positive_result": "Mode k=5 selects d5 exactly if k=5 is supplied as the readout channel.",
            "negative_result": "If the readout is unconstrained highest coherent power over modes k=1..6, k=6/Nyquist wins instead of d5.",
            "honest_limit": "This probe does not derive an internal source for choosing k=5 over k=6.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to carry a finite fifth-mode channel self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose k=5.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives fifth-mode channel selection as a strict source.",
            "Unconstrained highest-resonance selection is obstructed by the k=6 Nyquist winner.",
            "This is a symbolic no-go for unconstrained mode competition, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive an internal fifth-mode channel lock; otherwise keep k=5 as a conditional selector premise because unconstrained mode competition selects k=6.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian symbolic mode-competition no-go probe\n\n"
        "Status: unconstrained highest resonance selects Nyquist k=6, not d5.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- k=5 top: `{mode_certificates['5']['top_histogram']}`.\n"
        f"- k=6 top: `{mode_certificates['6']['top_histogram']}`.\n"
        f"- Global winner mode: `{global_winner_mode}` with power `{quadratic_text(global_winner_value)}`.\n"
        f"- Fifth-mode d5 power: `{quadratic_text(fifth_value)}`; global-minus-fifth `{quadratic_text(global_minus_fifth)}`.\n"
        f"- Unconstrained highest resonance selects d5: `{report['mode_competition_certificate']['unconstrained_highest_resonance_selects_d5']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: k=5 selects d5 only conditionally; unconstrained mode competition selects k=6/Nyquist.\n"
        "- No false pass: no strict fifth-mode source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
