#!/usr/bin/env python3
"""Scratch probe: distance-histogram proof of the Hebbian d5 margin.

The previous margin probe enumerated all 792 supports.  This probe compresses
that finite proof to the 35 possible pair-distance histograms for five active
nodes on Z_12 and replays the exact energy margin from those histograms.

Finite result:
- The d5 teacher orbit has unique top histogram h=(0,3,2,1,4,0).
- That histogram has exactly 12 supports, all the d5 translates.
- The closest non-d5 histogram is h=(1,3,2,1,3,0), with 24 supports.
- Binary/centered readouts have histogram gap 8; bipolar readouts have gap 32.
- Diagonal terms shift all five-active-support energies by a constant, so they
  do not change the histogram ordering or the margin.

No false pass: this is a finite conditional proof compression for supplied d5
teacher/readout premises, not a derivation of the teacher trace, not a Hebbian
law theorem, and not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_distance_histogram_margin_proof_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_distance_histogram_margin_proof_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
TARGET_STEP = 5
MINIMAL_REQUIRED_SUBGROUP = [1, 11]
D5_HISTOGRAM = (0, 3, 2, 1, 4, 0)
CLOSEST_NON_D5_HISTOGRAM = (1, 3, 2, 1, 3, 0)
VARIANTS = [
    ("binary_with_diagonal", "binary", False),
    ("binary_zero_self", "binary", True),
    ("centered_with_diagonal", "centered", False),
    ("centered_zero_self", "centered", True),
    ("bipolar_with_diagonal", "bipolar", False),
    ("bipolar_zero_self", "bipolar", True),
]

Support = tuple[int, ...]
Histogram = tuple[int, int, int, int, int, int]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def canonical_support(nodes: set[int] | tuple[int, ...]) -> Support:
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


def readout_vector(support: Support, vector_kind: str) -> list[Fraction]:
    active = set(support)
    if vector_kind == "binary":
        return [Fraction(1 if node in active else 0) for node in range(N)]
    if vector_kind == "centered":
        mean = Fraction(ACTIVE_COUNT, N)
        return [Fraction(1 if node in active else 0) - mean for node in range(N)]
    if vector_kind == "bipolar":
        return [Fraction(1 if node in active else -1) for node in range(N)]
    raise ValueError(f"unknown vector kind: {vector_kind}")


def first_row_weights_from_teacher(vector_kind: str, zero_self: bool) -> list[Fraction]:
    teacher = teacher_orbit(TARGET_STEP)
    weights = [Fraction(0) for _distance in range(N)]
    for support in teacher:
        vector = readout_vector(support, vector_kind)
        for distance in range(N):
            if zero_self and distance == 0:
                continue
            weights[distance] += vector[0] * vector[distance]
    return weights


def off_diagonal_histogram_energy(histogram: Histogram, first_row: list[Fraction]) -> Fraction:
    return 2 * sum(histogram[distance - 1] * first_row[distance] for distance in range(1, N // 2 + 1))


def diagonal_shift(first_row: list[Fraction], zero_self: bool) -> Fraction:
    if zero_self:
        return Fraction(0)
    return ACTIVE_COUNT * first_row[0]


def total_histogram_energy(histogram: Histogram, first_row: list[Fraction], zero_self: bool) -> Fraction:
    return diagonal_shift(first_row, zero_self) + off_diagonal_histogram_energy(histogram, first_row)


def unit_stabilizer_for_first_row(first_row: list[Fraction]) -> list[int]:
    def preserved_by(unit: int) -> bool:
        return all(first_row[(unit * distance) % N] == first_row[distance] for distance in range(N))

    return [unit for unit in AUT_UNITS if preserved_by(unit)]


def histogram_rows_for_variant(
    histogram_counter: Counter[Histogram],
    first_row: list[Fraction],
    zero_self: bool,
) -> list[dict[str, Any]]:
    rows = []
    for histogram, count in histogram_counter.items():
        rows.append(
            {
                "distance_histogram_d1_to_d6": list(histogram),
                "support_count": count,
                "off_diagonal_energy": fraction_text(off_diagonal_histogram_energy(histogram, first_row)),
                "total_energy": fraction_text(total_histogram_energy(histogram, first_row, zero_self)),
            }
        )
    return sorted(rows, key=lambda row: Fraction(row["total_energy"]), reverse=True)


def variant_certificate(
    variant_name: str,
    vector_kind: str,
    zero_self: bool,
    histogram_counter: Counter[Histogram],
) -> dict[str, Any]:
    first_row = first_row_weights_from_teacher(vector_kind, zero_self)
    rows = histogram_rows_for_variant(histogram_counter, first_row, zero_self)
    top = rows[0]
    second = rows[1]
    gap = Fraction(top["total_energy"]) - Fraction(second["total_energy"])
    d5_energy = total_histogram_energy(D5_HISTOGRAM, first_row, zero_self)
    closest_energy = total_histogram_energy(CLOSEST_NON_D5_HISTOGRAM, first_row, zero_self)
    return {
        "variant_name": variant_name,
        "vector_kind": vector_kind,
        "zero_self": zero_self,
        "unit_stabilizer": unit_stabilizer_for_first_row(first_row),
        "first_row_by_index_distance_0_to_11": [fraction_text(value) for value in first_row],
        "first_row_by_folded_distance_1_to_6": {
            str(distance): fraction_text(first_row[distance]) for distance in range(1, N // 2 + 1)
        },
        "diagonal_shift_for_five_active_supports": fraction_text(diagonal_shift(first_row, zero_self)),
        "top_histogram": top,
        "second_histogram": second,
        "d5_histogram_total_energy": fraction_text(d5_energy),
        "closest_non_d5_histogram_total_energy": fraction_text(closest_energy),
        "histogram_gap": fraction_text(gap),
        "top_histogram_is_d5": top["distance_histogram_d1_to_d6"] == list(D5_HISTOGRAM),
        "second_histogram_is_closest_non_d5": second["distance_histogram_d1_to_d6"] == list(CLOSEST_NON_D5_HISTOGRAM),
        "top_two_histogram_rows": rows[:2],
        "positive_histogram_margin_certificate": gap > 0
        and top["distance_histogram_d1_to_d6"] == list(D5_HISTOGRAM)
        and second["distance_histogram_d1_to_d6"] == list(CLOSEST_NON_D5_HISTOGRAM),
    }


def main() -> None:
    supports = all_supports()
    teacher = teacher_orbit(TARGET_STEP)
    teacher_set = set(teacher)
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    d5_supports_by_histogram = [support for support in supports if distance_histogram(support) == D5_HISTOGRAM]
    closest_supports_by_histogram = [support for support in supports if distance_histogram(support) == CLOSEST_NON_D5_HISTOGRAM]
    certificates = [
        variant_certificate(name, vector_kind, zero_self, histogram_counter)
        for name, vector_kind, zero_self in VARIANTS
    ]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_DISTANCE_HISTOGRAM_MARGIN_PROOF_PROBE__NOT_A_THEOREM",
        "status": "d5-energy-margin-reduced-to-distance-histogram-proof-not-origin-theorem",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "histogram_class_count": len(histogram_counter),
            "teacher_step": TARGET_STEP,
            "teacher_orbit_size": len(teacher),
            "automorphism_units": AUT_UNITS,
            "minimal_required_subgroup_from_previous_probe": MINIMAL_REQUIRED_SUBGROUP,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "histogram_proof_core": {
            "d5_histogram_d1_to_d6": list(D5_HISTOGRAM),
            "d5_histogram_support_count": len(d5_supports_by_histogram),
            "d5_histogram_supports_equal_teacher_orbit": set(d5_supports_by_histogram) == teacher_set,
            "closest_non_d5_histogram_d1_to_d6": list(CLOSEST_NON_D5_HISTOGRAM),
            "closest_non_d5_histogram_support_count": len(closest_supports_by_histogram),
            "closest_non_d5_examples": [list(support) for support in closest_supports_by_histogram[:8]],
            "support_enumeration_reduced_from": len(supports),
            "support_enumeration_reduced_to_histogram_classes": len(histogram_counter),
        },
        "variant_histogram_certificates": certificates,
        "histogram_margin_summary": {
            "tested_variant_count": len(certificates),
            "all_variants_have_positive_histogram_margin": all(row["positive_histogram_margin_certificate"] for row in certificates),
            "all_variants_have_required_stabilizer": all(row["unit_stabilizer"] == MINIMAL_REQUIRED_SUBGROUP for row in certificates),
            "gap_by_variant": {row["variant_name"]: row["histogram_gap"] for row in certificates},
            "diagonal_shift_by_variant": {row["variant_name"]: row["diagonal_shift_for_five_active_supports"] for row in certificates},
        },
        "candidate_source_interpretation": {
            "finite_gain": "The d5 energy-margin proof can be compressed from 792 support energies to 35 pair-distance histograms.",
            "conditional_positive_result": "For the supplied d5 teacher/readout premises, the unique top histogram is exactly the d5 histogram and the closest non-d5 histogram gives the previously observed gap.",
            "honest_limit": "The histogram proof does not derive the d5 teacher trace, the Hebbian update law, or the choice of distance-histogram representation from strict geometry.",
            "relation_to_previous_probe": "The previous margin probe found support-level gaps; this probe proves those gaps through distance-histogram classes and diagonal-shift replay.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to carry a finite pair-distance/Hebbian self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply the histogram representation.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the pair-distance histogram representation as a strict source.",
            "Histogram-margin replay across six tested readouts is not an exhaustive theorem over all learning laws.",
            "This is a conditional histogram proof compression, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the d5 teacher trace, Hebbian update, or pair-distance histogram representation internally; otherwise keep this as a conditional proof compression.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian distance-histogram margin proof probe\n\n"
        "Status: d5 margin reduced from support enumeration to distance-histogram proof; not an origin theorem.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- D5 histogram: `{D5_HISTOGRAM}` with support count `{len(d5_supports_by_histogram)}` and equals teacher orbit `{set(d5_supports_by_histogram) == teacher_set}`.\n"
        f"- Closest non-d5 histogram: `{CLOSEST_NON_D5_HISTOGRAM}` with support count `{len(closest_supports_by_histogram)}`.\n"
        f"- Gap by variant: `{report['histogram_margin_summary']['gap_by_variant']}`.\n"
        f"- Diagonal shift by variant: `{report['histogram_margin_summary']['diagonal_shift_by_variant']}`.\n"
        f"- All variants have positive histogram margin: `{report['histogram_margin_summary']['all_variants_have_positive_histogram_margin']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: histogram proof compresses the supplied Hebbian margin, but teacher trace, learning law, and histogram source remain supplied inputs.\n"
        "- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
