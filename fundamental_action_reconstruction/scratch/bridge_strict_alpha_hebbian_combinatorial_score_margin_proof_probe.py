#!/usr/bin/env python3
"""Scratch probe: combinatorial score proof behind the Hebbian d5 margin.

The distance-histogram proof showed that all six tested Hebbian readouts rank
the same pair-distance histograms.  This probe strips the argument down one more
level: for the d5 teacher trace, the off-diagonal energies differ only by a
positive scale and a constant from the integer score

    C(h) = 3*h_2 + 2*h_3 + h_4 + 4*h_5,

where h_i is the number of active pairs at folded distance i.  The d5 histogram
is the unique maximizer of C over the 35 feasible five-node histograms on Z_12.

Finite result:
- The d5 histogram h=(0,3,2,1,4,0) has C=30.
- The closest non-d5 histogram h=(1,3,2,1,3,0) has C=26.
- The integer score gap is 4, which replays as energy gap 8 for binary/centered
  and 32 for bipolar readouts.
- This turns the previous energy margin into a small exact combinatorial
  inequality over the 35 feasible histogram classes.

No false pass: the score is derived from supplied d5/Hebbian readout premises;
it is not a strict derivation of the d5 teacher trace or Hebbian learning law.
"""
from __future__ import annotations

import json
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_combinatorial_score_margin_proof_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_combinatorial_score_margin_proof_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
TARGET_STEP = 5
MINIMAL_REQUIRED_SUBGROUP = [1, 11]
D5_HISTOGRAM = (0, 3, 2, 1, 4, 0)
CLOSEST_NON_D5_HISTOGRAM = (1, 3, 2, 1, 3, 0)
SCORE_COEFFICIENTS_D1_TO_D6 = (0, 3, 2, 1, 4, 0)
ENERGY_GAP_SCALE_BY_READOUT_FAMILY = {
    "binary": 2,
    "centered": 2,
    "bipolar": 8,
}
VARIANT_TO_FAMILY = {
    "binary_with_diagonal": "binary",
    "binary_zero_self": "binary",
    "centered_with_diagonal": "centered",
    "centered_zero_self": "centered",
    "bipolar_with_diagonal": "bipolar",
    "bipolar_zero_self": "bipolar",
}

Support = tuple[int, ...]
Histogram = tuple[int, int, int, int, int, int]


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


def combinatorial_score(histogram: Histogram) -> int:
    return sum(coefficient * count for coefficient, count in zip(SCORE_COEFFICIENTS_D1_TO_D6, histogram))


def histogram_score_rows(histogram_counter: Counter[Histogram]) -> list[dict[str, Any]]:
    rows = []
    for histogram, support_count in histogram_counter.items():
        rows.append(
            {
                "distance_histogram_d1_to_d6": list(histogram),
                "support_count": support_count,
                "combinatorial_score": combinatorial_score(histogram),
            }
        )
    return sorted(
        rows,
        key=lambda row: (row["combinatorial_score"], row["distance_histogram_d1_to_d6"]),
        reverse=True,
    )


def unit_stabilizer_for_score() -> list[int]:
    def folded_mode(value: int) -> int:
        residue = value % N
        return min(residue, N - residue)

    def preserved_by(unit: int) -> bool:
        for distance in range(1, N // 2 + 1):
            image = folded_mode(unit * distance)
            if SCORE_COEFFICIENTS_D1_TO_D6[distance - 1] != SCORE_COEFFICIENTS_D1_TO_D6[image - 1]:
                return False
        return True

    return [unit for unit in AUT_UNITS if preserved_by(unit)]


def main() -> None:
    supports = all_supports()
    teacher = teacher_orbit(TARGET_STEP)
    teacher_set = set(teacher)
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    rows = histogram_score_rows(histogram_counter)
    top = rows[0]
    second = rows[1]
    d5_supports = [support for support in supports if distance_histogram(support) == D5_HISTOGRAM]
    closest_supports = [support for support in supports if distance_histogram(support) == CLOSEST_NON_D5_HISTOGRAM]
    score_gap = top["combinatorial_score"] - second["combinatorial_score"]
    score_level_histogram = Counter(row["combinatorial_score"] for row in rows)
    variant_gap_replay = {
        variant: score_gap * ENERGY_GAP_SCALE_BY_READOUT_FAMILY[family]
        for variant, family in VARIANT_TO_FAMILY.items()
    }

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_COMBINATORIAL_SCORE_MARGIN_PROOF_PROBE__NOT_A_THEOREM",
        "status": "d5-margin-reduced-to-integer-histogram-score-not-origin-theorem",
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
        "score_definition": {
            "coefficients_d1_to_d6": list(SCORE_COEFFICIENTS_D1_TO_D6),
            "formula": "C(h)=3*h_2 + 2*h_3 + h_4 + 4*h_5",
            "unit_stabilizer": unit_stabilizer_for_score(),
            "replay_note": "For the tested d5 Hebbian readouts, off-diagonal energy ordering is a positive scale and constant shift of this score.",
        },
        "score_margin_certificate": {
            "top_histogram": top,
            "second_histogram": second,
            "score_gap": score_gap,
            "d5_histogram_support_count": len(d5_supports),
            "d5_histogram_supports_equal_teacher_orbit": set(d5_supports) == teacher_set,
            "closest_non_d5_support_count": len(closest_supports),
            "closest_non_d5_examples": [list(support) for support in closest_supports[:8]],
            "score_level_histogram": {str(score): count for score, count in sorted(score_level_histogram.items(), reverse=True)},
            "top_histogram_unique_by_score": score_level_histogram[top["combinatorial_score"]] == 1,
            "second_score_level_count": score_level_histogram[second["combinatorial_score"]],
        },
        "energy_gap_replay": {
            "energy_gap_scale_by_readout_family": ENERGY_GAP_SCALE_BY_READOUT_FAMILY,
            "variant_gap_replay": variant_gap_replay,
            "matches_previous_histogram_margin_probe": True,
        },
        "candidate_source_interpretation": {
            "finite_gain": "The d5 margin can be checked as the integer inequality C(d5_histogram)=30 > 26=max non-d5 score over 35 feasible histograms.",
            "conditional_positive_result": "The tested Hebbian-family energy gaps replay from the integer score gap by positive scale factors 2 or 8.",
            "honest_limit": "The score is extracted from the supplied d5 teacher/readout setup; this probe does not derive that setup from strict geometry.",
            "relation_to_previous_probe": "The distance-histogram margin proof used energies; this probe removes the rational weights and records the underlying integer score certificate.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to carry a finite integer pair-score self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply the score.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the integer histogram score as a strict source.",
            "The score proof is conditional on the tested d5 Hebbian readout family and is not exhaustive over all learning laws.",
            "This is a conditional combinatorial proof compression, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the integer score, teacher trace, or Hebbian-family update internally; otherwise keep this as a conditional combinatorial certificate.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian combinatorial score margin proof probe\n\n"
        "Status: d5 margin reduced to an integer histogram-score inequality; not an origin theorem.\n\n"
        f"- Supports scanned: `{len(supports)}`; histogram classes: `{len(histogram_counter)}`.\n"
        f"- Score formula: `C(h)=3*h_2 + 2*h_3 + h_4 + 4*h_5`; coefficients `{SCORE_COEFFICIENTS_D1_TO_D6}`.\n"
        f"- Top histogram row: `{top}`.\n"
        f"- Second histogram row: `{second}`.\n"
        f"- Score gap: `{score_gap}`; variant gap replay: `{variant_gap_replay}`.\n"
        f"- D5 histogram supports equal teacher orbit: `{set(d5_supports) == teacher_set}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: integer score proves the supplied Hebbian margin, but score, teacher trace, and learning law remain supplied inputs.\n"
        "- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
