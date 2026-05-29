#!/usr/bin/env python3
"""Scratch probe: Hebbian rule-variant robustness audit for the d5 teacher trace.

The previous self-record audit used one centered zero-self Hebbian convention.
This probe checks whether that was a fragile convention artifact by replaying the
same finite d5 teacher trace under six exact Hebbian readouts:

- binary outer product, with and without diagonal,
- centered/covariance outer product, with and without diagonal,
- bipolar outer product, with and without diagonal.

Finite result:
- All six readouts have unit stabilizer {1,11} for the learned d5 weight record.
- All six readouts make the 12 d5 teacher translates the unique global maxima
  among all 792 five-active-node supports on Z_12.
- Therefore the previous d5 weight-stabilizer/energy maximum certificate is not
  an artifact of the centered zero-self convention alone.
- But the teacher trace and the choice to use a Hebbian-family update remain
  supplied inputs; this is robustness evidence, not a strict origin theorem.

No false pass: no Hebbian law is derived as strict-core dynamics, no strict d5
source theorem is claimed, and QW-2191 remains undischarged.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_rule_variant_robustness_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_rule_variant_robustness_audit_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
TARGET_STEP = 5
MINIMAL_REQUIRED_SUBGROUP = [1, 11]
VARIANTS = [
    ("binary_with_diagonal", "binary", False),
    ("binary_zero_self", "binary", True),
    ("centered_with_diagonal", "centered", False),
    ("centered_zero_self", "centered", True),
    ("bipolar_with_diagonal", "bipolar", False),
    ("bipolar_zero_self", "bipolar", True),
]

Support = tuple[int, ...]
Matrix = list[list[Fraction]]


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


def hebbian_weights(teacher: list[Support], vector_kind: str, zero_self: bool) -> Matrix:
    weights = [[Fraction(0) for _ in range(N)] for _ in range(N)]
    for support in teacher:
        vector = readout_vector(support, vector_kind)
        for left in range(N):
            for right in range(N):
                if zero_self and left == right:
                    continue
                weights[left][right] += vector[left] * vector[right]
    return weights


def learned_energy(support: Support, weights: Matrix) -> Fraction:
    return sum(weights[left][right] for left in support for right in support)


def matrix_preserved_by_unit(weights: Matrix, unit: int) -> bool:
    return all(
        weights[(unit * left) % N][(unit * right) % N] == weights[left][right]
        for left in range(N)
        for right in range(N)
    )


def unit_stabilizer(weights: Matrix) -> list[int]:
    return [unit for unit in AUT_UNITS if matrix_preserved_by_unit(weights, unit)]


def first_row_text(weights: Matrix) -> list[str]:
    return [fraction_text(weights[0][distance]) for distance in range(N)]


def folded_first_row_text(weights: Matrix) -> dict[str, str]:
    return {str(distance): fraction_text(weights[0][distance]) for distance in range(1, N // 2 + 1)}


def energy_landscape_certificate(
    supports: list[Support],
    teacher: list[Support],
    weights: Matrix,
) -> dict[str, Any]:
    energy_by_support = {support: learned_energy(support, weights) for support in supports}
    maximum = max(energy_by_support.values())
    maximizers = sorted(support for support, energy in energy_by_support.items() if energy == maximum)
    energy_histogram = Counter(energy_by_support.values())
    teacher_set = set(teacher)
    return {
        "maximum_energy": fraction_text(maximum),
        "maximizer_count": len(maximizers),
        "maximizers_equal_d5_teacher_orbit": set(maximizers) == teacher_set,
        "non_d5_global_maximizer_count": sum(1 for support in maximizers if support not in teacher_set),
        "energy_level_count": len(energy_histogram),
        "top_energy_levels": [
            {"energy": fraction_text(energy), "support_count": count}
            for energy, count in sorted(energy_histogram.items(), reverse=True)[:5]
        ],
    }


def variant_certificate(
    variant_name: str,
    vector_kind: str,
    zero_self: bool,
    supports: list[Support],
    teacher: list[Support],
) -> dict[str, Any]:
    weights = hebbian_weights(teacher, vector_kind, zero_self)
    stabilizer = unit_stabilizer(weights)
    energy_certificate = energy_landscape_certificate(supports, teacher, weights)
    return {
        "variant_name": variant_name,
        "vector_kind": vector_kind,
        "zero_self": zero_self,
        "diagonal_zero_check": all(weights[node][node] == 0 for node in range(N)),
        "unit_stabilizer": stabilizer,
        "stabilizer_equals_required_subgroup": stabilizer == MINIMAL_REQUIRED_SUBGROUP,
        "first_row_by_index_distance_0_to_11": first_row_text(weights),
        "first_row_by_folded_distance_1_to_6": folded_first_row_text(weights),
        "energy_landscape": energy_certificate,
        "passes_d5_robustness_gate": stabilizer == MINIMAL_REQUIRED_SUBGROUP
        and energy_certificate["maximizers_equal_d5_teacher_orbit"]
        and energy_certificate["non_d5_global_maximizer_count"] == 0,
    }


def main() -> None:
    supports = all_supports()
    teacher = teacher_orbit(TARGET_STEP)
    certificates = [
        variant_certificate(name, vector_kind, zero_self, supports, teacher)
        for name, vector_kind, zero_self in VARIANTS
    ]
    passing_variants = [row["variant_name"] for row in certificates if row["passes_d5_robustness_gate"]]
    failing_variants = [row["variant_name"] for row in certificates if not row["passes_d5_robustness_gate"]]
    maximum_energy_by_variant = {
        row["variant_name"]: row["energy_landscape"]["maximum_energy"]
        for row in certificates
    }

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_RULE_VARIANT_ROBUSTNESS_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "d5-hebbian-stabilizer-and-energy-maxima-robust-across-tested-readouts-not-origin-theorem",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "teacher_step": TARGET_STEP,
            "teacher_orbit_size": len(teacher),
            "automorphism_units": AUT_UNITS,
            "minimal_required_subgroup_from_previous_probe": MINIMAL_REQUIRED_SUBGROUP,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "variant_certificates": certificates,
        "robustness_summary": {
            "tested_variant_count": len(certificates),
            "passing_variant_count": len(passing_variants),
            "failing_variant_count": len(failing_variants),
            "passing_variants": passing_variants,
            "failing_variants": failing_variants,
            "all_tested_variants_have_required_stabilizer": all(row["stabilizer_equals_required_subgroup"] for row in certificates),
            "all_tested_variants_have_unique_d5_global_maxima": all(
                row["energy_landscape"]["maximizers_equal_d5_teacher_orbit"]
                and row["energy_landscape"]["non_d5_global_maximizer_count"] == 0
                for row in certificates
            ),
            "maximum_energy_by_variant": maximum_energy_by_variant,
        },
        "candidate_source_interpretation": {
            "finite_gain": "The d5 stabilizer/energy maximum certificate survives six exact Hebbian-family readouts, so the previous result is not merely a centered zero-self convention artifact.",
            "conditional_positive_result": "Given the d5 teacher/self-record trace and a Hebbian-family associative update, the learned record robustly carries the required subgroup and d5 maxima.",
            "honest_limit": "The scan does not derive the teacher trace, does not derive the Hebbian update law, and does not prove that these six readouts exhaust all plausible learning dynamics.",
            "relation_to_previous_probe": "The previous probe showed one learned d5 weight matrix has stabilizer {1,11}; this probe checks robustness under nearby Hebbian readout conventions.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may carry a finite Hebbian-family self-record pattern.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose the teacher trace or learning law.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "Robustness across six tested Hebbian-family readouts is not an exhaustive theorem over all learning laws.",
            "This is a conditional robustness audit, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the d5 teacher/self-record trace or Hebbian-family update internally; otherwise keep this as robustness evidence for a supplied learning premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian rule-variant robustness audit probe\n\n"
        "Status: d5 Hebbian stabilizer and energy maxima are robust across tested readouts; not an origin theorem.\n\n"
        f"- Tested variants: `{len(certificates)}`.\n"
        f"- Passing variants: `{passing_variants}`.\n"
        f"- Failing variants: `{failing_variants}`.\n"
        f"- All tested variants have required stabilizer: `{report['robustness_summary']['all_tested_variants_have_required_stabilizer']}`.\n"
        f"- All tested variants have unique d5 global maxima: `{report['robustness_summary']['all_tested_variants_have_unique_d5_global_maxima']}`.\n"
        f"- Maximum energy by variant: `{maximum_energy_by_variant}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: robust inside tested Hebbian-family readouts, but teacher trace and learning law remain supplied inputs.\n"
        "- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
