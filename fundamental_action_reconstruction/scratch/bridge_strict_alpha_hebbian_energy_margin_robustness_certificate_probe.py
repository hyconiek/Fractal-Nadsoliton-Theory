#!/usr/bin/env python3
"""Scratch probe: exact Hebbian energy-margin robustness certificate for d5.

The rule-variant robustness probe showed that six Hebbian-family readouts all
make the d5 teacher orbit the unique global maximum.  This probe adds the next
finite proof layer: compute the exact positive gap between d5 maxima and the
best non-d5 supports, plus simple perturbation radii that preserve the winner.

Finite result:
- Binary and centered readouts have exact non-d5 energy gap 8.
- Bipolar readouts have exact non-d5 energy gap 32.
- The closest non-d5 competitor count is 24 in every tested readout.
- Therefore the d5 maximum is isolated by a positive finite margin under the
  supplied d5 teacher trace and tested Hebbian-family readouts.

No false pass: these are conditional margins for supplied teacher/update
premises, not a derivation of the teacher trace, not a Hebbian-law theorem, and
not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_energy_margin_robustness_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_energy_margin_robustness_certificate_report.md"

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


def margin_certificate(
    variant_name: str,
    vector_kind: str,
    zero_self: bool,
    supports: list[Support],
    teacher: list[Support],
) -> dict[str, Any]:
    teacher_set = set(teacher)
    weights = hebbian_weights(teacher, vector_kind, zero_self)
    energy_by_support = {support: learned_energy(support, weights) for support in supports}
    maximum = max(energy_by_support.values())
    d5_maximizers = sorted(support for support in teacher if energy_by_support[support] == maximum)
    non_d5_rows = [(energy, support) for support, energy in energy_by_support.items() if support not in teacher_set]
    best_non_d5_energy = max(energy for energy, _support in non_d5_rows)
    best_non_d5_supports = sorted(support for energy, support in non_d5_rows if energy == best_non_d5_energy)
    gap = maximum - best_non_d5_energy
    gap_histogram = Counter(maximum - energy for energy, _support in non_d5_rows)
    support_energy_perturbation_radius = gap / 2
    entrywise_weight_perturbation_bound = gap / (2 * ACTIVE_COUNT * ACTIVE_COUNT)
    return {
        "variant_name": variant_name,
        "vector_kind": vector_kind,
        "zero_self": zero_self,
        "unit_stabilizer": unit_stabilizer(weights),
        "maximum_energy": fraction_text(maximum),
        "d5_maximizer_count": len(d5_maximizers),
        "d5_maximizers_equal_teacher_orbit": set(d5_maximizers) == teacher_set,
        "best_non_d5_energy": fraction_text(best_non_d5_energy),
        "energy_gap_to_best_non_d5": fraction_text(gap),
        "best_non_d5_competitor_count": len(best_non_d5_supports),
        "best_non_d5_competitor_examples": [list(support) for support in best_non_d5_supports[:8]],
        "non_d5_gap_histogram": {
            fraction_text(delta): count for delta, count in sorted(gap_histogram.items())
        },
        "support_energy_perturbation_safe_radius_strict": fraction_text(support_energy_perturbation_radius),
        "entrywise_weight_perturbation_sufficient_bound_strict": fraction_text(entrywise_weight_perturbation_bound),
        "positive_margin_certificate": gap > 0
        and set(d5_maximizers) == teacher_set
        and len(best_non_d5_supports) > 0,
    }


def main() -> None:
    supports = all_supports()
    teacher = teacher_orbit(TARGET_STEP)
    certificates = [
        margin_certificate(name, vector_kind, zero_self, supports, teacher)
        for name, vector_kind, zero_self in VARIANTS
    ]
    gap_by_variant = {row["variant_name"]: row["energy_gap_to_best_non_d5"] for row in certificates}
    best_non_d5_count_by_variant = {row["variant_name"]: row["best_non_d5_competitor_count"] for row in certificates}
    support_radius_by_variant = {row["variant_name"]: row["support_energy_perturbation_safe_radius_strict"] for row in certificates}
    entrywise_bound_by_variant = {row["variant_name"]: row["entrywise_weight_perturbation_sufficient_bound_strict"] for row in certificates}

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_ENERGY_MARGIN_ROBUSTNESS_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "status": "positive-d5-energy-margin-across-tested-hebbian-readouts-not-origin-theorem",
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
        "margin_certificates": certificates,
        "margin_summary": {
            "tested_variant_count": len(certificates),
            "all_tested_variants_have_positive_margin": all(row["positive_margin_certificate"] for row in certificates),
            "all_tested_variants_have_required_stabilizer": all(row["unit_stabilizer"] == MINIMAL_REQUIRED_SUBGROUP for row in certificates),
            "gap_by_variant": gap_by_variant,
            "minimum_gap": fraction_text(min(Fraction(row["energy_gap_to_best_non_d5"]) for row in certificates)),
            "best_non_d5_count_by_variant": best_non_d5_count_by_variant,
            "support_energy_perturbation_safe_radius_by_variant_strict": support_radius_by_variant,
            "entrywise_weight_perturbation_sufficient_bound_by_variant_strict": entrywise_bound_by_variant,
        },
        "candidate_source_interpretation": {
            "finite_gain": "The d5 maximum is not merely unique; it is separated from the best non-d5 supports by an exact positive gap in every tested Hebbian-family readout.",
            "conditional_positive_result": "Given the supplied d5 teacher/self-record trace and one of the tested Hebbian readouts, small enough support-energy or entrywise-weight perturbations cannot change the d5 winner.",
            "honest_limit": "The margin protects a supplied readout; it does not derive the teacher trace, the Hebbian learning law, or the perturbation model from strict geometry.",
            "relation_to_previous_probe": "The previous rule-variant audit checked uniqueness; this probe records exact gaps and sufficient perturbation bounds.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may carry a finite Hebbian-family self-record with a positive energy margin.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to choose the teacher trace, learning law, or perturbation norm.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "No theorem derives the perturbation model or robustness norm from strict geometry.",
            "Positive margins across six tested readouts are not an exhaustive theorem over all learning laws or all perturbations.",
            "This is a conditional margin certificate, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the teacher trace, Hebbian-family update, or perturbation norm internally; otherwise keep the positive margin as a conditional stability certificate.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian energy-margin robustness certificate probe\n\n"
        "Status: positive d5 energy margins across tested Hebbian readouts; not an origin theorem.\n\n"
        f"- Tested variants: `{len(certificates)}`.\n"
        f"- Gap by variant: `{gap_by_variant}`.\n"
        f"- Minimum gap: `{report['margin_summary']['minimum_gap']}`.\n"
        f"- Best non-d5 competitor count by variant: `{best_non_d5_count_by_variant}`.\n"
        f"- Support-energy perturbation safe radius by variant: `{support_radius_by_variant}`.\n"
        f"- Entrywise-weight perturbation sufficient bound by variant: `{entrywise_bound_by_variant}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: exact positive margins protect the supplied Hebbian readouts, but the teacher trace, learning law, and perturbation norm remain supplied inputs.\n"
        "- No false pass: no strict d5-origin theorem, no Hebbian-law theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
