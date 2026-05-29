#!/usr/bin/env python3
"""Scratch probe: Hebbian weight stabilizer self-record audit for d5.

Previous probes showed that the d5 selector becomes invariant after reducing
Aut(Z_12) to {1,11}, and that an externally supplied cycle/locality metric has
exactly that stabilizer.  This probe asks a narrower, more internal finite
question:

    if the Hebbian teacher trace is already supplied, does the learned weight
    matrix itself carry the same {1,11} subgroup as a self-record?

Finite result:
- The centered zero-self Hebbian matrix learned from the d5 teacher orbit has
  exact unit stabilizer {1,11} under multiplicative units of Z_12.
- The obstruction is explicit: units 5 and 7 map first-row distance 1 to
  distance 5, but the learned weights at those distances are -25/12 and 23/12.
- This is useful because the learned Hebbian state can preserve the same
  symmetry reduction required by the d5 selector without adding a separate
  informational layer underneath the nadsoliton.
- It is not a strict source theorem: the d5 teacher/self-record trace is still
  supplied as input, and the same {1,11} stabilizer also appears for the
  contiguous unit-related teacher orbit.

No false pass: this is a conditional self-record/stabilizer audit, not a d5
origin theorem, not strict selector closure, and not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_weight_stabilizer_self_record_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_weight_stabilizer_self_record_audit_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
TARGET_D5_STEP = 5
CANONICAL_STEPS = {
    "contiguous_step_1_11": 1,
    "parity_minus_one_step_2_10": 2,
    "third_step_3_9_degenerate": 3,
    "fourth_step_4_8_degenerate": 4,
    "fifth_step_d5_step_5_7": 5,
    "nyquist_step_6": 6,
}
MINIMAL_REQUIRED_SUBGROUP = [1, 11]

Support = tuple[int, ...]
Matrix = list[list[Fraction]]


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def canonical_support(nodes: set[int]) -> Support:
    return tuple(sorted(node % N for node in nodes))


def teacher_support(step: int, translate: int) -> Support:
    return canonical_support({translate + index * step for index in range(ACTIVE_COUNT)})


def teacher_orbit(step: int) -> list[Support]:
    return sorted({teacher_support(step, translate) for translate in range(N)})


def centered_vector(support: Support) -> list[Fraction]:
    active = set(support)
    mean = Fraction(ACTIVE_COUNT, N)
    return [Fraction(1 if node in active else 0) - mean for node in range(N)]


def hebbian_weights(teacher: list[Support]) -> Matrix:
    weights = [[Fraction(0) for _ in range(N)] for _ in range(N)]
    for support in teacher:
        vector = centered_vector(support)
        for left in range(N):
            for right in range(N):
                if left != right:
                    weights[left][right] += vector[left] * vector[right]
    return weights


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


def folded_distance(distance: int) -> int:
    residue = distance % N
    return min(residue, N - residue)


def first_row_by_folded_distance(weights: Matrix) -> dict[str, str]:
    return {str(distance): fraction_text(weights[0][distance]) for distance in range(1, N // 2 + 1)}


def unit_distance_action(unit: int) -> dict[str, int]:
    return {str(distance): folded_distance(unit * distance) for distance in range(1, N // 2 + 1)}


def unit_obstruction_rows(weights: Matrix) -> list[dict[str, Any]]:
    rows = []
    for unit in AUT_UNITS:
        mismatches = []
        for distance in range(1, N // 2 + 1):
            image_distance = folded_distance(unit * distance)
            source_weight = weights[0][distance]
            image_weight = weights[0][image_distance]
            if source_weight != image_weight:
                mismatches.append(
                    {
                        "source_folded_distance": distance,
                        "image_folded_distance": image_distance,
                        "source_weight": fraction_text(source_weight),
                        "image_weight": fraction_text(image_weight),
                    }
                )
        rows.append(
            {
                "unit": unit,
                "preserves_weight_matrix": matrix_preserved_by_unit(weights, unit),
                "folded_distance_action": unit_distance_action(unit),
                "mismatch_count_on_first_row_folded_distances": len(mismatches),
                "mismatch_examples": mismatches[:6],
            }
        )
    return rows


def class_certificate(class_name: str, step: int) -> dict[str, Any]:
    teacher = teacher_orbit(step)
    weights = hebbian_weights(teacher)
    stabilizer = unit_stabilizer(weights)
    return {
        "class_name": class_name,
        "step": step,
        "teacher_orbit_size": len(teacher),
        "teacher_orbit_examples": [list(support) for support in teacher[:4]],
        "first_row_by_index_distance_0_to_11": first_row_text(weights),
        "first_row_by_folded_distance_1_to_6": first_row_by_folded_distance(weights),
        "unit_stabilizer": stabilizer,
        "equals_minimal_required_subgroup": stabilizer == MINIMAL_REQUIRED_SUBGROUP,
        "unit_obstruction_rows": unit_obstruction_rows(weights),
    }


def main() -> None:
    class_certificates = {
        class_name: class_certificate(class_name, step)
        for class_name, step in CANONICAL_STEPS.items()
    }
    d5_certificate = class_certificates["fifth_step_d5_step_5_7"]
    contiguous_certificate = class_certificates["contiguous_step_1_11"]
    full_stabilizer_classes = sorted(
        name for name, certificate in class_certificates.items()
        if certificate["unit_stabilizer"] == AUT_UNITS
    )
    minimal_stabilizer_classes = sorted(
        name for name, certificate in class_certificates.items()
        if certificate["unit_stabilizer"] == MINIMAL_REQUIRED_SUBGROUP
    )

    d5_weights = hebbian_weights(teacher_orbit(TARGET_D5_STEP))
    unit5_row = next(row for row in unit_obstruction_rows(d5_weights) if row["unit"] == 5)
    unit7_row = next(row for row in unit_obstruction_rows(d5_weights) if row["unit"] == 7)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_WEIGHT_STABILIZER_SELF_RECORD_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "d5-hebbian-weight-carries-required-subgroup-conditionally-not-origin-theorem",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "automorphism_units": AUT_UNITS,
            "canonical_step_classes": CANONICAL_STEPS,
            "minimal_required_subgroup_from_previous_probe": MINIMAL_REQUIRED_SUBGROUP,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "d5_weight_self_record_certificate": {
            "unit_stabilizer": d5_certificate["unit_stabilizer"],
            "equals_minimal_required_subgroup": d5_certificate["equals_minimal_required_subgroup"],
            "first_row_by_folded_distance_1_to_6": d5_certificate["first_row_by_folded_distance_1_to_6"],
            "unit_5_obstruction_example": unit5_row["mismatch_examples"][0],
            "unit_7_obstruction_example": unit7_row["mismatch_examples"][0],
            "readout": "The learned d5 Hebbian matrix itself preserves exactly units {1,11}; units 5 and 7 are rejected by exact rational first-row weight mismatches.",
        },
        "cross_teacher_stabilizer_audit": {
            "class_certificates": class_certificates,
            "classes_with_minimal_required_stabilizer": minimal_stabilizer_classes,
            "classes_with_full_aut_stabilizer": full_stabilizer_classes,
            "contiguous_has_same_stabilizer_as_d5": contiguous_certificate["unit_stabilizer"] == d5_certificate["unit_stabilizer"],
            "nonuniqueness_warning": "The {1,11} stabilizer is not unique to d5: the contiguous unit-related teacher orbit has the same stabilizer, so this audit does not derive d5 as origin.",
        },
        "candidate_source_interpretation": {
            "finite_gain": "If a d5 teacher/self-record trace is supplied, its centered zero-self Hebbian weight matrix carries the same subgroup needed by the selector stack.",
            "conditional_positive_result": "The subgroup can be read from the learned state itself rather than imposed as an extra informational layer below the nadsoliton.",
            "honest_limit": "The supplied teacher/self-record trace is still the source input; this probe does not derive why the trace is d5 rather than a unit-related or other teacher orbit.",
            "relation_to_cycle_metric_probe": "The previous cycle-metric probe supplied {1,11} from a locality record; this probe shows the d5 Hebbian weight record also has {1,11} as its exact unit stabilizer.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may carry a finite Hebbian weight/self-record pattern.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply the weight record.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a Hebbian learning law as strict-core dynamics.",
            "The same {1,11} stabilizer also appears for the contiguous teacher orbit, so stabilizer evidence alone is not a d5 origin theorem.",
            "This is a conditional self-record/stabilizer audit, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, cycle metric source, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the d5 teacher/self-record trace or Hebbian update law internally; otherwise keep the Hebbian stabilizer as a conditional downstream record only.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian weight stabilizer self-record audit probe\n\n"
        "Status: d5 Hebbian weight carries the required {1,11} subgroup conditionally; not a d5 origin theorem.\n\n"
        f"- D5 unit stabilizer: `{d5_certificate['unit_stabilizer']}`.\n"
        f"- D5 stabilizer equals required subgroup: `{d5_certificate['equals_minimal_required_subgroup']}`.\n"
        f"- D5 folded-distance weights: `{d5_certificate['first_row_by_folded_distance_1_to_6']}`.\n"
        f"- Classes with required {MINIMAL_REQUIRED_SUBGROUP} stabilizer: `{minimal_stabilizer_classes}`.\n"
        f"- Classes with full Aut stabilizer: `{full_stabilizer_classes}`.\n"
        f"- Contiguous has same stabilizer as d5: `{contiguous_certificate['unit_stabilizer'] == d5_certificate['unit_stabilizer']}`.\n"
        f"- Unit 5 obstruction example: `{unit5_row['mismatch_examples'][0]}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: learned d5 Hebbian weight can carry the subgroup record, but the d5 teacher trace and Hebbian law remain supplied inputs.\n"
        "- No false pass: no strict d5-origin theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
