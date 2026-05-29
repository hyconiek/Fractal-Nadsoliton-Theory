#!/usr/bin/env python3
"""Scratch probe: Hebbian teacher-orbit nonuniqueness proof audit.

The previous finite landscape probe proved that, after the d5 teacher trace is
supplied, the centered zero-self Hebbian energy has exactly the d5 orbit as its
global-maximum family.  This probe checks the missing source question more
honestly: is that property special to d5, or does Hebbian energy simply make
whatever translated teacher orbit was supplied into its own attractor/maximizer?

Finite result on Z_12 five-node supports:
- the admissible step-generated translated teacher families are step classes
  {1,11}, {2,10}, and {5,7};
- each class has 12 teacher supports;
- for each class, exact centered zero-self Hebbian energy has exactly those 12
  teacher supports as its global maximizers;
- therefore the Hebbian energy-landscape theorem is an associative-memory
  theorem, not a derivation of the fifth-step/d5 teacher source.

No false pass: d5 remains selected only if the d5 teacher/self-record trace is
already supplied.  This file does not discharge QW-2191 or close ToE.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_teacher_orbit_nonuniqueness_proof_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_teacher_orbit_nonuniqueness_proof_report.md"

N = 12
ACTIVE_COUNT = 5
RHO = Fraction(ACTIVE_COUNT, N)
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
STEP_CLASSES = {
    "contiguous_step_1_11": [1, 11],
    "parity_minus_one_step_2_10": [2, 10],
    "fifth_step_d5_step_5_7": [5, 7],
}

Support = tuple[int, ...]
DistanceHistogram = tuple[int, ...]
Matrix = list[list[Fraction]]


def all_supports() -> list[Support]:
    return list(combinations(range(N), ACTIVE_COUNT))


def canonical_support(nodes: list[int] | tuple[int, ...]) -> Support:
    return tuple(sorted(int(node) % N for node in nodes))


def step_teacher_orbit(step: int) -> list[Support]:
    return sorted(
        {
            canonical_support([(start + step * index) % N for index in range(ACTIVE_COUNT)])
            for start in range(N)
            if len({(start + step * index) % N for index in range(ACTIVE_COUNT)}) == ACTIVE_COUNT
        }
    )


def centered_activity(support: Support) -> list[Fraction]:
    active = set(support)
    return [Fraction(1) - RHO if node in active else -RHO for node in range(N)]


def build_centered_zero_self_hebbian(teacher: list[Support]) -> Matrix:
    weights = [[Fraction(0) for _col in range(N)] for _row in range(N)]
    for support in teacher:
        activity = centered_activity(support)
        for row in range(N):
            for col in range(N):
                weights[row][col] += activity[row] * activity[col]
    for index in range(N):
        weights[index][index] = Fraction(0)
    return weights


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def learned_energy(support: Support, weights: Matrix) -> Fraction:
    return sum(weights[row][col] for row in support for col in support)


def distance_histogram(support: Support) -> DistanceHistogram:
    counts = Counter()
    for left, right in combinations(support, 2):
        distance = (right - left) % N
        folded = min(distance, N - distance)
        counts[folded] += 1
    return tuple(counts[distance] for distance in range(1, N // 2 + 1))


def histogram_energy(histogram: DistanceHistogram, distance_weights: dict[int, Fraction]) -> Fraction:
    return 2 * sum(histogram[distance - 1] * distance_weights[distance] for distance in range(1, N // 2 + 1))


def orbit_class_certificate(name: str, steps: list[int], supports: list[Support]) -> dict[str, Any]:
    teacher = step_teacher_orbit(steps[0])
    step_orbits_equal = all(step_teacher_orbit(step) == teacher for step in steps)
    weights = build_centered_zero_self_hebbian(teacher)
    teacher_set = set(teacher)
    energy_by_support = {support: learned_energy(support, weights) for support in supports}
    maximum = max(energy_by_support.values())
    maximizers = sorted(support for support, energy in energy_by_support.items() if energy == maximum)
    energy_histogram = Counter(energy_by_support.values())
    histogram_counter = Counter(distance_histogram(support) for support in supports)
    distance_weights = {distance: weights[0][distance] for distance in range(1, N // 2 + 1)}
    top_histogram = max(
        histogram_counter,
        key=lambda histogram: (histogram_energy(histogram, distance_weights), histogram),
    )
    replay_failures = [
        support
        for support in supports
        if learned_energy(support, weights) != histogram_energy(distance_histogram(support), distance_weights)
    ]
    return {
        "class_name": name,
        "steps": steps,
        "step_orbits_equal": step_orbits_equal,
        "teacher_support_count": len(teacher),
        "teacher_supports": [list(row) for row in teacher],
        "first_row_by_cyclic_distance": [fraction_text(value) for value in weights[0]],
        "distance_weights_d1_to_d6": {str(distance): fraction_text(weight) for distance, weight in distance_weights.items()},
        "energy_level_count": len(energy_histogram),
        "energy_level_histogram": {fraction_text(energy): count for energy, count in sorted(energy_histogram.items())},
        "maximum_energy": fraction_text(maximum),
        "global_maximizer_count": len(maximizers),
        "global_maximizers_equal_teacher_orbit": set(maximizers) == teacher_set,
        "non_teacher_global_maximizer_count": sum(1 for support in maximizers if support not in teacher_set),
        "top_distance_histogram": {
            "distance_histogram_d1_to_d6": list(top_histogram),
            "support_count": histogram_counter[top_histogram],
            "energy": fraction_text(histogram_energy(top_histogram, distance_weights)),
        },
        "histogram_class_count": len(histogram_counter),
        "direct_energy_equals_histogram_replay_for_all_supports": not replay_failures,
        "replay_failure_count": len(replay_failures),
    }


def main() -> None:
    supports = all_supports()
    certificates = {
        name: orbit_class_certificate(name, steps, supports)
        for name, steps in STEP_CLASSES.items()
    }
    all_classes_self_maximize = all(
        certificate["global_maximizers_equal_teacher_orbit"]
        for certificate in certificates.values()
    )
    d5_certificate = certificates["fifth_step_d5_step_5_7"]
    non_d5_self_maximizing_classes = [
        name for name, certificate in certificates.items()
        if name != "fifth_step_d5_step_5_7" and certificate["global_maximizers_equal_teacher_orbit"]
    ]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_TEACHER_ORBIT_NONUNIQUENESS_PROOF_PROBE__NOT_A_THEOREM",
        "status": "finite-hebbian-associative-memory-nonuniqueness-not-d5-source",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "teacher_rule": "translated step-generated five-node supports with five distinct nodes",
            "learning_rule": "centered zero-self Hebbian W from the supplied teacher orbit",
            "energy": "E(S)=sum_{i,j in S} W_ij = x_S^T W x_S, diag(W)=0",
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "orbit_class_certificates": certificates,
        "cross_orbit_nonuniqueness_certificate": {
            "admissible_step_classes": STEP_CLASSES,
            "all_classes_self_maximize": all_classes_self_maximize,
            "non_d5_self_maximizing_classes": non_d5_self_maximizing_classes,
            "d5_maximum_energy": d5_certificate["maximum_energy"],
            "d5_global_maximizer_count": d5_certificate["global_maximizer_count"],
            "honest_obstruction": "Hebbian energy makes the supplied translated teacher orbit its own global maximum; it does not decide which teacher orbit should be supplied.",
        },
        "proof_reading": {
            "finite_theorem": "For each admissible translated step-teacher class tested here, centered zero-self Hebbian energy has exactly that teacher orbit as its global maximizer family.",
            "negative_consequence_for_d5_source": "The d5 landscape theorem is not source-selective: contiguous and parity-minus-one teacher records pass the same self-maximization test.",
            "relation_to_previous_probe": "The previous d5 global-maximum proof is valid conditionally, but this cross-orbit replay proves it is an associative-memory certificate rather than a d5-source theorem.",
            "remaining_gap": "A strict-side source is still needed to select the fifth-step/d5 teacher trace before Hebbian learning begins.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain one finite self-recorded teacher trace and internally stabilize it Hebbian-style.",
            "forbidden_reading": "The family of possible teacher traces is not a separate informational layer underneath the nadsoliton.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives the centered zero-self Hebbian energy from strict geometry.",
            "Hebbian self-maximization is nonunique across teacher orbits and therefore is not a d5-source theorem.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict internal reason that the supplied self-record teacher trace must be the fifth-step/d5 class rather than another Hebbian-self-maximizing teacher class.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian teacher-orbit nonuniqueness proof probe\n\n"
        "Status: finite Hebbian associative-memory nonuniqueness; not a d5 source theorem.\n\n"
        f"- Supports scanned: `{len(supports)}` five-active-node states on `Z_12`.\n"
        f"- Admissible step classes: `{STEP_CLASSES}`.\n"
        f"- All tested teacher classes self-maximize: `{all_classes_self_maximize}`.\n"
        f"- Non-d5 self-maximizing classes: `{non_d5_self_maximizing_classes}`.\n"
        f"- d5 class maximum energy: `{d5_certificate['maximum_energy']}` with maximizer count `{d5_certificate['global_maximizer_count']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: Hebbian energy stabilizes whichever translated teacher orbit is supplied; it does not choose d5 by itself.\n"
        "- No false pass: no d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
