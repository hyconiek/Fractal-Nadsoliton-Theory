#!/usr/bin/env python3
"""Scratch probe: Hebbian teacher resonance-label selector audit.

The previous teacher-orbit nonuniqueness probe showed that centered Hebbian
energy self-maximizes every supplied translated teacher orbit tested here.  This
probe asks the next selector question:

    Can the "nadsoliton seeks highest resonance" idea distinguish which
    Hebbian teacher orbit should be supplied?

Finite answer on five-node supports of Z_12:
- Hebbian self-maximization alone leaves three translated teacher classes:
  contiguous {1,11}, parity-minus-one {2,10}, and fifth-step/d5 {5,7}.
- Each class carries a learned Fourier resonance label: contiguous -> k=1,
  parity-minus-one -> k=6, fifth-step/d5 -> k=5.
- Unqualified highest learned label/power selects the parity/Nyquist class k=6,
  not d5.
- If one adds an explicit non-Nyquist unit-channel / primitive-resonance premise,
  the highest remaining learned label is k=5 and d5 is selected over contiguous.

No false pass: the non-Nyquist unit-channel premise is an extra selector premise,
not derived here from strict nadsoliton geometry and not a QW-2191 discharge.
"""
from __future__ import annotations

import json
import math
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_teacher_resonance_label_selector_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_teacher_resonance_label_selector_audit_report.md"

N = 12
ACTIVE_COUNT = 5
RHO = Fraction(ACTIVE_COUNT, N)
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
INDEPENDENT_MODES = [1, 2, 3, 4, 5, 6]
STEP_CLASSES = {
    "contiguous_step_1_11": [1, 11],
    "parity_minus_one_step_2_10": [2, 10],
    "fifth_step_d5_step_5_7": [5, 7],
}
TARGET_D5_CLASS = "fifth_step_d5_step_5_7"
FLOAT_TOL = 1e-9

Support = tuple[int, ...]
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


def learned_energy(support: Support, weights: Matrix) -> Fraction:
    return sum(weights[row][col] for row in support for col in support)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def coherent_power(support: Support, mode: int) -> float:
    real = sum(math.cos(2.0 * math.pi * mode * node / N) for node in support)
    imag = sum(math.sin(2.0 * math.pi * mode * node / N) for node in support)
    return real * real + imag * imag


def rounded(value: float) -> float:
    return round(value, 12)


def averaged_mode_scores(teacher: list[Support]) -> dict[int, float]:
    return {
        mode: sum(coherent_power(support, mode) for support in teacher) / len(teacher)
        for mode in INDEPENDENT_MODES
    }


def leading_modes(scores: dict[int, float]) -> list[int]:
    maximum = max(scores.values())
    return [mode for mode in INDEPENDENT_MODES if abs(scores[mode] - maximum) <= FLOAT_TOL]


def orbit_certificate(name: str, steps: list[int], supports: list[Support]) -> dict[str, Any]:
    teacher = step_teacher_orbit(steps[0])
    weights = build_centered_zero_self_hebbian(teacher)
    energy_by_support = {support: learned_energy(support, weights) for support in supports}
    maximum = max(energy_by_support.values())
    maximizers = sorted(support for support, energy in energy_by_support.items() if energy == maximum)
    scores = averaged_mode_scores(teacher)
    leaders = leading_modes(scores)
    return {
        "class_name": name,
        "steps": steps,
        "step_orbits_equal": all(step_teacher_orbit(step) == teacher for step in steps),
        "teacher_support_count": len(teacher),
        "hebbian_global_maximizers_equal_teacher_orbit": set(maximizers) == set(teacher),
        "hebbian_maximum_energy": fraction_text(maximum),
        "hebbian_global_maximizer_count": len(maximizers),
        "learned_average_mode_scores": {str(mode): rounded(scores[mode]) for mode in INDEPENDENT_MODES},
        "learned_leading_modes": leaders,
        "learned_leading_score": rounded(max(scores.values())),
        "is_unit_step_class": all(math.gcd(step, N) == 1 for step in steps),
        "is_nyquist_labeled": 6 in leaders,
        "teacher_supports_sample": [list(row) for row in teacher[:12]],
    }


def class_names_with_max_label(certificates: dict[str, dict[str, Any]], names: list[str]) -> list[str]:
    maximum_label = max(max(certificates[name]["learned_leading_modes"]) for name in names)
    return [name for name in names if max(certificates[name]["learned_leading_modes"]) == maximum_label]


def main() -> None:
    supports = all_supports()
    certificates = {
        name: orbit_certificate(name, steps, supports)
        for name, steps in STEP_CLASSES.items()
    }
    hebbian_self_max_classes = [
        name for name, certificate in certificates.items()
        if certificate["hebbian_global_maximizers_equal_teacher_orbit"]
    ]
    unit_classes = [name for name in hebbian_self_max_classes if certificates[name]["is_unit_step_class"]]
    non_nyquist_unit_classes = [name for name in unit_classes if not certificates[name]["is_nyquist_labeled"]]
    unqualified_highest_label_winners = class_names_with_max_label(certificates, hebbian_self_max_classes)
    unit_highest_label_winners = class_names_with_max_label(certificates, unit_classes)
    fifth_mode_winners = [name for name in hebbian_self_max_classes if 5 in certificates[name]["learned_leading_modes"]]

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_TEACHER_RESONANCE_LABEL_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "conditional-nonnyquist-unit-resonance-label-selects-d5-not-strict-source",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "teacher_classes": STEP_CLASSES,
            "learning_rule": "centered zero-self Hebbian W from supplied translated teacher orbit",
            "resonance_label": "leading averaged non-DC Fourier mode over the teacher orbit",
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "orbit_class_certificates": certificates,
        "selector_stack_audit": {
            "hebbian_self_max_only": {
                "survivors": hebbian_self_max_classes,
                "survivor_count": len(hebbian_self_max_classes),
                "selects_d5_uniquely": hebbian_self_max_classes == [TARGET_D5_CLASS],
                "verdict": "nonunique: Hebbian self-maximization is associative memory, not a source selector",
            },
            "unit_step_filter": {
                "survivors": unit_classes,
                "survivor_count": len(unit_classes),
                "selects_d5_uniquely": unit_classes == [TARGET_D5_CLASS],
                "verdict": "still nonunique: contiguous and d5 are both unit-step classes",
            },
            "unqualified_highest_learned_label": {
                "winners": unqualified_highest_label_winners,
                "winning_labels": [certificates[name]["learned_leading_modes"] for name in unqualified_highest_label_winners],
                "selects_d5_uniquely": unqualified_highest_label_winners == [TARGET_D5_CLASS],
                "verdict": "negative for d5: unqualified highest learned label selects the k=6 parity/Nyquist class",
            },
            "non_nyquist_unit_highest_label": {
                "premise": "exclude non-unit/Nyquist channel, then choose the highest learned label among unit translated teacher classes",
                "eligible_classes": non_nyquist_unit_classes,
                "winners": unit_highest_label_winners,
                "winning_labels": [certificates[name]["learned_leading_modes"] for name in unit_highest_label_winners],
                "selects_d5_uniquely": unit_highest_label_winners == [TARGET_D5_CLASS],
                "verdict": "conditional positive: with this explicit primitive non-Nyquist resonance premise, d5 beats contiguous",
            },
            "explicit_fifth_mode_lock": {
                "winners": fifth_mode_winners,
                "selects_d5_uniquely": fifth_mode_winners == [TARGET_D5_CLASS],
                "verdict": "positive but tautological/source-premised: a fifth-mode lock selects the fifth-step teacher class",
            },
        },
        "proof_reading": {
            "finite_gain": "The teacher orbit carries a computable resonance label; d5 is exactly the unit teacher class labeled by k=5 rather than k=1.",
            "negative_result": "Hebbian self-maximization plus unqualified highest resonance is insufficient and points to k=6 parity/Nyquist instead.",
            "conditional_positive_result": "Hebbian self-maximization plus an explicit primitive non-Nyquist/unit-channel highest-label premise selects d5 uniquely among translated step-teacher classes.",
            "remaining_gap": "The primitive non-Nyquist/unit-channel premise is not derived here from strict nadsoliton geometry.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain one self-recorded resonance label for its finite teacher trace.",
            "forbidden_reading": "The resonance-label selector is not a separate informational layer underneath the nadsoliton.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives the primitive non-Nyquist/unit-channel premise from strict geometry.",
            "Unqualified highest resonance remains insufficient: it selects the k=6 parity/Nyquist teacher class, not d5.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the primitive non-Nyquist/unit-channel premise internally; otherwise keep it marked as an explicit non-strict selector assumption.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian teacher resonance-label selector audit probe\n\n"
        "Status: non-Nyquist unit highest-label selects d5 conditionally; not a strict source theorem.\n\n"
        f"- Supports scanned: `{len(supports)}` five-active-node states on `Z_12`.\n"
        f"- Hebbian self-max-only survivors: `{hebbian_self_max_classes}`.\n"
        f"- Unit-step survivors: `{unit_classes}`.\n"
        f"- Unqualified highest learned-label winners: `{unqualified_highest_label_winners}`.\n"
        f"- Non-Nyquist unit highest-label winners: `{unit_highest_label_winners}`.\n"
        f"- Explicit fifth-mode-lock winners: `{fifth_mode_winners}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: d5 is selected only after adding a primitive non-Nyquist/unit-channel or fifth-mode premise.\n"
        "- No false pass: no d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
