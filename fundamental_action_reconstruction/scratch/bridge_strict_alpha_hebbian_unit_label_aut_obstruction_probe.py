#!/usr/bin/env python3
"""Scratch probe: Aut(Z_12) obstruction for unit resonance-label selection.

The previous resonance-label selector audit found a conditional positive stack:
Hebbian self-maximization + primitive non-Nyquist/unit-channel + highest learned
label selects d5 over the contiguous teacher class.  This probe checks whether
that last numeric-label preference can be strict-core invariant.

Finite result:
- The unit teacher classes contiguous {1,11} and fifth-step/d5 {5,7} are in the
  same Aut(Z_12) orbit.  Multiplication by unit 5 sends step 1 to step 5, and
  multiplication by unit 5 sends learned label k=1 to k=5.
- Therefore the numeric rule "among unit classes choose the larger label, k=5
  over k=1" is not Aut(Z_12)-invariant.  It imports a generator/label-order
  choice.
- The parity/Nyquist class is separate because step 2 is non-unit and label k=6
  is fixed under all units.

No false pass: the d5 choice can be recovered from a unit highest-label premise,
but this probe shows that premise is symmetry-breaking/non-strict unless a strict
internal orientation/generator/label-order source is supplied.
"""
from __future__ import annotations

import json
import math
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_unit_label_aut_obstruction_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_unit_label_aut_obstruction_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
INDEPENDENT_MODES = [1, 2, 3, 4, 5, 6]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
FLOAT_TOL = 1e-9
STEP_CLASSES = {
    "contiguous_step_1_11": [1, 11],
    "parity_minus_one_step_2_10": [2, 10],
    "fifth_step_d5_step_5_7": [5, 7],
}
TARGET_D5_CLASS = "fifth_step_d5_step_5_7"

Support = tuple[int, ...]


def canonical_support(nodes: list[int] | tuple[int, ...]) -> Support:
    return tuple(sorted(int(node) % N for node in nodes))


def all_supports() -> list[Support]:
    return list(combinations(range(N), ACTIVE_COUNT))


def step_teacher_orbit(step: int) -> list[Support]:
    return sorted(
        {
            canonical_support([(start + step * index) % N for index in range(ACTIVE_COUNT)])
            for start in range(N)
            if len({(start + step * index) % N for index in range(ACTIVE_COUNT)}) == ACTIVE_COUNT
        }
    )


def multiply_support(support: Support, unit: int) -> Support:
    return canonical_support([unit * node for node in support])


def multiply_orbit(orbit: list[Support], unit: int) -> list[Support]:
    return sorted({multiply_support(support, unit) for support in orbit})


def folded_mode(mode: int) -> int:
    residue = mode % N
    if residue == 0:
        return 0
    return min(residue, N - residue)


def transformed_mode_label(mode: int, unit: int) -> int:
    return folded_mode(mode * unit)


def coherent_power(support: Support, mode: int) -> float:
    real = sum(math.cos(2.0 * math.pi * mode * node / N) for node in support)
    imag = sum(math.sin(2.0 * math.pi * mode * node / N) for node in support)
    return real * real + imag * imag


def averaged_mode_scores(orbit: list[Support]) -> dict[int, float]:
    return {
        mode: sum(coherent_power(support, mode) for support in orbit) / len(orbit)
        for mode in INDEPENDENT_MODES
    }


def leading_modes(scores: dict[int, float]) -> list[int]:
    maximum = max(scores.values())
    return [mode for mode in INDEPENDENT_MODES if abs(scores[mode] - maximum) <= FLOAT_TOL]


def rounded(value: float) -> float:
    return round(value, 12)


def class_name_for_orbit(orbit: list[Support], class_orbits: dict[str, list[Support]]) -> str | None:
    for name, candidate in class_orbits.items():
        if orbit == candidate:
            return name
    return None


def main() -> None:
    supports = all_supports()
    class_orbits = {name: step_teacher_orbit(steps[0]) for name, steps in STEP_CLASSES.items()}
    class_labels = {name: leading_modes(averaged_mode_scores(orbit)) for name, orbit in class_orbits.items()}

    aut_action_rows: list[dict[str, Any]] = []
    for name, orbit in class_orbits.items():
        for unit in UNITS:
            image_orbit = multiply_orbit(orbit, unit)
            image_class = class_name_for_orbit(image_orbit, class_orbits)
            image_labels_from_action = sorted({transformed_mode_label(label, unit) for label in class_labels[name]})
            aut_action_rows.append(
                {
                    "source_class": name,
                    "unit": unit,
                    "image_class": image_class,
                    "source_leading_labels": class_labels[name],
                    "transformed_leading_labels": image_labels_from_action,
                    "image_class_leading_labels": class_labels.get(image_class),
                    "label_equivariance_holds": image_class is not None
                    and image_labels_from_action == class_labels[image_class],
                }
            )

    unit_class_names = ["contiguous_step_1_11", TARGET_D5_CLASS]
    numeric_unit_highest_winner = max(unit_class_names, key=lambda name: max(class_labels[name]))
    mapped_by_unit_5 = class_name_for_orbit(multiply_orbit(class_orbits["contiguous_step_1_11"], 5), class_orbits)
    label_1_by_unit_5 = transformed_mode_label(1, 5)

    aut_invariant_functions_on_unit_pair = {
        "unit_pair_orbit": unit_class_names,
        "aut_connecting_unit": 5,
        "contiguous_maps_to": mapped_by_unit_5,
        "label_1_maps_to": label_1_by_unit_5,
        "d5_numeric_label": class_labels[TARGET_D5_CLASS],
        "contiguous_numeric_label": class_labels["contiguous_step_1_11"],
        "numeric_highest_label_winner": numeric_unit_highest_winner,
        "is_numeric_highest_label_aut_invariant": False,
        "reason": "An Aut(Z_12) unit maps the contiguous unit class and label k=1 to the d5 unit class and label k=5, so choosing k=5 over k=1 uses a non-invariant numeric label order.",
    }

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_UNIT_LABEL_AUT_OBSTRUCTION_PROBE__NOT_A_THEOREM",
        "status": "unit-highest-label-d5-selector-is-aut-breaking-not-strict",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": len(supports),
            "automorphism_units": UNITS,
            "teacher_classes": STEP_CLASSES,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "class_label_replay": {
            name: {
                "teacher_support_count": len(orbit),
                "learned_leading_labels": class_labels[name],
                "learned_average_mode_scores": {
                    str(mode): rounded(score) for mode, score in averaged_mode_scores(orbit).items()
                },
            }
            for name, orbit in class_orbits.items()
        },
        "aut_z12_action_audit": {
            "rows": aut_action_rows,
            "all_class_images_identified": all(row["image_class"] is not None for row in aut_action_rows),
            "all_label_equivariance_checks_pass": all(row["label_equivariance_holds"] for row in aut_action_rows),
            "unit_classes_in_same_aut_orbit": mapped_by_unit_5 == TARGET_D5_CLASS,
            "parity_class_separate_under_units": all(
                row["image_class"] == "parity_minus_one_step_2_10"
                for row in aut_action_rows
                if row["source_class"] == "parity_minus_one_step_2_10"
            ),
        },
        "numeric_label_selector_obstruction": aut_invariant_functions_on_unit_pair,
        "proof_reading": {
            "finite_gain": "The Aut(Z_12) action on teacher classes and learned labels is computed explicitly and label-equivariantly.",
            "negative_result": "The conditional rule 'choose the larger non-Nyquist unit label' is not Aut-invariant; it distinguishes two Aut-equivalent unit classes.",
            "relation_to_previous_probe": "The previous non-Nyquist/unit highest-label selector selects d5, but this probe marks it as a symmetry-breaking premise rather than strict-core output.",
            "remaining_gap": "A strict internal orientation, generator, or label-order source is still needed to make k=5 preferable to k=1.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may contain an internal orientation/generator record if such a strict source is later proved.",
            "forbidden_reading": "The orientation/generator record is not introduced here as a separate informational substrate beneath the nadsoliton.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a strict orientation/generator/label-order source that prefers k=5 over k=1.",
            "The non-Nyquist unit highest-label selector is explicitly classified as Aut(Z_12)-breaking unless such a source is added.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Search for a strict internal orientation/generator source; without it, keep the k=5-over-k=1 label order as an explicit non-strict selector premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian unit-label Aut obstruction probe\n\n"
        "Status: non-Nyquist unit highest-label d5 selector is Aut-breaking unless an orientation/generator source is added.\n\n"
        f"- Supports scanned: `{len(supports)}` five-active-node states on `Z_12`.\n"
        f"- Automorphism units checked: `{UNITS}`.\n"
        f"- Class labels: `{class_labels}`.\n"
        f"- Unit `5` maps contiguous class to: `{mapped_by_unit_5}` and label `1` to `{label_1_by_unit_5}`.\n"
        f"- All label equivariance checks pass: `{all(row['label_equivariance_holds'] for row in aut_action_rows)}`.\n"
        f"- Numeric highest-label unit winner: `{numeric_unit_highest_winner}`; Aut-invariant: `False`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: selecting k=5 over k=1 needs a symmetry-breaking orientation/generator premise.\n"
        "- No false pass: no strict d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
