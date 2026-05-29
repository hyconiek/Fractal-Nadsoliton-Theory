#!/usr/bin/env python3
"""Scratch probe: Aut(Z_12)-invariant selector no-go for Hebbian unit classes.

The previous unit-label audit showed that the rule "choose k=5 over k=1" is
Aut-breaking because unit 5 maps the contiguous unit teacher class to the d5
unit teacher class.  This probe turns that observation into a finite selector
no-go certificate.

Finite result:
- On the three Hebbian self-maximizing translated teacher classes C={1,11},
  P={2,10}, D={5,7}, Aut(Z_12) has two orbits: {C,D} and {P}.
- Therefore the only Aut-invariant subsets are: empty, {P}, {C,D}, and all
  three classes.
- No Aut-invariant subset is the singleton {D}; no Aut-invariant selector can
  choose d5 from the unit pair {C,D}.
- The previous non-Nyquist/unit highest-label selector is explicitly identified
  as a non-invariant singleton selector.

No false pass: this is a finite strict-symmetry obstruction for this candidate
selector stack, not a theorem about all possible future strict selector sources.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_aut_invariant_selector_no_go_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_aut_invariant_selector_no_go_report.md"

N = 12
ACTIVE_COUNT = 5
UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
CLASSES = {
    "contiguous_step_1_11": [1, 11],
    "parity_minus_one_step_2_10": [2, 10],
    "fifth_step_d5_step_5_7": [5, 7],
}
CLASS_LABELS = {
    "contiguous_step_1_11": [1],
    "parity_minus_one_step_2_10": [6],
    "fifth_step_d5_step_5_7": [5],
}
TARGET_D5_CLASS = "fifth_step_d5_step_5_7"
CONTIGUOUS_CLASS = "contiguous_step_1_11"
PARITY_CLASS = "parity_minus_one_step_2_10"


def folded_step(step: int) -> list[int]:
    residue = step % N
    return sorted({residue, (-residue) % N})


def class_for_step(step: int) -> str:
    folded = folded_step(step)
    for name, steps in CLASSES.items():
        if sorted(step % N for step in steps) == folded:
            return name
    raise ValueError(f"No class for step {step}")


def unit_action_on_class(class_name: str, unit: int) -> str:
    representative_step = CLASSES[class_name][0]
    return class_for_step(representative_step * unit)


def folded_mode(mode: int) -> int:
    residue = mode % N
    if residue == 0:
        return 0
    return min(residue, N - residue)


def unit_action_on_label(label: int, unit: int) -> int:
    return folded_mode(label * unit)


def permutation_for_unit(unit: int) -> dict[str, str]:
    return {class_name: unit_action_on_class(class_name, unit) for class_name in CLASSES}


def orbit_of_class(start: str, action_table: dict[str, dict[str, str]]) -> list[str]:
    return sorted({action_table[str(unit)][start] for unit in UNITS})


def all_subsets(items: list[str]) -> list[list[str]]:
    subsets: list[list[str]] = []
    for size in range(len(items) + 1):
        subsets.extend([list(row) for row in combinations(items, size)])
    return subsets


def image_subset(subset: list[str], unit: int, action_table: dict[str, dict[str, str]]) -> list[str]:
    return sorted({action_table[str(unit)][item] for item in subset})


def is_invariant_subset(subset: list[str], action_table: dict[str, dict[str, str]]) -> bool:
    normalized = sorted(subset)
    return all(image_subset(normalized, unit, action_table) == normalized for unit in UNITS)


def main() -> None:
    class_names = list(CLASSES)
    support_count = len(list(combinations(range(N), ACTIVE_COUNT)))
    action_table = {str(unit): permutation_for_unit(unit) for unit in UNITS}
    label_action_table = {
        str(unit): {
            name: [unit_action_on_label(label, unit) for label in labels]
            for name, labels in CLASS_LABELS.items()
        }
        for unit in UNITS
    }
    orbits = []
    seen: set[str] = set()
    for class_name in class_names:
        if class_name in seen:
            continue
        orbit = orbit_of_class(class_name, action_table)
        seen.update(orbit)
        orbits.append(orbit)

    subset_rows = []
    invariant_subsets = []
    for subset in all_subsets(class_names):
        invariant = is_invariant_subset(subset, action_table)
        if invariant:
            invariant_subsets.append(subset)
        subset_rows.append(
            {
                "subset": subset,
                "is_aut_invariant": invariant,
                "images_by_unit": {str(unit): image_subset(subset, unit, action_table) for unit in UNITS},
            }
        )

    d5_singleton = [TARGET_D5_CLASS]
    unit_pair = sorted([CONTIGUOUS_CLASS, TARGET_D5_CLASS])
    previous_selector_subset = d5_singleton
    previous_selector_images = {
        str(unit): image_subset(previous_selector_subset, unit, action_table) for unit in UNITS
    }

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_AUT_INVARIANT_SELECTOR_NO_GO_PROBE__NOT_A_THEOREM",
        "status": "finite-aut-invariant-selector-no-go-for-d5-unit-pair",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": support_count,
            "automorphism_units": UNITS,
            "teacher_classes": CLASSES,
            "learned_labels": CLASS_LABELS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "aut_action_table": action_table,
        "label_action_table": label_action_table,
        "orbit_partition": {
            "orbits": orbits,
            "unit_pair_orbit": unit_pair,
            "parity_orbit": [PARITY_CLASS],
            "orbit_count": len(orbits),
        },
        "invariant_subset_enumeration": {
            "candidate_subset_count": len(subset_rows),
            "invariant_subset_count": len(invariant_subsets),
            "invariant_subsets": invariant_subsets,
            "subset_rows": subset_rows,
        },
        "d5_selector_no_go_certificate": {
            "d5_singleton": d5_singleton,
            "d5_singleton_is_aut_invariant": is_invariant_subset(d5_singleton, action_table),
            "unit_pair_is_aut_invariant": is_invariant_subset(unit_pair, action_table),
            "aut_invariant_singletons": [subset for subset in invariant_subsets if len(subset) == 1],
            "can_aut_invariant_selector_pick_d5_from_unit_pair": False,
            "reason": "The d5 singleton is not closed under unit 5/7; the smallest Aut-invariant subset containing d5 also contains the contiguous unit class.",
        },
        "previous_selector_classification": {
            "previous_rule": "non-Nyquist/unit highest learned label chooses k=5 over k=1",
            "selected_subset": previous_selector_subset,
            "selected_subset_images_by_unit": previous_selector_images,
            "is_aut_invariant": is_invariant_subset(previous_selector_subset, action_table),
            "classification": "symmetry-breaking selector premise, not strict Aut-invariant closure",
        },
        "proof_reading": {
            "finite_theorem": "For the Hebbian self-maximizing teacher-class set {C,P,D}, every Aut(Z_12)-invariant subset is a union of the orbits {C,D} and {P}.",
            "negative_result": "No Aut-invariant selector can choose the singleton d5 class from the Aut-equivalent unit pair {contiguous,d5}.",
            "relation_to_previous_probe": "The previous k=5-over-k=1 label preference is now certified as Aut-breaking rather than merely suspicious.",
            "remaining_gap": "A strict internal orientation/generator/label-order source is still required before d5 can be selected strictly from the unit pair.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later supply an internal orientation/generator record if a strict source is proved.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply that record here.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives a strict orientation/generator/label-order source that prefers k=5 over k=1.",
            "No Aut(Z_12)-invariant selector on this finite candidate set can pick d5 alone from the unit pair.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Either derive an internal Aut-breaking orientation/generator source, or keep any k=5-over-k=1 selector explicitly non-strict.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian Aut-invariant selector no-go probe\n\n"
        "Status: finite Aut-invariant selector no-go for the d5 unit pair.\n\n"
        f"- Supports scanned: `{support_count}` five-active-node states on `Z_12`.\n"
        f"- Aut action table: `{action_table}`.\n"
        f"- Orbit partition: `{orbits}`.\n"
        f"- Invariant subsets: `{invariant_subsets}`.\n"
        f"- d5 singleton Aut-invariant: `{is_invariant_subset(d5_singleton, action_table)}`.\n"
        f"- Unit pair Aut-invariant: `{is_invariant_subset(unit_pair, action_table)}`.\n"
        f"- Previous selector classification: `symmetry-breaking selector premise`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: strict Aut-invariant logic can keep {contiguous,d5} together but cannot choose d5 alone.\n"
        "- No false pass: no strict d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
