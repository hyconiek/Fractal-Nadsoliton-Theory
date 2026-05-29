#!/usr/bin/env python3
"""Scratch probe: minimal Aut-breaking needed for the Hebbian d5 selector.

The previous Aut-invariant no-go proved that full Aut(Z_12) cannot select the
singleton d5 teacher class from the unit pair {contiguous,d5}.  This probe asks
the next finite question:

    What exact subgroup reduction is sufficient for the previous numeric
    highest-label rule (k=5 over k=1) to become invariant?

Finite result:
- Aut(Z_12) = {1,5,7,11} has five subgroups.
- The full group and the two swap subgroups {1,5}, {1,7} do not preserve the d5
  singleton.
- The subgroup {1,11} and the trivial subgroup preserve both unit singletons.
- Under {1,11}, the learned labels k=1 and k=5 are stable, so the non-Nyquist
  unit highest-label rule selects d5 invariantly inside that reduced symmetry.

No false pass: this only computes the minimal symmetry reduction required by the
candidate selector.  It does not derive the subgroup {1,11}, an orientation
record, or a generator/label-order source from strict nadsoliton geometry.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_minimal_aut_breaking_selector_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_minimal_aut_breaking_selector_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
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
UNIT_PAIR = [CONTIGUOUS_CLASS, TARGET_D5_CLASS]


def mul_mod(left: int, right: int) -> int:
    return (left * right) % N


def is_subgroup(candidate: list[int]) -> bool:
    values = set(candidate)
    return 1 in values and all(mul_mod(a, b) in values for a in values for b in values)


def all_subgroups() -> list[list[int]]:
    subgroups: list[list[int]] = []
    for size in range(1, len(AUT_UNITS) + 1):
        for row in combinations(AUT_UNITS, size):
            candidate = sorted(row)
            if is_subgroup(candidate):
                subgroups.append(candidate)
    return sorted(subgroups, key=lambda row: (len(row), row))


def folded_step(step: int) -> list[int]:
    residue = step % N
    return sorted({residue, (-residue) % N})


def class_for_step(step: int) -> str:
    folded = folded_step(step)
    for name, steps in CLASSES.items():
        if sorted(value % N for value in steps) == folded:
            return name
    raise ValueError(f"No class for step {step}")


def unit_action_on_class(class_name: str, unit: int) -> str:
    return class_for_step(CLASSES[class_name][0] * unit)


def folded_mode(mode: int) -> int:
    residue = mode % N
    if residue == 0:
        return 0
    return min(residue, N - residue)


def unit_action_on_label(label: int, unit: int) -> int:
    return folded_mode(label * unit)


def orbit_under_subgroup(class_name: str, subgroup: list[int]) -> list[str]:
    return sorted({unit_action_on_class(class_name, unit) for unit in subgroup})


def label_orbit_under_subgroup(label: int, subgroup: list[int]) -> list[int]:
    return sorted({unit_action_on_label(label, unit) for unit in subgroup})


def subset_is_invariant(subset: list[str], subgroup: list[int]) -> bool:
    normalized = sorted(subset)
    for unit in subgroup:
        image = sorted({unit_action_on_class(class_name, unit) for class_name in subset})
        if image != normalized:
            return False
    return True


def highest_label_winner(classes: list[str]) -> str:
    return max(classes, key=lambda name: max(CLASS_LABELS[name]))


def subgroup_certificate(subgroup: list[int]) -> dict[str, Any]:
    d5_orbit = orbit_under_subgroup(TARGET_D5_CLASS, subgroup)
    contiguous_orbit = orbit_under_subgroup(CONTIGUOUS_CLASS, subgroup)
    unit_pair_invariant = subset_is_invariant(UNIT_PAIR, subgroup)
    d5_singleton_invariant = subset_is_invariant([TARGET_D5_CLASS], subgroup)
    contiguous_singleton_invariant = subset_is_invariant([CONTIGUOUS_CLASS], subgroup)
    label_orbits = {
        name: label_orbit_under_subgroup(CLASS_LABELS[name][0], subgroup)
        for name in UNIT_PAIR
    }
    winner = highest_label_winner(UNIT_PAIR)
    winner_singleton_invariant = subset_is_invariant([winner], subgroup)
    return {
        "subgroup": subgroup,
        "order": len(subgroup),
        "d5_orbit": d5_orbit,
        "contiguous_orbit": contiguous_orbit,
        "unit_pair_invariant": unit_pair_invariant,
        "d5_singleton_invariant": d5_singleton_invariant,
        "contiguous_singleton_invariant": contiguous_singleton_invariant,
        "unit_label_orbits": label_orbits,
        "highest_label_winner_on_unit_pair": winner,
        "highest_label_winner_singleton_invariant": winner_singleton_invariant,
        "highest_label_rule_is_well_defined_in_reduced_symmetry": unit_pair_invariant
        and d5_singleton_invariant
        and contiguous_singleton_invariant
        and winner == TARGET_D5_CLASS
        and winner_singleton_invariant,
    }


def main() -> None:
    support_count = len(list(combinations(range(N), ACTIVE_COUNT)))
    subgroups = all_subgroups()
    subgroup_rows = [subgroup_certificate(subgroup) for subgroup in subgroups]
    sufficient_rows = [row for row in subgroup_rows if row["highest_label_rule_is_well_defined_in_reduced_symmetry"]]
    nontrivial_sufficient_rows = [row for row in sufficient_rows if row["order"] > 1]
    minimal_nontrivial = min(nontrivial_sufficient_rows, key=lambda row: row["order"], default=None)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_MINIMAL_AUT_BREAKING_SELECTOR_PROBE__NOT_A_THEOREM",
        "status": "minimal-subgroup-reduction-for-d5-label-selector-not-strict-source",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": support_count,
            "automorphism_units": AUT_UNITS,
            "teacher_classes": CLASSES,
            "learned_labels": CLASS_LABELS,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "aut_subgroup_enumeration": {
            "subgroup_count": len(subgroups),
            "subgroups": subgroups,
            "rows": subgroup_rows,
        },
        "minimal_breaking_certificate": {
            "full_aut_group": AUT_UNITS,
            "full_aut_d5_singleton_invariant": next(row for row in subgroup_rows if row["subgroup"] == AUT_UNITS)["d5_singleton_invariant"],
            "sufficient_subgroups": [row["subgroup"] for row in sufficient_rows],
            "nontrivial_sufficient_subgroups": [row["subgroup"] for row in nontrivial_sufficient_rows],
            "minimal_nontrivial_sufficient_subgroup": minimal_nontrivial["subgroup"] if minimal_nontrivial else None,
            "minimal_nontrivial_sufficient_subgroup_order": minimal_nontrivial["order"] if minimal_nontrivial else None,
            "excluded_swap_units": [unit for unit in AUT_UNITS if minimal_nontrivial and unit not in minimal_nontrivial["subgroup"]],
            "interpretation": "The candidate d5 label selector becomes invariant only after excluding the swap units 5 and 7, i.e. after reducing Aut(Z_12) to {1,11} or further to {1}.",
        },
        "selector_stack_status": {
            "candidate_rule": "non-Nyquist unit highest learned label chooses k=5 over k=1",
            "works_under_full_aut": False,
            "works_under_minimal_nontrivial_subgroup": minimal_nontrivial is not None,
            "selected_class_under_reduced_symmetry": TARGET_D5_CLASS if minimal_nontrivial else None,
            "still_requires_extra_source": "A strict reason for preserving only {1,11} rather than full Aut(Z_12) is not derived here.",
        },
        "proof_reading": {
            "finite_gain": "The exact subgroup reduction needed by the candidate d5 label selector is enumerated: {1,11} is the minimal nontrivial sufficient subgroup.",
            "negative_result": "Full Aut(Z_12), {1,5}, and {1,7} all fail because they swap d5 with contiguous.",
            "conditional_positive_result": "If a strict internal source reduces the relevant symmetry to {1,11}, then k=5 and k=1 labels are stable and the highest-label rule selects d5.",
            "remaining_gap": "This probe does not derive the {1,11} subgroup, orientation record, or generator/label-order source.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to carry an internal symmetry-reduction record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply this symmetry reduction.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives the subgroup {1,11}, orientation record, generator choice, or label-order source from strict geometry.",
            "This is a minimal symmetry-breaking requirement certificate, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive an internal reason for the {1,11} subgroup reduction; without it, keep the d5 label selector explicitly non-strict.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian minimal Aut-breaking selector probe\n\n"
        "Status: minimal subgroup reduction for d5 label selector; not a strict source theorem.\n\n"
        f"- Supports scanned: `{support_count}` five-active-node states on `Z_12`.\n"
        f"- Aut subgroups: `{subgroups}`.\n"
        f"- Sufficient subgroups for d5 highest-label selector: `{[row['subgroup'] for row in sufficient_rows]}`.\n"
        f"- Minimal nontrivial sufficient subgroup: `{minimal_nontrivial['subgroup'] if minimal_nontrivial else None}`.\n"
        f"- Excluded swap units: `{[unit for unit in AUT_UNITS if minimal_nontrivial and unit not in minimal_nontrivial['subgroup']]}`.\n"
        f"- Full Aut d5 singleton invariant: `{next(row for row in subgroup_rows if row['subgroup'] == AUT_UNITS)['d5_singleton_invariant']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: d5 selection needs at least the symmetry reduction Aut(Z_12) -> {1,11}; that reduction is not derived here.\n"
        "- No false pass: no strict d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
