#!/usr/bin/env python3
"""Scratch probe: cycle-metric subgroup audit for the Hebbian d5 selector.

The previous minimal Aut-breaking probe proved that the candidate d5 selector
needs the symmetry reduction Aut(Z_12) -> {1,11}.  This probe tests one concrete
finite candidate for such a reduction:

    add a cyclic nearest-neighbor metric/locality record on Z_12.

Finite result:
- Multiplicative units preserving the undirected cycle edge distance d=1 are
  exactly {1,11}.
- The same subgroup is the minimal nontrivial subgroup previously required for
  the non-Nyquist unit highest-label selector to choose d5 over contiguous.
- Therefore a cycle-metric/locality record would be sufficient to make the d5
  selector invariant inside the reduced symmetry.
- But this probe does not derive that metric/locality record from strict
  nadsoliton geometry; it only verifies that if supplied, it has exactly the
  required subgroup effect.

No false pass: this is a conditional source-candidate audit, not strict selector
closure and not a QW-2191 discharge.
"""
from __future__ import annotations

import json
from itertools import combinations
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_strict_alpha_hebbian_cycle_metric_subgroup_selector_audit_report.json"
OUT_MD = HERE / "bridge_strict_alpha_hebbian_cycle_metric_subgroup_selector_audit_report.md"

N = 12
ACTIVE_COUNT = 5
AUT_UNITS = [1, 5, 7, 11]
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
TARGET_D5_CLASS = "fifth_step_d5_step_5_7"
CONTIGUOUS_CLASS = "contiguous_step_1_11"
PARITY_CLASS = "parity_minus_one_step_2_10"
UNIT_PAIR = [CONTIGUOUS_CLASS, TARGET_D5_CLASS]
CLASS_LABELS = {
    CONTIGUOUS_CLASS: [1],
    PARITY_CLASS: [6],
    TARGET_D5_CLASS: [5],
}
CLASSES = {
    CONTIGUOUS_CLASS: [1, 11],
    PARITY_CLASS: [2, 10],
    TARGET_D5_CLASS: [5, 7],
}
MINIMAL_REQUIRED_SUBGROUP = [1, 11]


def folded_distance(distance: int) -> int:
    residue = distance % N
    return min(residue, N - residue)


def unit_preserves_distance(unit: int, distance: int) -> bool:
    return folded_distance(unit * distance) == folded_distance(distance)


def stabilizer_for_distance(distance: int) -> list[int]:
    return [unit for unit in AUT_UNITS if unit_preserves_distance(unit, distance)]


def cycle_edges(distance: int) -> list[list[int]]:
    edges = {
        tuple(sorted((node, (node + distance) % N)))
        for node in range(N)
    }
    return [list(edge) for edge in sorted(edges)]


def image_edges(edges: list[list[int]], unit: int) -> list[list[int]]:
    return sorted([sorted([(unit * left) % N, (unit * right) % N]) for left, right in edges])


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


def subset_is_invariant(subset: list[str], subgroup: list[int]) -> bool:
    normalized = sorted(subset)
    for unit in subgroup:
        image = sorted({unit_action_on_class(class_name, unit) for class_name in subset})
        if image != normalized:
            return False
    return True


def highest_label_winner(classes: list[str]) -> str:
    return max(classes, key=lambda name: max(CLASS_LABELS[name]))


def selector_certificate_for_subgroup(subgroup: list[int]) -> dict[str, Any]:
    winner = highest_label_winner(UNIT_PAIR)
    label_orbits = {
        class_name: sorted({unit_action_on_label(CLASS_LABELS[class_name][0], unit) for unit in subgroup})
        for class_name in UNIT_PAIR
    }
    return {
        "subgroup": subgroup,
        "unit_pair_invariant": subset_is_invariant(UNIT_PAIR, subgroup),
        "d5_singleton_invariant": subset_is_invariant([TARGET_D5_CLASS], subgroup),
        "contiguous_singleton_invariant": subset_is_invariant([CONTIGUOUS_CLASS], subgroup),
        "label_orbits_under_subgroup": label_orbits,
        "highest_label_winner_on_unit_pair": winner,
        "selects_d5_in_reduced_symmetry": winner == TARGET_D5_CLASS
        and subset_is_invariant([TARGET_D5_CLASS], subgroup)
        and subset_is_invariant([CONTIGUOUS_CLASS], subgroup),
    }


def main() -> None:
    support_count = len(list(combinations(range(N), ACTIVE_COUNT)))
    distance_stabilizers = {str(distance): stabilizer_for_distance(distance) for distance in range(1, N // 2 + 1)}
    distance_rows = []
    for distance in range(1, N // 2 + 1):
        edges = cycle_edges(distance)
        stabilizer = distance_stabilizers[str(distance)]
        distance_rows.append(
            {
                "folded_distance": distance,
                "edge_count": len(edges),
                "stabilizer_units": stabilizer,
                "equals_minimal_required_subgroup": stabilizer == MINIMAL_REQUIRED_SUBGROUP,
                "unit_edge_preservation": {
                    str(unit): image_edges(edges, unit) == edges for unit in AUT_UNITS
                },
            }
        )

    cycle_metric_subgroup = distance_stabilizers["1"]
    d5_step_metric_subgroup = distance_stabilizers["5"]
    selector_under_cycle_metric = selector_certificate_for_subgroup(cycle_metric_subgroup)
    full_aut_selector = selector_certificate_for_subgroup(AUT_UNITS)

    report: dict[str, Any] = {
        "result_kind": "SCRATCH_STRICT_ALPHA_HEBBIAN_CYCLE_METRIC_SUBGROUP_SELECTOR_AUDIT_PROBE__NOT_A_THEOREM",
        "status": "cycle-metric-record-supplies-required-subgroup-conditionally-not-strict-source",
        "finite_model": {
            "ring": "Z_12",
            "active_count": ACTIVE_COUNT,
            "support_count": support_count,
            "automorphism_units": AUT_UNITS,
            "minimal_required_subgroup_from_previous_probe": MINIMAL_REQUIRED_SUBGROUP,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
        },
        "distance_stabilizer_audit": {
            "distance_stabilizers": distance_stabilizers,
            "rows": distance_rows,
            "nearest_neighbor_distance_1_stabilizer": cycle_metric_subgroup,
            "d5_step_distance_5_stabilizer": d5_step_metric_subgroup,
            "distance_1_matches_minimal_required_subgroup": cycle_metric_subgroup == MINIMAL_REQUIRED_SUBGROUP,
            "distance_5_matches_minimal_required_subgroup": d5_step_metric_subgroup == MINIMAL_REQUIRED_SUBGROUP,
        },
        "selector_replay": {
            "full_aut_selector": full_aut_selector,
            "cycle_metric_reduced_selector": selector_under_cycle_metric,
            "cycle_metric_makes_d5_selector_invariant": selector_under_cycle_metric["selects_d5_in_reduced_symmetry"],
            "full_aut_makes_d5_selector_invariant": full_aut_selector["selects_d5_in_reduced_symmetry"],
        },
        "candidate_source_interpretation": {
            "finite_gain": "A nearest-neighbor cycle metric/locality record has exactly the subgroup effect required by the previous minimal Aut-breaking audit: it reduces multiplicative units to {1,11}.",
            "conditional_positive_result": "If that metric/locality record is admitted as an internal nadsoliton self-record, the non-Nyquist unit highest-label rule selects d5 invariantly inside the reduced symmetry.",
            "honest_limit": "The metric/locality record itself is not derived here from strict nadsoliton geometry, so this remains a conditional selector-source candidate.",
            "relation_to_previous_probe": "This supplies a concrete finite candidate for the previously missing subgroup {1,11}, but does not promote it to strict-core source data.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself may later be shown to carry a finite locality/metric self-record.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced to supply the metric record here.",
            "preferred_order_preserved": "nadsoliton -> light -> matter -> emergent observer",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is asserted or used.",
            "No legacy physical-role claim is transferred onto K_strict_gate.",
            "No theorem derives the d5 teacher/self-record trace from strict nadsoliton geometry.",
            "No theorem derives the cycle metric/locality record from strict geometry.",
            "This is a conditional subgroup-source audit, not strict-core selector closure.",
            "No endpoint, arrow orientation, ledger selector, positive lambda action, or fifth-mode source theorem is claimed.",
            "No QW-2191 discharge and no strict-core selector closure are claimed.",
            "No ToE closure is claimed.",
        ],
        "next_honest_step": "Try to derive the cycle metric/locality self-record internally; otherwise keep the {1,11} subgroup reduction explicitly conditional.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha Hebbian cycle-metric subgroup selector audit probe\n\n"
        "Status: cycle metric supplies the required {1,11} subgroup conditionally; not a strict source theorem.\n\n"
        f"- Supports scanned: `{support_count}` five-active-node states on `Z_12`.\n"
        f"- Distance stabilizers: `{distance_stabilizers}`.\n"
        f"- Distance `1` stabilizer equals required subgroup: `{cycle_metric_subgroup == MINIMAL_REQUIRED_SUBGROUP}`.\n"
        f"- Distance `5` stabilizer equals required subgroup: `{d5_step_metric_subgroup == MINIMAL_REQUIRED_SUBGROUP}`.\n"
        f"- Full Aut d5 selector invariant: `{full_aut_selector['selects_d5_in_reduced_symmetry']}`.\n"
        f"- Cycle-metric reduced d5 selector invariant: `{selector_under_cycle_metric['selects_d5_in_reduced_symmetry']}`.\n"
        f"- Target replay: `q^5={TARGET_Q_POWER}`, eta `{TARGET_ETA}`.\n"
        "- Honest read: a cycle metric/locality record would supply the needed subgroup, but that record is not derived here.\n"
        "- No false pass: no strict d5-source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
