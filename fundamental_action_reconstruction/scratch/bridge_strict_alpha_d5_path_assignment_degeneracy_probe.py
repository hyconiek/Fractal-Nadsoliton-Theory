#!/usr/bin/env python3
"""Scratch probe: distance-5 path assignment degeneracy after support selection.

The distance-5 band-pass probe gives a clean finite support selector: maximizing
h_5(S) forces the 5-node support into the fifth-step orbit.  This packet takes
the next honest step and asks what remains after that support is selected.

On any selected support, the internal distance-5 edges form a path of length 4.
For the balanced ledger (2,2,2,1,1), the simplest path-smooth assignment action

    E_path(a)=sum_{r=0..3} (a_{r+1}-a_r)^2

is exactly the number of 1/2 transitions along the ordered distance-5 path.  It
selects only the two endpoint-block assignments

    22211  and  11222

(up to the chosen path orientation), not a unique source phase, orientation, or
assignment.  Thus the support selector plus convex ledger selector still leaves
a real finite degeneracy unless an additional boundary/source-phase/orientation
premise is supplied.
"""
from __future__ import annotations

import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BANDPASS = HERE / "bridge_strict_alpha_bandpass_d5_support_selector_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_d5_path_assignment_degeneracy_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d5_path_assignment_degeneracy_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
BALANCED_LEDGER = (2, 2, 2, 1, 1)
REFERENCE_ORDERED_SUPPORT = (0, 7, 2, 9, 4)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def cyclic_distance(left: int, right: int) -> int:
    raw = abs(left - right) % Z12_NODE_COUNT
    return min(raw, Z12_NODE_COUNT - raw)


def ordered_d5_support(start: int, orientation: int) -> tuple[int, ...]:
    if orientation not in {-1, 1}:
        raise ValueError("orientation must be +/-1")
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((start + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def canonical_support(ordered_support: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(ordered_support))


def unique_balanced_assignments() -> list[tuple[int, ...]]:
    return sorted(set(itertools.permutations(BALANCED_LEDGER)))


def path_transition_count(assignment: tuple[int, ...]) -> int:
    return sum(1 for left, right in zip(assignment, assignment[1:]) if left != right)


def path_variation_energy(assignment: tuple[int, ...]) -> int:
    return sum((right - left) ** 2 for left, right in zip(assignment, assignment[1:]))


def high_node_set(ordered_support: tuple[int, ...], assignment: tuple[int, ...]) -> list[int]:
    return sorted(node for node, value in zip(ordered_support, assignment) if value == 2)


def assignment_rows_for_ordered_support(ordered_support: tuple[int, ...]) -> list[dict[str, Any]]:
    rows = []
    for assignment in unique_balanced_assignments():
        rows.append(
            {
                "ordered_support": list(ordered_support),
                "assignment": list(assignment),
                "path_variation_energy": path_variation_energy(assignment),
                "transition_count": path_transition_count(assignment),
                "high_value_nodes": high_node_set(ordered_support, assignment),
            }
        )
    return rows


def orbit_assignment_scan() -> dict[str, Any]:
    rows = []
    seen_ordered = set()
    for start in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            ordered_support = ordered_d5_support(start, orientation)
            if ordered_support in seen_ordered:
                continue
            seen_ordered.add(ordered_support)
            for row in assignment_rows_for_ordered_support(ordered_support):
                row["start"] = start
                row["orientation"] = orientation
                row["canonical_support"] = list(canonical_support(ordered_support))
                rows.append(row)
    min_energy = min(row["path_variation_energy"] for row in rows)
    max_energy = max(row["path_variation_energy"] for row in rows)
    minimizers = [row for row in rows if row["path_variation_energy"] == min_energy]
    maximizers = [row for row in rows if row["path_variation_energy"] == max_energy]
    unique_min_high_sets = sorted({tuple(row["high_value_nodes"]) for row in minimizers})
    return {
        "ordered_support_count": len(seen_ordered),
        "assignments_per_ordered_support": len(unique_balanced_assignments()),
        "total_assignment_rows": len(rows),
        "min_path_variation_energy": min_energy,
        "max_path_variation_energy": max_energy,
        "minimizer_count": len(minimizers),
        "maximizer_count": len(maximizers),
        "unique_minimizer_high_node_sets_count": len(unique_min_high_sets),
        "unique_minimizer_high_node_sets": [list(nodes) for nodes in unique_min_high_sets],
        "reference_ordered_support_rows": assignment_rows_for_ordered_support(REFERENCE_ORDERED_SUPPORT),
        "reference_minimizers": [
            row for row in assignment_rows_for_ordered_support(REFERENCE_ORDERED_SUPPORT)
            if row["path_variation_energy"] == min_energy
        ],
        "reference_maximizers": [
            row for row in assignment_rows_for_ordered_support(REFERENCE_ORDERED_SUPPORT)
            if row["path_variation_energy"] == max_energy
        ],
    }


def transition_histogram() -> dict[str, int]:
    counts: dict[str, int] = {}
    for assignment in unique_balanced_assignments():
        key = str(path_transition_count(assignment))
        counts[key] = counts.get(key, 0) + 1
    return counts


def main() -> None:
    bandpass_report = load_json(BANDPASS)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    scan = orbit_assignment_scan()
    reference_rows = scan["reference_ordered_support_rows"]

    report = {
        "status": "OPEN_STRICT_ALPHA_D5_PATH_ASSIGNMENT_DEGENERACY_NO_ASSIGNMENT_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_D5_PATH_ASSIGNMENT_DEGENERACY_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "bandpass_d5_support_selector_certificate": str(BANDPASS.relative_to(ROOT)),
        },
        "previous_bandpass_replay": {
            "result_kind": bandpass_report["result_kind"],
            "maximizers_equal_fifth_step_orbit": bandpass_report["bandpass_selector_scan"]["maximizers_equal_fifth_step_orbit"],
            "max_h5": bandpass_report["bandpass_selector_scan"]["max_h5"],
            "maximizer_count": bandpass_report["bandpass_selector_scan"]["maximizer_count"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "balanced_ledger": list(BALANCED_LEDGER),
        },
        "path_assignment_model": {
            "ordered_support_example": list(REFERENCE_ORDERED_SUPPORT),
            "all_adjacent_distances_are_5": all(cyclic_distance(left, right) == DISTANCE_SELECTED for left, right in zip(REFERENCE_ORDERED_SUPPORT, REFERENCE_ORDERED_SUPPORT[1:])),
            "energy": "E_path(a)=sum_{r=0..3}(a_{r+1}-a_r)^2",
            "energy_equals_transition_count_for_values_1_and_2": True,
            "transition_histogram_for_one_ordered_support": transition_histogram(),
            "reference_rows": reference_rows,
        },
        "orbit_assignment_scan": scan,
        "conditional_assignment_schema": {
            "support_premise": "distance-5 band-pass support selector has already chosen the fifth-step support orbit",
            "ledger_premise": "convex/Schur-convex ledger selector has already chosen the multiset (2,2,2,1,1)",
            "path_smooth_assignment_premise": "minimize E_path along the induced distance-5 path",
            "conditional_output": "assignments reduce to endpoint-block patterns 22211 or 11222 on each oriented path",
            "remaining_obstruction": "endpoint side, translate/source phase, path orientation, and physical origin of E_path remain unselected",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": True,
            "content": "Path-smooth assignment is a useful finite refinement after d5 support selection, but it leaves a controlled degeneracy rather than a unique selector.",
            "why_this_is_more_proof_like": "The assignment problem is reduced to a 10-row exact transition-count table per ordered support and a 240-row orbit scan.",
            "why_this_is_not_enough": "No strict theorem derives the path-smooth assignment action, endpoint choice, translate/source phase, orientation, or final balanced-ledger assignment.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The path-smooth assignment selector is a conditional model-internal premise, not a derived strict nadsoliton action theorem.",
            "Distance-5 support plus path-smooth assignment still leaves endpoint/translate/orientation degeneracy.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, and assignment premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Look for a strict boundary/source-phase/orientation term that chooses one endpoint block and one translate; otherwise keep assignment as an explicit non-strict premise.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha d5 path assignment degeneracy probe\n\n"
        "Status: assignment degeneracy certificate after distance-5 support selection; no assignment selector theorem.\n\n"
        f"- Reference ordered support `{list(REFERENCE_ORDERED_SUPPORT)}` has all adjacent cyclic distances equal to `{DISTANCE_SELECTED}`.\n"
        f"- One-support assignment table: `{len(unique_balanced_assignments())}` balanced assignments with transition histogram `{transition_histogram()}`.\n"
        f"- Orbit scan: `{scan['ordered_support_count']}` ordered supports and `{scan['total_assignment_rows']}` assignment rows.\n"
        f"- Path-smooth minimizers: energy `{scan['min_path_variation_energy']}`, count `{scan['minimizer_count']}`, unique high-node sets `{scan['unique_minimizer_high_node_sets_count']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: path smoothing narrows assignments to endpoint blocks, but does not choose endpoint, translate, orientation, or source phase.\n"
        "- No false pass: no strict assignment theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
