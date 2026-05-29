#!/usr/bin/env python3
"""Scratch probe: anchored endpoint assignment selector after d5 path degeneracy.

The previous d5 path assignment probe showed a controlled but real degeneracy:
after distance-5 support selection and the balanced ledger (2,2,2,1,1), path
smoothing leaves the two endpoint-block assignments 22211 and 11222 on every
oriented distance-5 path.

This packet checks the smallest finite selector that would remove that last
assignment ambiguity once a source endpoint and orientation have already been
supplied.  For an ordered path a=(a_0,...,a_4), define the lexicographic selector

    first minimize  E_path(a)=sum_r (a_{r+1}-a_r)^2,
    then minimize   M_source(a)=sum_r r*a_r.

The second term puts larger ledger entries closer to the source endpoint r=0.
For the balanced ledger this uniquely selects 22211 on each ordered path.  The
opposite endpoint convention would select 11222.  Therefore source endpoint +
orientation + path-smooth endpoint bias can close the finite assignment problem,
but all three are still extra selector premises; this is not a strict-core
source theorem and does not discharge QW-2191.
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
PATH_DEGENERACY = HERE / "bridge_strict_alpha_d5_path_assignment_degeneracy_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_d5_anchored_endpoint_assignment_selector_report.json"
OUT_MD = HERE / "bridge_strict_alpha_d5_anchored_endpoint_assignment_selector_report.md"

Z12_NODE_COUNT = 12
SUPPORT_SIZE = 5
DISTANCE_SELECTED = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
BALANCED_LEDGER = (2, 2, 2, 1, 1)
REFERENCE_SOURCE = 0
REFERENCE_ORIENTATION = -1


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def cyclic_distance(left: int, right: int) -> int:
    raw = abs(left - right) % Z12_NODE_COUNT
    return min(raw, Z12_NODE_COUNT - raw)


def ordered_d5_support(source: int, orientation: int) -> tuple[int, ...]:
    if orientation not in {-1, 1}:
        raise ValueError("orientation must be +/-1")
    step = (orientation * DISTANCE_SELECTED) % Z12_NODE_COUNT
    return tuple((source + index * step) % Z12_NODE_COUNT for index in range(SUPPORT_SIZE))


def unique_balanced_assignments() -> list[tuple[int, ...]]:
    return sorted(set(itertools.permutations(BALANCED_LEDGER)))


def path_energy(assignment: tuple[int, ...]) -> int:
    return sum((right - left) ** 2 for left, right in zip(assignment, assignment[1:]))


def source_moment(assignment: tuple[int, ...]) -> int:
    return sum(index * value for index, value in enumerate(assignment))


def reverse_source_moment(assignment: tuple[int, ...]) -> int:
    return sum((SUPPORT_SIZE - 1 - index) * value for index, value in enumerate(assignment))


def high_nodes(ordered_support: tuple[int, ...], assignment: tuple[int, ...]) -> list[int]:
    return sorted(node for node, value in zip(ordered_support, assignment) if value == 2)


def assignment_row(ordered_support: tuple[int, ...], assignment: tuple[int, ...]) -> dict[str, Any]:
    return {
        "ordered_support": list(ordered_support),
        "assignment": list(assignment),
        "path_energy": path_energy(assignment),
        "source_moment": source_moment(assignment),
        "reverse_source_moment": reverse_source_moment(assignment),
        "lex_forward_key": [path_energy(assignment), source_moment(assignment)],
        "lex_reverse_key": [path_energy(assignment), reverse_source_moment(assignment)],
        "high_value_nodes": high_nodes(ordered_support, assignment),
    }


def rows_for_ordered_support(ordered_support: tuple[int, ...]) -> list[dict[str, Any]]:
    return [assignment_row(ordered_support, assignment) for assignment in unique_balanced_assignments()]


def selector_for_ordered_support(ordered_support: tuple[int, ...]) -> dict[str, Any]:
    rows = rows_for_ordered_support(ordered_support)
    forward = min(rows, key=lambda row: tuple(row["lex_forward_key"]))
    reverse = min(rows, key=lambda row: tuple(row["lex_reverse_key"]))
    path_min_energy = min(row["path_energy"] for row in rows)
    path_minimizers = [row for row in rows if row["path_energy"] == path_min_energy]
    return {
        "ordered_support": list(ordered_support),
        "all_adjacent_distances_are_5": all(cyclic_distance(left, right) == DISTANCE_SELECTED for left, right in zip(ordered_support, ordered_support[1:])),
        "row_count": len(rows),
        "rows": rows,
        "path_min_energy": path_min_energy,
        "path_minimizer_count_before_endpoint_bias": len(path_minimizers),
        "path_minimizers_before_endpoint_bias": path_minimizers,
        "forward_endpoint_selected_row": forward,
        "reverse_endpoint_selected_row": reverse,
        "forward_assignment_unique": sum(row["lex_forward_key"] == forward["lex_forward_key"] for row in rows) == 1,
        "reverse_assignment_unique": sum(row["lex_reverse_key"] == reverse["lex_reverse_key"] for row in rows) == 1,
    }


def orbit_scan() -> dict[str, Any]:
    rows = []
    selected_forward = []
    selected_reverse = []
    for source in range(Z12_NODE_COUNT):
        for orientation in (-1, 1):
            ordered_support = ordered_d5_support(source, orientation)
            result = selector_for_ordered_support(ordered_support)
            result["source"] = source
            result["orientation"] = orientation
            rows.append(result)
            selected_forward.append(tuple(result["forward_endpoint_selected_row"]["assignment"]))
            selected_reverse.append(tuple(result["reverse_endpoint_selected_row"]["assignment"]))
    return {
        "source_count": Z12_NODE_COUNT,
        "orientation_count_per_source": 2,
        "ordered_support_count": len(rows),
        "assignments_per_ordered_support": len(unique_balanced_assignments()),
        "all_forward_assignments_unique": all(row["forward_assignment_unique"] for row in rows),
        "all_reverse_assignments_unique": all(row["reverse_assignment_unique"] for row in rows),
        "forward_selected_assignment_set": [list(item) for item in sorted(set(selected_forward))],
        "reverse_selected_assignment_set": [list(item) for item in sorted(set(selected_reverse))],
        "rows_summary": [
            {
                "source": row["source"],
                "orientation": row["orientation"],
                "ordered_support": row["ordered_support"],
                "forward_assignment": row["forward_endpoint_selected_row"]["assignment"],
                "reverse_assignment": row["reverse_endpoint_selected_row"]["assignment"],
                "forward_high_value_nodes": row["forward_endpoint_selected_row"]["high_value_nodes"],
                "reverse_high_value_nodes": row["reverse_endpoint_selected_row"]["high_value_nodes"],
            }
            for row in rows
        ],
    }


def main() -> None:
    path_report = load_json(PATH_DEGENERACY)
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** SUPPORT_SIZE)
    reference_support = ordered_d5_support(REFERENCE_SOURCE, REFERENCE_ORIENTATION)
    reference = selector_for_ordered_support(reference_support)
    scan = orbit_scan()

    report = {
        "status": "OPEN_STRICT_ALPHA_D5_ANCHORED_ENDPOINT_ASSIGNMENT_SELECTOR_NO_STRICT_SOURCE_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_D5_ANCHORED_ENDPOINT_ASSIGNMENT_SELECTOR_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "d5_path_assignment_degeneracy": str(PATH_DEGENERACY.relative_to(ROOT)),
        },
        "previous_path_degeneracy_replay": {
            "result_kind": path_report["result_kind"],
            "minimizer_count": path_report["orbit_assignment_scan"]["minimizer_count"],
            "unique_minimizer_high_node_sets_count": path_report["orbit_assignment_scan"]["unique_minimizer_high_node_sets_count"],
        },
        "target_identity_replay": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, SUPPORT_SIZE) - STRICT_TARGET_ETA,
            "balanced_ledger": list(BALANCED_LEDGER),
        },
        "anchored_selector_definition": {
            "support_rule": "ordered support is (source + r*orientation*5) mod 12 for r=0..4",
            "path_energy": "E_path(a)=sum_r(a_{r+1}-a_r)^2",
            "source_moment": "M_source(a)=sum_r r*a_r",
            "forward_endpoint_selector": "lexicographically minimize (E_path, M_source)",
            "reverse_endpoint_selector": "lexicographically minimize (E_path, M_reverse) with M_reverse=sum_r(4-r)*a_r",
        },
        "reference_source_orientation_certificate": {
            "source": REFERENCE_SOURCE,
            "orientation": REFERENCE_ORIENTATION,
            **reference,
        },
        "orbit_scan": scan,
        "conditional_selector_schema": {
            "support_premise": "distance-5 band-pass support selector plus explicit source and orientation chooses an ordered support",
            "ledger_premise": "convex/Schur-convex ledger selector supplies the multiset (2,2,2,1,1)",
            "endpoint_bias_premise": "lexicographic endpoint bias chooses high values nearest the declared source endpoint",
            "conditional_output": "for each source and orientation, the forward endpoint selector uniquely chooses assignment 22211",
            "remaining_obstruction": "the strict origin of source, orientation, and endpoint-bias action is not derived",
        },
        "candidate_interpretation": {
            "supported_by_this_probe": bool(scan["all_forward_assignments_unique"] and scan["forward_selected_assignment_set"] == [[2, 2, 2, 1, 1]]),
            "content": "An explicit source endpoint and orientation remove the final finite assignment degeneracy for the d5 path model.",
            "why_this_is_more_proof_like": "The selector is exact and lexicographic: it first preserves path smoothing, then uniquely minimizes a source moment across every ordered support in the 24-case orbit scan.",
            "why_this_is_not_enough": "The source endpoint, orientation, and endpoint-bias action remain premises; no strict theorem exports them from nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "The anchored endpoint assignment selector is conditional on explicit source, orientation, and endpoint-bias premises.",
            "No theorem derives the endpoint-bias/source-moment action from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, or binary-rescale quotient.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged; the source/orientation premises are exactly extra selector data.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, support premise, ledger selector, source/orientation premise, and assignment premise.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Try to derive the source/orientation/endpoint-bias terms from strict nadsoliton geometry; otherwise mark them as non-strict selector axioms.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha d5 anchored endpoint assignment selector probe\n\n"
        "Status: anchored endpoint assignment selector; no strict source theorem.\n\n"
        f"- Reference ordered support `{list(reference_support)}` from source `{REFERENCE_SOURCE}`, orientation `{REFERENCE_ORIENTATION}`.\n"
        f"- Before endpoint bias: path minimizer count `{reference['path_minimizer_count_before_endpoint_bias']}`.\n"
        f"- Forward endpoint selector chooses `{reference['forward_endpoint_selected_row']['assignment']}` uniquely: `{reference['forward_assignment_unique']}`.\n"
        f"- Orbit scan: `{scan['ordered_support_count']}` ordered supports; all forward assignments unique: `{scan['all_forward_assignments_unique']}`.\n"
        f"- Target replay: `q^5={product.numerator}/{product.denominator}`, eta residual `{eta_from_product(product, SUPPORT_SIZE)-STRICT_TARGET_ETA:.3e}`.\n"
        "- Honest read: source+orientation+endpoint bias closes the finite assignment table, but those are extra selector premises.\n"
        "- No false pass: no strict source/orientation/action theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
