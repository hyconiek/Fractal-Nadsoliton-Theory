#!/usr/bin/env python3
"""Scratch probe: phase-sign node-support interval-boundary certificate.

The previous GF(2) certificates reconstruct the audited node bits and edge bits
on the finite path 0--1--...--11.  This probe records the next finite
combinatorial fact: the 1-node support is the disjoint union of maximal path
intervals

    [2,5] union [8,9],

and the edge flip cochain is exactly the GF(2) boundary of those intervals.  In
matrix form, edge_bits = interval_boundary_matrix * all_one_interval_vector over
GF(2).  Because neither interval touches an endpoint, each component contributes
two boundary edges and the four-edge flip support is explained as 2 * number of
negative support components.

This is finite path support bookkeeping only.  It does not derive phase zeros,
omega/phi, damping, or transport from strict nadsoliton dynamics, and it does not
claim a bridge theorem, selector discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_node_support_interval_boundary_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
EDGE_SUPPORT_REPORT = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"
REDUCED_COBOUNDARY_REPORT = HERE / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_report.json"
CYCLE_CLOSURE_REPORT = HERE / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EXPECTED_INTERVALS = [{"start": 2, "end": 5}, {"start": 8, "end": 9}]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def gf2_rref(matrix: list[list[int]]) -> dict[str, Any]:
    work = [row[:] for row in matrix]
    row_count = len(work)
    col_count = len(work[0]) if work else 0
    pivot_row = 0
    pivots = []
    for col in range(col_count):
        pivot = None
        for candidate in range(pivot_row, row_count):
            if work[candidate][col]:
                pivot = candidate
                break
        if pivot is None:
            continue
        if pivot != pivot_row:
            work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        for target in range(row_count):
            if target != pivot_row and work[target][col]:
                work[target] = [a ^ b for a, b in zip(work[target], work[pivot_row])]
        pivots.append({"row": pivot_row, "pivot_col": col})
        pivot_row += 1
        if pivot_row == row_count:
            break
    return {
        "rank": len(pivots),
        "nullity": col_count - len(pivots),
        "pivot_rows": pivots,
        "rref_matrix": work,
    }


def support(bits: list[int], labels: list[str]) -> list[str]:
    return [label for label, bit in zip(labels, bits) if bit]


def maximal_one_intervals(node_bits: list[int]) -> list[dict[str, int]]:
    intervals = []
    cursor = 0
    while cursor < len(node_bits):
        if node_bits[cursor] == 0:
            cursor += 1
            continue
        start = cursor
        while cursor + 1 < len(node_bits) and node_bits[cursor + 1] == 1:
            cursor += 1
        intervals.append({"start": start, "end": cursor})
        cursor += 1
    return intervals


def interval_node_vector(interval: dict[str, int]) -> list[int]:
    return [1 if interval["start"] <= node <= interval["end"] else 0 for node in range(NODE_COUNT)]


def interval_boundary_edge_vector(interval: dict[str, int]) -> list[int]:
    vector = [0] * EDGE_COUNT
    start = interval["start"]
    end = interval["end"]
    if start > 0:
        vector[start - 1] = 1
    if end < NODE_COUNT - 1:
        vector[end] ^= 1
    return vector


def interval_rows(intervals: list[dict[str, int]]) -> list[dict[str, Any]]:
    rows = []
    for index, interval in enumerate(intervals):
        node_vector = interval_node_vector(interval)
        boundary_vector = interval_boundary_edge_vector(interval)
        touches_left_endpoint = interval["start"] == 0
        touches_right_endpoint = interval["end"] == NODE_COUNT - 1
        rows.append({
            "interval_index": index,
            "interval": interval,
            "node_vector": node_vector,
            "node_labels": [node for node, bit in enumerate(node_vector) if bit],
            "boundary_edge_vector": boundary_vector,
            "boundary_edges": support(boundary_vector, EDGE_LABELS),
            "touches_left_endpoint": touches_left_endpoint,
            "touches_right_endpoint": touches_right_endpoint,
            "endpoint_touch_count": int(touches_left_endpoint) + int(touches_right_endpoint),
            "boundary_weight": sum(boundary_vector),
        })
    return rows


def transpose_columns_to_rows(columns: list[list[int]]) -> list[list[int]]:
    return [[columns[col][row] for col in range(len(columns))] for row in range(len(columns[0]))]


def edge_rows(boundary_matrix: list[list[int]], expected_edge_bits: list[int]) -> list[dict[str, Any]]:
    return [
        {
            "edge": edge,
            "interval_boundary_row": row,
            "computed_edge_bit_from_interval_boundaries": sum(row) % 2,
            "audited_edge_bit": expected,
            "is_flip_edge": expected == 1,
            "matches_audited_edge_bit": (sum(row) % 2) == expected,
        }
        for edge, row, expected in zip(EDGE_LABELS, boundary_matrix, expected_edge_bits)
    ]


def component_count_formula(intervals: list[dict[str, int]], edge_bits: list[int]) -> dict[str, Any]:
    endpoint_touch_count = sum(1 for interval in intervals if interval["start"] == 0) + sum(
        1 for interval in intervals if interval["end"] == NODE_COUNT - 1
    )
    return {
        "component_count": len(intervals),
        "endpoint_touch_count": endpoint_touch_count,
        "predicted_boundary_weight_2c_minus_endpoint_touches": 2 * len(intervals) - endpoint_touch_count,
        "actual_edge_bit_hamming_weight": sum(edge_bits),
        "formula_matches": 2 * len(intervals) - endpoint_touch_count == sum(edge_bits),
    }


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    edge_support = load_json(EDGE_SUPPORT_REPORT)
    reduced_coboundary = load_json(REDUCED_COBOUNDARY_REPORT)
    cycle_closure = load_json(CYCLE_CLOSURE_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    intervals = maximal_one_intervals(node_bits)
    rows = interval_rows(intervals)
    boundary_columns = [row["boundary_edge_vector"] for row in rows]
    boundary_matrix = transpose_columns_to_rows(boundary_columns)
    interval_vector = [1] * len(intervals)
    edge_bits_from_intervals = matvec_gf2(boundary_matrix, interval_vector)
    node_bits_from_intervals = [sum(row["node_vector"][node] for row in rows) % 2 for node in range(NODE_COUNT)]
    rank = gf2_rref(boundary_matrix)
    component_formula = component_count_formula(intervals, edge_bits)
    edge_boundary_rows = edge_rows(boundary_matrix, edge_bits)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_NODE_SUPPORT_INTERVAL_BOUNDARY_CERTIFICATE__GF2_COMPONENT_BOUNDARY",
        "status": "node-support-maximal-interval-boundary-recovers-edge-flip-cochain",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_edge_support_minimality_certificate": str(EDGE_SUPPORT_REPORT.relative_to(ROOT)),
            "phase_sign_reduced_coboundary_inverse_certificate": str(REDUCED_COBOUNDARY_REPORT.relative_to(ROOT)),
            "phase_sign_cycle_closure_obstruction_certificate": str(CYCLE_CLOSURE_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "negative node support",
                "node support interval",
                "support interval",
                "component boundary",
                "boundary of negative",
                "maximal negative",
                "interval boundary",
            ],
            "finding": "No existing strict-completion scratch report exported the maximal node-support interval decomposition with a GF(2) interval-boundary matrix before this file.",
        },
        "node_support_interval_definition": {
            "field": "GF(2)",
            "graph": "finite path 0--1--...--11",
            "node_bits": node_bits,
            "one_node_support": [node for node, bit in enumerate(node_bits) if bit],
            "maximal_one_intervals": intervals,
            "interval_multiplicity_vector": interval_vector,
            "interval_boundary_matrix_edges_by_intervals": boundary_matrix,
            "equation": "edge_bits = interval_boundary_matrix * interval_multiplicity_vector mod 2",
        },
        "interval_rows": rows,
        "edge_boundary_rows": edge_boundary_rows,
        "rank_certificate": {
            "interval_boundary_column_rank_over_gf2": rank["rank"],
            "interval_boundary_column_nullity_over_gf2": rank["nullity"],
            "full_column_rank_on_interval_boundary_subspace": rank["rank"] == len(intervals),
            "pivot_rows": rank["pivot_rows"],
            "rref_matrix": rank["rref_matrix"],
        },
        "component_boundary_formula": component_formula,
        "reconstruction_certificate": {
            "node_bits_from_interval_support_union": node_bits_from_intervals,
            "edge_bits_from_interval_boundaries": edge_bits_from_intervals,
            "flip_edges_from_interval_boundaries": support(edge_bits_from_intervals, EDGE_LABELS),
            "audited_edge_bits": edge_bits,
            "audited_flip_edges": support(edge_bits, EDGE_LABELS),
        },
        "node_support_interval_boundary_summary": {
            "maximal_interval_count": len(intervals),
            "maximal_intervals": intervals,
            "matches_expected_intervals": intervals == EXPECTED_INTERVALS,
            "interval_boundary_rank": rank["rank"],
            "interval_boundary_nullity": rank["nullity"],
            "interval_boundary_full_column_rank": rank["rank"] == len(intervals),
            "node_bits_recovered_from_interval_union": node_bits_from_intervals == node_bits,
            "edge_bits_recovered_from_interval_boundaries": edge_bits_from_intervals == edge_bits,
            "all_edge_boundary_rows_match": all(row["matches_audited_edge_bit"] for row in edge_boundary_rows),
            "boundary_weight_formula_matches": component_formula["formula_matches"],
            "matches_expected_node_bits": node_bits == EXPECTED_NODE_BITS,
            "matches_expected_edge_bits": edge_bits == EXPECTED_EDGE_BITS,
            "matches_expected_flip_edges": support(edge_bits_from_intervals, EDGE_LABELS) == EXPECTED_FLIP_EDGES,
            "matches_edge_support_minimality_flip_edges": support(edge_bits_from_intervals, EDGE_LABELS) == edge_support["edge_support_minimality_summary"]["derived_phase_sign_flip_edges"],
            "matches_reduced_coboundary_edge_bits": edge_bits_from_intervals == reduced_coboundary["reconstruction_certificate"]["edge_bits_from_R_tail_node_bits"],
            "contrasts_cycle_odd_closing_obstruction": cycle_closure["cycle_closure_summary"]["odd_closing_obstructed_by_cycle_parity"],
        },
        "blocker_context": {
            "resolved_locally": [
                "The negative/one node support is partitioned into maximal path intervals [2,5] and [8,9].",
                "The GF(2) boundary of those intervals is exactly the audited four-edge flip cochain.",
                "The boundary matrix has full interval-column rank 2.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "support_step": "The audited node-bit 1-support is exactly the disjoint maximal interval union [2,5] U [8,9].",
            "boundary_step": "The GF(2) interval-boundary matrix maps the all-one interval vector to edge bits [0,1,0,0,0,1,0,1,0,1,0].",
            "rank_step": "The two interval-boundary columns have rank 2 over GF(2), so neither component boundary is redundant in the interval-boundary subspace.",
            "component_count_step": "Because neither interval touches a path endpoint, boundary_weight = 2*component_count = 4, matching the four flip edges.",
            "chain_step": "The boundary result matches Z2, edge-support minimality, reduced-coboundary reconstruction, and the cycle-closure parity sanity check.",
            "theoretical_limit": "This is finite node-support boundary bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    summary = payload["node_support_interval_boundary_summary"]
    formula = payload["component_boundary_formula"]
    lines = [
        "# Phase-sign node-support interval-boundary certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate decomposes the audited 1-node support into maximal path",
        "intervals and verifies that their GF(2) boundary is exactly the edge-flip",
        "cochain.",
        "",
        "## Summary",
        "",
        f"- Maximal intervals: `{summary['maximal_intervals']}`",
        f"- Interval-boundary rank: `{summary['interval_boundary_rank']}`",
        f"- Interval-boundary nullity: `{summary['interval_boundary_nullity']}`",
        f"- Node bits recovered from interval union: `{summary['node_bits_recovered_from_interval_union']}`",
        f"- Edge bits recovered from interval boundaries: `{summary['edge_bits_recovered_from_interval_boundaries']}`",
        f"- Boundary weight formula: `2*{formula['component_count']} - {formula['endpoint_touch_count']} = {formula['actual_edge_bit_hamming_weight']}`",
        f"- Flip edges: `{payload['reconstruction_certificate']['flip_edges_from_interval_boundaries']}`",
        "",
        "## Hard limits",
        "",
    ]
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload["node_support_interval_boundary_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
