#!/usr/bin/env python3
"""Scratch probe: phase-sign component-quotient lift matrix certificate.

The component-quotient adjacency certificate collapses maximal constant-sign runs
of the path to five quotient vertices.  This probe records the finite GF(2)
linear maps that lift that quotient data back to the original 12-node path:

    node_bits = S * quotient_component_bits,
    edge_bits = E * quotient_edge_bits,
    B_path * S = E * B_quotient.

Here S is the 12x5 component-expansion matrix, B_quotient is the 4x5 quotient
path coboundary, E embeds quotient edges into the original 11 path edges, and
B_path is the original 11x12 path coboundary.  This is finite matrix
bookkeeping only; it does not derive phase zeros, omega/phi, damping, or
transport from strict nadsoliton dynamics, and it does not claim a kernel bridge,
selector discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
QUOTIENT_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
EXPECTED_COMPONENT_BITS = [0, 1, 0, 1, 0]
EXPECTED_QUOTIENT_EDGE_BITS = [1, 1, 1, 1]
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def matmul_gf2(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    if not left or not right:
        return []
    return [
        [sum(left[row][k] & right[k][col] for k in range(len(right))) % 2 for col in range(len(right[0]))]
        for row in range(len(left))
    ]


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


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [[matrix[row][col] for row in range(len(matrix))] for col in range(len(matrix[0]))]


def path_coboundary_matrix() -> list[list[int]]:
    rows = []
    for edge_index in range(EDGE_COUNT):
        row = [0] * NODE_COUNT
        row[edge_index] = 1
        row[edge_index + 1] = 1
        rows.append(row)
    return rows


def quotient_coboundary_matrix(component_count: int) -> list[list[int]]:
    rows = []
    for index in range(component_count - 1):
        row = [0] * component_count
        row[index] = 1
        row[index + 1] = 1
        rows.append(row)
    return rows


def component_expansion_matrix(components: list[dict[str, Any]]) -> list[list[int]]:
    matrix = [[0] * len(components) for _ in range(NODE_COUNT)]
    for component in components:
        for node in component["nodes"]:
            matrix[node][component["component_index"]] = 1
    return matrix


def quotient_edge_embedding_matrix(quotient_edges: list[dict[str, Any]]) -> list[list[int]]:
    matrix = [[0] * len(quotient_edges) for _ in range(EDGE_COUNT)]
    for quotient_edge in quotient_edges:
        matrix[quotient_edge["path_edge_index"]][quotient_edge["quotient_edge_index"]] = 1
    return matrix


def nonzero_row_labels(matrix: list[list[int]], labels: list[str]) -> list[str]:
    return [label for label, row in zip(labels, matrix) if any(row)]


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    quotient = load_json(QUOTIENT_REPORT)
    gf2 = load_json(GF2_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    components = quotient["component_quotient_definition"]["component_rows"]
    quotient_edges = quotient["component_quotient_definition"]["quotient_edge_rows"]
    quotient_component_bits = [component["bit"] for component in components]
    s_matrix = component_expansion_matrix(components)
    b_path = path_coboundary_matrix()
    b_quotient = quotient_coboundary_matrix(len(components))
    e_matrix = quotient_edge_embedding_matrix(quotient_edges)
    quotient_edge_bits = matvec_gf2(b_quotient, quotient_component_bits)
    lifted_node_bits = matvec_gf2(s_matrix, quotient_component_bits)
    embedded_edge_bits = matvec_gf2(e_matrix, quotient_edge_bits)
    b_path_times_s = matmul_gf2(b_path, s_matrix)
    e_times_b_quotient = matmul_gf2(e_matrix, b_quotient)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_LIFT_MATRIX_CERTIFICATE__COMMUTING_LIFT_DIAGRAM",
        "status": "component-quotient-lift-and-path-coboundary-square-commute-over-gf2",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_adjacency_certificate": str(QUOTIENT_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "quotient lift",
                "lift matrix",
                "component lift",
                "Bq",
                "component expansion",
                "quotient edge embedding",
                "B*S",
                "E*Bq",
            ],
            "finding": "Existing reports export the component quotient and path GF(2) reconstruction, but not the explicit component-expansion/edge-embedding matrices S and E or the commuting square B_path*S = E*B_quotient before this file.",
        },
        "matrix_definition": {
            "field": "GF(2)",
            "component_bits": quotient_component_bits,
            "quotient_edge_bits": quotient_edge_bits,
            "S_component_expansion_matrix_nodes_by_components": s_matrix,
            "B_path_coboundary_matrix_edges_by_nodes": b_path,
            "B_quotient_coboundary_matrix_quotient_edges_by_components": b_quotient,
            "E_quotient_edge_embedding_matrix_path_edges_by_quotient_edges": e_matrix,
            "B_path_times_S": b_path_times_s,
            "E_times_B_quotient": e_times_b_quotient,
            "commuting_square_equation": "B_path * S = E * B_quotient over GF(2)",
        },
        "rank_certificate": {
            "rank_S_component_expansion": gf2_rref(s_matrix)["rank"],
            "rank_E_quotient_edge_embedding": gf2_rref(e_matrix)["rank"],
            "rank_B_quotient": gf2_rref(b_quotient)["rank"],
            "rank_B_path_times_S": gf2_rref(b_path_times_s)["rank"],
            "S_has_full_column_rank": gf2_rref(s_matrix)["rank"] == len(components),
            "E_has_full_column_rank": gf2_rref(e_matrix)["rank"] == len(quotient_edges),
            "B_quotient_has_path_rank_component_count_minus_one": gf2_rref(b_quotient)["rank"] == len(components) - 1,
        },
        "reconstruction_certificate": {
            "node_bits_from_S_component_bits": lifted_node_bits,
            "edge_bits_from_E_quotient_edge_bits": embedded_edge_bits,
            "edge_bits_from_B_path_node_bits": matvec_gf2(b_path, node_bits),
            "nonzero_embedding_rows": nonzero_row_labels(e_matrix, EDGE_LABELS),
            "gf2_linear_system_edge_bits": gf2["gf2_linear_system_summary"]["solution_edge_bit_pattern"],
        },
        "component_quotient_lift_matrix_summary": {
            "matches_expected_component_bits": quotient_component_bits == EXPECTED_COMPONENT_BITS,
            "matches_expected_quotient_edge_bits": quotient_edge_bits == EXPECTED_QUOTIENT_EDGE_BITS,
            "S_lifts_component_bits_to_node_bits": lifted_node_bits == node_bits == EXPECTED_NODE_BITS,
            "E_embeds_quotient_edge_bits_to_path_edge_bits": embedded_edge_bits == edge_bits == EXPECTED_EDGE_BITS,
            "path_boundary_of_lifted_nodes_matches_edge_bits": matvec_gf2(b_path, lifted_node_bits) == edge_bits,
            "commuting_square_BS_equals_EBq": b_path_times_s == e_times_b_quotient,
            "embedding_rows_are_flip_edges": nonzero_row_labels(e_matrix, EDGE_LABELS) == [edge["path_edge"] for edge in quotient_edges],
            "matches_gf2_linear_system_edge_bits": embedded_edge_bits == gf2["gf2_linear_system_summary"]["solution_edge_bit_pattern"],
            "rank_S_full": gf2_rref(s_matrix)["rank"] == len(components),
            "rank_E_full": gf2_rref(e_matrix)["rank"] == len(quotient_edges),
            "rank_B_quotient_full_path_rank": gf2_rref(b_quotient)["rank"] == len(components) - 1,
        },
        "blocker_context": {
            "resolved_locally": [
                "The component-expansion matrix S lifts quotient component bits to all 12 audited node bits.",
                "The quotient-edge embedding matrix E embeds the four quotient flips into the 11 path edges.",
                "The square B_path*S = E*B_quotient commutes over GF(2).",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "lift_step": "S expands quotient component bits [0,1,0,1,0] to the audited 12 node bits.",
            "edge_embedding_step": "E embeds quotient edge bits [1,1,1,1] into path edges 1->2, 5->6, 7->8, and 9->10.",
            "commuting_square_step": "The finite GF(2) matrices satisfy B_path*S = E*B_quotient, so taking a path boundary after lifting equals embedding the quotient boundary.",
            "rank_step": "S has full column rank 5, E has full column rank 4, and B_quotient has path rank 4.",
            "theoretical_limit": "This is finite quotient-lift matrix bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["component_quotient_lift_matrix_summary"]
    rank = payload["rank_certificate"]
    lines = [
        "# Phase-sign component-quotient lift matrix certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate exports the component-expansion and quotient-edge embedding",
        "matrices and checks the GF(2) commuting square B_path*S = E*B_quotient.",
        "",
        "## Summary",
        "",
        f"- S lifts component bits to node bits: `{summary['S_lifts_component_bits_to_node_bits']}`",
        f"- E embeds quotient edge bits to path edge bits: `{summary['E_embeds_quotient_edge_bits_to_path_edge_bits']}`",
        f"- Commuting square passes: `{summary['commuting_square_BS_equals_EBq']}`",
        f"- Rank S: `{rank['rank_S_component_expansion']}`",
        f"- Rank E: `{rank['rank_E_quotient_edge_embedding']}`",
        f"- Rank B_quotient: `{rank['rank_B_quotient']}`",
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
    print(json.dumps(payload["component_quotient_lift_matrix_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
