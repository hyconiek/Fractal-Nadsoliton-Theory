#!/usr/bin/env python3
"""Scratch probe: phase-sign component-quotient projection certificate.

The component-quotient lift certificate records how quotient bits expand back to
path-node bits.  This companion certificate records the reverse finite GF(2)
projection/collapse map:

    q = Q * node_bits,
    Q * S = I_5,
    S * Q * node_bits = node_bits,
    G * B_path * S = B_quotient,
    H * B_path * S = 0.

Here Q selects one representative node from each maximal constant-sign component,
S is the component expansion matrix, G selects the quotient-boundary path edges,
and H selects interior path edges.  The result certifies that the quotient is a
valid finite collapse of the audited path data.  It is matrix bookkeeping only;
it does not derive phase zeros, omega/phi, damping, transport, a kernel bridge,
selector discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
QUOTIENT_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_adjacency_certificate_report.json"
LIFT_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
EXPECTED_COMPONENT_BITS = [0, 1, 0, 1, 0]
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_BOUNDARY_EDGE_INDICES = [1, 5, 7, 9]
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


def identity_matrix(size: int) -> list[list[int]]:
    return [[1 if row == col else 0 for col in range(size)] for row in range(size)]


def zero_matrix(rows: int, cols: int) -> list[list[int]]:
    return [[0 for _ in range(cols)] for _ in range(rows)]


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


def representative_projection_matrix(components: list[dict[str, Any]]) -> list[list[int]]:
    matrix = [[0] * NODE_COUNT for _ in components]
    for component in components:
        matrix[component["component_index"]][component["start"]] = 1
    return matrix


def edge_selector_matrix(edge_indices: list[int]) -> list[list[int]]:
    matrix = [[0] * EDGE_COUNT for _ in edge_indices]
    for row, edge_index in enumerate(edge_indices):
        matrix[row][edge_index] = 1
    return matrix


def projector_trace(matrix: list[list[int]]) -> int:
    return sum(matrix[index][index] for index in range(min(len(matrix), len(matrix[0]))))


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    quotient = load_json(QUOTIENT_REPORT)
    lift = load_json(LIFT_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    components = quotient["component_quotient_definition"]["component_rows"]
    quotient_edges = quotient["component_quotient_definition"]["quotient_edge_rows"]
    component_bits = [component["bit"] for component in components]
    boundary_edge_indices = [edge["path_edge_index"] for edge in quotient_edges]
    interior_edge_indices = [index for index in range(EDGE_COUNT) if index not in boundary_edge_indices]

    s_matrix = component_expansion_matrix(components)
    q_matrix = representative_projection_matrix(components)
    b_path = path_coboundary_matrix()
    b_quotient = quotient_coboundary_matrix(len(components))
    g_boundary = edge_selector_matrix(boundary_edge_indices)
    h_interior = edge_selector_matrix(interior_edge_indices)
    q_times_s = matmul_gf2(q_matrix, s_matrix)
    s_times_q = matmul_gf2(s_matrix, q_matrix)
    path_boundary_after_lift = matmul_gf2(b_path, s_matrix)
    boundary_restriction = matmul_gf2(g_boundary, path_boundary_after_lift)
    interior_restriction = matmul_gf2(h_interior, path_boundary_after_lift)
    projected_component_bits = matvec_gf2(q_matrix, node_bits)
    projected_then_lifted_node_bits = matvec_gf2(s_matrix, projected_component_bits)
    quotient_edge_bits_from_projection = matvec_gf2(b_quotient, projected_component_bits)
    boundary_edge_bits_from_path = matvec_gf2(g_boundary, edge_bits)
    interior_edge_bits_from_path = matvec_gf2(h_interior, edge_bits)
    projector_squared = matmul_gf2(s_times_q, s_times_q)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_PROJECTION_CERTIFICATE__COLLAPSE_SECTION_PROJECTOR",
        "status": "component-quotient-projection-section-and-boundary-restriction-certified-over-gf2",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_adjacency_certificate": str(QUOTIENT_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_lift_matrix_certificate": str(LIFT_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "component quotient projection",
                "representative projection",
                "Q*S",
                "S*Q",
                "component projector",
                "boundary selector",
                "interior edge selector",
                "G*B_path*S",
                "H*B_path*S",
            ],
            "finding": "Existing reports export the quotient, lift matrices S/E, and path boundary data, but not the reverse representative projection Q, the section identity Q*S=I, the projector S*Q, or the boundary/interior selector factorization G*B_path*S=B_quotient and H*B_path*S=0 before this file.",
        },
        "matrix_definition": {
            "field": "GF(2)",
            "component_bits": component_bits,
            "boundary_edge_indices": boundary_edge_indices,
            "interior_edge_indices": interior_edge_indices,
            "boundary_edge_labels": [EDGE_LABELS[index] for index in boundary_edge_indices],
            "interior_edge_labels": [EDGE_LABELS[index] for index in interior_edge_indices],
            "Q_representative_projection_matrix_components_by_nodes": q_matrix,
            "S_component_expansion_matrix_nodes_by_components": s_matrix,
            "S_times_Q_projector_matrix_nodes_by_nodes": s_times_q,
            "B_path_coboundary_matrix_edges_by_nodes": b_path,
            "B_quotient_coboundary_matrix_quotient_edges_by_components": b_quotient,
            "G_boundary_edge_selector_matrix_quotient_edges_by_path_edges": g_boundary,
            "H_interior_edge_selector_matrix_interior_edges_by_path_edges": h_interior,
            "Q_times_S": q_times_s,
            "G_times_B_path_times_S": boundary_restriction,
            "H_times_B_path_times_S": interior_restriction,
            "section_equation": "Q * S = I_5 over GF(2)",
            "projector_equation_on_audited_nodes": "S * Q * node_bits = node_bits",
            "boundary_restriction_equation": "G * B_path * S = B_quotient over GF(2)",
            "interior_vanish_equation": "H * B_path * S = 0 over GF(2)",
        },
        "rank_certificate": {
            "rank_Q_representative_projection": gf2_rref(q_matrix)["rank"],
            "rank_S_component_expansion": gf2_rref(s_matrix)["rank"],
            "rank_S_times_Q_projector": gf2_rref(s_times_q)["rank"],
            "projector_trace": projector_trace(s_times_q),
            "rank_boundary_restriction": gf2_rref(boundary_restriction)["rank"],
            "Q_has_full_row_rank": gf2_rref(q_matrix)["rank"] == len(components),
            "S_has_full_column_rank": gf2_rref(s_matrix)["rank"] == len(components),
            "S_times_Q_has_component_subspace_rank": gf2_rref(s_times_q)["rank"] == len(components),
            "S_times_Q_is_idempotent": projector_squared == s_times_q,
        },
        "reconstruction_certificate": {
            "component_bits_from_Q_node_bits": projected_component_bits,
            "node_bits_from_S_Q_node_bits": projected_then_lifted_node_bits,
            "quotient_edge_bits_from_B_quotient_Q_node_bits": quotient_edge_bits_from_projection,
            "boundary_edge_bits_from_G_path_edge_bits": boundary_edge_bits_from_path,
            "interior_edge_bits_from_H_path_edge_bits": interior_edge_bits_from_path,
            "lift_report_node_bits_from_S_component_bits": lift["reconstruction_certificate"]["node_bits_from_S_component_bits"],
            "lift_report_edge_bits_from_E_quotient_edge_bits": lift["reconstruction_certificate"]["edge_bits_from_E_quotient_edge_bits"],
        },
        "component_quotient_projection_summary": {
            "matches_expected_component_bits": component_bits == EXPECTED_COMPONENT_BITS,
            "Q_projects_node_bits_to_component_bits": projected_component_bits == component_bits == EXPECTED_COMPONENT_BITS,
            "Q_times_S_is_identity": q_times_s == identity_matrix(len(components)),
            "S_times_Q_is_idempotent_projector": projector_squared == s_times_q,
            "S_Q_preserves_audited_node_bits": projected_then_lifted_node_bits == node_bits == EXPECTED_NODE_BITS,
            "boundary_selector_edges_are_flip_edges": boundary_edge_indices == EXPECTED_BOUNDARY_EDGE_INDICES,
            "boundary_restriction_equals_quotient_coboundary": boundary_restriction == b_quotient,
            "interior_restriction_vanishes": interior_restriction == zero_matrix(len(interior_edge_indices), len(components)),
            "quotient_boundary_matches_selected_path_edges": quotient_edge_bits_from_projection == boundary_edge_bits_from_path == [1, 1, 1, 1],
            "interior_path_edges_have_zero_bits": interior_edge_bits_from_path == [0] * len(interior_edge_indices),
            "lift_report_round_trip_matches": projected_then_lifted_node_bits == lift["reconstruction_certificate"]["node_bits_from_S_component_bits"] and edge_bits == lift["reconstruction_certificate"]["edge_bits_from_E_quotient_edge_bits"] == EXPECTED_EDGE_BITS,
            "rank_Q_full": gf2_rref(q_matrix)["rank"] == len(components),
            "projector_rank_is_component_count": gf2_rref(s_times_q)["rank"] == len(components),
        },
        "blocker_context": {
            "resolved_locally": [
                "The representative projection Q collapses audited node bits to quotient component bits.",
                "The expansion S is a section of Q because Q*S=I_5.",
                "Boundary selectors recover B_quotient while interior selectors vanish on the component subspace.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "projection_step": "Q selects component representatives and maps the audited 12 node bits to quotient bits [0,1,0,1,0].",
            "section_step": "The finite GF(2) matrices satisfy Q*S=I_5, so S is an injective section of the quotient projection Q.",
            "projector_step": "S*Q is idempotent with rank 5 and fixes the audited node-bit vector because it is constant on quotient components.",
            "boundary_selector_step": "The selectors satisfy G*B_path*S=B_quotient on flip edges 1->2, 5->6, 7->8, 9->10 and H*B_path*S=0 on interior edges.",
            "theoretical_limit": "This is finite quotient projection bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["component_quotient_projection_summary"]
    rank = payload["rank_certificate"]
    lines = [
        "# Phase-sign component-quotient projection certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate exports the representative projection Q and boundary/interior",
        "selector matrices and checks Q*S=I_5, S*Q projector idempotence,",
        "G*B_path*S=B_quotient, and H*B_path*S=0 over GF(2).",
        "",
        "## Summary",
        "",
        f"- Q projects node bits to component bits: `{summary['Q_projects_node_bits_to_component_bits']}`",
        f"- Q*S is identity: `{summary['Q_times_S_is_identity']}`",
        f"- S*Q fixes audited node bits: `{summary['S_Q_preserves_audited_node_bits']}`",
        f"- Boundary restriction equals quotient coboundary: `{summary['boundary_restriction_equals_quotient_coboundary']}`",
        f"- Interior restriction vanishes: `{summary['interior_restriction_vanishes']}`",
        f"- Rank Q: `{rank['rank_Q_representative_projection']}`",
        f"- Rank S*Q: `{rank['rank_S_times_Q_projector']}`",
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
    print(json.dumps(payload["component_quotient_projection_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
