#!/usr/bin/env python3
"""Scratch probe: component-quotient exact-sequence certificate over GF(2).

The quotient projection certificate verifies a representative projection Q and
section S on the audited vector.  This probe records the ambient linear-algebra
fact behind that collapse:

    0 -> C^0(quotient) --S--> C^0(path) --F=H*B_path--> C^1(interior) -> 0

is exact for the audited component partition.  Equivalently, the image of S is
precisely the kernel of the interior-edge coboundary F: nodes are in the quotient
subspace exactly when every interior edge has zero coboundary.  This is finite
GF(2) bookkeeping only; it does not derive phase zeros, omega/phi, damping,
transport, a kernel bridge, selector discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_exact_sequence_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_exact_sequence_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
PROJECTION_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_report.json"
LIFT_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_lift_matrix_certificate_report.json"

NODE_COUNT = 12
INTERIOR_EDGE_COUNT = 7
COMPONENT_COUNT = 5
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_COMPONENT_BITS = [0, 1, 0, 1, 0]


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


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [[matrix[row][col] for row in range(len(matrix))] for col in range(len(matrix[0]))]


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
        "pivot_columns": [pivot["pivot_col"] for pivot in pivots],
        "rref_matrix": work,
    }


def rank_gf2(matrix: list[list[int]]) -> int:
    return gf2_rref(matrix)["rank"]


def nullspace_basis_gf2(matrix: list[list[int]]) -> list[list[int]]:
    rref = gf2_rref(matrix)
    col_count = len(matrix[0]) if matrix else 0
    pivot_columns = rref["pivot_columns"]
    free_columns = [col for col in range(col_count) if col not in pivot_columns]
    basis = []
    for free_col in free_columns:
        vector = [0] * col_count
        vector[free_col] = 1
        for pivot in reversed(rref["pivot_rows"]):
            row = pivot["row"]
            pivot_col = pivot["pivot_col"]
            vector[pivot_col] = sum(rref["rref_matrix"][row][col] & vector[col] for col in free_columns) % 2
        basis.append(vector)
    return basis


def columns_to_matrix(columns: list[list[int]]) -> list[list[int]]:
    return [[columns[col][row] for col in range(len(columns))] for row in range(len(columns[0]))]


def column_span_rank(columns: list[list[int]]) -> int:
    return rank_gf2(transpose(columns)) if columns else 0


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    projection = load_json(PROJECTION_REPORT)
    lift = load_json(LIFT_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    matrix_definition = projection["matrix_definition"]
    q_matrix = matrix_definition["Q_representative_projection_matrix_components_by_nodes"]
    s_matrix = matrix_definition["S_component_expansion_matrix_nodes_by_components"]
    b_path = matrix_definition["B_path_coboundary_matrix_edges_by_nodes"]
    h_interior = matrix_definition["H_interior_edge_selector_matrix_interior_edges_by_path_edges"]
    f_matrix = matmul_gf2(h_interior, b_path)
    f_times_s = matmul_gf2(f_matrix, s_matrix)
    sq_projector = matrix_definition["S_times_Q_projector_matrix_nodes_by_nodes"]
    kernel_basis_rows = nullspace_basis_gf2(f_matrix)
    kernel_basis_as_columns = columns_to_matrix(kernel_basis_rows)
    s_columns = transpose(s_matrix)
    combined_columns = s_columns + kernel_basis_rows
    kernel_basis_projected = [matvec_gf2(sq_projector, basis_vector) for basis_vector in kernel_basis_rows]
    projected_component_bits = matvec_gf2(q_matrix, node_bits)
    projected_then_lifted = matvec_gf2(s_matrix, projected_component_bits)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_EXACT_SEQUENCE_CERTIFICATE__IMAGE_EQUALS_INTERIOR_KERNEL",
        "status": "component-quotient-exact-sequence-imS-equals-ker-HBpath-certified-over-gf2",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_projection_certificate": str(PROJECTION_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_lift_matrix_certificate": str(LIFT_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "exact sequence",
                "im(S)=ker",
                "image equals kernel",
                "H*B_path",
                "interior kernel",
                "component subspace",
                "quotient subspace",
                "ker(HB)",
            ],
            "finding": "Existing quotient projection/lift reports export Q, S, and selector equations, but not the exact-sequence certificate proving im(S)=ker(H*B_path) by rank-nullity and kernel-basis projection before this file.",
        },
        "exact_sequence_definition": {
            "field": "GF(2)",
            "sequence": "0 -> GF(2)^5 --S--> GF(2)^12 --F=H*B_path--> GF(2)^7 -> 0",
            "S_component_expansion_matrix_nodes_by_components": s_matrix,
            "Q_representative_projection_matrix_components_by_nodes": q_matrix,
            "H_interior_edge_selector_matrix_interior_edges_by_path_edges": h_interior,
            "B_path_coboundary_matrix_edges_by_nodes": b_path,
            "F_interior_coboundary_matrix_equals_H_times_B_path": f_matrix,
            "F_times_S": f_times_s,
            "S_times_Q_projector": sq_projector,
            "kernel_basis_rows_for_F": kernel_basis_rows,
            "kernel_basis_as_node_by_basis_matrix": kernel_basis_as_columns,
            "kernel_basis_projected_by_SQ": kernel_basis_projected,
        },
        "rank_nullity_certificate": {
            "rank_F_interior_coboundary": rank_gf2(f_matrix),
            "nullity_F_interior_coboundary": NODE_COUNT - rank_gf2(f_matrix),
            "rank_S_component_expansion": rank_gf2(s_matrix),
            "rank_kernel_basis": column_span_rank(kernel_basis_rows),
            "rank_combined_S_columns_and_kernel_basis": column_span_rank(combined_columns),
            "F_has_full_row_rank": rank_gf2(f_matrix) == INTERIOR_EDGE_COUNT,
            "rank_nullity_matches_component_count": NODE_COUNT - rank_gf2(f_matrix) == COMPONENT_COUNT,
            "S_rank_matches_kernel_dimension": rank_gf2(s_matrix) == NODE_COUNT - rank_gf2(f_matrix),
            "combined_span_has_no_extra_dimension": column_span_rank(combined_columns) == COMPONENT_COUNT,
        },
        "reconstruction_certificate": {
            "F_node_bits": matvec_gf2(f_matrix, node_bits),
            "Q_node_bits": projected_component_bits,
            "S_Q_node_bits": projected_then_lifted,
            "lift_report_node_bits_from_S_component_bits": lift["reconstruction_certificate"]["node_bits_from_S_component_bits"],
            "projection_report_component_bits_from_Q_node_bits": projection["reconstruction_certificate"]["component_bits_from_Q_node_bits"],
            "projection_report_node_bits_from_S_Q_node_bits": projection["reconstruction_certificate"]["node_bits_from_S_Q_node_bits"],
        },
        "component_quotient_exact_sequence_summary": {
            "F_equals_H_times_B_path_rows_by_nodes": len(f_matrix) == INTERIOR_EDGE_COUNT and len(f_matrix[0]) == NODE_COUNT,
            "F_times_S_is_zero": f_times_s == zero_matrix(INTERIOR_EDGE_COUNT, COMPONENT_COUNT),
            "rank_F_is_interior_edge_count": rank_gf2(f_matrix) == INTERIOR_EDGE_COUNT,
            "nullity_F_is_component_count": NODE_COUNT - rank_gf2(f_matrix) == COMPONENT_COUNT,
            "S_columns_span_kernel_F": rank_gf2(s_matrix) == COMPONENT_COUNT and column_span_rank(combined_columns) == COMPONENT_COUNT,
            "S_Q_projects_kernel_basis_to_itself": kernel_basis_projected == kernel_basis_rows,
            "audited_node_bits_lie_in_kernel_F": matvec_gf2(f_matrix, node_bits) == [0] * INTERIOR_EDGE_COUNT,
            "Q_then_S_recovers_audited_node_bits": projected_then_lifted == node_bits == EXPECTED_NODE_BITS,
            "Q_recovers_expected_component_bits": projected_component_bits == EXPECTED_COMPONENT_BITS,
            "matches_projection_report_round_trip": projected_component_bits == projection["reconstruction_certificate"]["component_bits_from_Q_node_bits"] and projected_then_lifted == projection["reconstruction_certificate"]["node_bits_from_S_Q_node_bits"],
            "matches_lift_report_node_bits": projected_then_lifted == lift["reconstruction_certificate"]["node_bits_from_S_component_bits"],
        },
        "blocker_context": {
            "resolved_locally": [
                "The interior-edge coboundary F=H*B_path has rank 7 and nullity 5.",
                "The component expansion image im(S) has dimension 5 and lies in ker(F).",
                "Rank-nullity plus the combined-span check proves im(S)=ker(F).",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "exactness_step": "F=H*B_path has rank 7 and nullity 5, while S has rank 5 and F*S=0; therefore im(S)=ker(F).",
            "projector_step": "The projector S*Q fixes every exported kernel-basis vector, so the quotient projector is the identity on ker(F).",
            "audited_vector_step": "The audited node-bit vector has F*b=0 and is recovered by S*Q*b with quotient bits [0,1,0,1,0].",
            "theoretical_limit": "This is finite exact-sequence bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["component_quotient_exact_sequence_summary"]
    rank = payload["rank_nullity_certificate"]
    lines = [
        "# Phase-sign component-quotient exact-sequence certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate proves the finite GF(2) exactness statement im(S)=ker(F)",
        "for F=H*B_path, where H selects interior component edges.",
        "",
        "## Summary",
        "",
        f"- F*S is zero: `{summary['F_times_S_is_zero']}`",
        f"- rank(F): `{rank['rank_F_interior_coboundary']}`",
        f"- nullity(F): `{rank['nullity_F_interior_coboundary']}`",
        f"- S columns span ker(F): `{summary['S_columns_span_kernel_F']}`",
        f"- S*Q fixes kernel basis: `{summary['S_Q_projects_kernel_basis_to_itself']}`",
        f"- Audited node bits lie in ker(F): `{summary['audited_node_bits_lie_in_kernel_F']}`",
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
    print(json.dumps(payload["component_quotient_exact_sequence_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
