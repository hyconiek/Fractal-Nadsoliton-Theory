#!/usr/bin/env python3
"""Scratch probe: component-quotient coordinate-isomorphism certificate.

The complement-inverse certificate proves a direct-sum splitting.  This probe
packages that splitting as one explicit GF(2) coordinate isomorphism:

    T = [Q ; F] : C^0(path) -> C^0(quotient) plus C^1(interior),
    U = [S | N*(F*N)^(-1)] : quotient plus interior coordinates -> C^0(path),
    T*U = I_12 and U*T = I_12.

It certifies that path node bits are equivalently quotient component bits plus
interior residual coboundary bits.  This is finite linear algebra only; it does
not derive phase zeros, omega/phi, damping, transport, a kernel bridge, selector
discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.md"
COMPLEMENT_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_complement_inverse_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"

NODE_COUNT = 12
COMPONENT_COUNT = 5
RESIDUAL_DIMENSION = 7
EXPECTED_COMPONENT_BITS = [0, 1, 0, 1, 0]
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]


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


def vertical_stack(top: list[list[int]], bottom: list[list[int]]) -> list[list[int]]:
    return [row[:] for row in top] + [row[:] for row in bottom]


def horizontal_stack(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    return [[*left[row], *right[row]] for row in range(len(left))]


def split_coordinate_vector(vector: list[int]) -> dict[str, list[int]]:
    return {
        "quotient_component_coordinates": vector[:COMPONENT_COUNT],
        "interior_residual_coordinates": vector[COMPONENT_COUNT:],
    }


def build_payload() -> dict[str, Any]:
    complement = load_json(COMPLEMENT_REPORT)
    z2 = load_json(Z2_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    complement_definition = complement["complement_definition"]
    q_matrix = complement_definition["Q_representative_projection_matrix_components_by_nodes"]
    s_matrix = complement_definition["S_component_expansion_matrix_nodes_by_components"]
    n_matrix = complement_definition["N_residual_basis_matrix_nodes_by_residuals"]
    f_matrix = complement_definition["F_interior_coboundary_matrix"]
    fn_inverse = complement_definition["F_times_N_inverse"]
    coordinate_matrix_t = vertical_stack(q_matrix, f_matrix)
    residual_lift_matrix = matmul_gf2(n_matrix, fn_inverse)
    inverse_matrix_u = horizontal_stack(s_matrix, residual_lift_matrix)
    t_times_u = matmul_gf2(coordinate_matrix_t, inverse_matrix_u)
    u_times_t = matmul_gf2(inverse_matrix_u, coordinate_matrix_t)
    coordinate_vector = matvec_gf2(coordinate_matrix_t, node_bits)
    reconstructed_node_bits = matvec_gf2(inverse_matrix_u, coordinate_vector)
    split_coordinates = split_coordinate_vector(coordinate_vector)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_ISOMORPHISM_CERTIFICATE__TWO_SIDED_BLOCK_INVERSE",
        "status": "component-quotient-coordinate-isomorphism-two-sided-inverse-certified-over-gf2",
        "source_reports": {
            "phase_sign_component_quotient_complement_inverse_certificate": str(COMPLEMENT_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "coordinate isomorphism",
                "T*U",
                "U*T",
                "[Q;F]",
                "[S | N",
                "block inverse",
                "quotient plus interior coordinates",
            ],
            "finding": "Existing complement reports prove [S N] and F*N separately, but not the single coordinate isomorphism T=[Q;F] with explicit two-sided inverse U=[S | N*(F*N)^(-1)] before this file.",
        },
        "coordinate_isomorphism_definition": {
            "field": "GF(2)",
            "coordinate_map_T_rows_quotient_then_interior": coordinate_matrix_t,
            "inverse_map_U_columns_quotient_then_interior": inverse_matrix_u,
            "residual_lift_matrix_N_times_FN_inverse": residual_lift_matrix,
            "T_times_U": t_times_u,
            "U_times_T": u_times_t,
            "coordinate_map_statement": "T=[Q;F] sends node bits to quotient coordinates followed by interior residual coordinates.",
            "inverse_map_statement": "U=[S | N*(F*N)^(-1)] reconstructs node bits from quotient plus interior coordinates.",
        },
        "rank_inverse_certificate": {
            "rank_T_coordinate_map": rank_gf2(coordinate_matrix_t),
            "rank_U_inverse_map": rank_gf2(inverse_matrix_u),
            "T_has_full_rank": rank_gf2(coordinate_matrix_t) == NODE_COUNT,
            "U_has_full_rank": rank_gf2(inverse_matrix_u) == NODE_COUNT,
            "T_times_U_is_identity": t_times_u == identity_matrix(NODE_COUNT),
            "U_times_T_is_identity": u_times_t == identity_matrix(NODE_COUNT),
        },
        "reconstruction_certificate": {
            "node_bits": node_bits,
            "coordinate_vector_T_node_bits": coordinate_vector,
            "quotient_component_coordinates": split_coordinates["quotient_component_coordinates"],
            "interior_residual_coordinates": split_coordinates["interior_residual_coordinates"],
            "node_bits_from_U_T_node_bits": reconstructed_node_bits,
            "complement_report_component_bits": complement["reconstruction_certificate"]["component_bits_Q_node_bits"],
            "complement_report_residual_coordinates": complement["reconstruction_certificate"]["residual_coordinates_from_FN_inverse_F_node_bits"],
        },
        "component_quotient_coordinate_isomorphism_summary": {
            "T_rank_full_12": rank_gf2(coordinate_matrix_t) == NODE_COUNT,
            "U_rank_full_12": rank_gf2(inverse_matrix_u) == NODE_COUNT,
            "T_U_identity": t_times_u == identity_matrix(NODE_COUNT),
            "U_T_identity": u_times_t == identity_matrix(NODE_COUNT),
            "coordinate_vector_splits_expected": split_coordinates["quotient_component_coordinates"] == EXPECTED_COMPONENT_BITS and split_coordinates["interior_residual_coordinates"] == [0] * RESIDUAL_DIMENSION,
            "U_reconstructs_audited_node_bits": reconstructed_node_bits == node_bits == EXPECTED_NODE_BITS,
            "matches_complement_report_coordinates": split_coordinates["quotient_component_coordinates"] == complement["reconstruction_certificate"]["component_bits_Q_node_bits"] and split_coordinates["interior_residual_coordinates"] == complement["reconstruction_certificate"]["residual_coordinates_from_FN_inverse_F_node_bits"],
        },
        "blocker_context": {
            "resolved_locally": [
                "The block coordinate map T=[Q;F] has rank 12.",
                "The explicit inverse U=[S | N*(F*N)^(-1)] is verified on both sides.",
                "The audited vector has quotient coordinates [0,1,0,1,0] and zero interior residual coordinates.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "coordinate_step": "T=[Q;F] maps C^0(path) to quotient coordinates plus interior residual coordinates and has rank 12.",
            "inverse_step": "U=[S | N*(F*N)^(-1)] satisfies T*U=I_12 and U*T=I_12 over GF(2).",
            "audited_vector_step": "The audited node bits map to coordinates [0,1,0,1,0 | 0,0,0,0,0,0,0] and reconstruct exactly under U.",
            "theoretical_limit": "This is finite coordinate-isomorphism bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["component_quotient_coordinate_isomorphism_summary"]
    rank = payload["rank_inverse_certificate"]
    lines = [
        "# Phase-sign component-quotient coordinate-isomorphism certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate exports T=[Q;F] and U=[S | N*(F*N)^(-1)] and checks",
        "the two-sided identities T*U=I_12 and U*T=I_12 over GF(2).",
        "",
        "## Summary",
        "",
        f"- rank(T): `{rank['rank_T_coordinate_map']}`",
        f"- rank(U): `{rank['rank_U_inverse_map']}`",
        f"- T*U identity: `{summary['T_U_identity']}`",
        f"- U*T identity: `{summary['U_T_identity']}`",
        f"- audited coordinates expected: `{summary['coordinate_vector_splits_expected']}`",
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
    print(json.dumps(payload["component_quotient_coordinate_isomorphism_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
