#!/usr/bin/env python3
"""Scratch probe: component-quotient complement inverse certificate.

The exact-sequence certificate proves im(S)=ker(F) for F=H*B_path.  This probe
adds the complementary finite GF(2) splitting: choose a residual basis N for
ker(Q), prove

    C^0(path) = im(S) direct_sum im(N),
    Q*S = I_5, Q*N = 0,
    F*S = 0, F*N is invertible.

Thus every node vector has a quotient component plus an interior residual, and
the residual is uniquely determined by its interior-edge coboundary.  This is
finite linear algebra only; it does not derive phase zeros, omega/phi, damping,
transport, a kernel bridge, selector discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_complement_inverse_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_complement_inverse_certificate_report.md"
PROJECTION_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_projection_certificate_report.json"
EXACT_SEQUENCE_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_exact_sequence_certificate_report.json"
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


def transpose(matrix: list[list[int]]) -> list[list[int]]:
    return [[matrix[row][col] for row in range(len(matrix))] for col in range(len(matrix[0]))]


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
        "pivot_columns": [pivot["pivot_col"] for pivot in pivots],
        "rref_matrix": work,
    }


def rank_gf2(matrix: list[list[int]]) -> int:
    return gf2_rref(matrix)["rank"]


def invert_square_gf2(matrix: list[list[int]]) -> dict[str, Any]:
    size = len(matrix)
    work = [row[:] + identity for row, identity in zip(matrix, identity_matrix(size))]
    pivot_row = 0
    for col in range(size):
        pivot = None
        for candidate in range(pivot_row, size):
            if work[candidate][col]:
                pivot = candidate
                break
        if pivot is None:
            raise ValueError("matrix is singular over GF(2)")
        if pivot != pivot_row:
            work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        for target in range(size):
            if target != pivot_row and work[target][col]:
                work[target] = [a ^ b for a, b in zip(work[target], work[pivot_row])]
        pivot_row += 1
    inverse = [row[size:] for row in work]
    return {
        "inverse": inverse,
        "left_inverse_check": matmul_gf2(inverse, matrix),
        "right_inverse_check": matmul_gf2(matrix, inverse),
    }


def residual_basis_matrix(
    q_matrix: list[list[int]],
    s_matrix: list[list[int]],
) -> tuple[list[list[int]], list[dict[str, int]]]:
    representatives = [row.index(1) for row in q_matrix]
    residual_columns = []
    rows = []
    for node, s_row in enumerate(s_matrix):
        component_index = s_row.index(1)
        representative = representatives[component_index]
        if node == representative:
            continue
        column = [0] * NODE_COUNT
        column[node] = 1
        residual_columns.append(column)
        rows.append({"component_index": component_index, "representative_node": representative, "residual_node": node})
    return [[column[row] for column in residual_columns] for row in range(NODE_COUNT)], rows


def build_payload() -> dict[str, Any]:
    projection = load_json(PROJECTION_REPORT)
    exact_sequence = load_json(EXACT_SEQUENCE_REPORT)
    z2 = load_json(Z2_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    matrix_definition = projection["matrix_definition"]
    q_matrix = matrix_definition["Q_representative_projection_matrix_components_by_nodes"]
    s_matrix = matrix_definition["S_component_expansion_matrix_nodes_by_components"]
    f_matrix = exact_sequence["exact_sequence_definition"]["F_interior_coboundary_matrix_equals_H_times_B_path"]
    n_matrix, residual_rows = residual_basis_matrix(q_matrix, s_matrix)
    qs = matmul_gf2(q_matrix, s_matrix)
    qn = matmul_gf2(q_matrix, n_matrix)
    fs = matmul_gf2(f_matrix, s_matrix)
    fn = matmul_gf2(f_matrix, n_matrix)
    fn_inverse_data = invert_square_gf2(fn)
    block_basis = [[*s_matrix[row], *n_matrix[row]] for row in range(NODE_COUNT)]
    component_bits = matvec_gf2(q_matrix, node_bits)
    quotient_part = matvec_gf2(s_matrix, component_bits)
    residual_part = [left ^ right for left, right in zip(node_bits, quotient_part)]
    residual_coordinates_from_F = matvec_gf2(fn_inverse_data["inverse"], matvec_gf2(f_matrix, node_bits))
    residual_from_coordinates = matvec_gf2(n_matrix, residual_coordinates_from_F)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COMPLEMENT_INVERSE_CERTIFICATE__DIRECT_SUM_AND_INTERIOR_INVERSE",
        "status": "component-quotient-complement-direct-sum-and-FN-inverse-certified-over-gf2",
        "source_reports": {
            "phase_sign_component_quotient_projection_certificate": str(PROJECTION_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_exact_sequence_certificate": str(EXACT_SEQUENCE_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "component quotient complement",
                "residual basis",
                "direct sum",
                "F*N",
                "FN inverse",
                "ker(Q)",
                "interior residual inverse",
                "block basis [S N]",
            ],
            "finding": "Existing quotient reports prove lift/projection/exactness, but not a complement basis N for ker(Q), the direct-sum matrix [S N], or the invertibility of F*N on the residual subspace before this file.",
        },
        "complement_definition": {
            "field": "GF(2)",
            "Q_representative_projection_matrix_components_by_nodes": q_matrix,
            "S_component_expansion_matrix_nodes_by_components": s_matrix,
            "N_residual_basis_matrix_nodes_by_residuals": n_matrix,
            "residual_basis_rows": residual_rows,
            "F_interior_coboundary_matrix": f_matrix,
            "Q_times_S": qs,
            "Q_times_N": qn,
            "F_times_S": fs,
            "F_times_N": fn,
            "F_times_N_inverse": fn_inverse_data["inverse"],
            "FN_inverse_times_FN": fn_inverse_data["left_inverse_check"],
            "FN_times_FN_inverse": fn_inverse_data["right_inverse_check"],
            "block_basis_matrix_S_then_N": block_basis,
        },
        "rank_certificate": {
            "rank_S": rank_gf2(s_matrix),
            "rank_N": rank_gf2(n_matrix),
            "rank_block_basis_S_then_N": rank_gf2(block_basis),
            "rank_F_times_N": rank_gf2(fn),
            "S_and_N_span_all_nodes": rank_gf2(block_basis) == NODE_COUNT,
            "N_has_residual_rank": rank_gf2(n_matrix) == RESIDUAL_DIMENSION,
            "F_times_N_invertible": rank_gf2(fn) == RESIDUAL_DIMENSION,
        },
        "reconstruction_certificate": {
            "component_bits_Q_node_bits": component_bits,
            "quotient_part_S_Q_node_bits": quotient_part,
            "residual_part_node_bits_plus_SQnode_bits": residual_part,
            "F_node_bits": matvec_gf2(f_matrix, node_bits),
            "residual_coordinates_from_FN_inverse_F_node_bits": residual_coordinates_from_F,
            "residual_from_coordinates": residual_from_coordinates,
            "exact_sequence_F_node_bits": exact_sequence["reconstruction_certificate"]["F_node_bits"],
        },
        "component_quotient_complement_inverse_summary": {
            "Q_times_S_is_identity": qs == identity_matrix(COMPONENT_COUNT),
            "Q_times_N_is_zero": qn == zero_matrix(COMPONENT_COUNT, RESIDUAL_DIMENSION),
            "F_times_S_is_zero": fs == zero_matrix(RESIDUAL_DIMENSION, COMPONENT_COUNT),
            "F_times_N_is_invertible": rank_gf2(fn) == RESIDUAL_DIMENSION and fn_inverse_data["left_inverse_check"] == identity_matrix(RESIDUAL_DIMENSION) and fn_inverse_data["right_inverse_check"] == identity_matrix(RESIDUAL_DIMENSION),
            "block_basis_is_full_rank_12": rank_gf2(block_basis) == NODE_COUNT,
            "direct_sum_dimensions_add_to_12": rank_gf2(s_matrix) + rank_gf2(n_matrix) == NODE_COUNT,
            "audited_node_bits_have_zero_residual_part": residual_part == [0] * NODE_COUNT,
            "audited_residual_coordinates_zero": residual_coordinates_from_F == [0] * RESIDUAL_DIMENSION,
            "audited_quotient_part_matches_node_bits": quotient_part == node_bits == EXPECTED_NODE_BITS,
            "component_bits_match_expected": component_bits == EXPECTED_COMPONENT_BITS,
            "matches_exact_sequence_F_node_bits": matvec_gf2(f_matrix, node_bits) == exact_sequence["reconstruction_certificate"]["F_node_bits"],
        },
        "blocker_context": {
            "resolved_locally": [
                "A residual basis N for ker(Q) is exported explicitly.",
                "The block matrix [S N] has full rank 12, proving a direct-sum decomposition.",
                "F*N is invertible, so interior-edge coboundaries uniquely coordinate residuals.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "direct_sum_step": "The block matrix [S N] has rank 12, with rank(S)=5 and rank(N)=7, so C^0(path)=im(S) direct_sum im(N).",
            "annihilation_step": "Q*N=0 and F*S=0, while Q*S=I_5, so N is a residual complement to quotient constants.",
            "inverse_step": "F*N has rank 7 and a two-sided inverse, so F identifies the residual complement with interior-edge cochains.",
            "audited_vector_step": "For the audited node bits, Q*b=[0,1,0,1,0], S*Q*b=b, and the residual coordinates are all zero.",
            "theoretical_limit": "This is finite complement-inverse bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["component_quotient_complement_inverse_summary"]
    rank = payload["rank_certificate"]
    lines = [
        "# Phase-sign component-quotient complement inverse certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate exports a residual basis N and proves [S N] is a full",
        "GF(2) basis with F*N invertible on the residual complement.",
        "",
        "## Summary",
        "",
        f"- Q*N is zero: `{summary['Q_times_N_is_zero']}`",
        f"- F*S is zero: `{summary['F_times_S_is_zero']}`",
        f"- F*N invertible: `{summary['F_times_N_is_invertible']}`",
        f"- rank [S N]: `{rank['rank_block_basis_S_then_N']}`",
        f"- audited residual coordinates zero: `{summary['audited_residual_coordinates_zero']}`",
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
    print(json.dumps(payload["component_quotient_complement_inverse_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
