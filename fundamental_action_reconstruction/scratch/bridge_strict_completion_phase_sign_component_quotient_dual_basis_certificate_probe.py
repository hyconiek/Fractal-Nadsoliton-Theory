#!/usr/bin/env python3
"""Scratch probe: component-quotient dual-basis certificate over GF(2).

The coordinate-isomorphism certificate proves T=[Q;F] and
U=[S | N*(F*N)^(-1)] are two-sided inverses.  This probe unfolds that block
inverse as a finite dual-basis/decomposition certificate:

    <T_i, U_j> = delta_ij,
    b = sum_i (T_i b) U_i.

For the audited phase-sign node vector, only quotient coordinates 1 and 3 are
active and all seven interior residual coordinates are zero.  This is finite
linear algebra only; it does not derive phase zeros, omega/phi, damping,
transport, a kernel bridge, selector discharge, or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_dual_basis_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_dual_basis_certificate_report.md"
COORDINATE_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"

NODE_COUNT = 12
COMPONENT_COUNT = 5
RESIDUAL_DIMENSION = 7
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_COORDINATE_VECTOR = [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]


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


def dot_gf2(row: list[int], column: list[int]) -> int:
    return sum(left & right for left, right in zip(row, column)) % 2


def add_vectors_gf2(left: list[int], right: list[int]) -> list[int]:
    return [a ^ b for a, b in zip(left, right)]


def scale_vector_gf2(bit: int, vector: list[int]) -> list[int]:
    return [bit & entry for entry in vector]


def coordinate_label(index: int) -> str:
    if index < COMPONENT_COUNT:
        return f"quotient_component_{index}"
    return f"interior_residual_{index - COMPONENT_COUNT}"


def active_decomposition_rows(coordinates: list[int], u_columns: list[list[int]]) -> list[dict[str, Any]]:
    rows = []
    partial = [0] * NODE_COUNT
    for index, bit in enumerate(coordinates):
        contribution = scale_vector_gf2(bit, u_columns[index])
        partial = add_vectors_gf2(partial, contribution)
        rows.append({
            "coordinate_index": index,
            "coordinate_label": coordinate_label(index),
            "coordinate_bit": bit,
            "basis_column": u_columns[index],
            "contribution": contribution,
            "partial_sum_after_row": partial[:],
        })
    return rows


def build_payload() -> dict[str, Any]:
    coordinate = load_json(COORDINATE_REPORT)
    z2 = load_json(Z2_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    definition = coordinate["coordinate_isomorphism_definition"]
    t_matrix = definition["coordinate_map_T_rows_quotient_then_interior"]
    u_matrix = definition["inverse_map_U_columns_quotient_then_interior"]
    u_columns = transpose(u_matrix)
    pairing_matrix = [[dot_gf2(t_row, u_column) for u_column in u_columns] for t_row in t_matrix]
    coordinate_vector = matvec_gf2(t_matrix, node_bits)
    decomposition_rows = active_decomposition_rows(coordinate_vector, u_columns)
    reconstructed_from_dual_basis = decomposition_rows[-1]["partial_sum_after_row"]
    active_coordinate_labels = [row["coordinate_label"] for row in decomposition_rows if row["coordinate_bit"]]
    active_contribution_columns = [row["basis_column"] for row in decomposition_rows if row["coordinate_bit"]]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_DUAL_BASIS_CERTIFICATE__PAIRING_AND_DECOMPOSITION",
        "status": "component-quotient-dual-basis-pairing-and-node-decomposition-certified-over-gf2",
        "source_reports": {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate": str(COORDINATE_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "dual basis",
                "Kronecker pairing",
                "basis decomposition",
                "active coordinate",
                "pairing matrix",
                "sum_i",
                "T_i b",
                "U_i",
            ],
            "finding": "Existing coordinate-isomorphism reports prove T and U are inverse matrices, but not the expanded dual-basis pairing rows <T_i,U_j>=delta_ij or the audited vector decomposition b=sum_i(T_i b)U_i before this file.",
        },
        "dual_basis_definition": {
            "field": "GF(2)",
            "coordinate_functional_rows_T_i": t_matrix,
            "basis_columns_U_i": u_columns,
            "pairing_matrix_T_i_on_U_j": pairing_matrix,
            "coordinate_vector_T_b": coordinate_vector,
            "decomposition_rows": decomposition_rows,
            "active_coordinate_labels": active_coordinate_labels,
            "active_contribution_columns": active_contribution_columns,
            "dual_pairing_statement": "<T_i,U_j>=delta_ij over GF(2)",
            "decomposition_statement": "b=sum_i (T_i b) U_i over GF(2)",
        },
        "rank_pairing_certificate": {
            "rank_T_functionals": rank_gf2(t_matrix),
            "rank_U_basis": rank_gf2(u_matrix),
            "rank_pairing_matrix": rank_gf2(pairing_matrix),
            "pairing_matrix_is_identity": pairing_matrix == identity_matrix(NODE_COUNT),
            "T_functionals_full_rank": rank_gf2(t_matrix) == NODE_COUNT,
            "U_basis_full_rank": rank_gf2(u_matrix) == NODE_COUNT,
        },
        "reconstruction_certificate": {
            "node_bits": node_bits,
            "coordinate_vector_T_b": coordinate_vector,
            "reconstructed_from_dual_basis": reconstructed_from_dual_basis,
            "coordinate_report_coordinate_vector": coordinate["reconstruction_certificate"]["coordinate_vector_T_node_bits"],
            "coordinate_report_reconstructed_node_bits": coordinate["reconstruction_certificate"]["node_bits_from_U_T_node_bits"],
        },
        "component_quotient_dual_basis_summary": {
            "pairing_matrix_is_identity": pairing_matrix == identity_matrix(NODE_COUNT),
            "T_and_U_full_rank": rank_gf2(t_matrix) == NODE_COUNT and rank_gf2(u_matrix) == NODE_COUNT,
            "coordinate_vector_matches_expected": coordinate_vector == EXPECTED_COORDINATE_VECTOR,
            "active_coordinates_are_quotient_components_1_and_3": active_coordinate_labels == ["quotient_component_1", "quotient_component_3"],
            "all_interior_residual_coordinates_zero": coordinate_vector[COMPONENT_COUNT:] == [0] * RESIDUAL_DIMENSION,
            "dual_basis_reconstructs_node_bits": reconstructed_from_dual_basis == node_bits == EXPECTED_NODE_BITS,
            "matches_coordinate_isomorphism_report": coordinate_vector == coordinate["reconstruction_certificate"]["coordinate_vector_T_node_bits"] and reconstructed_from_dual_basis == coordinate["reconstruction_certificate"]["node_bits_from_U_T_node_bits"],
        },
        "blocker_context": {
            "resolved_locally": [
                "The rows T_i and columns U_j are exported as a Kronecker dual basis over GF(2).",
                "The audited vector decomposes as b=sum_i(T_i b)U_i.",
                "Only quotient_component_1 and quotient_component_3 are active; all interior residual coordinates vanish.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "pairing_step": "The exported rows T_i and columns U_j satisfy <T_i,U_j>=delta_ij, so they are dual bases over GF(2).",
            "decomposition_step": "The audited node vector satisfies b=sum_i(T_i b)U_i, with active coordinates quotient_component_1 and quotient_component_3 only.",
            "residual_step": "All seven interior residual coordinates vanish for the audited phase-sign vector.",
            "theoretical_limit": "This is finite dual-basis bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["component_quotient_dual_basis_summary"]
    rank = payload["rank_pairing_certificate"]
    lines = [
        "# Phase-sign component-quotient dual-basis certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate exports the rows T_i and columns U_j as a Kronecker",
        "dual basis and reconstructs the audited node vector from its coordinates.",
        "",
        "## Summary",
        "",
        f"- pairing matrix identity: `{summary['pairing_matrix_is_identity']}`",
        f"- rank pairing matrix: `{rank['rank_pairing_matrix']}`",
        f"- expected coordinates: `{summary['coordinate_vector_matches_expected']}`",
        f"- active quotient coordinates only: `{summary['active_coordinates_are_quotient_components_1_and_3']}`",
        f"- reconstructs node bits: `{summary['dual_basis_reconstructs_node_bits']}`",
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
    print(json.dumps(payload["component_quotient_dual_basis_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
