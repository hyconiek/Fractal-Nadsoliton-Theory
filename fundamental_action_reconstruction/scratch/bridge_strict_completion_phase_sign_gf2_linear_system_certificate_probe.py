#!/usr/bin/env python3
"""Scratch probe: phase-sign GF(2) linear-system certificate.

The Z2 coboundary and edge-support minimality reports certify the phase signs as
finite path data.  This probe adds a linear-algebra certificate over GF(2): the
prefix reconstruction equations form a unit lower-triangular 11x11 system

    b(d) xor b(0) = sum_{i=0}^{d-1} e(i,i+1)  mod 2,  d=1..11.

The matrix is full-rank with determinant 1 over GF(2), has an explicit inverse
(first differences), and its unique solution is exactly the four-edge phase-flip
support.  This is finite algebraic bookkeeping only; it is not a new phase fit,
not a bridge theorem, and not a strict dynamical derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
EDGE_SUPPORT_REPORT = HERE / "bridge_strict_completion_phase_sign_edge_support_minimality_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = NODE_COUNT - 1
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def prefix_matrix() -> list[list[int]]:
    return [[1 if col <= row else 0 for col in range(EDGE_COUNT)] for row in range(EDGE_COUNT)]


def first_difference_inverse() -> list[list[int]]:
    rows = []
    for row in range(EDGE_COUNT):
        values = [0] * EDGE_COUNT
        values[row] = 1
        if row > 0:
            values[row - 1] = 1
        rows.append(values)
    return rows


def matmul_gf2(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    width = len(right[0])
    inner = len(right)
    return [
        [sum(left_row[k] & right[k][col] for k in range(inner)) % 2 for col in range(width)]
        for left_row in left
    ]


def identity(size: int) -> list[list[int]]:
    return [[1 if row == col else 0 for col in range(size)] for row in range(size)]


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def gauss_jordan_gf2(matrix: list[list[int]], rhs: list[int]) -> dict[str, Any]:
    augmented = [row[:] + [rhs_value] for row, rhs_value in zip(matrix, rhs)]
    row_count = len(augmented)
    col_count = len(matrix[0])
    pivot_rows = []
    operations = []
    pivot_row = 0

    for col in range(col_count):
        pivot = None
        for candidate in range(pivot_row, row_count):
            if augmented[candidate][col] == 1:
                pivot = candidate
                break
        if pivot is None:
            continue
        if pivot != pivot_row:
            augmented[pivot_row], augmented[pivot] = augmented[pivot], augmented[pivot_row]
            operations.append({
                "operation": "swap_rows",
                "row_a": pivot_row,
                "row_b": pivot,
                "pivot_col": col,
            })
        operations.append({
            "operation": "pivot",
            "row": pivot_row,
            "pivot_col": col,
        })
        for target in range(row_count):
            if target != pivot_row and augmented[target][col] == 1:
                augmented[target] = [a ^ b for a, b in zip(augmented[target], augmented[pivot_row])]
                operations.append({
                    "operation": "row_xor",
                    "target_row": target,
                    "source_row": pivot_row,
                    "pivot_col": col,
                })
        pivot_rows.append({"row": pivot_row, "pivot_col": col})
        pivot_row += 1
        if pivot_row == row_count:
            break

    return {
        "rref_augmented_matrix": augmented,
        "rank": len(pivot_rows),
        "nullity": col_count - len(pivot_rows),
        "pivot_rows": pivot_rows,
        "row_operations": operations,
        "solution": [row[-1] for row in augmented],
        "left_block_is_identity": [row[:-1] for row in augmented] == identity(col_count),
    }


def support(edge_bits: list[int]) -> list[str]:
    return [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit == 1]


def equation_rows(matrix: list[list[int]], rhs: list[int], solution: list[int]) -> list[dict[str, Any]]:
    evaluated = matvec_gf2(matrix, solution)
    rows = []
    for row_index, row in enumerate(matrix):
        prefix_edges = [EDGE_LABELS[col] for col, value in enumerate(row) if value]
        rows.append({
            "node_d": row_index + 1,
            "prefix_edges": prefix_edges,
            "rhs_node_bit_xor_anchor": rhs[row_index],
            "evaluated_prefix_sum_mod2": evaluated[row_index],
            "equation_passes": evaluated[row_index] == rhs[row_index],
        })
    return rows


def solution_rows(solution: list[int], z2_edge_bits: list[int]) -> list[dict[str, Any]]:
    return [
        {
            "edge": edge,
            "solution_edge_bit": bit,
            "z2_edge_bit": z2_bit,
            "is_flip_edge": bit == 1,
            "matches_z2_edge_bit": bit == z2_bit,
        }
        for edge, bit, z2_bit in zip(EDGE_LABELS, solution, z2_edge_bits)
    ]


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    edge_support = load_json(EDGE_SUPPORT_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    z2_edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    anchor = node_bits[0]
    rhs = [bit ^ anchor for bit in node_bits[1:]]
    matrix = prefix_matrix()
    inverse = first_difference_inverse()
    elimination = gauss_jordan_gf2(matrix, rhs)
    solution = elimination["solution"]
    inverse_solution = matvec_gf2(inverse, rhs)
    matrix_times_inverse = matmul_gf2(matrix, inverse)
    inverse_times_matrix = matmul_gf2(inverse, matrix)
    identity_11 = identity(EDGE_COUNT)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_GF2_LINEAR_SYSTEM_CERTIFICATE__FULL_RANK_PREFIX_RECONSTRUCTION",
        "status": "gf2-prefix-system-full-rank-with-unique-four-edge-phase-flip-solution",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_edge_support_minimality_certificate": str(EDGE_SUPPORT_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "incidence rank",
                "GF(2) linear system",
                "row reduction",
                "prefix matrix",
                "phase_sign_gf2",
            ],
            "conclusion": "No existing strict-completion probe exported a GF(2) full-rank prefix linear-system certificate before this file.",
        },
        "linear_system_definition": {
            "field": "GF(2)",
            "unknowns": EDGE_LABELS,
            "equations": "b(d) xor b(0) = sum_{i=0}^{d-1} e(i,i+1) mod 2 for d=1..11",
            "anchor_bit_b0": anchor,
            "prefix_matrix_L": matrix,
            "rhs_node_bits_minus_anchor": rhs,
            "explicit_inverse_first_difference_matrix": inverse,
        },
        "gf2_row_reduction_certificate": elimination,
        "equation_check_rows": equation_rows(matrix, rhs, solution),
        "solution_edge_rows": solution_rows(solution, z2_edge_bits),
        "inverse_checks": {
            "L_times_inverse_is_identity": matrix_times_inverse == identity_11,
            "inverse_times_L_is_identity": inverse_times_matrix == identity_11,
            "inverse_solution_edge_bits": inverse_solution,
            "inverse_solution_matches_row_reduction_solution": inverse_solution == solution,
        },
        "gf2_linear_system_summary": {
            "node_bit_pattern": node_bits,
            "solution_edge_bit_pattern": solution,
            "solution_flip_edges": support(solution),
            "rank": elimination["rank"],
            "nullity": elimination["nullity"],
            "determinant_mod2": 1 if elimination["rank"] == EDGE_COUNT else 0,
            "left_block_reduces_to_identity": elimination["left_block_is_identity"],
            "unique_solution": elimination["rank"] == EDGE_COUNT and elimination["left_block_is_identity"],
            "solution_hamming_weight": sum(solution),
            "all_equations_pass": all(row["equation_passes"] for row in equation_rows(matrix, rhs, solution)),
            "all_solution_edges_match_z2": all(row["matches_z2_edge_bit"] for row in solution_rows(solution, z2_edge_bits)),
            "matches_expected_node_bits": node_bits == EXPECTED_NODE_BITS,
            "matches_expected_edge_bits": solution == EXPECTED_EDGE_BITS,
            "matches_expected_flip_edges": support(solution) == EXPECTED_FLIP_EDGES,
            "matches_edge_support_minimality_solution": solution == edge_support["matching_edge_assignment_rows"][0]["edge_bits"],
            "inherits_edge_support_minimality": edge_support["edge_support_minimality_summary"]["unique_matching_assignment"] and edge_support["edge_support_minimality_summary"]["all_lower_support_assignments_fail"],
        },
        "blocker_context": {
            "what_this_certifies": "full-rank GF(2) linear algebra for finite prefix reconstruction of the already exported phase-sign edge bits",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "strict_damping_parameter_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "QW-2191_selector_obstruction",
            ],
        },
        "proof_certificate": {
            "system_step": "The prefix reconstruction laws are written as L e = b_tail xor b0 over GF(2).",
            "rank_step": "Gaussian elimination over GF(2) reduces L to the 11x11 identity, so rank(L)=11 and nullity(L)=0.",
            "inverse_step": "The first-difference matrix is verified as a two-sided inverse of L over GF(2).",
            "solution_step": "The unique solution is the Z2 edge-bit pattern with support 1->2, 5->6, 7->8, and 9->10.",
            "theoretical_limit": "This proves finite linear-algebra uniqueness only; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["gf2_linear_system_summary"]
    lines = [
        "# Phase-sign GF(2) linear-system certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate writes finite prefix reconstruction as a GF(2) linear",
        "system `L e = b_tail xor b0`.  The prefix matrix is full-rank, has an",
        "explicit first-difference inverse, and has exactly the audited four-flip",
        "solution.",
        "",
        "## Summary",
        "",
        f"- Rank: `{summary['rank']}`",
        f"- Nullity: `{summary['nullity']}`",
        f"- Determinant mod 2: `{summary['determinant_mod2']}`",
        f"- Unique solution: `{summary['unique_solution']}`",
        f"- Solution Hamming weight: `{summary['solution_hamming_weight']}`",
        f"- Flip edges: `{summary['solution_flip_edges']}`",
        f"- All equations pass: `{summary['all_equations_pass']}`",
        f"- Inherits edge-support minimality: `{summary['inherits_edge_support_minimality']}`",
        "",
        "## Equation rows",
        "",
        "| node d | rhs | evaluated | passes |",
        "| ---: | ---: | ---: | :---: |",
    ]
    for row in payload["equation_check_rows"]:
        lines.append(
            f"| {row['node_d']} | {row['rhs_node_bit_xor_anchor']} | {row['evaluated_prefix_sum_mod2']} | `{row['equation_passes']}` |"
        )
    lines.extend([
        "",
        "## Hard limits",
        "",
    ])
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
