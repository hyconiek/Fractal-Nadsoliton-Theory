#!/usr/bin/env python3
"""Scratch probe: phase-sign anchored reduced coboundary inverse certificate.

The previous path-cohomology certificate proves that the GF(2) coboundary on the
finite path 0--1--...--11 has only the constant-node kernel and no H^1 cycle
ambiguity.  This probe records the anchored version as an explicit square matrix
certificate: delete the anchor column b(0), obtain the reduced coboundary

    R: C^0_{b(0)=0} -> C^1,

and verify that R is invertible over GF(2).  Its two-sided inverse is exactly the
prefix matrix P, so tail_node_bits = P * edge_bits and R * tail_node_bits =
edge_bits.

This is finite graph linear algebra only.  It does not derive the phase zeros,
omega/phi, damping, or transport from strict nadsoliton dynamics, and it does not
claim a bridge theorem or ToE closure.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_reduced_coboundary_inverse_certificate_report.md"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"
PATH_COHOMOLOGY_REPORT = HERE / "bridge_strict_completion_phase_sign_path_cohomology_triviality_certificate_report.json"
CYCLE_CLOSURE_REPORT = HERE / "bridge_strict_completion_phase_sign_cycle_closure_obstruction_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
ANCHOR_NODE = 0
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]
TAIL_NODE_LABELS = [str(d) for d in range(1, NODE_COUNT)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def identity_matrix(size: int) -> list[list[int]]:
    return [[1 if row == col else 0 for col in range(size)] for row in range(size)]


def matmul_gf2(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    if not left or not right:
        return []
    right_cols = len(right[0])
    return [
        [sum(left[row][k] & right[k][col] for k in range(len(right))) % 2 for col in range(right_cols)]
        for row in range(len(left))
    ]


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def coboundary_matrix() -> list[list[int]]:
    rows = []
    for edge_index in range(EDGE_COUNT):
        row = [0] * NODE_COUNT
        row[edge_index] = 1
        row[edge_index + 1] = 1
        rows.append(row)
    return rows


def reduced_coboundary_matrix(full_delta: list[list[int]]) -> list[list[int]]:
    return [row[1:] for row in full_delta]


def prefix_inverse_matrix() -> list[list[int]]:
    return [[1 if col <= row else 0 for col in range(EDGE_COUNT)] for row in range(EDGE_COUNT)]


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
        "left_block_is_identity": work == identity_matrix(row_count) if row_count == col_count else False,
        "determinant_mod2": 1 if row_count == col_count and len(pivots) == row_count else 0,
    }


def support(edge_bits: list[int]) -> list[str]:
    return [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit]


def reduced_equation_rows(reduced_delta: list[list[int]], tail_node_bits: list[int], edge_bits: list[int]) -> list[dict[str, Any]]:
    computed_edges = matvec_gf2(reduced_delta, tail_node_bits)
    rows = []
    for edge_index, (edge, computed, expected) in enumerate(zip(EDGE_LABELS, computed_edges, edge_bits)):
        left_node = edge_index
        right_node = edge_index + 1
        rows.append({
            "edge": edge,
            "equation": f"b{left_node} xor b{right_node}",
            "anchor_column_removed": left_node == ANCHOR_NODE,
            "reduced_row": reduced_delta[edge_index],
            "computed_edge_bit_from_R_tail": computed,
            "audited_edge_bit": expected,
            "matches_audited_edge_bit": computed == expected,
        })
    return rows


def prefix_inverse_rows(prefix_inverse: list[list[int]], edge_bits: list[int], tail_node_bits: list[int]) -> list[dict[str, Any]]:
    computed_tail = matvec_gf2(prefix_inverse, edge_bits)
    rows = []
    for tail_index, (node_label, computed, expected) in enumerate(zip(TAIL_NODE_LABELS, computed_tail, tail_node_bits)):
        rows.append({
            "tail_node": node_label,
            "prefix_edges": EDGE_LABELS[: tail_index + 1],
            "prefix_inverse_row": prefix_inverse[tail_index],
            "computed_tail_node_bit": computed,
            "audited_tail_node_bit": expected,
            "matches_audited_tail_node_bit": computed == expected,
        })
    return rows


def build_payload() -> dict[str, Any]:
    z2 = load_json(Z2_REPORT)
    gf2 = load_json(GF2_REPORT)
    path_cohomology = load_json(PATH_COHOMOLOGY_REPORT)
    cycle_closure = load_json(CYCLE_CLOSURE_REPORT)

    node_bits = [row["node_bit"] for row in z2["node_bit_rows"]]
    edge_bits = [row["edge_bit"] for row in z2["edge_bit_rows"]]
    tail_node_bits = node_bits[1:]
    full_delta = coboundary_matrix()
    reduced_delta = reduced_coboundary_matrix(full_delta)
    prefix_inverse = prefix_inverse_matrix()
    identity = identity_matrix(EDGE_COUNT)
    r_times_p = matmul_gf2(reduced_delta, prefix_inverse)
    p_times_r = matmul_gf2(prefix_inverse, reduced_delta)
    rref = gf2_rref(reduced_delta)
    tail_from_edges = matvec_gf2(prefix_inverse, edge_bits)
    edges_from_tail = matvec_gf2(reduced_delta, tail_node_bits)
    full_nodes_from_reduced_inverse = [ANCHOR_NODE] + tail_from_edges
    equation_rows = reduced_equation_rows(reduced_delta, tail_node_bits, edge_bits)
    inverse_rows = prefix_inverse_rows(prefix_inverse, edge_bits, tail_node_bits)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_REDUCED_COBOUNDARY_INVERSE_CERTIFICATE__ANCHORED_GF2_ISOMORPHISM",
        "status": "anchored-reduced-coboundary-invertible-over-gf2-prefix-inverse-certified",
        "source_reports": {
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
            "phase_sign_path_cohomology_triviality_certificate": str(PATH_COHOMOLOGY_REPORT.relative_to(ROOT)),
            "phase_sign_cycle_closure_obstruction_certificate": str(CYCLE_CLOSURE_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "reduced coboundary",
                "anchored reduced",
                "relative cohomology",
                "delta_reduced",
                "anchored inverse",
            ],
            "finding": "Existing reports certify prefix reconstruction, path H1=0, and artificial cycle closure, but did not export an anchored reduced coboundary R with an explicit two-sided prefix inverse P.",
        },
        "anchored_reduced_complex_definition": {
            "field": "GF(2)",
            "graph": "finite path 0--1--...--11",
            "anchor_node": ANCHOR_NODE,
            "tail_node_labels": TAIL_NODE_LABELS,
            "edge_labels": EDGE_LABELS,
            "full_coboundary_delta_matrix": full_delta,
            "reduced_coboundary_R_matrix_anchor_column_removed": reduced_delta,
            "prefix_inverse_P_matrix": prefix_inverse,
            "equation": "R * tail_node_bits = edge_bits and tail_node_bits = P * edge_bits",
        },
        "rank_certificate": rref,
        "two_sided_inverse_certificate": {
            "R_times_P": r_times_p,
            "P_times_R": p_times_r,
            "identity_11": identity,
            "R_times_P_is_identity": r_times_p == identity,
            "P_times_R_is_identity": p_times_r == identity,
        },
        "reconstruction_certificate": {
            "audited_node_bits": node_bits,
            "audited_tail_node_bits": tail_node_bits,
            "audited_edge_bits": edge_bits,
            "tail_node_bits_from_P_edge_bits": tail_from_edges,
            "edge_bits_from_R_tail_node_bits": edges_from_tail,
            "full_node_bits_from_anchor_and_reduced_inverse": full_nodes_from_reduced_inverse,
            "flip_edges_from_edge_bits": support(edge_bits),
        },
        "reduced_equation_rows": equation_rows,
        "prefix_inverse_rows": inverse_rows,
        "reduced_coboundary_inverse_summary": {
            "rank_R": rref["rank"],
            "nullity_R": rref["nullity"],
            "determinant_mod2_R": rref["determinant_mod2"],
            "anchored_map_is_isomorphism": rref["rank"] == EDGE_COUNT and rref["nullity"] == 0,
            "two_sided_prefix_inverse_verified": r_times_p == identity and p_times_r == identity,
            "tail_node_bits_match_z2": tail_from_edges == tail_node_bits,
            "full_node_bits_match_z2": full_nodes_from_reduced_inverse == node_bits,
            "edge_bits_match_z2": edges_from_tail == edge_bits,
            "flip_edges": support(edge_bits),
            "matches_expected_node_bits": node_bits == EXPECTED_NODE_BITS,
            "matches_expected_edge_bits": edge_bits == EXPECTED_EDGE_BITS,
            "matches_expected_flip_edges": support(edge_bits) == EXPECTED_FLIP_EDGES,
            "matches_gf2_linear_system_prefix_inverse": prefix_inverse == gf2["linear_system_definition"]["prefix_matrix_L"],
            "matches_gf2_linear_system_solution": edge_bits == gf2["gf2_linear_system_summary"]["solution_edge_bit_pattern"],
            "matches_path_cohomology_anchor_reconstruction": full_nodes_from_reduced_inverse == path_cohomology["path_cohomology_summary"]["anchored_reconstructed_node_bits"],
            "inherits_path_h1_zero": path_cohomology["path_cohomology_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 0,
            "contrasts_cycle_h1_one_boundary_check": cycle_closure["cycle_closure_summary"]["h1_dimension_dim_C1_minus_rank_delta"] == 1,
            "all_reduced_equations_pass": all(row["matches_audited_edge_bit"] for row in equation_rows),
            "all_prefix_inverse_rows_pass": all(row["matches_audited_tail_node_bit"] for row in inverse_rows),
        },
        "blocker_context": {
            "resolved_locally": [
                "The anchor b(0)=0 removes the constant-node kernel on the finite path.",
                "The reduced coboundary R is a square invertible 11x11 GF(2) matrix.",
                "The prefix matrix P is certified as the exact two-sided inverse of R.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "anchor_step": "Delete the anchor column b(0) from delta:C0->C1 to form the reduced coboundary R on tail node bits b(1)..b(11).",
            "rank_step": "R has rank(R)=11, nullity(R)=0, and determinant 1 over GF(2).",
            "inverse_step": "The prefix matrix P satisfies R*P=I_11 and P*R=I_11, so P is a two-sided inverse for the anchored coboundary.",
            "reconstruction_step": "P*edge_bits recovers the audited tail node bits, and R*tail_node_bits recovers the audited edge bits with flip edges 1->2, 5->6, 7->8, and 9->10.",
            "cohomology_boundary_step": "This refines the path H1=0 result after choosing b(0)=0 and contrasts with the artificial closed-cycle H1=1 sanity check.",
            "theoretical_limit": "This is finite anchored GF(2) bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["reduced_coboundary_inverse_summary"]
    lines = [
        "# Phase-sign anchored reduced coboundary inverse certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate turns the path cohomology anchor into an explicit square",
        "GF(2) matrix statement.  Removing the anchor column b(0) gives an",
        "11x11 reduced coboundary R, and the prefix matrix P is verified as its",
        "two-sided inverse.",
        "",
        "## Summary",
        "",
        f"- Rank of R: `{summary['rank_R']}`",
        f"- Nullity of R: `{summary['nullity_R']}`",
        f"- Determinant mod 2 of R: `{summary['determinant_mod2_R']}`",
        f"- Two-sided prefix inverse verified: `{summary['two_sided_prefix_inverse_verified']}`",
        f"- Node bits recovered: `{payload['reconstruction_certificate']['full_node_bits_from_anchor_and_reduced_inverse']}`",
        f"- Edge bits recovered: `{payload['reconstruction_certificate']['edge_bits_from_R_tail_node_bits']}`",
        f"- Flip edges: `{summary['flip_edges']}`",
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
    print(json.dumps(payload["reduced_coboundary_inverse_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
