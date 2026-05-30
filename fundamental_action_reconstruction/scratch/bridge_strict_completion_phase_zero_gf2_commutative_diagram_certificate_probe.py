#!/usr/bin/env python3
"""Scratch probe: phase-zero GF(2) commutative diagram certificate.

The recent finite certificates define the same phase-sign bookkeeping through
several maps: zero-carriers to edges, edges to nodes, and zero-carriers directly
to nodes.  This probe audits the whole diagram over GF(2):

    C_tail = L*M,
    D*C_tail = M,
    B*C_full = M,
    C_full*1 = node_bits,
    M*1 = edge_bits,
    B*node_bits = edge_bits.

Here M is carrier/edge incidence, L is prefix summation, D is first difference,
B is adjacent node boundary, and C is the carrier-prefix node matrix.  This is a
finite commutative-diagram proof check only; it is not a new zero-location proof,
not a bridge theorem, and not a strict dynamical derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_gf2_commutative_diagram_certificate_report.md"
CARRIER_EDGE_REPORT = HERE / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_report.json"
CARRIER_PREFIX_REPORT = HERE / "bridge_strict_completion_phase_zero_carrier_prefix_node_matrix_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"

NODE_COUNT = 12
EDGE_COUNT = 11
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(EDGE_COUNT)]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def matmul_gf2(left: list[list[int]], right: list[list[int]]) -> list[list[int]]:
    if not left or not right:
        return []
    width = len(right[0])
    inner = len(right)
    return [
        [sum(left_row[k] & right[k][col] for k in range(inner)) % 2 for col in range(width)]
        for left_row in left
    ]


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def identity(size: int) -> list[list[int]]:
    return [[1 if row == col else 0 for col in range(size)] for row in range(size)]


def adjacent_boundary_matrix() -> list[list[int]]:
    rows = []
    for edge_index in range(EDGE_COUNT):
        row = [0] * NODE_COUNT
        row[edge_index] = 1
        row[edge_index + 1] = 1
        rows.append(row)
    return rows


def support(edge_bits: list[int]) -> list[str]:
    return [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit]


def diagram_check_rows(checks: dict[str, bool]) -> list[dict[str, Any]]:
    descriptions = {
        "C_tail_equals_L_times_M": "carrier-prefix tail matrix equals prefix matrix times carrier/edge incidence",
        "D_times_C_tail_equals_M": "first difference of carrier-prefix tail recovers carrier/edge incidence",
        "B_times_C_full_equals_M": "adjacent node boundary of full carrier-prefix matrix recovers carrier/edge incidence",
        "C_full_times_one_equals_node_bits": "carrier-prefix full matrix maps all-one carrier vector to node bits",
        "M_times_one_equals_edge_bits": "carrier/edge incidence maps all-one carrier vector to edge bits",
        "B_times_node_bits_equals_edge_bits": "adjacent node boundary maps node bits to edge bits",
        "D_times_L_is_identity": "first-difference matrix is a left inverse of prefix matrix",
        "L_times_D_is_identity": "first-difference matrix is a right inverse of prefix matrix",
    }
    return [
        {"check": key, "description": descriptions[key], "passes": value}
        for key, value in checks.items()
    ]


def edge_rows(edge_bits_from_boundary: list[int], expected_edge_bits: list[int]) -> list[dict[str, Any]]:
    return [
        {
            "edge": edge,
            "boundary_edge_bit": bit,
            "expected_edge_bit": expected,
            "is_flip_edge": bit == 1,
            "matches_expected_edge_bit": bit == expected,
        }
        for edge, bit, expected in zip(EDGE_LABELS, edge_bits_from_boundary, expected_edge_bits)
    ]


def build_payload() -> dict[str, Any]:
    carrier_edge = load_json(CARRIER_EDGE_REPORT)
    carrier_prefix = load_json(CARRIER_PREFIX_REPORT)
    gf2 = load_json(GF2_REPORT)
    z2 = load_json(Z2_REPORT)

    m_matrix = carrier_edge["carrier_edge_incidence_matrix"]
    carrier_vector = carrier_edge["incidence_definition"]["carrier_multiplicity_vector"]
    l_matrix = gf2["linear_system_definition"]["prefix_matrix_L"]
    d_matrix = gf2["linear_system_definition"]["explicit_inverse_first_difference_matrix"]
    c_tail = carrier_prefix["matrix_definition"]["carrier_prefix_tail_matrix_C_tail_equals_LM"]
    c_full = carrier_prefix["matrix_definition"]["carrier_prefix_full_node_matrix_C"]
    b_matrix = adjacent_boundary_matrix()
    node_bits = carrier_prefix["carrier_prefix_node_matrix_summary"]["node_bit_pattern_from_carrier_prefix_matrix"]
    edge_bits = carrier_edge["carrier_edge_incidence_summary"]["edge_bit_pattern_from_carrier_incidence"]

    l_times_m = matmul_gf2(l_matrix, m_matrix)
    d_times_c_tail = matmul_gf2(d_matrix, c_tail)
    b_times_c_full = matmul_gf2(b_matrix, c_full)
    c_full_times_one = matvec_gf2(c_full, carrier_vector)
    m_times_one = matvec_gf2(m_matrix, carrier_vector)
    b_times_node_bits = matvec_gf2(b_matrix, node_bits)
    d_times_l = matmul_gf2(d_matrix, l_matrix)
    l_times_d = matmul_gf2(l_matrix, d_matrix)

    checks = {
        "C_tail_equals_L_times_M": c_tail == l_times_m,
        "D_times_C_tail_equals_M": d_times_c_tail == m_matrix,
        "B_times_C_full_equals_M": b_times_c_full == m_matrix,
        "C_full_times_one_equals_node_bits": c_full_times_one == node_bits,
        "M_times_one_equals_edge_bits": m_times_one == edge_bits,
        "B_times_node_bits_equals_edge_bits": b_times_node_bits == edge_bits,
        "D_times_L_is_identity": d_times_l == identity(EDGE_COUNT),
        "L_times_D_is_identity": l_times_d == identity(EDGE_COUNT),
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_GF2_COMMUTATIVE_DIAGRAM_CERTIFICATE__CARRIER_EDGE_NODE_MAPS",
        "status": "carrier-edge-node-gf2-diagram-commutes-for-audited-phase-sign-data",
        "source_reports": {
            "carrier_edge_incidence_certificate": str(CARRIER_EDGE_REPORT.relative_to(ROOT)),
            "carrier_prefix_node_matrix_certificate": str(CARRIER_PREFIX_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "commutative diagram",
                "exact sequence",
                "chain complex",
                "B*C",
                "D*C",
                "D*L",
                "carrier edge node diagram",
            ],
            "conclusion": "No strict-completion report exported the full GF(2) carrier-edge-node commutative diagram before this file.",
        },
        "diagram_definition": {
            "field": "GF(2)",
            "M_carrier_edge_incidence": m_matrix,
            "L_prefix_matrix": l_matrix,
            "D_first_difference_matrix": d_matrix,
            "B_adjacent_node_boundary_matrix": b_matrix,
            "C_tail_carrier_prefix_matrix": c_tail,
            "C_full_node_carrier_prefix_matrix": c_full,
            "carrier_multiplicity_vector": carrier_vector,
            "equations_checked": list(checks),
        },
        "diagram_check_rows": diagram_check_rows(checks),
        "edge_boundary_rows": edge_rows(b_times_node_bits, edge_bits),
        "diagram_summary": {
            "all_diagram_checks_pass": all(checks.values()),
            "edge_bits_from_boundary": b_times_node_bits,
            "node_bits_from_carrier_prefix": node_bits,
            "flip_edges_from_boundary": support(b_times_node_bits),
            "matches_expected_edge_bits": b_times_node_bits == EXPECTED_EDGE_BITS,
            "matches_expected_node_bits": node_bits == EXPECTED_NODE_BITS,
            "matches_expected_flip_edges": support(b_times_node_bits) == EXPECTED_FLIP_EDGES,
            "matches_z2_node_bits": node_bits == [row["node_bit"] for row in z2["node_bit_rows"]],
            "matches_z2_edge_bits": b_times_node_bits == [row["edge_bit"] for row in z2["edge_bit_rows"]],
            "matches_carrier_edge_incidence_edge_bits": b_times_node_bits == edge_bits,
            "matches_gf2_solution_edge_bits": b_times_node_bits == gf2["gf2_linear_system_summary"]["solution_edge_bit_pattern"],
            "inherits_carrier_prefix_rank_4": carrier_prefix["rank_certificate"]["carrier_prefix_rank_over_gf2"] == 4,
            "inherits_carrier_edge_rank_4": carrier_edge["rank_certificate"]["column_rank_over_gf2"] == 4,
        },
        "blocker_context": {
            "what_this_certifies": "finite GF(2) commutative diagram tying carrier, edge, and node phase-sign bookkeeping together",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "strict_damping_parameter_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "QW-2191_selector_obstruction",
            ],
        },
        "proof_certificate": {
            "prefix_composition_step": "C_tail=L*M verifies that carrier-prefix rows are prefix sums of carrier/edge incidence rows.",
            "difference_step": "D*C_tail=M and B*C_full=M verify that first differences recover the carrier/edge incidence matrix.",
            "vector_step": "C_full*1 gives node bits, M*1 gives edge bits, and B*node_bits gives the same edge bits.",
            "inverse_step": "D*L and L*D both equal the 11x11 identity over GF(2).",
            "theoretical_limit": "This proves finite diagram commutativity only; it does not derive zero carriers, omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["diagram_summary"]
    lines = [
        "# Phase-zero GF(2) commutative diagram certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate checks the finite carrier-edge-node maps as one GF(2)",
        "commutative diagram: `C_tail=L*M`, first differences recover `M`, and",
        "all-one carrier/node vectors recover the audited node and edge bits.",
        "",
        "## Summary",
        "",
        f"- All diagram checks pass: `{summary['all_diagram_checks_pass']}`",
        f"- Node bits from carrier prefix: `{summary['node_bits_from_carrier_prefix']}`",
        f"- Edge bits from boundary: `{summary['edge_bits_from_boundary']}`",
        f"- Flip edges from boundary: `{summary['flip_edges_from_boundary']}`",
        f"- Matches Z2 node bits: `{summary['matches_z2_node_bits']}`",
        f"- Matches Z2 edge bits: `{summary['matches_z2_edge_bits']}`",
        f"- Inherits carrier-prefix rank 4: `{summary['inherits_carrier_prefix_rank_4']}`",
        f"- Inherits carrier-edge rank 4: `{summary['inherits_carrier_edge_rank_4']}`",
        "",
        "## Diagram checks",
        "",
        "| check | passes |",
        "| --- | :---: |",
    ]
    for row in payload["diagram_check_rows"]:
        lines.append(f"| `{row['check']}` | `{row['passes']}` |")
    lines.extend(["", "## Hard limits", ""])
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
