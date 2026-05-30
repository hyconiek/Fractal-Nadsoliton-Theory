#!/usr/bin/env python3
"""Scratch probe: phase-zero carrier-prefix node-matrix certificate.

The carrier/edge incidence certificate maps ordered zero-carriers to edge bits,
and the GF(2) prefix linear-system certificate maps edge bits to node bits.  This
probe composes those two finite maps:

    node_tail_bits = prefix_matrix_L * carrier_edge_incidence_M * 1  over GF(2).

Equivalently, each node row records exactly which rational zero-carriers lie
strictly to the left of that integer node.  The resulting carrier-prefix matrix
recovers the cell-sign node-bit pattern and has rank four on the carrier
subspace.  This is finite matrix bookkeeping only; it is not a new zero-location
proof, not a bridge theorem, and not a strict dynamical derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_carrier_prefix_node_matrix_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_carrier_prefix_node_matrix_certificate_report.md"
CARRIER_EDGE_REPORT = HERE / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_report.json"
CELL_SIGN_REPORT = HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"

NODE_LABELS = list(range(12))
EDGE_LABELS = [f"{d}->{d + 1}" for d in range(11)]
EXPECTED_CARRIER_ORDER = ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"]
EXPECTED_NODE_BITS = [0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0]
EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def prefix_matrix_tail() -> list[list[int]]:
    return [[1 if col < d else 0 for col in range(11)] for d in range(1, 12)]


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


def gf2_rank(matrix: list[list[int]]) -> tuple[int, list[dict[str, int]]]:
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
    return len(pivots), pivots


def signs_from_bits(bits: list[int]) -> list[int]:
    return [1 if bit == 0 else -1 for bit in bits]


def node_rows(carrier_order: list[str], full_prefix_matrix: list[list[int]], node_bits: list[int], cell_sign_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    cell_by_d = {row["d"]: row for row in cell_sign_rows}
    for d, row, bit in zip(NODE_LABELS, full_prefix_matrix, node_bits):
        carriers_left = [carrier for carrier, value in zip(carrier_order, row) if value]
        cell_row = cell_by_d[d]
        rows.append({
            "d": d,
            "carrier_prefix_row": row,
            "carriers_left_by_matrix": carriers_left,
            "carrier_count_left_mod2": sum(row) % 2,
            "node_bit_from_matrix": bit,
            "phase_sign_from_matrix": 1 if bit == 0 else -1,
            "cell_sign_zero_carriers_left": cell_row["zero_carriers_left_of_node"],
            "matches_cell_sign_left_carriers": carriers_left == cell_row["zero_carriers_left_of_node"],
            "matches_cell_sign_node_bit": bit == (0 if cell_row["derived_phase_transport_sign"] == 1 else 1),
        })
    return rows


def edge_difference_rows(full_prefix_matrix: list[list[int]], edge_matrix: list[list[int]]) -> list[dict[str, Any]]:
    rows = []
    for edge_index, edge in enumerate(EDGE_LABELS):
        difference = [a ^ b for a, b in zip(full_prefix_matrix[edge_index + 1], full_prefix_matrix[edge_index])]
        rows.append({
            "edge": edge,
            "node_prefix_difference_row": difference,
            "carrier_edge_incidence_row": edge_matrix[edge_index],
            "difference_recovers_edge_incidence": difference == edge_matrix[edge_index],
        })
    return rows


def build_payload() -> dict[str, Any]:
    carrier_edge = load_json(CARRIER_EDGE_REPORT)
    cell_sign = load_json(CELL_SIGN_REPORT)
    z2 = load_json(Z2_REPORT)
    gf2 = load_json(GF2_REPORT)

    carrier_order = carrier_edge["incidence_definition"]["carrier_order"]
    carrier_vector = carrier_edge["incidence_definition"]["carrier_multiplicity_vector"]
    edge_matrix = carrier_edge["carrier_edge_incidence_matrix"]
    prefix_tail = prefix_matrix_tail()
    carrier_prefix_tail = matmul_gf2(prefix_tail, edge_matrix)
    full_prefix_matrix = [[0 for _carrier in carrier_order]] + carrier_prefix_tail
    tail_bits = matvec_gf2(carrier_prefix_tail, carrier_vector)
    node_bits = [0] + tail_bits
    rank, pivots = gf2_rank(full_prefix_matrix)
    rows = node_rows(carrier_order, full_prefix_matrix, node_bits, cell_sign["node_sign_rows"])
    edge_rows = edge_difference_rows(full_prefix_matrix, edge_matrix)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CARRIER_PREFIX_NODE_MATRIX_CERTIFICATE__GF2_CARRIER_TO_NODE_PARITY",
        "status": "carrier-prefix-node-matrix-recovers-cell-sign-and-z2-node-bits",
        "source_reports": {
            "carrier_edge_incidence_certificate": str(CARRIER_EDGE_REPORT.relative_to(ROOT)),
            "cell_sign_certificate": str(CELL_SIGN_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "carrier prefix",
                "prefix carrier",
                "carrier-prefix",
                "carrier node matrix",
                "L * M",
                "carrier left of node",
            ],
            "conclusion": "Existing cell-sign rows list carriers left of nodes, but no strict-completion report exported the composed carrier-prefix node matrix C=L*M before this file.",
        },
        "matrix_definition": {
            "field": "GF(2)",
            "carrier_order": carrier_order,
            "node_order": NODE_LABELS,
            "edge_order": EDGE_LABELS,
            "prefix_tail_matrix_L": prefix_tail,
            "carrier_edge_incidence_matrix_M": edge_matrix,
            "carrier_prefix_tail_matrix_C_tail_equals_LM": carrier_prefix_tail,
            "carrier_prefix_full_node_matrix_C": full_prefix_matrix,
            "carrier_multiplicity_vector": carrier_vector,
            "node_bit_rule": "node_bits = [0] + C_tail * carrier_multiplicity_vector mod 2",
        },
        "node_prefix_rows": rows,
        "edge_difference_rows": edge_rows,
        "rank_certificate": {
            "carrier_prefix_rank_over_gf2": rank,
            "pivot_rows": pivots,
            "carrier_count": len(carrier_order),
            "full_column_rank_on_carrier_prefix_subspace": rank == len(carrier_order),
        },
        "carrier_prefix_node_matrix_summary": {
            "carrier_order": carrier_order,
            "node_bit_pattern_from_carrier_prefix_matrix": node_bits,
            "phase_sign_pattern_from_carrier_prefix_matrix": signs_from_bits(node_bits),
            "carrier_prefix_rank_over_gf2": rank,
            "all_node_rows_match_cell_sign_left_carriers": all(row["matches_cell_sign_left_carriers"] for row in rows),
            "all_node_bits_match_cell_sign": all(row["matches_cell_sign_node_bit"] for row in rows),
            "all_edge_differences_recover_carrier_edge_incidence": all(row["difference_recovers_edge_incidence"] for row in edge_rows),
            "matches_expected_carrier_order": carrier_order == EXPECTED_CARRIER_ORDER,
            "matches_expected_node_bits": node_bits == EXPECTED_NODE_BITS,
            "matches_expected_sign_pattern": signs_from_bits(node_bits) == EXPECTED_SIGN_PATTERN,
            "matches_z2_node_bits": node_bits == [row["node_bit"] for row in z2["node_bit_rows"]],
            "matches_gf2_prefix_solution": node_bits[1:] == gf2["linear_system_definition"]["rhs_node_bits_minus_anchor"],
            "inherits_carrier_edge_incidence_rank": carrier_edge["rank_certificate"]["full_column_rank_on_carrier_subspace"],
        },
        "blocker_context": {
            "what_this_certifies": "finite GF(2) composition from rational zero-carrier incidence to integer node sign bits",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "strict_damping_parameter_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "QW-2191_selector_obstruction",
            ],
        },
        "proof_certificate": {
            "composition_step": "Compose the prefix matrix L with the carrier/edge incidence matrix M over GF(2), giving C=L*M.",
            "node_step": "Multiplying C by the all-one carrier vector gives exactly the audited node-bit pattern from d=1..11, with b(0)=0 prepended.",
            "cell_sign_step": "Every row of C lists the same carriers-left-of-node set as the rational cell-sign certificate.",
            "edge_difference_step": "Adjacent node-prefix row differences recover the carrier/edge incidence rows, so edge flips remain the same four carrier crossings.",
            "rank_step": "The carrier-prefix node matrix has GF(2) column rank 4 on the four-carrier subspace.",
            "theoretical_limit": "This proves finite matrix composition only; it does not derive zero carriers, omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["carrier_prefix_node_matrix_summary"]
    lines = [
        "# Phase-zero carrier-prefix node-matrix certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate composes the GF(2) prefix matrix with the carrier/edge",
        "incidence matrix, producing a carrier-prefix node matrix `C=L*M` whose",
        "rows list the zero-carriers left of each integer node.",
        "",
        "## Summary",
        "",
        f"- Carrier order: `{summary['carrier_order']}`",
        f"- Node bits from carrier-prefix matrix: `{summary['node_bit_pattern_from_carrier_prefix_matrix']}`",
        f"- Phase signs from carrier-prefix matrix: `{summary['phase_sign_pattern_from_carrier_prefix_matrix']}`",
        f"- Carrier-prefix rank over GF(2): `{summary['carrier_prefix_rank_over_gf2']}`",
        f"- All node rows match cell-sign left carriers: `{summary['all_node_rows_match_cell_sign_left_carriers']}`",
        f"- All node bits match cell-sign: `{summary['all_node_bits_match_cell_sign']}`",
        f"- Edge differences recover carrier/edge incidence: `{summary['all_edge_differences_recover_carrier_edge_incidence']}`",
        f"- Matches Z2 node bits: `{summary['matches_z2_node_bits']}`",
        "",
        "## Node prefix rows",
        "",
        "| d | carriers left by matrix | node bit | sign | matches cell-sign |",
        "| ---: | --- | ---: | ---: | :---: |",
    ]
    for row in payload["node_prefix_rows"]:
        lines.append(
            f"| {row['d']} | {row['carriers_left_by_matrix']} | {row['node_bit_from_matrix']} | {row['phase_sign_from_matrix']} | `{row['matches_cell_sign_left_carriers'] and row['matches_cell_sign_node_bit']}` |"
        )
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
