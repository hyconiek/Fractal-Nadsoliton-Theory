#!/usr/bin/env python3
"""Scratch probe: phase-zero carrier/edge incidence certificate.

The cell-partition and cell-sign reports already say which zero-carriers are
crossed by each integer edge.  This probe records the same fact as a finite
GF(2) incidence map from ordered zero-carriers to integer edges:

    edge_bit = carrier_edge_incidence_matrix * carrier_multiplicity_vector mod 2.

Each ordered carrier is strictly contained in one open integer edge by rational
boundary-clearance data.  The resulting 11x4 incidence matrix has four standard
basis columns at the four flip edges, so its carrier-column rank is four and its
parity image is exactly the audited Z2 edge-bit vector.  This is finite carrier
bookkeeping only; it is not a new zero-location proof, not a bridge theorem, and
not a strict dynamical derivation.
"""
from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_zero_carrier_edge_incidence_certificate_report.md"
CELL_PARTITION_REPORT = HERE / "bridge_strict_completion_phase_zero_cell_partition_certificate_report.json"
CELL_SIGN_REPORT = HERE / "bridge_strict_completion_phase_zero_cell_sign_certificate_report.json"
Z2_REPORT = HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json"
GF2_REPORT = HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json"

EDGE_LABELS = [f"{d}->{d + 1}" for d in range(11)]
EXPECTED_CARRIER_ORDER = ["legacy_z0", "legacy_z1", "strict_z0", "legacy_z2"]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EXPECTED_EDGE_BITS = [0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 0]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def fraction_from_payload(payload: dict[str, Any]) -> Fraction:
    return Fraction(payload["numerator"], payload["denominator"])


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "decimal": float(value),
        "text": f"{value.numerator}/{value.denominator}",
    }


def edge_index(edge: str) -> int:
    return EDGE_LABELS.index(edge)


def carrier_containment_rows(carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for carrier in carriers:
        edge = carrier["edge_or_domain_location"]
        left_node, right_node = [int(part) for part in edge.split("->")]
        lower = fraction_from_payload(carrier["lower"])
        upper = fraction_from_payload(carrier["upper"])
        left_clearance = lower - Fraction(left_node, 1)
        right_clearance = Fraction(right_node, 1) - upper
        rows.append({
            "carrier_label": carrier["label"],
            "source": carrier["source"],
            "edge": edge,
            "lower": carrier["lower"],
            "upper": carrier["upper"],
            "left_boundary_clearance": fraction_payload(left_clearance),
            "right_boundary_clearance": fraction_payload(right_clearance),
            "strictly_inside_open_edge_by_rational_bounds": left_clearance > 0 and right_clearance > 0,
        })
    return rows


def incidence_matrix(carriers: list[dict[str, Any]]) -> list[list[int]]:
    matrix = [[0 for _carrier in carriers] for _edge in EDGE_LABELS]
    for col, carrier in enumerate(carriers):
        matrix[edge_index(carrier["edge_or_domain_location"])][col] = 1
    return matrix


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def gf2_rank(matrix: list[list[int]]) -> tuple[int, list[dict[str, Any]]]:
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


def edge_incidence_rows(matrix: list[list[int]], carriers: list[dict[str, Any]], edge_bits: list[int]) -> list[dict[str, Any]]:
    rows = []
    for edge, row, bit in zip(EDGE_LABELS, matrix, edge_bits):
        incident = [carrier["label"] for carrier, value in zip(carriers, row) if value]
        rows.append({
            "edge": edge,
            "incident_zero_carriers": incident,
            "incident_carrier_count": len(incident),
            "incidence_row": row,
            "edge_bit_from_incidence": bit,
            "odd_carrier_parity": len(incident) % 2,
            "edge_bit_equals_odd_carrier_parity": bit == len(incident) % 2,
        })
    return rows


def carrier_column_rows(matrix: list[list[int]], carriers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = []
    for col, carrier in enumerate(carriers):
        column = [row[col] for row in matrix]
        support_edges = [edge for edge, value in zip(EDGE_LABELS, column) if value]
        rows.append({
            "carrier_label": carrier["label"],
            "column": column,
            "support_edges": support_edges,
            "column_weight": sum(column),
            "is_single_edge_standard_basis_column": sum(column) == 1,
            "matches_declared_carrier_edge": support_edges == [carrier["edge_or_domain_location"]],
        })
    return rows


def build_payload() -> dict[str, Any]:
    cell_partition = load_json(CELL_PARTITION_REPORT)
    cell_sign = load_json(CELL_SIGN_REPORT)
    z2 = load_json(Z2_REPORT)
    gf2 = load_json(GF2_REPORT)

    carriers = cell_partition["domain_zero_carriers_ordered"]
    carrier_vector = [1 for _carrier in carriers]
    matrix = incidence_matrix(carriers)
    edge_bits = matvec_gf2(matrix, carrier_vector)
    rank, pivots = gf2_rank(matrix)
    edge_rows = edge_incidence_rows(matrix, carriers, edge_bits)
    column_rows = carrier_column_rows(matrix, carriers)
    containment_rows = carrier_containment_rows(carriers)
    flip_edges = [edge for edge, bit in zip(EDGE_LABELS, edge_bits) if bit]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_ZERO_CARRIER_EDGE_INCIDENCE_CERTIFICATE__GF2_CARRIER_TO_EDGE_PARITY",
        "status": "zero-carrier-to-edge-incidence-map-reconstructs-phase-z2-edge-bits",
        "source_reports": {
            "cell_partition_certificate": str(CELL_PARTITION_REPORT.relative_to(ROOT)),
            "cell_sign_certificate": str(CELL_SIGN_REPORT.relative_to(ROOT)),
            "phase_sign_z2_coboundary_certificate": str(Z2_REPORT.relative_to(ROOT)),
            "phase_sign_gf2_linear_system_certificate": str(GF2_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "carrier edge incidence",
                "zero-carrier edge incidence",
                "carrier_edge_incidence_matrix",
                "GF(2) carrier",
                "crossing matrix",
            ],
            "conclusion": "Existing reports list crossed carriers, but no strict-completion report exported the carrier-to-edge GF(2) incidence matrix and rank certificate before this file.",
        },
        "incidence_definition": {
            "field": "GF(2)",
            "edge_order": EDGE_LABELS,
            "carrier_order": [carrier["label"] for carrier in carriers],
            "carrier_multiplicity_vector": carrier_vector,
            "incidence_rule": "M[edge,carrier]=1 iff the rational zero-carrier interval is strictly contained in that open integer edge",
            "edge_bit_rule": "edge_bits = M * carrier_multiplicity_vector mod 2",
        },
        "carrier_rational_containment_rows": containment_rows,
        "carrier_edge_incidence_matrix": matrix,
        "carrier_column_rows": column_rows,
        "edge_incidence_rows": edge_rows,
        "rank_certificate": {
            "column_rank_over_gf2": rank,
            "pivot_rows": pivots,
            "carrier_count": len(carriers),
            "full_column_rank_on_carrier_subspace": rank == len(carriers),
        },
        "carrier_edge_incidence_summary": {
            "carrier_order": [carrier["label"] for carrier in carriers],
            "edge_bit_pattern_from_carrier_incidence": edge_bits,
            "derived_flip_edges_from_carrier_incidence": flip_edges,
            "column_rank_over_gf2": rank,
            "all_carriers_strictly_inside_open_edges_by_rational_bounds": all(row["strictly_inside_open_edge_by_rational_bounds"] for row in containment_rows),
            "all_carrier_columns_are_single_edge_standard_basis_columns": all(row["is_single_edge_standard_basis_column"] for row in column_rows),
            "all_edges_have_at_most_one_incident_carrier": max(row["incident_carrier_count"] for row in edge_rows) <= 1,
            "all_edge_bits_equal_odd_carrier_parity": all(row["edge_bit_equals_odd_carrier_parity"] for row in edge_rows),
            "matches_expected_carrier_order": [carrier["label"] for carrier in carriers] == EXPECTED_CARRIER_ORDER,
            "matches_expected_edge_bits": edge_bits == EXPECTED_EDGE_BITS,
            "matches_expected_flip_edges": flip_edges == EXPECTED_FLIP_EDGES,
            "matches_cell_partition_flip_edges": flip_edges == cell_partition["cell_partition_summary"]["derived_phase_sign_flip_edges"],
            "matches_cell_sign_flip_edges": flip_edges == cell_sign["cell_sign_summary"]["derived_phase_sign_flip_edges"],
            "matches_z2_edge_bits": edge_bits == [row["edge_bit"] for row in z2["edge_bit_rows"]],
            "matches_gf2_solution_edge_bits": edge_bits == gf2["gf2_linear_system_summary"]["solution_edge_bit_pattern"],
        },
        "blocker_context": {
            "what_this_certifies": "finite GF(2) incidence from rational zero-carriers to the audited phase edge bits",
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_transport_derivation_from_nadsoliton_dynamics",
                "strict_damping_parameter_derivation_from_nadsoliton_dynamics",
                "orientation_chi11_source",
                "QW-2191_selector_obstruction",
            ],
        },
        "proof_certificate": {
            "containment_step": "Each ordered zero-carrier has positive rational clearance from both endpoints of its containing integer edge.",
            "incidence_step": "The carrier/edge incidence matrix has one 1 in each carrier column, at the containing edge.",
            "parity_step": "Multiplying by the all-one carrier multiplicity vector over GF(2) gives the audited edge-bit vector.",
            "rank_step": "The four carrier columns are distinct standard basis columns, so the carrier subspace has GF(2) column rank 4.",
            "theoretical_limit": "This proves finite carrier-to-edge parity only; it does not derive zero carriers, omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["carrier_edge_incidence_summary"]
    lines = [
        "# Phase-zero carrier/edge incidence certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate writes rational zero-carrier crossings as a GF(2)",
        "carrier-to-edge incidence matrix.  Multiplying by the all-one carrier",
        "multiplicity vector recovers the audited phase edge bits.",
        "",
        "## Summary",
        "",
        f"- Carrier order: `{summary['carrier_order']}`",
        f"- Edge bits from incidence: `{summary['edge_bit_pattern_from_carrier_incidence']}`",
        f"- Flip edges from incidence: `{summary['derived_flip_edges_from_carrier_incidence']}`",
        f"- Column rank over GF(2): `{summary['column_rank_over_gf2']}`",
        f"- All carriers strictly inside open edges: `{summary['all_carriers_strictly_inside_open_edges_by_rational_bounds']}`",
        f"- All carrier columns single-edge standard basis: `{summary['all_carrier_columns_are_single_edge_standard_basis_columns']}`",
        f"- All edges have at most one incident carrier: `{summary['all_edges_have_at_most_one_incident_carrier']}`",
        f"- Matches Z2 edge bits: `{summary['matches_z2_edge_bits']}`",
        f"- Matches GF(2) solution edge bits: `{summary['matches_gf2_solution_edge_bits']}`",
        "",
        "## Edge incidence rows",
        "",
        "| edge | incident carriers | edge bit | parity pass |",
        "| --- | --- | ---: | :---: |",
    ]
    for row in payload["edge_incidence_rows"]:
        lines.append(
            f"| {row['edge']} | {row['incident_zero_carriers']} | {row['edge_bit_from_incidence']} | `{row['edge_bit_equals_odd_carrier_parity']}` |"
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
