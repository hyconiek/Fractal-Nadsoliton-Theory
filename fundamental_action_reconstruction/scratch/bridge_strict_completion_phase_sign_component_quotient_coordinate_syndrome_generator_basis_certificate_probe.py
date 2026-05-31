#!/usr/bin/env python3
"""Scratch probe: component-quotient coordinate syndrome-generator basis certificate.

The syndrome-decoder certificate verifies the global table c+T*r(c)=c_target.
This probe records the local generator-level reason behind that table: flipping
coordinate i changes the residual syndrome by column U_i, and decoding that
residual generator by T returns the unit coordinate vector e_i.

The probe checks all 12 coordinate generators and all 4096*12 hypercube edges:

    r(c+e_i)+r(c) = U_i,      T*U_i = e_i    over GF(2).

This is finite GF(2) bookkeeping only; it does not derive phase zeros,
omega/phi, damping, transport, a kernel bridge, selector discharge, or ToE
closure.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_generator_basis_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_generator_basis_certificate_report.md"
COORDINATE_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json"
DECODER_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_decoder_certificate_report.json"

NODE_COUNT = 12
COORDINATE_COUNT = 12


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def xor_vectors(left: list[int], right: list[int]) -> list[int]:
    return [a ^ b for a, b in zip(left, right)]


def unit_vector(index: int, length: int) -> list[int]:
    return [1 if col == index else 0 for col in range(length)]


def column(matrix: list[list[int]], index: int) -> list[int]:
    return [row[index] for row in matrix]


def rank_gf2(rows: list[list[int]]) -> int:
    matrix = [row[:] for row in rows if any(row)]
    if not matrix:
        return 0
    rank = 0
    col = 0
    width = len(matrix[0])
    while rank < len(matrix) and col < width:
        pivot = next((row for row in range(rank, len(matrix)) if matrix[row][col]), None)
        if pivot is None:
            col += 1
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for row in range(len(matrix)):
            if row != rank and matrix[row][col]:
                matrix[row] = [a ^ b for a, b in zip(matrix[row], matrix[rank])]
        rank += 1
        col += 1
    return rank


def build_generator_rows(coordinate_map_t: list[list[int]], inverse_matrix_u: list[list[int]]) -> list[dict[str, Any]]:
    rows = []
    for index in range(COORDINATE_COUNT):
        residual_generator = column(inverse_matrix_u, index)
        decoded_coordinate_generator = matvec_gf2(coordinate_map_t, residual_generator)
        rows.append({
            "coordinate_index": index,
            "coordinate_unit_e_i": unit_vector(index, COORDINATE_COUNT),
            "residual_generator_U_i": residual_generator,
            "residual_generator_weight": sum(residual_generator),
            "decoded_T_U_i": decoded_coordinate_generator,
            "decoded_matches_unit": decoded_coordinate_generator == unit_vector(index, COORDINATE_COUNT),
        })
    return rows


def check_hypercube_edges(inverse_matrix_u: list[list[int]], target_node_bits: list[int]) -> dict[str, Any]:
    edge_count = 0
    failures = []
    examples = []
    generator_columns = [column(inverse_matrix_u, index) for index in range(COORDINATE_COUNT)]
    for bits in itertools.product([0, 1], repeat=COORDINATE_COUNT):
        coordinate_vector = list(bits)
        residual = xor_vectors(matvec_gf2(inverse_matrix_u, coordinate_vector), target_node_bits)
        for index in range(COORDINATE_COUNT):
            flipped = coordinate_vector[:]
            flipped[index] ^= 1
            flipped_residual = xor_vectors(matvec_gf2(inverse_matrix_u, flipped), target_node_bits)
            residual_delta = xor_vectors(flipped_residual, residual)
            expected_delta = generator_columns[index]
            edge_count += 1
            row = {
                "coordinate_vector": coordinate_vector,
                "flipped_index": index,
                "residual_delta": residual_delta,
                "expected_residual_generator_U_i": expected_delta,
            }
            if residual_delta != expected_delta:
                failures.append(row)
            elif len(examples) < 12:
                examples.append(row)
    return {
        "hypercube_edge_checks": edge_count,
        "hypercube_edge_failures": failures,
        "hypercube_edge_examples": examples,
    }


def build_payload() -> dict[str, Any]:
    coordinate = load_json(COORDINATE_REPORT)
    decoder = load_json(DECODER_REPORT)

    coordinate_map_t = coordinate["coordinate_isomorphism_definition"]["coordinate_map_T_rows_quotient_then_interior"]
    inverse_matrix_u = coordinate["coordinate_isomorphism_definition"]["inverse_map_U_columns_quotient_then_interior"]
    target_node_bits = decoder["coordinate_syndrome_decoder_definition"]["target_node_bits"]
    generator_rows = build_generator_rows(coordinate_map_t, inverse_matrix_u)
    edge_certificate = check_hypercube_edges(inverse_matrix_u, target_node_bits)
    residual_generators = [row["residual_generator_U_i"] for row in generator_rows]
    decoded_generators = [row["decoded_T_U_i"] for row in generator_rows]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_SYNDROME_GENERATOR_BASIS_CERTIFICATE__HYPERCUBE_EDGE_CHECKS",
        "status": "component-quotient-coordinate-syndrome-generator-basis-certified",
        "source_reports": {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate": str(COORDINATE_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_coordinate_syndrome_decoder_certificate": str(DECODER_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "coordinate syndrome generator basis",
                "hypercube edge syndrome",
                "T U_i e_i",
                "residual generator U_i",
                "coordinate generator basis",
                "syndrome generator columns",
            ],
            "finding": "Existing syndrome-decoder reports verify the full correction table, but do not export the local generator basis U_i and all hypercube edge checks before this file.",
        },
        "generator_basis_definition": {
            "field": "GF(2)",
            "coordinate_count": COORDINATE_COUNT,
            "node_count": NODE_COUNT,
            "generator_formula": "r(c+e_i)+r(c)=U_i over GF(2)",
            "decoder_generator_formula": "T*U_i=e_i",
            "coordinate_map_T_rows_quotient_then_interior": coordinate_map_t,
            "inverse_map_U_columns_quotient_then_interior": inverse_matrix_u,
        },
        "generator_basis_certificate": {
            "generator_rows": generator_rows,
            "residual_generator_column_rank": rank_gf2(residual_generators),
            "decoded_generator_row_rank": rank_gf2(decoded_generators),
            "residual_generator_weights": [row["residual_generator_weight"] for row in generator_rows],
        },
        "hypercube_edge_certificate": edge_certificate,
        "coordinate_syndrome_generator_basis_summary": {
            "all_12_generators_checked": len(generator_rows) == COORDINATE_COUNT,
            "all_generators_decode_to_coordinate_units": all(row["decoded_matches_unit"] for row in generator_rows),
            "residual_generators_have_full_rank_12": rank_gf2(residual_generators) == COORDINATE_COUNT,
            "decoded_generators_have_full_rank_12": rank_gf2(decoded_generators) == COORDINATE_COUNT,
            "checked_all_4096_times_12_hypercube_edges": edge_certificate["hypercube_edge_checks"] == 4096 * COORDINATE_COUNT,
            "all_hypercube_edge_deltas_match_generators": edge_certificate["hypercube_edge_failures"] == [],
            "matches_syndrome_decoder_report": decoder["coordinate_syndrome_decoder_summary"]["all_coordinate_vectors_decode_to_audited_coordinate"] and decoder["coordinate_syndrome_decoder_summary"]["all_residual_syndromes_reencode_correctly"],
        },
        "blocker_context": {
            "resolved_locally": [
                "All 12 coordinate generator residual columns U_i were exported.",
                "Each residual generator decodes by T to its coordinate unit e_i.",
                "All 4096*12 coordinate hypercube edges have residual delta equal to the corresponding U_i.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "generator_step": "For each coordinate unit e_i, export U_i and verify T*U_i=e_i.",
            "rank_step": "The 12 residual generators U_i have GF(2) rank 12, so they form a node-syndrome basis.",
            "edge_step": "All 4096*12 hypercube edges satisfy r(c+e_i)+r(c)=U_i.",
            "decoder_link_step": "The local generator checks explain the full decoder relation c+T*r(c)=c_target exported by the syndrome-decoder report.",
            "theoretical_limit": "This is finite syndrome-generator bookkeeping; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["coordinate_syndrome_generator_basis_summary"]
    edge_certificate = payload["hypercube_edge_certificate"]
    lines = [
        "# Phase-sign component-quotient coordinate syndrome-generator basis certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate exports the local residual-syndrome generators and checks",
        "all coordinate hypercube edge deltas against them.",
        "",
        "## Summary",
        "",
        f"- Generators checked: `{payload['generator_basis_definition']['coordinate_count']}`",
        f"- Hypercube edge checks: `{edge_certificate['hypercube_edge_checks']}`",
        f"- Generators decode to coordinate units: `{summary['all_generators_decode_to_coordinate_units']}`",
        f"- Edge deltas match generators: `{summary['all_hypercube_edge_deltas_match_generators']}`",
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
    print(json.dumps(payload["coordinate_syndrome_generator_basis_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
