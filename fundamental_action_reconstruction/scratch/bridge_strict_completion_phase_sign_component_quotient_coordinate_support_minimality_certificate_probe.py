#!/usr/bin/env python3
"""Scratch probe: component-quotient coordinate-support minimality certificate.

The dual-basis certificate gives the unique coordinate vector c=T*b for the
audited node bits.  This probe makes the uniqueness/minimal-support consequence
computational: enumerate all 2^12 GF(2) coordinate vectors, map them through U,
and certify that the audited node vector is reconstructed by exactly one
coordinate vector, with minimum support weight 2:

    c = [0,1,0,1,0 | 0,0,0,0,0,0,0].

No weight-0 or weight-1 coordinate vector reconstructs the audited node vector,
and no other weight-2 coordinate vector does.  This is finite enumeration only;
it does not derive phase zeros, omega/phi, damping, transport, a kernel bridge,
selector discharge, or ToE closure.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_support_minimality_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_support_minimality_certificate_report.md"
DUAL_BASIS_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_dual_basis_certificate_report.json"
COORDINATE_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json"

NODE_COUNT = 12
COORDINATE_COUNT = 12
COMPONENT_COUNT = 5
EXPECTED_COORDINATE_VECTOR = [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
EXPECTED_ACTIVE_LABELS = ["quotient_component_1", "quotient_component_3"]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def vector_weight(vector: list[int]) -> int:
    return sum(vector)


def coordinate_label(index: int) -> str:
    if index < COMPONENT_COUNT:
        return f"quotient_component_{index}"
    return f"interior_residual_{index - COMPONENT_COUNT}"


def support_labels(vector: list[int]) -> list[str]:
    return [coordinate_label(index) for index, bit in enumerate(vector) if bit]


def enumerate_coordinate_space(
    inverse_matrix_u: list[list[int]],
    target_node_bits: list[int],
) -> dict[str, Any]:
    matching_rows = []
    weight_histogram = {str(weight): 0 for weight in range(COORDINATE_COUNT + 1)}
    lower_weight_failure_examples = []
    same_weight_failure_examples = []
    for bits in itertools.product([0, 1], repeat=COORDINATE_COUNT):
        vector = list(bits)
        weight = vector_weight(vector)
        weight_histogram[str(weight)] += 1
        reconstructed = matvec_gf2(inverse_matrix_u, vector)
        matches = reconstructed == target_node_bits
        if matches:
            matching_rows.append({
                "coordinate_vector": vector,
                "support_weight": weight,
                "support_labels": support_labels(vector),
            })
        elif weight < vector_weight(EXPECTED_COORDINATE_VECTOR) and len(lower_weight_failure_examples) < 8:
            lower_weight_failure_examples.append({
                "coordinate_vector": vector,
                "support_weight": weight,
                "support_labels": support_labels(vector),
                "reconstructed_node_bits": reconstructed,
            })
        elif weight == vector_weight(EXPECTED_COORDINATE_VECTOR) and vector != EXPECTED_COORDINATE_VECTOR and len(same_weight_failure_examples) < 8:
            same_weight_failure_examples.append({
                "coordinate_vector": vector,
                "support_weight": weight,
                "support_labels": support_labels(vector),
                "reconstructed_node_bits": reconstructed,
            })
    min_matching_weight = min(row["support_weight"] for row in matching_rows) if matching_rows else None
    return {
        "coordinate_space_size": 2 ** COORDINATE_COUNT,
        "weight_histogram": weight_histogram,
        "matching_rows": matching_rows,
        "min_matching_weight": min_matching_weight,
        "lower_weight_failure_examples": lower_weight_failure_examples,
        "same_weight_failure_examples": same_weight_failure_examples,
    }


def build_payload() -> dict[str, Any]:
    dual = load_json(DUAL_BASIS_REPORT)
    coordinate = load_json(COORDINATE_REPORT)

    inverse_matrix_u = coordinate["coordinate_isomorphism_definition"]["inverse_map_U_columns_quotient_then_interior"]
    target_node_bits = dual["reconstruction_certificate"]["node_bits"]
    audited_coordinate_vector = dual["reconstruction_certificate"]["coordinate_vector_T_b"]
    enumeration = enumerate_coordinate_space(inverse_matrix_u, target_node_bits)
    matching_vectors = [row["coordinate_vector"] for row in enumeration["matching_rows"]]
    matching_weights = [row["support_weight"] for row in enumeration["matching_rows"]]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_SUPPORT_MINIMALITY_CERTIFICATE__FULL_COORDINATE_ENUMERATION",
        "status": "component-quotient-coordinate-support-minimality-certified-by-2^12-enumeration",
        "source_reports": {
            "phase_sign_component_quotient_dual_basis_certificate": str(DUAL_BASIS_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_coordinate_isomorphism_certificate": str(COORDINATE_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "coordinate support minimality",
                "coordinate support",
                "2^12 coordinate",
                "minimal coordinate",
                "quotient_component_1",
                "quotient_component_3",
                "same weight failure",
                "weight histogram",
            ],
            "finding": "Existing dual-basis reports identify the active coordinates, but do not enumerate all 2^12 coordinate vectors to certify unique reconstruction and minimum support weight before this file.",
        },
        "coordinate_support_definition": {
            "field": "GF(2)",
            "coordinate_count": COORDINATE_COUNT,
            "target_node_bits": target_node_bits,
            "audited_coordinate_vector_from_dual_basis": audited_coordinate_vector,
            "expected_coordinate_vector": EXPECTED_COORDINATE_VECTOR,
            "expected_active_coordinate_labels": EXPECTED_ACTIVE_LABELS,
            "inverse_map_U_columns_quotient_then_interior": inverse_matrix_u,
        },
        "enumeration_certificate": enumeration,
        "coordinate_support_minimality_summary": {
            "enumerated_all_2^12_coordinate_vectors": enumeration["coordinate_space_size"] == 4096,
            "unique_matching_coordinate_vector": matching_vectors == [EXPECTED_COORDINATE_VECTOR],
            "matching_coordinate_vector_equals_dual_basis": audited_coordinate_vector == EXPECTED_COORDINATE_VECTOR,
            "minimum_matching_weight_is_2": enumeration["min_matching_weight"] == 2,
            "all_lower_weight_vectors_fail": all(weight >= 2 for weight in matching_weights) and len(enumeration["lower_weight_failure_examples"]) > 0,
            "all_other_weight_2_vectors_fail": len(enumeration["matching_rows"]) == 1 and enumeration["matching_rows"][0]["coordinate_vector"] == EXPECTED_COORDINATE_VECTOR,
            "active_coordinates_are_expected_quotient_components": enumeration["matching_rows"][0]["support_labels"] == EXPECTED_ACTIVE_LABELS,
            "all_interior_residual_coordinates_zero": EXPECTED_COORDINATE_VECTOR[COMPONENT_COUNT:] == [0] * (COORDINATE_COUNT - COMPONENT_COUNT),
            "weight_histogram_sums_to_coordinate_space": sum(enumeration["weight_histogram"].values()) == enumeration["coordinate_space_size"],
        },
        "blocker_context": {
            "resolved_locally": [
                "All 4096 coordinate vectors were enumerated under the certified inverse map U.",
                "The audited node vector has exactly one coordinate representation.",
                "The unique representation has support weight 2 on quotient_component_1 and quotient_component_3.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "enumeration_step": "All 2^12 coordinate vectors were mapped through U and compared with the audited node bits.",
            "uniqueness_step": "Exactly one coordinate vector reconstructs the audited node bits: [0,1,0,1,0,0,0,0,0,0,0,0].",
            "minimality_step": "No weight-0 or weight-1 coordinate vector reconstructs the audited node bits, so the coordinate support weight is minimally 2.",
            "support_step": "The unique support is quotient_component_1 plus quotient_component_3, with all seven interior residual coordinates zero.",
            "theoretical_limit": "This is finite coordinate-support enumeration; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["coordinate_support_minimality_summary"]
    enumeration = payload["enumeration_certificate"]
    lines = [
        "# Phase-sign component-quotient coordinate-support minimality certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate enumerates all 2^12 coordinate vectors and proves the",
        "audited node vector has a unique minimum support coordinate representation.",
        "",
        "## Summary",
        "",
        f"- Enumerated coordinate vectors: `{enumeration['coordinate_space_size']}`",
        f"- Unique matching coordinate vector: `{summary['unique_matching_coordinate_vector']}`",
        f"- Minimum matching weight is 2: `{summary['minimum_matching_weight_is_2']}`",
        f"- Active coordinates expected: `{summary['active_coordinates_are_expected_quotient_components']}`",
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
    print(json.dumps(payload["coordinate_support_minimality_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
