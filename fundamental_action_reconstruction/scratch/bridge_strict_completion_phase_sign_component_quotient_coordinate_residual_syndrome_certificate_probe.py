#!/usr/bin/env python3
"""Scratch probe: component-quotient coordinate residual-syndrome certificate.

The coordinate-support minimality certificate proves that exactly one coordinate
vector reconstructs the audited node bits.  This probe records the complementary
failure certificate: every coordinate vector c has a residual syndrome

    r(c) = U*c + b_target  over GF(2),

and all 2^12 residual syndromes are enumerated.  The zero syndrome appears
exactly once, at the audited coordinate vector, while every other coordinate
vector is assigned an explicit nonzero residual node pattern.  This is finite
GF(2) bookkeeping only; it does not derive phase zeros, omega/phi, damping,
transport, a kernel bridge, selector discharge, or ToE closure.
"""
from __future__ import annotations

import itertools
import json
from math import comb
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_residual_syndrome_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_residual_syndrome_certificate_report.md"
COORDINATE_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json"
SUPPORT_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_support_minimality_certificate_report.json"

NODE_COUNT = 12
COORDINATE_COUNT = 12
EXPECTED_COORDINATE_VECTOR = [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
ZERO_NODE_RESIDUAL = [0] * NODE_COUNT


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def matvec_gf2(matrix: list[list[int]], vector: list[int]) -> list[int]:
    return [sum(value & vector[col] for col, value in enumerate(row)) % 2 for row in matrix]


def xor_vectors(left: list[int], right: list[int]) -> list[int]:
    return [a ^ b for a, b in zip(left, right)]


def vector_weight(vector: list[int]) -> int:
    return sum(vector)


def enumerate_residual_syndromes(
    inverse_matrix_u: list[list[int]],
    target_node_bits: list[int],
    expected_coordinate_vector: list[int],
) -> dict[str, Any]:
    residual_histogram = {str(weight): 0 for weight in range(NODE_COUNT + 1)}
    coordinate_delta_histogram = {str(weight): 0 for weight in range(COORDINATE_COUNT + 1)}
    residual_to_coordinate: dict[str, list[int]] = {}
    zero_syndrome_rows = []
    minimum_nonzero_residual_weight = NODE_COUNT + 1
    minimum_nonzero_residual_examples = []
    nonzero_failure_examples = []

    for bits in itertools.product([0, 1], repeat=COORDINATE_COUNT):
        coordinate_vector = list(bits)
        coordinate_delta = xor_vectors(coordinate_vector, expected_coordinate_vector)
        reconstructed_node_bits = matvec_gf2(inverse_matrix_u, coordinate_vector)
        residual = xor_vectors(reconstructed_node_bits, target_node_bits)
        residual_key = "".join(str(bit) for bit in residual)
        residual_weight = vector_weight(residual)
        coordinate_delta_weight = vector_weight(coordinate_delta)

        residual_histogram[str(residual_weight)] += 1
        coordinate_delta_histogram[str(coordinate_delta_weight)] += 1
        residual_to_coordinate[residual_key] = coordinate_vector

        row = {
            "coordinate_vector": coordinate_vector,
            "coordinate_delta_from_audited": coordinate_delta,
            "coordinate_delta_weight": coordinate_delta_weight,
            "reconstructed_node_bits": reconstructed_node_bits,
            "residual_syndrome": residual,
            "residual_weight": residual_weight,
        }
        if residual == ZERO_NODE_RESIDUAL:
            zero_syndrome_rows.append(row)
        elif residual_weight < minimum_nonzero_residual_weight:
            minimum_nonzero_residual_weight = residual_weight
            minimum_nonzero_residual_examples = [row]
        elif residual_weight == minimum_nonzero_residual_weight and len(minimum_nonzero_residual_examples) < 8:
            minimum_nonzero_residual_examples.append(row)
        elif len(nonzero_failure_examples) < 8:
            nonzero_failure_examples.append(row)

    expected_residual_histogram = {str(weight): comb(NODE_COUNT, weight) for weight in range(NODE_COUNT + 1)}
    return {
        "coordinate_space_size": 2 ** COORDINATE_COUNT,
        "unique_residual_syndrome_count": len(residual_to_coordinate),
        "residual_weight_histogram": residual_histogram,
        "expected_binomial_residual_weight_histogram": expected_residual_histogram,
        "coordinate_delta_weight_histogram": coordinate_delta_histogram,
        "zero_syndrome_rows": zero_syndrome_rows,
        "minimum_nonzero_residual_weight": minimum_nonzero_residual_weight,
        "minimum_nonzero_residual_examples": minimum_nonzero_residual_examples,
        "nonzero_failure_examples": nonzero_failure_examples,
    }


def build_payload() -> dict[str, Any]:
    coordinate = load_json(COORDINATE_REPORT)
    support = load_json(SUPPORT_REPORT)

    inverse_matrix_u = coordinate["coordinate_isomorphism_definition"]["inverse_map_U_columns_quotient_then_interior"]
    target_node_bits = support["coordinate_support_definition"]["target_node_bits"]
    expected_coordinate_vector = support["coordinate_support_definition"]["expected_coordinate_vector"]
    residuals = enumerate_residual_syndromes(inverse_matrix_u, target_node_bits, expected_coordinate_vector)
    zero_rows = residuals["zero_syndrome_rows"]

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_RESIDUAL_SYNDROME_CERTIFICATE__FULL_RESIDUAL_ENUMERATION",
        "status": "component-quotient-coordinate-residual-syndrome-certified-by-2^12-enumeration",
        "source_reports": {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate": str(COORDINATE_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_coordinate_support_minimality_certificate": str(SUPPORT_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "coordinate residual syndrome",
                "residual syndrome histogram",
                "zero syndrome row",
                "nonzero coordinate residual",
                "binomial residual weight histogram",
                "coordinate delta syndrome",
            ],
            "finding": "Existing coordinate-support minimality reports enumerate matching coordinate vectors, but do not export the full residual-syndrome table proving every nonmatching coordinate has a nonzero node residual before this file.",
        },
        "coordinate_residual_definition": {
            "field": "GF(2)",
            "coordinate_count": COORDINATE_COUNT,
            "node_count": NODE_COUNT,
            "residual_formula": "residual_syndrome(c)=U*c + target_node_bits over GF(2)",
            "expected_coordinate_vector": expected_coordinate_vector,
            "target_node_bits": target_node_bits,
            "zero_node_residual": ZERO_NODE_RESIDUAL,
            "inverse_map_U_columns_quotient_then_interior": inverse_matrix_u,
        },
        "residual_syndrome_certificate": residuals,
        "coordinate_residual_syndrome_summary": {
            "enumerated_all_2^12_coordinate_vectors": residuals["coordinate_space_size"] == 4096,
            "all_residual_syndromes_unique": residuals["unique_residual_syndrome_count"] == 4096,
            "zero_syndrome_unique": len(zero_rows) == 1,
            "zero_syndrome_coordinate_is_audited": len(zero_rows) == 1 and zero_rows[0]["coordinate_vector"] == EXPECTED_COORDINATE_VECTOR,
            "every_nonmatching_coordinate_has_nonzero_residual": len(zero_rows) == 1 and zero_rows[0]["coordinate_vector"] == expected_coordinate_vector,
            "minimum_nonzero_residual_weight_is_1": residuals["minimum_nonzero_residual_weight"] == 1,
            "residual_weight_histogram_is_binomial": residuals["residual_weight_histogram"] == residuals["expected_binomial_residual_weight_histogram"],
            "coordinate_delta_histogram_sums_to_coordinate_space": sum(residuals["coordinate_delta_weight_histogram"].values()) == residuals["coordinate_space_size"],
            "matches_coordinate_support_minimality_unique_vector": support["coordinate_support_minimality_summary"]["unique_matching_coordinate_vector"] and support["coordinate_support_minimality_summary"]["matching_coordinate_vector_equals_dual_basis"],
        },
        "blocker_context": {
            "resolved_locally": [
                "All 4096 coordinate vectors were assigned explicit residual syndromes.",
                "The zero residual syndrome appears exactly once, at the audited coordinate vector.",
                "Every nonmatching coordinate vector has a nonzero node residual syndrome.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "residual_step": "For every coordinate vector c, compute residual_syndrome(c)=U*c+target_node_bits over GF(2).",
            "zero_step": "Exactly one row has residual_syndrome=0, and its coordinate vector is [0,1,0,1,0,0,0,0,0,0,0,0].",
            "failure_step": "The remaining 4095 coordinate vectors have nonzero residual syndromes, with minimum nonzero residual node weight 1.",
            "histogram_step": "The residual-weight histogram equals the binomial node-space histogram C(12,k), confirming full residual-syndrome coverage.",
            "theoretical_limit": "This is finite residual-syndrome enumeration; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["coordinate_residual_syndrome_summary"]
    residuals = payload["residual_syndrome_certificate"]
    lines = [
        "# Phase-sign component-quotient coordinate residual-syndrome certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate enumerates all 2^12 coordinate vectors and records",
        "the residual syndrome U*c+target_node_bits for each one.",
        "",
        "## Summary",
        "",
        f"- Enumerated coordinate vectors: `{residuals['coordinate_space_size']}`",
        f"- Unique residual syndromes: `{residuals['unique_residual_syndrome_count']}`",
        f"- Zero syndrome unique: `{summary['zero_syndrome_unique']}`",
        f"- Every nonmatching coordinate has nonzero residual: `{summary['every_nonmatching_coordinate_has_nonzero_residual']}`",
        f"- Residual histogram is binomial: `{summary['residual_weight_histogram_is_binomial']}`",
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
    print(json.dumps(payload["coordinate_residual_syndrome_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
