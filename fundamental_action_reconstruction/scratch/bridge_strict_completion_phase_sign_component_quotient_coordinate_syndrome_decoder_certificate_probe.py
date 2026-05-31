#!/usr/bin/env python3
"""Scratch probe: component-quotient coordinate syndrome-decoder certificate.

The residual-syndrome certificate assigns every coordinate vector c a residual
node syndrome r(c)=U*c+b_target.  This probe records the finite decoder that
turns any such residual syndrome back into the unique coordinate correction:

    delta(c) = T*r(c),     c + delta(c) = c_target    over GF(2).

Equivalently, for every node residual s in GF(2)^12, the coordinate correction
T*s satisfies U*T*s=s.  This is a full 2^12-table decoding certificate only; it
does not derive phase zeros, omega/phi, damping, transport, a kernel bridge,
selector discharge, or ToE closure.
"""
from __future__ import annotations

import itertools
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_decoder_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_syndrome_decoder_certificate_report.md"
COORDINATE_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_isomorphism_certificate_report.json"
RESIDUAL_REPORT = HERE / "bridge_strict_completion_phase_sign_component_quotient_coordinate_residual_syndrome_certificate_report.json"

NODE_COUNT = 12
COORDINATE_COUNT = 12
EXPECTED_COORDINATE_VECTOR = [0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
ZERO_VECTOR = [0] * NODE_COUNT


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


def enumerate_decoder(
    coordinate_map_t: list[list[int]],
    inverse_matrix_u: list[list[int]],
    target_node_bits: list[int],
    expected_coordinate_vector: list[int],
) -> dict[str, Any]:
    coordinate_rows_checked = 0
    decoded_coordinate_rows = []
    correction_weight_histogram = {str(weight): 0 for weight in range(COORDINATE_COUNT + 1)}
    residual_weight_to_correction_weight_examples = []
    coordinate_decode_failures = []

    for bits in itertools.product([0, 1], repeat=COORDINATE_COUNT):
        coordinate_vector = list(bits)
        reconstructed_node_bits = matvec_gf2(inverse_matrix_u, coordinate_vector)
        residual_syndrome = xor_vectors(reconstructed_node_bits, target_node_bits)
        correction_delta = matvec_gf2(coordinate_map_t, residual_syndrome)
        corrected_coordinate_vector = xor_vectors(coordinate_vector, correction_delta)
        correction_weight = vector_weight(correction_delta)
        correction_weight_histogram[str(correction_weight)] += 1
        coordinate_rows_checked += 1
        row = {
            "coordinate_vector": coordinate_vector,
            "residual_syndrome": residual_syndrome,
            "residual_weight": vector_weight(residual_syndrome),
            "correction_delta_T_residual": correction_delta,
            "correction_weight": correction_weight,
            "corrected_coordinate_vector": corrected_coordinate_vector,
        }
        if corrected_coordinate_vector == expected_coordinate_vector:
            if len(decoded_coordinate_rows) < 8:
                decoded_coordinate_rows.append(row)
            elif correction_weight <= 1 and len(decoded_coordinate_rows) < 16:
                decoded_coordinate_rows.append(row)
        else:
            coordinate_decode_failures.append(row)
        if len(residual_weight_to_correction_weight_examples) < 16:
            residual_weight_to_correction_weight_examples.append({
                "residual_weight": row["residual_weight"],
                "correction_weight": correction_weight,
                "residual_syndrome": residual_syndrome,
                "correction_delta_T_residual": correction_delta,
            })

    residual_rows_checked = 0
    residual_decode_failures = []
    zero_residual_decoder_rows = []
    for bits in itertools.product([0, 1], repeat=NODE_COUNT):
        residual_syndrome = list(bits)
        correction_delta = matvec_gf2(coordinate_map_t, residual_syndrome)
        reencoded_residual = matvec_gf2(inverse_matrix_u, correction_delta)
        residual_rows_checked += 1
        row = {
            "residual_syndrome": residual_syndrome,
            "correction_delta_T_residual": correction_delta,
            "reencoded_residual_U_T_residual": reencoded_residual,
            "residual_weight": vector_weight(residual_syndrome),
            "correction_weight": vector_weight(correction_delta),
        }
        if residual_syndrome == ZERO_VECTOR:
            zero_residual_decoder_rows.append(row)
        if reencoded_residual != residual_syndrome:
            residual_decode_failures.append(row)

    return {
        "coordinate_rows_checked": coordinate_rows_checked,
        "residual_rows_checked": residual_rows_checked,
        "decoded_coordinate_examples": decoded_coordinate_rows,
        "residual_weight_to_correction_weight_examples": residual_weight_to_correction_weight_examples,
        "correction_weight_histogram": correction_weight_histogram,
        "coordinate_decode_failures": coordinate_decode_failures,
        "residual_decode_failures": residual_decode_failures,
        "zero_residual_decoder_rows": zero_residual_decoder_rows,
    }


def build_payload() -> dict[str, Any]:
    coordinate = load_json(COORDINATE_REPORT)
    residual = load_json(RESIDUAL_REPORT)

    coordinate_map_t = coordinate["coordinate_isomorphism_definition"]["coordinate_map_T_rows_quotient_then_interior"]
    inverse_matrix_u = coordinate["coordinate_isomorphism_definition"]["inverse_map_U_columns_quotient_then_interior"]
    target_node_bits = residual["coordinate_residual_definition"]["target_node_bits"]
    expected_coordinate_vector = residual["coordinate_residual_definition"]["expected_coordinate_vector"]
    decoder = enumerate_decoder(coordinate_map_t, inverse_matrix_u, target_node_bits, expected_coordinate_vector)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_PHASE_SIGN_COMPONENT_QUOTIENT_COORDINATE_SYNDROME_DECODER_CERTIFICATE__FULL_DECODER_TABLE",
        "status": "component-quotient-coordinate-syndrome-decoder-certified-by-2^12-enumeration",
        "source_reports": {
            "phase_sign_component_quotient_coordinate_isomorphism_certificate": str(COORDINATE_REPORT.relative_to(ROOT)),
            "phase_sign_component_quotient_coordinate_residual_syndrome_certificate": str(RESIDUAL_REPORT.relative_to(ROOT)),
        },
        "grep_disambiguation": {
            "searched_terms": [
                "coordinate syndrome decoder",
                "T residual correction",
                "coordinate correction decoder",
                "U T residual decoder",
                "decoded coordinate vector",
                "syndrome correction table",
            ],
            "finding": "Existing coordinate residual-syndrome reports enumerate residuals, but do not export the full T-based syndrome decoder proving c+T*r(c)=c_target before this file.",
        },
        "coordinate_syndrome_decoder_definition": {
            "field": "GF(2)",
            "coordinate_count": COORDINATE_COUNT,
            "node_count": NODE_COUNT,
            "residual_formula": "r(c)=U*c+target_node_bits over GF(2)",
            "decoder_formula": "delta(c)=T*r(c); corrected_coordinate=c+delta(c)",
            "residual_decoder_formula": "U*T*s=s for every residual syndrome s",
            "expected_coordinate_vector": expected_coordinate_vector,
            "target_node_bits": target_node_bits,
            "coordinate_map_T_rows_quotient_then_interior": coordinate_map_t,
            "inverse_map_U_columns_quotient_then_interior": inverse_matrix_u,
        },
        "decoder_certificate": decoder,
        "coordinate_syndrome_decoder_summary": {
            "enumerated_all_2^12_coordinate_vectors": decoder["coordinate_rows_checked"] == 4096,
            "enumerated_all_2^12_residual_syndromes": decoder["residual_rows_checked"] == 4096,
            "all_coordinate_vectors_decode_to_audited_coordinate": len(decoder["coordinate_decode_failures"]) == 0,
            "all_residual_syndromes_reencode_correctly": len(decoder["residual_decode_failures"]) == 0,
            "zero_residual_decodes_to_zero_correction": len(decoder["zero_residual_decoder_rows"]) == 1 and decoder["zero_residual_decoder_rows"][0]["correction_delta_T_residual"] == ZERO_VECTOR,
            "correction_weight_histogram_sums_to_coordinate_space": sum(decoder["correction_weight_histogram"].values()) == decoder["coordinate_rows_checked"],
            "matches_residual_syndrome_certificate": residual["coordinate_residual_syndrome_summary"]["all_residual_syndromes_unique"] and residual["coordinate_residual_syndrome_summary"]["zero_syndrome_coordinate_is_audited"],
        },
        "blocker_context": {
            "resolved_locally": [
                "All 4096 coordinate vectors decode back to the audited coordinate vector by c+T*r(c).",
                "All 4096 residual syndromes satisfy U*T*s=s.",
                "The zero residual syndrome decodes to the zero coordinate correction.",
            ],
            "still_open": [
                "strict_phase_frequency_derivation_from_nadsoliton_dynamics",
                "strict_damping_transport_derivation_from_nadsoliton_dynamics",
                "QW-2191_selector_obstruction",
                "strict ToE closure",
            ],
        },
        "proof_certificate": {
            "coordinate_decoder_step": "For every coordinate vector c, compute r(c)=U*c+target_node_bits and verify c+T*r(c)=c_target.",
            "residual_decoder_step": "For every residual syndrome s in GF(2)^12, verify U*T*s=s.",
            "zero_decoder_step": "The zero residual syndrome has zero coordinate correction, so the audited coordinate vector is fixed.",
            "histogram_step": "The correction-weight histogram covers all 4096 coordinate vectors and records the finite decoder table.",
            "theoretical_limit": "This is finite syndrome-decoder enumeration; it does not derive omega/phi, beta/eta, damping, or transport from strict nadsoliton dynamics.",
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
    summary = payload["coordinate_syndrome_decoder_summary"]
    decoder = payload["decoder_certificate"]
    lines = [
        "# Phase-sign component-quotient coordinate syndrome-decoder certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This certificate enumerates all coordinate vectors and all residual",
        "syndromes to verify the finite T-based decoder.",
        "",
        "## Summary",
        "",
        f"- Coordinate vectors checked: `{decoder['coordinate_rows_checked']}`",
        f"- Residual syndromes checked: `{decoder['residual_rows_checked']}`",
        f"- All coordinates decode to target: `{summary['all_coordinate_vectors_decode_to_audited_coordinate']}`",
        f"- All residual syndromes reencode correctly: `{summary['all_residual_syndromes_reencode_correctly']}`",
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
    print(json.dumps(payload["coordinate_syndrome_decoder_summary"], indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
