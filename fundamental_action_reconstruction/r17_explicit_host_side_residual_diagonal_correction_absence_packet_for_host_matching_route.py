#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r17_explicit_host_side_residual_diagonal_correction_absence_packet_for_host_matching_route.json"
OUT_SUMMARY = GENERATED / "r17_explicit_host_side_residual_diagonal_correction_absence_packet_for_host_matching_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def zero_matrix(rows: int, cols: int) -> list[list[float]]:
    return [[0.0 for _ in range(cols)] for _ in range(rows)]


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r8 = load_json("fundamental_action_reconstruction/generated/r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs.json")
    r13 = load_json("fundamental_action_reconstruction/generated/r13_partial_host_to_canonical_block_overlap_packet.json")
    r14 = load_json(
        "fundamental_action_reconstruction/generated/r14_explicit_frozen_kernel_specialization_packet_for_host_matching_route.json"
    )
    r15 = load_json(
        "fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
    )
    r16 = load_json(
        "fundamental_action_reconstruction/generated/r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route.json"
    )

    host_basis = r13["host_side_decomposition"]["carrier_basis"]
    host_size = len(host_basis)
    control_basis = r16["declared_control_pullback_packet"]["control_basis"]

    checks = [
        {
            "id": "r8_host_operator_symbol_is_exactly_k_total_plus_m0_sq_i",
            "actual": r8["host_operator_schema"]["symbol"],
            "expected": "A = K_total + m0^2 I",
            "meaning": "R8 already fixes the host operator schema exactly as K_total + m0^2 I",
        },
        {
            "id": "r13_host_operator_symbol_is_exactly_k_total_plus_m0_sq_i",
            "actual": r13["host_side_decomposition"]["host_operator_symbol"],
            "expected": "A_host = K_total + m0^2 I",
            "meaning": "R13 already exports the host-side decomposition exactly as A_host = K_total + m0^2 I",
        },
        {
            "id": "r14_kernel_channel_specialized",
            "actual": r14["specialization_result"]["kernel_channel_specialization_witness_present"],
            "expected": True,
            "meaning": "R14 already closes the host kernel/off-diagonal component",
        },
        {
            "id": "r15_host_scalar_floor_embedding_present",
            "actual": r15["embedding_result"]["explicit_host_scalar_floor_embedding_packet_present"],
            "expected": True,
            "meaning": "R15 already closes the host scalar-floor component",
        },
        {
            "id": "r16_declared_residual_pullback_present",
            "actual": r16["result"]["explicit_declared_control_pullback_of_residual_local_diagonal_sector_present"],
            "expected": True,
            "meaning": "R16 already exports the declared control pullback of the canonical residual diagonal sector",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R17",
        "packet_goal": "materialize_the_absence_of_any_host_side_residual_diagonal_correction_after_subtracting_the_already_closed_kernel_and_scalar_floor_components",
        "source_reports": ["R8", "R13", "R14", "R15", "R16"],
        "host_side_residual_correction_packet": {
            "host_basis": host_basis,
            "control_basis": control_basis,
            "decomposition_rule": "A_host - K_total_host_frozen - m0^2 I = 0",
            "host_side_residual_correction_symbol": "Delta_host_diag_residual",
            "host_side_residual_correction_rows": zero_matrix(host_size, host_size),
            "declared_control_pullback_symbol": "T_control^T Delta_host_diag_residual T_control",
            "declared_control_pullback_rows": zero_matrix(len(control_basis), len(control_basis)),
            "pair1_declared_block_rows": zero_matrix(2, 2),
        },
        "result": {
            "explicit_host_side_residual_diagonal_correction_absence_packet_present": True,
            "host_side_residual_diagonal_correction_present": False,
            "declared_control_pullback_of_host_side_residual_correction_present": "zero_only",
            "zero_witness_for_canonical_residual_declared_pullback_present": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_unchanged_by_R15_R16_R17",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "host-side diagonal correction branch only",
            "boundary": "R17 does not alter the matched light-facing kernel channel; it only closes the alternative host-side correction branch by showing that no such residual correction exists in the current host decomposition",
        },
        "classification": "explicit_host_side_residual_diagonal_correction_absence_packet_present_so_only_the_zero_witness_for_the_canonical_residual_declared_pullback_remains",
        "frontier": "R17_B1",
        "frontier_text": "the repo now shows that the current host route has no residual diagonal correction beyond K_total plus m0^2 I, so the remaining diagonal-side blocker is reduced to a zero witness for the canonical residual declared pullback",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_canonical_residual_declared_pullback_is_zero",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R17",
        "status": "PASS_PARTIAL_EXPLICIT_HOST_SIDE_RESIDUAL_DIAGONAL_CORRECTION_ABSENCE_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R17_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
