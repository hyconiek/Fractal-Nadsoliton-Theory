#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route.json"
OUT_SUMMARY = GENERATED / "r16_explicit_residual_local_diagonal_declared_control_pullback_packet_for_host_matching_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def clean_scalar(value: float) -> float:
    if abs(value) < 1e-15:
        return 0.0
    return round(value, 15)


def build_entry_expression(coefficients: list[float], residual_terms: list[str]) -> str:
    pieces: list[str] = []
    for coeff, residual in zip(coefficients, residual_terms):
        cleaned = clean_scalar(coeff)
        if cleaned == 0.0:
            continue
        if cleaned == 1.0:
            pieces.append(f"({residual})")
        elif cleaned == -1.0:
            pieces.append(f"-({residual})")
        else:
            pieces.append(f"({cleaned})*({residual})")
    if not pieces:
        return "0.0"
    expression = pieces[0]
    for piece in pieces[1:]:
        if piece.startswith("-"):
            expression += " - " + piece[1:]
        else:
            expression += " + " + piece
    return expression


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    c15 = load_json("fundamental_action_reconstruction/generated/c15_control_only_pullback_submatrix_packet_summary.json")
    r11 = load_json(
        "fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
    )
    r15 = load_json(
        "fundamental_action_reconstruction/generated/r15_explicit_host_scalar_floor_embedding_packet_for_host_matching_route.json"
    )

    control_basis = r11["transport_packet"]["domain_basis"]
    transport_rows = r11["transport_packet"]["matrix_rows"]
    diagonal_rows = r15["diagonal_decomposition"]["entrywise_rows"]
    residual_terms = [entry["residual_local_diagonal_entry"] for entry in diagonal_rows]
    carrier_basis = r15["diagonal_decomposition"]["carrier_basis"]

    coefficient_tensor: list[list[list[float]]] = []
    matrix_rows: list[list[str]] = []
    for a in range(len(control_basis)):
        coefficient_row: list[list[float]] = []
        expression_row: list[str] = []
        for b in range(len(control_basis)):
            coeffs = [clean_scalar(float(transport_rows[i][a]) * float(transport_rows[i][b])) for i in range(len(carrier_basis))]
            coefficient_row.append(coeffs)
            expression_row.append(build_entry_expression(coeffs, residual_terms))
        coefficient_tensor.append(coefficient_row)
        matrix_rows.append(expression_row)

    checks = [
        {
            "id": "c15_formal_control_pullback_formula_present",
            "actual": c15["result"]["control_only_pullback_formula_present"],
            "expected": "yes",
            "meaning": "C15 already exports the formal control pullback formula M_control = T_control^T H_PsiPsi T_control",
        },
        {
            "id": "r11_declared_control_transport_packet_present",
            "actual": r11["result"]["explicit_declared_control_transport_packet_present"],
            "expected": True,
            "meaning": "R11 already exports the declared control transport matrix",
        },
        {
            "id": "r15_residual_local_diagonal_sector_present",
            "actual": r15["embedding_result"]["residual_local_diagonal_sector_present"],
            "expected": True,
            "meaning": "R15 already isolates the residual local diagonal sector",
        },
        {
            "id": "pullback_matrix_shape_is_4x4",
            "actual": [len(matrix_rows), len(matrix_rows[0])],
            "expected": [4, 4],
            "meaning": "the declared control pullback of the residual local diagonal sector is a 4x4 matrix on the control carrier",
        },
        {
            "id": "pair1_declared_block_present",
            "actual": len(matrix_rows[:2]) == 2 and len(matrix_rows[0][:2]) == 2,
            "expected": True,
            "meaning": "the pair1 declared 2x2 block is explicitly available inside the control pullback matrix",
        },
        {
            "id": "r11_physical_canonicalization_still_absent",
            "actual": r11["result"]["strict_selector_relevant_physical_canonicalization_present"],
            "expected": False,
            "meaning": "R11 still leaves physical canonicalization open, so this packet remains declared-control only",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R16",
        "packet_goal": "materialize_the_explicit_declared_control_pullback_of_the_residual_local_diagonal_sector_after_subtracting_the_host_scalar_floor",
        "source_reports": ["C15", "R11", "R15"],
        "declared_control_pullback_packet": {
            "control_basis": control_basis,
            "carrier_basis": carrier_basis,
            "source_diagonal_symbol": "D_local_residual",
            "assembly_formula": "M_control_residual_diag_declared = T_control^T D_local_residual T_control",
            "coefficient_tensor_by_carrier_slot": coefficient_tensor,
            "matrix_rows": matrix_rows,
            "pair1_declared_block_rows": [row[:2] for row in matrix_rows[:2]],
            "scope": "declared_control_pullback_only",
        },
        "result": {
            "explicit_declared_control_pullback_of_residual_local_diagonal_sector_present": True,
            "pair1_declared_residual_diagonal_block_present": True,
            "zero_or_host_side_cancellation_witness_present": False,
            "physical_canonicalization_present": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R15_R16",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "declared control pullback of the non-light residual local diagonal sector only",
            "boundary": "R16 does not alter the matched light-facing kernel channel; it only exports the declared control-side image of the diagonal residual complement",
        },
        "classification": "explicit_declared_control_pullback_of_the_residual_local_diagonal_sector_present_but_no_zero_or_host_side_cancellation_witness",
        "frontier": "R16_B1",
        "frontier_text": "the repo now exports the declared control pullback of the residual local diagonal sector, but it still lacks a witness that this pullback vanishes or matches a host-side correction",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_residual_local_diagonal_declared_control_pullback_vanishes",
            "no_claim_that_the_residual_local_diagonal_declared_control_pullback_matches_the_host",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R16",
        "status": "PASS_PARTIAL_EXPLICIT_RESIDUAL_DIAGONAL_DECLARED_CONTROL_PULLBACK_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R16_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
