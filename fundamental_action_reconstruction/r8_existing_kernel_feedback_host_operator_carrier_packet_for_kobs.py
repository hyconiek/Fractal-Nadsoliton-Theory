#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs.json"
OUT_SUMMARY = GENERATED / "r8_existing_kernel_feedback_host_operator_carrier_packet_for_kobs_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2186 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2186_ktotal_spectral_stability_margin_gate.json"
    )
    c9 = load_json("fundamental_action_reconstruction/generated/c9_action_origin_host_carrier_audit_summary.json")
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")
    c14 = load_json("fundamental_action_reconstruction/generated/c14_control_mode_to_psi_transport_schema_summary.json")
    p11 = load_json(
        "fundamental_action_reconstruction/generated/p11_existing_kernel_feedback_to_explicit_h3_chain_factorization_map_probe.json"
    )

    psi_basis = c14["control_transport"]["psi_basis_indices"]
    matrix_metrics = q2186["matrix"]

    checks = [
        {
            "id": "p11_originally_missing_explicit_legacy_carrier",
            "actual": "explicit_operator_level_existing_kernel_feedback_carrier_with_declared_basis_or_finite_state_space"
            in p11["missing_upstream_objects"],
            "expected": True,
            "meaning": "P11 indeed localized the explicit operator-level legacy carrier as one of the four missing factorization subobjects",
        },
        {
            "id": "qw2186_strict_branch_scope_verdict",
            "actual": q2186["verdict"],
            "expected": "KTOTAL_SPECTRAL_STABILITY_MARGIN_GATE_PASS_STRICT_BRANCH_SCOPE",
            "meaning": "QW-2186 exports a strict branch-scope certified host operator",
        },
        {
            "id": "qw2186_has_finite_12_octave_host_carrier",
            "actual": int(matrix_metrics["n_octaves"]),
            "expected": 12,
            "meaning": "QW-2186 exports a finite host carrier with 12 octave slots",
        },
        {
            "id": "c9_host_operator_present_branch_scope",
            "actual": c9["carrier_schema"]["host_operator"],
            "expected": "QW-2186 exports a branch-scope positive host operator A = K_total + m0^2 I.",
            "meaning": "C9 already identifies the QW-2186 object as a real host operator",
        },
        {
            "id": "c10_psi_sector_carrier_schema_present",
            "actual": c10["result"]["psi_sector_carrier_schema_present"],
            "expected": "yes",
            "meaning": "C10 already places the host inside the canonical Psi-sector carrier family",
        },
        {
            "id": "c10_host_to_concrete_block_still_not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "expected": "not_shown",
            "meaning": "the host carrier is present, but concrete block identification is still absent",
        },
        {
            "id": "c14_declared_psi_basis_available",
            "actual": psi_basis,
            "expected": [
                "psi0",
                "psi1",
                "psi2",
                "psi3",
                "psi4",
                "psi5",
                "psi6",
                "psi7",
                "psi8",
                "psi9",
                "psi10",
                "psi11",
            ],
            "meaning": "C14 already exports a declared canonical Psi-basis order on the 12-slot carrier",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R8",
        "operator_scope": "existing_kernel_feedback_host_operator_carrier_packet",
        "source_reports": ["QW-2186", "C9", "C10", "C14", "P11"],
        "host_operator_schema": {
            "symbol": "A = K_total + m0^2 I",
            "existing_feedback_role": "the existing kernel feedback lives inside the K_total component of the certified host operator",
            "scope": "host_scope_branch_level",
            "certificate_verdict": q2186["verdict"],
        },
        "declared_host_carrier": {
            "carrier_label": "Psi_host_12",
            "finite_state_space_size": int(matrix_metrics["n_octaves"]),
            "basis_order": psi_basis,
            "basis_status": "declared_carrier_basis_only_not_selector_sector_projection",
            "carrier_family_statement": c10["carrier_schema"]["canonical_hessian"],
        },
        "spectral_certificate": {
            "required_shift_ge": float(matrix_metrics["required_shift_ge"]),
            "broken_floor_used": float(matrix_metrics["broken_floor_used"]),
            "lambda_min_A": float(matrix_metrics["lambda_min_A"]),
            "lambda_max_A": float(matrix_metrics["lambda_max_A"]),
            "condition_number_2": float(matrix_metrics["condition_number_2"]),
        },
        "construction_rule": {
            "type": "host_scope_packetization_of_existing_feedback_carrier",
            "rule": "promote the already certified QW-2186 host operator into an explicit factorization-lane carrier packet by attaching its finite host size and declared canonical Psi-basis order, while keeping the packet strictly at host scope",
        },
        "consistency_checks": checks,
        "projection_to_explicit_h3_chain_present": False,
        "selector_sector_reduction_present": False,
        "intertwiner_or_equality_witness_present": False,
        "host_scope_only": True,
        "factorization_status": "host_scope_operator_level_existing_kernel_feedback_carrier_present_but_not_yet_projected_to_the_explicit_h3_chain",
        "classification": "explicit_operator_level_existing_kernel_feedback_host_carrier_present_but_host_scope_only_and_not_yet_projected_to_the_explicit_h3_chain",
        "frontier": "R8_B1",
        "no_false_pass": True,
    }

    summary = {
        "stage": "R8",
        "status": "PASS_PARTIAL_EXPLICIT_OPERATOR_LEVEL_EXISTING_KERNEL_FEEDBACK_HOST_CARRIER_PACKET_READY",
        "result": artifact["classification"],
        "frontier": [
            "R8_B1",
            "H14_B1",
            "H15_B1",
            "H16_B1",
            "C10_B1",
            "C9_B2",
        ],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
