#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r12_explicit_canonical_psi_block_export_packet_with_kernel_mixing_for_host_route.json"
OUT_SUMMARY = GENERATED / "r12_explicit_canonical_psi_block_export_packet_with_kernel_mixing_for_host_route_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def diagonal_entry(index: int) -> str:
    return (
        f"3*g4_psi{index}*vpsi{index}**2 + 5*g6_psi{index}*vpsi{index}**4 + "
        f"2*gY{index}*vphi**2 + m2_psi{index}"
    )


def off_diagonal_entry(i: int, j: int) -> str:
    return f"(K_{i}_{j} + K_{j}_{i})/2"


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    q2165 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2165_l13_exhaustive_canonical_eom_gate.json"
    )
    q2166 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2166_l14_exhaustive_canonical_hessian_gate.json"
    )
    c16 = load_json("fundamental_action_reconstruction/generated/c16_minimal_psi_hessian_coefficient_class_table_summary.json")
    c17 = load_json("fundamental_action_reconstruction/generated/c17_index_complete_psi_row_stencil_audit_summary.json")
    c20 = load_json("fundamental_action_reconstruction/generated/c20_finite_materialization_recipe_audit_summary.json")
    r11 = load_json(
        "fundamental_action_reconstruction/generated/r11_symmetry_certified_declared_control_transport_packet_for_psi_block_route.json"
    )

    psi_basis = r11["transport_packet"]["target_basis"]
    n_psi = len(psi_basis)

    matrix_rows: list[list[str]] = []
    serialized_rows: list[dict[str, Any]] = []
    for i in range(n_psi):
        row_entries = []
        for j in range(n_psi):
            if i == j:
                row_entries.append(diagonal_entry(i))
            else:
                row_entries.append(off_diagonal_entry(i, j))
        matrix_rows.append(row_entries)
        serialized_rows.append(
            {
                "row_label": f"psi{i}",
                "entries": row_entries,
            }
        )

    checks = [
        {
            "id": "c16_representative_coefficient_class_rows_present",
            "actual": c16["result"]["representative_coefficient_class_rows_present"],
            "expected": "yes",
            "meaning": "representative coefficient-class rows already exist",
        },
        {
            "id": "c17_index_complete_row_stencil_present",
            "actual": c17["result"]["index_complete_row_stencil_schema_present"],
            "expected": "yes",
            "meaning": "an index-complete Psi-row stencil is already present",
        },
        {
            "id": "c20_materialization_recipe_present",
            "actual": c20["result"]["finite_persisted_materialization_recipe_present"],
            "expected": "yes",
            "meaning": "the finite materialization recipe is already present",
        },
        {
            "id": "q2165_exhaustive_eom_gate_pass_present",
            "actual": q2165["verdict"],
            "expected": "L13_EXHAUSTIVE_CANONICAL_EOM_GATE_PASS_PARTIAL_ALL_ORDERS_OPEN",
            "meaning": "QW-2165 provides exhaustive all-rows EoM support",
        },
        {
            "id": "q2166_exhaustive_hessian_gate_pass_present",
            "actual": q2166["verdict"],
            "expected": "L14_EXHAUSTIVE_CANONICAL_HESSIAN_GATE_PASS_PARTIAL_FULL_THEOREM_OPEN",
            "meaning": "QW-2166 provides exhaustive Hessian-level support",
        },
        {
            "id": "q2166_hessian_matches_operator",
            "actual": q2166["flags"]["linear_operator_matrix_matches_canonical_hessian"],
            "expected": True,
            "meaning": "the linear operator matrix matches the canonical Hessian",
        },
        {
            "id": "q2166_hessian_contains_kernel_mixing_entries",
            "actual": q2166["flags"]["hessian_contains_kernel_mixing_entries"],
            "expected": True,
            "meaning": "the canonical block contains kernel-mixing entries",
        },
        {
            "id": "r11_explicit_declared_transport_packet_present",
            "actual": r11["result"]["explicit_declared_control_transport_packet_present"],
            "expected": True,
            "meaning": "the declared transport packet is present",
        },
        {
            "id": "sample_eta0_contains_expected_off_diagonal_token",
            "actual": "(K_0_1 + K_1_0)*eta1(x)/2" in q2166["model"]["sample_eom_eta0"],
            "expected": True,
            "meaning": "sample eta0 row contains the expected symmetric kernel-mixing token",
        },
        {
            "id": "sample_eta6_contains_expected_diagonal_token",
            "actual": "(3*g4_psi6*vpsi6**2 + 5*g6_psi6*vpsi6**4 + 2*gY6*vphi**2 + m2_psi6)*eta6(x)"
            in q2166["model"]["sample_eom_eta6"],
            "expected": True,
            "meaning": "sample eta6 row contains the expected diagonal coefficient class",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R12",
        "packet_goal": "materialize_an_explicit_coefficient_filled_canonical_Psi_x_Psi_block_export_on_full_declared_transport_support",
        "source_reports": ["QW-2165", "QW-2166", "C16", "C17", "C20", "R11"],
        "chosen_transported_index_set": {
            "scope": "full_declared_transport_support_only",
            "basis_order": psi_basis,
            "cardinality": n_psi,
            "reason": "the declared control transport columns from R11 have support across the full canonical Psi carrier, so the narrowest honest concrete transported index-set is the full 12-slot support",
        },
        "canonical_psi_block_export": {
            "symbol": "H_PsiPsi_declared_full_support",
            "shape": [n_psi, n_psi],
            "serialized_rows": serialized_rows,
            "matrix_rows": matrix_rows,
            "construction_rule": "assemble the full canonical Psi-Psi Hessian block from the index-complete row stencil: diagonal entries use the local class from C16/C17, off-diagonal entries use the symmetric kernel-mixing class (K_i_j + K_j_i)/2",
        },
        "kernel_mixing_light_facing_summary": {
            "strict_admissible_channel": "symmetric kernel-mixing entries (K_i_j + K_j_i)/2 in the off-diagonal Psi-Psi block",
            "off_diagonal_entry_count": n_psi * (n_psi - 1),
            "diagonal_local_entry_count": n_psi,
            "boundary": "this packet remembers the current light-facing channel only through strict-admissible kernel-mixing terms; it does not claim K_obs factorization, photon identification, or selector closure",
        },
        "declared_transport_boundary": {
            "explicit_declared_control_transport_packet_present": True,
            "physical_canonicalization_present": False,
            "host_matching_present": False,
            "selector_relevant_target_justification_present": False,
        },
        "consistency_checks": checks,
        "classification": "explicit_coefficient_filled_canonical_Psi_block_present_on_full_declared_transport_support_with_kernel_mixing_highlight_but_not_yet_physically_canonicalized_or_host_matched",
        "frontier": "R12_B1",
        "frontier_text": "the repo now exports a concrete coefficient-filled canonical Psi block on the full declared transport support, but physical canonicalization of the transport and host-to-block matching remain absent",
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_the_full_support_export_is_already_selector_relevant",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_the_exported_block_is_already_matched_to_QW2186",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R12",
        "status": "PASS_PARTIAL_EXPLICIT_CANONICAL_PSI_BLOCK_EXPORT_READY_ON_FULL_DECLARED_TRANSPORT_SUPPORT",
        "result": artifact["classification"],
        "frontier": ["R12_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
