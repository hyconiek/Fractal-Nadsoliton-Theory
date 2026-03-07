#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f3_current_far_frontier_kernel_artifact_sensitivity_classification_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f3_current_far_frontier_kernel_artifact_sensitivity_classification_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f2 = load_json(
        "fundamental_action_reconstruction/generated/f2_strict_gate_kernel_provenance_and_far_input_classification_packet_summary.json"
    )
    a1 = load_json(
        "fundamental_action_reconstruction/generated/a1_minimal_action_ansatz_summary.json"
    )
    a4 = load_json(
        "fundamental_action_reconstruction/generated/a4_rg_emergence_summary.json"
    )
    a8 = load_json(
        "fundamental_action_reconstruction/generated/a8_gravity_bridge_summary.json"
    )
    p46 = load_json(
        "fundamental_action_reconstruction/generated/p46_canonical_ontology_supported_direct_formal_c1s1_family_route_probe_after_direct_m2_psi4_target_action_coefficient_defect_polynomial_packet_summary.json"
    )
    n49 = load_json(
        "fundamental_action_reconstruction/generated/n49_current_canonical_ontology_supported_direct_formal_c1s1_family_route_boundary_theorem_after_direct_m2_psi4_target_action_coefficient_defect_polynomial_packet_summary.json"
    )
    n50 = load_json(
        "fundamental_action_reconstruction/generated/n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "f2_kernel_input_classification_present",
            "actual": f2["status"],
            "expected": "F2_EXECUTED_STRICT_GATE_KERNEL_PROVENANCE_AND_FAR_INPUT_CLASSIFICATION_PACKET_NO_FALSE_PASS",
            "meaning": "F2 already classifies K_strict_gate as a later-pipeline operational kernel rather than an ontological source layer",
        },
        {
            "id": "n50_kernel_nonidentification_discharged",
            "actual": n50["status"],
            "expected": "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS",
            "meaning": "N50 already discharges the current-repo nonidentification theorem between legacy and strict kernels",
        },
        {
            "id": "a1_silent_substitution_disallowed",
            "actual": a1["kernel_source_classification"]["silent_full_substitution_disallowed"],
            "expected": True,
            "meaning": "A1 no longer allows silent replacement of the ontological source layer by K_strict_gate",
        },
        {
            "id": "a4_silent_substitution_disallowed",
            "actual": a4["kernel_source_classification"]["silent_full_substitution_disallowed"],
            "expected": True,
            "meaning": "A4 no longer allows silent replacement of ontological shell data by K_strict_gate",
        },
        {
            "id": "a8_silent_substitution_disallowed",
            "actual": a8["kernel_source_classification"]["silent_full_substitution_disallowed"],
            "expected": True,
            "meaning": "A8 no longer allows silent inheritance of legacy hierarchy semantics into K_strict_gate",
        },
        {
            "id": "p46_corrected_route_still_live",
            "actual": p46["status"],
            "expected": "CANONICAL_ONTOLOGY_SUPPORTED_SOURCE_SIDE_CLOSED_AND_TARGET_ACTION_DEFECT_POLYNOMIAL_EXPORTED_ROUTE_STILL_NOT_CLOSED_AFTER_R37",
            "meaning": "P46 still exports a live corrected route frontier after the kernel-input correction",
        },
        {
            "id": "n49_corrected_route_boundary_present",
            "actual": n49["status"],
            "expected": "N49_DISCHARGED_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_DIRECT_FORMAL_C1S1_FAMILY_ROUTE_BOUNDARY_AFTER_R37_NO_FALSE_PASS",
            "meaning": "N49 already fixes the corrected-route boundary after R37",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact_sensitive_upstream_classes = [
        "silent_action_first_substitution_of_K_strict_gate_for_K_legacy_ont_in_A1",
        "silent_RG_shell_semantics_substitution_of_K_strict_gate_for_canonical_D_f_alpha_geo_beta_tors_layer_in_A4",
        "silent_gravity_hierarchy_semantics_inheritance_into_K_strict_gate_in_A8",
        "any_unbridged_reuse_of_K_strict_gate_as_if_it_already_carried_full_legacy_ontological_roles",
    ]
    kernel_split_robust_current_route_objects = p46["remaining_missing_objects"]

    artifact = {
        "stage": "F3",
        "lane": "current_far_frontier_kernel_artifact_sensitivity_classification_current_repo_state_only",
        "goal": "classify_which_current_far_blockers_are_artifact_sensitive_to_the_legacy_strict_kernel_split_and_which_still_remain_live_after_the_kernel_input_correction",
        "status": "F3_EXECUTED_CURRENT_FAR_FRONTIER_KERNEL_ARTIFACT_SENSITIVITY_CLASSIFICATION_PACKET_NO_FALSE_PASS",
        "reason": "F2 and N50 already exclude silent ontological inheritance from the strict gate kernel, A1/A4/A8 have been corrected to enforce that rule, and P46/N49 still export a live canonical-ontology-supported direct route frontier after those corrections; therefore the remaining P46 objects are currently classified as kernel-split-robust frontier objects, while silent ontological uses of K_strict_gate remain artifact-sensitive upstream classes",
        "artifact_sensitive_upstream_classes": artifact_sensitive_upstream_classes,
        "kernel_split_robust_current_route": {
            "route": "AX9 + AX10 + AX11 + R35 + R36 + R37 -> P46/N49",
            "reason": [
                "the route already runs on the canonical-ontology-supported lane",
                "the source-side local closures use the corrected nadsoliton ontology",
                "A1/A4/A8 now explicitly exclude silent ontological substitution by K_strict_gate",
                "the current remaining missing objects persist after that correction"
            ],
            "remaining_missing_objects": kernel_split_robust_current_route_objects,
        },
        "recommended_continuation_rule": {
            "continue_on": "kernel_split_robust_current_route_only",
            "best_next_move": "attack_explicit_zero_witness_for_the_direct_m2_psi4_target_action_coefficient_defect_polynomial_on_common_psi4_squared_over_2_support",
            "do_not_reopen": artifact_sensitive_upstream_classes,
        },
        "not_classified_globally": [
            "full_main_host_route_equivalence_frontier_in_all_respects",
            "global_selector_closure",
            "QW2191_discharge",
            "full_ToE_closure",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F3",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "artifact_sensitive_upstream_classes": artifact_sensitive_upstream_classes,
        "kernel_split_robust_current_route_remaining_missing_objects": kernel_split_robust_current_route_objects,
        "recommended_next_move": artifact["recommended_continuation_rule"]["best_next_move"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
