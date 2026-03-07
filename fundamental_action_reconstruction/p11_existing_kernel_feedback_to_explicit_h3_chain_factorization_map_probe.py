#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p11_existing_kernel_feedback_to_explicit_h3_chain_factorization_map_probe.json"
OUT_SUMMARY = GENERATED / "p11_existing_kernel_feedback_to_explicit_h3_chain_factorization_map_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    sources = {
        "H14": load_json("fundamental_action_reconstruction/generated/h14_existing_kernel_feedback_vs_new_kobs_separation_audit_summary.json"),
        "H15": load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json"),
        "H16": load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json"),
        "R2": load_json("fundamental_action_reconstruction/generated/r2_existing_internal_feedback_parameter_packet_for_kobs_summary.json"),
        "R7": load_json("fundamental_action_reconstruction/generated/r7_shared_frozen_kernel_provenance_packet_for_kobs_summary.json"),
        "P10": load_json("fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"),
    }

    route_checks = [
        {
            "id": "existing_feedback_recognized_but_not_identified",
            "pass": sources["H14"]["status"] == "PASS_PARTIAL_EXISTING_FEEDBACK_RECOGNIZED_BUT_NOT_IDENTIFIED_WITH_KOBS",
            "expected": "PASS_PARTIAL_EXISTING_FEEDBACK_RECOGNIZED_BUT_NOT_IDENTIFIED_WITH_KOBS",
            "actual": sources["H14"]["status"],
            "meaning": "existing kernel feedback is recognized, but no identification map exists yet",
        },
        {
            "id": "existing_feedback_has_no_selector_sector_export",
            "pass": sources["H15"]["result"] == "existing_feedback_not_identified_with_kobs",
            "expected": "existing_feedback_not_identified_with_kobs",
            "actual": sources["H15"]["result"],
            "meaning": "existing kernel feedback still has no exported selector-sector reduction of its own",
        },
        {
            "id": "operator_origin_has_only_partial_witnesses",
            "pass": sources["H16"]["result"] == "admissible_values_have_asymmetric_partial_witnesses_only",
            "expected": "admissible_values_have_asymmetric_partial_witnesses_only",
            "actual": sources["H16"]["result"],
            "meaning": "the repository still has only partial operator-origin witnesses and no provenance-valid factorization",
        },
        {
            "id": "r2_parameter_packet_present",
            "pass": sources["R2"]["result"]
            == "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "expected": "existing_internal_feedback_parameter_packet_present_but_not_yet_an_operator_level_kobs_factorization",
            "actual": sources["R2"]["result"],
            "meaning": "the old light/matter/observer data is packetized but not factorized",
        },
        {
            "id": "r7_shared_provenance_present",
            "pass": sources["R7"]["result"] == "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
            "expected": "shared_frozen_kernel_provenance_present_but_not_operator_factorization",
            "actual": sources["R7"]["result"],
            "meaning": "the explicit chain and old QW family share frozen-kernel provenance",
        },
        {
            "id": "p10_full_current_pair_block_present_but_nonidentified",
            "pass": sources["P10"]["status"] == "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "actual": sources["P10"]["status"],
            "meaning": "the explicit current-pair H3 block exists, but identification with existing kernel feedback is still absent",
        },
    ]

    route_state = {
        "shared_frozen_kernel_provenance_present": True,
        "existing_internal_feedback_parameter_packet_present": True,
        "explicit_current_pair_h3_chain_present": True,
        "explicit_current_pair_h3_block_present": True,
        "explicit_operator_level_existing_kernel_feedback_carrier_present": False,
        "typed_projection_from_existing_feedback_to_explicit_chain_present": False,
        "existing_feedback_selector_sector_reduction_present": False,
        "intertwiner_or_equality_witness_present": False,
        "equivalence_factorization_map_present": False,
    }

    missing_upstream_objects = [
        "explicit_operator_level_existing_kernel_feedback_carrier_with_declared_basis_or_finite_state_space",
        "typed_projection_or_pushforward_map_from_existing_kernel_feedback_into_the_explicit_H3_slot_chain_or_current_pair_block",
        "selector_sector_reduction_of_existing_kernel_feedback_onto_pair1_or_an_equivalent_actual_pair_target",
        "intertwiner_or_equality_witness_identifying_the_reduced_existing_feedback_object_with_the_computed_current_pair_H3_block",
    ]

    report = {
        "stage": "P11",
        "goal": "compute_or_fail_existing_kernel_feedback_to_explicit_h3_chain_factorization_map",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_TO_EXPLICIT_CHAIN_FACTORIZATION_ROUTE",
        "reason": "the repo now has shared frozen-kernel provenance and a fully explicit current-pair H3 chain with a computed selector-facing block, but it still lacks an operator-level existing-kernel-feedback carrier, a typed projection into that chain, a selector-sector reduction on the legacy side, and an intertwiner/equality witness identifying the resulting reduced legacy object with the computed block",
        "lane": "existing_kernel_feedback_to_explicit_h3_chain_factorization_route_after_R7",
        "route_under_test": [
            "existing_kernel_feedback_inside_K_total",
            "shared_frozen_kernel_provenance_packet",
            "R2_parameter_packet",
            "explicit_current_pair_H3_chain",
            "computed_current_pair_H3_block",
            "equivalence_or_factorization_map",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "existing_kernel_feedback_inside_K_total",
            "R2_internal_feedback_parameter_packet",
            "R7_shared_frozen_kernel_provenance_packet",
            "explicit_current_pair_H3_chain_packets",
            "computed_current_pair_H3_block",
        ],
        "decomposition_of_P10_missing_object": {
            "from_P10": "equivalence_or_factorization_map_from_existing_kernel_feedback_and_R2_parameter_packet_to_H3_operator_chain",
            "into_current_blockers": missing_upstream_objects,
        },
        "missing_upstream_objects": missing_upstream_objects,
        "blocking_frontier": {
            "R7_B1": "shared frozen-kernel provenance is present, but no operator-level factorization object is exported from that provenance",
            "H14_B1": sources["H14"]["frontier"]["H14_B1"],
            "H15_B1": "existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback",
            "H16_B1": "both admissible operator-origin values have only asymmetric partial witnesses and still no provenance-valid factorization route",
        },
        "required_next_step": "IMPLEMENT_ONE_OF_THE_FOUR_FACTORIZATION_SUBOBJECTS_OR_FORMALIZE_A_ROUTE_SPECIFIC_NONFACTORIZATION_THEOREM_BEFORE_CLAIMING_EXISTING_KERNEL_FEEDBACK_ALREADY_CONTAINS_SELECTOR_FACING_KOBS",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": report["stage"],
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P10_missing_object": report["decomposition_of_P10_missing_object"],
        "missing_upstream_objects": report["missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
