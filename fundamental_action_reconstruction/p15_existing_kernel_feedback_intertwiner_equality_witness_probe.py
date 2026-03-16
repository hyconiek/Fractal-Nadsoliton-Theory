#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p15_existing_kernel_feedback_intertwiner_equality_witness_probe.json"
OUT_SUMMARY = GENERATED / "p15_existing_kernel_feedback_intertwiner_equality_witness_probe_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p10 = load_json(
        "fundamental_action_reconstruction/generated/p10_existing_kernel_feedback_to_kobs_rerun_after_explicit_current_pair_chain_summary.json"
    )
    p14 = load_json(
        "fundamental_action_reconstruction/generated/p14_existing_kernel_feedback_factorization_rerun_after_current_pair_chart_reduction_packet.json"
    )
    r10 = load_json(
        "fundamental_action_reconstruction/generated/r10_legacy_control_to_current_pair_chart_reduction_packet_for_kobs_summary.json"
    )
    c10 = load_json("fundamental_action_reconstruction/generated/c10_psi_sector_host_identification_audit_summary.json")
    c15 = load_json("fundamental_action_reconstruction/generated/c15_control_only_pullback_submatrix_packet_summary.json")
    h15 = load_json("fundamental_action_reconstruction/generated/h15_existing_feedback_selector_sector_reduction_audit_summary.json")
    h16 = load_json("fundamental_action_reconstruction/generated/h16_operator_origin_partial_witness_audit_summary.json")
    h18 = load_json("fundamental_action_reconstruction/generated/h18_composite_route_a_provenance_binding_instance_summary.json")
    h18_artifact = load_json("fundamental_action_reconstruction/generated/route_a_provenance_valid_binding_instance.json")
    o2 = load_json("fundamental_action_reconstruction/generated/o2_exported_composite_a1_ext_instance.json")

    p14_missing = p14["remaining_missing_upstream_objects"]
    remaining_missing = [
        "explicit_coefficient_filled_legacy_chart_reduced_operator_object_on_the_chosen_current_pair_chart_pair1_or_equivalent_actual_target",
        "intertwiner_or_equality_witness_identifying_that_legacy_chart_reduced_operator_object_with_the_computed_current_pair_H3_block",
    ]

    route_checks = [
        {
            "id": "p14_originally_missing_intertwiner_or_equality_witness",
            "pass": "intertwiner_or_equality_witness_identifying_the_chart_reduced_legacy_object_with_the_computed_current_pair_H3_block"
            in p14_missing,
            "expected": True,
            "actual": "intertwiner_or_equality_witness_identifying_the_chart_reduced_legacy_object_with_the_computed_current_pair_H3_block"
            in p14_missing,
            "meaning": "P14 indeed reduced the factorization route to a single nominal witness object",
        },
        {
            "id": "r10_current_pair_chart_reduction_present",
            "pass": r10["result"]
            == "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "expected": "explicit_current_pair_chart_reduction_present_but_chart_scoped_only_and_not_a_strict_selector_target_justification",
            "actual": r10["result"],
            "meaning": "the legacy route reaches the chosen current-pair chart only at map level",
        },
        {
            "id": "p10_current_pair_h3_block_present",
            "pass": p10["status"] == "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "expected": "CURRENT_PAIR_H3_BLOCK_COMPUTABLE_BUT_NOT_IDENTIFIED_WITH_EXISTING_KERNEL_FEEDBACK",
            "actual": p10["status"],
            "meaning": "the explicit current-pair H3 block is present on the target side",
        },
        {
            "id": "c10_host_to_concrete_psi_block_identification_not_shown",
            "pass": c10["result"]["host_to_concrete_psi_block_identification"] == "not_shown",
            "expected": "not_shown",
            "actual": c10["result"]["host_to_concrete_psi_block_identification"],
            "meaning": "the legacy host operator is still not identified with a concrete Psi-sector block",
        },
        {
            "id": "c15_coefficient_filled_control_pullback_present",
            "pass": c15["result"]["coefficient_filled_M_control_present"] == "yes",
            "expected": "yes",
            "actual": c15["result"]["coefficient_filled_M_control_present"],
            "meaning": "a coefficient-filled declared control pullback M_control is exported (declared scope; see P476 pointer in C15)",
        },
        {
            "id": "h15_existing_feedback_not_identified_with_kobs",
            "pass": h15["result"] == "existing_feedback_not_identified_with_kobs",
            "expected": "existing_feedback_not_identified_with_kobs",
            "actual": h15["result"],
            "meaning": "existing kernel feedback still has no selector-facing export identified as its own reduction",
        },
        {
            "id": "h16_operator_origin_route_still_partial_only",
            "pass": h16["result"] == "admissible_values_have_asymmetric_partial_witnesses_only",
            "expected": "admissible_values_have_asymmetric_partial_witnesses_only",
            "actual": h16["result"],
            "meaning": "the operator-origin route still has only asymmetric partial witnesses at the audit level",
        },
        {
            "id": "h18_extension_lane_route_a_witness_exists",
            "pass": h18["result"] == "route_a_binding_completed_but_coefficients_unevaluated"
            and bool(h18_artifact["provenance_valid_route_a_instance"]),
            "expected": True,
            "actual": h18["result"] == "route_a_binding_completed_but_coefficients_unevaluated"
            and bool(h18_artifact["provenance_valid_route_a_instance"]),
            "meaning": "the strongest extension-lane composite witness exists, but only as a provenance-valid unevaluated Route A instance",
        },
        {
            "id": "o2_exported_composite_matrix_entries_still_unresolved",
            "pass": (o2["computable_now"] is False)
            and o2["coefficient_state"] == {
                "a_1": "UNRESOLVED_VALUE",
                "b_1": "UNRESOLVED_VALUE",
                "d_1": "UNRESOLVED_VALUE",
            },
            "expected": True,
            "actual": (o2["computable_now"] is False)
            and o2["coefficient_state"] == {
                "a_1": "UNRESOLVED_VALUE",
                "b_1": "UNRESOLVED_VALUE",
                "d_1": "UNRESOLVED_VALUE",
            },
            "meaning": "the strongest exported composite object is still not matrix-evaluated on pair1",
        },
    ]

    route_state = {
        "chosen_current_pair_chart_reduction_present": True,
        "computed_current_pair_h3_block_present": True,
        "provenance_valid_extension_lane_route_a_witness_present": True,
        "evaluated_exported_composite_a1_instance_present": False,
        "coefficient_filled_legacy_chart_reduced_operator_present": False,
        "intertwiner_or_equality_witness_present": False,
        "existing_kernel_feedback_identified_with_computed_current_pair_block": False,
    }

    report = {
        "stage": "P15",
        "goal": "compute_or_fail_existing_kernel_feedback_intertwiner_or_equality_witness_after_current_pair_chart_reduction",
        "status": "NOT_COMPUTABLE_FROM_CURRENT_EXISTING_KERNEL_FEEDBACK_INTERTWINER_EQUALITY_WITNESS_ROUTE",
        "reason": "the repo now reaches the chosen current-pair chart and contains a computed current-pair H3 block; it also exports a coefficient-filled declared control pullback M_control (C15 via P476), but it still exports neither a coefficient-filled legacy-side chart-reduced operator object identified as existing kernel feedback on that chart nor an intertwiner/equality witness identifying such a legacy object with the computed block; the strongest extension-lane composite witness also remains unevaluated",
        "lane": "existing_kernel_feedback_intertwiner_equality_witness_route_after_R10",
        "route_under_test": [
            "existing_kernel_feedback_inside_K_total",
            "host_scope_operator_level_legacy_carrier",
            "typed_host_to_control_pushforward",
            "chosen_current_pair_chart_reduction",
            "explicit_coefficient_filled_legacy_chart_reduced_operator_object",
            "computed_current_pair_H3_block",
            "intertwiner_or_equality_witness",
        ],
        "route_checks": route_checks,
        "route_state": route_state,
        "supporting_present_but_insufficient_objects": [
            "R10_current_pair_chart_reduction_packet",
            "P10_computed_current_pair_H3_block",
            "C15_coefficient_filled_declared_control_pullback_M_control (declared; not host-matched)",
            "H18_provenance_valid_extension_lane_route_A_witness",
            "O2_exported_composite_A_1_instance",
        ],
        "decomposition_of_P14_missing_object": {
            "from_P14": "intertwiner_or_equality_witness_identifying_the_chart_reduced_legacy_object_with_the_computed_current_pair_H3_block",
            "into_current_blockers": remaining_missing,
        },
        "remaining_missing_upstream_objects": remaining_missing,
        "computed_current_pair_block": {
            "a_1": p10["computed_a_1"],
            "b_1": p10["computed_b_1"],
            "d_1": p10["computed_d_1"],
        },
        "blocking_frontier": {
            "C10_B1": c10["residual_blockers"]["C10_B1"],
            "H15_B1": "existing kernel feedback has no explicit residual-selector-sector reduction or projected selector-block export in the current repository, so K_obs remains a distinct extension hypothesis rather than an identified reformulation of existing kernel feedback",
            "H16_B1": "both admissible operator-origin values still have only asymmetric partial witnesses and no provenance-valid equality witness to the computed current-pair block",
            "H18_B1": "a provenance-valid Route A witness exists on the hypothesis-extension lane for pair1, but no evaluated coefficient triple is exported from it",
            "O2_B1": "a persisted exported_composite_A_1 instance exists on pair1, but its coefficient entries remain unresolved and unevaluated",
        },
        "required_next_step": "EITHER_EXPORT_THE_COEFFICIENT_FILLED_LEGACY_CHART_REDUCED_OPERATOR_OBJECT_ON_PAIR1_FIRST_OR_KEEP_THE_FACTORIZATION_ROUTE_NEGATIVE_WITHOUT_CLAIMING_ANY_INTERTWINER_OR_EQUALITY_WITNESS",
        "strict_core_promotion": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P15",
        "status": report["status"],
        "reason": report["reason"],
        "decomposition_of_P14_missing_object": report["decomposition_of_P14_missing_object"],
        "remaining_missing_upstream_objects": report["remaining_missing_upstream_objects"],
        "required_next_step": report["required_next_step"],
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(json.dumps(summary, ensure_ascii=True))


if __name__ == "__main__":
    main()
