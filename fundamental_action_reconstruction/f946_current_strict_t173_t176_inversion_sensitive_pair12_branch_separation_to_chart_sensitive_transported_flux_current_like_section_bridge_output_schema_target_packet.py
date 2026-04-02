#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F945 = GENERATED / "f945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_target_packet.json"
IN_P946_CANDIDATE = GENERATED / "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_candidate_nonexport_audit_probe.json"
IN_P946_CLASS = GENERATED / "p946_current_strict_t173_t176_output_schema_class_candidate_supported_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_blocked.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"
IN_T218 = ROOT / "T218_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_SPEC.md"

OUT = GENERATED / "f946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_packet.json"
OUT_SUMMARY = GENERATED / "f946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_packet_summary.json"

F945_STATUS = (
    "F945_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
)
P946_CANDIDATE_STATUS = (
    "PASS_T173_T176_INVERSION_SENSITIVE_PAIR12_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_CANDIDATE_NONEXPORT_AUDITED"
)
P946_CLASS_STATUS = (
    "P946_CURRENT_STRICT_T173_T176_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_BLOCKED"
)
TARGET_OBJECT_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_v1"
)
EXPECTED_SUPPORT_INTERFACE_NAME = (
    "chart_sensitive_pair12_typed_descent_from_Sigma_sel_src_target_v1_to_the_surviving_F301_pair12_carrier_without_Q_basis_sel_v1_terminal_collapse_and_without_projector_only_atlas_collapse"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F945, IN_P946_CANDIDATE, IN_P946_CLASS, IN_P763, IN_P764, IN_T218]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F946",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f945 = load_json(IN_F945)
    p946_candidate = load_json(IN_P946_CANDIDATE)
    p946_class = load_json(IN_P946_CLASS)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)
    t218_text = load_text(IN_T218)

    p946_candidate_conclusion = p946_candidate.get("audit_conclusion") or {}
    p946_class_theorem = p946_class.get("theorem_result") or {}
    p946_class_support_stack = p946_class.get("output_schema_class_support_stack") or {}

    preconditions = {
        "f945_target_packet_executed": f945.get("status") == F945_STATUS,
        "p946_candidate_probe_passed": p946_candidate.get("status") == P946_CANDIDATE_STATUS,
        "p946_class_probe_passed": p946_class.get("status") == P946_CLASS_STATUS,
        "no_actual_bridge_output_schema_supplier_exported": (
            p946_candidate_conclusion.get("current_repo_already_exports_actual_bridge_output_schema_supplier")
            is False
        ),
        "support_class_candidate_grade_present": (
            p946_class_theorem.get("bridge_output_schema_class_candidate_supported_now") is True
        ),
        "exact_bridge_output_schema_still_unexported": (
            p946_class_theorem.get("bridge_output_schema_exported_now") is False
        ),
        "sharp_blocker_is_bridge_output_schema": (
            p946_class_theorem.get("sharp_blocker_field") == "bridge_output_schema"
        ),
        "future_only_support_interface_target_candidate_present": (
            p946_candidate_conclusion.get("current_repo_already_exports_relevant_future_only_interface_target_candidate")
            is True
        ),
        "support_interface_name_matches_p763_boundary": (
            p946_candidate_conclusion.get("current_best_support_interface_name") == EXPECTED_SUPPORT_INTERFACE_NAME
            and p763.get("exact_named_missing_interface") == EXPECTED_SUPPORT_INTERFACE_NAME
        ),
        "p763_interface_still_unexported": (
            p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported") is True
        ),
        "p764_target_future_only_and_below_t176": (
            p764.get("t218_target_exported_on_current_repo_state") is True
            and p764.get("current_t218_target_is_future_route_only") is True
            and p764.get("current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge")
            is True
        ),
        "t218_scope_still_route_local_not_full_c_v1_bridge": (
            "Sigma_sel_src_target_v1" in t218_text
            and "target_ends_in_surviving_pair12_residual_datum_carrier_lane := yes" in t218_text
            and "future_route_only := yes" in t218_text
            and "target_remains_below_global_t176_discharge := yes" in t218_text
        ),
    }
    all_preconditions_pass = all(preconditions.values())

    status = (
        "F946_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
        if all_preconditions_pass
        else "F946_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "F946",
        "packet_name": "CurrentStrictT173T176InversionSensitivePair12BranchSeparationToChartSensitiveTransportedFluxCurrentLikeSectionBridgeOutputSchemaTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f945_bridge_target_packet": rel(IN_F945),
            "p946_supplier_vs_support_candidate_audit": rel(IN_P946_CANDIDATE),
            "p946_output_schema_class_support_audit": rel(IN_P946_CLASS),
            "p763_exact_missing_interface_boundary": rel(IN_P763),
            "p764_future_only_interface_target": rel(IN_P764),
            "t218_interface_target_spec": rel(IN_T218),
        },
        "preconditions": preconditions,
        "target_object": {
            "object_id": TARGET_OBJECT_ID,
            "goal": "Freeze the exact bridge_output_schema object still missing for the F945 bridge target after class-level support and one future-only support interface target are already present.",
            "required_fields": [
                {
                    "name": "inversion_sensitive_pair12_to_chart_sensitive_bridge_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F945 bridge target on full C_v1."
                },
                {
                    "name": "output_schema_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve that support exists only at candidate/class grade around the lane."
                },
                {
                    "name": "future_only_support_interface_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exported T218 future-only support target without promoting it to an actual bridge supplier."
                },
                {
                    "name": "exact_missing_support_interface_name",
                    "required": True,
                    "hard_limit": "Must preserve the exact P763/T218 interface name and keep it unexported."
                },
                {
                    "name": "support_candidate_scope_gap_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that the current support candidate remains route-local and below full-C_v1 T176 discharge."
                },
                {
                    "name": "exact_bridge_output_schema_statement",
                    "required": True,
                    "hard_limit": "Must state the exact missing bridge_output_schema object without claiming that it already exists."
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny T176 discharge, T177 discharge, T185 discharge, QW-2191 discharge, and ToE closure."
                }
            ]
        },
        "target_refs": {
            "inversion_sensitive_pair12_to_chart_sensitive_bridge_target_ref": (
                (f945.get("target_object") or {}).get("object_id")
            ),
            "output_schema_class_candidate_support_refs": p946_class_support_stack.get("candidate_support_refs"),
            "future_only_support_interface_target_ref": p946_candidate_conclusion.get(
                "current_best_future_only_interface_target_candidate"
            ),
            "exact_missing_support_interface_name": p946_candidate_conclusion.get(
                "current_best_support_interface_name"
            ),
            "support_candidate_scope_gap_ref": {
                "starts_at": "Sigma_sel_src_target_v1",
                "ends_at": "surviving_F301_pair12_residual_datum_carrier_lane",
                "future_route_only": p764.get("current_t218_target_is_future_route_only"),
                "below_actual_interface_export_and_below_t176_discharge": p764.get(
                    "current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge"
                ),
            },
            "exact_bridge_output_schema_statement": {
                "object_id": TARGET_OBJECT_ID,
                "for_bridge_target": (f945.get("target_object") or {}).get("object_id"),
                "required_minimum_properties": (f945.get("bridge_output_schema") or {}).get("minimum_properties"),
            },
            "supplier_gap_ref": {
                "current_repo_already_exports_actual_bridge_output_schema_supplier": p946_candidate_conclusion.get(
                    "current_repo_already_exports_actual_bridge_output_schema_supplier"
                ),
                "bridge_output_schema_exported_now": p946_class_theorem.get("bridge_output_schema_exported_now"),
            },
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class support around the active F945 bridge lane at candidate grade.",
            "The repo also already names one sharp future-only support interface target through P763/P764/T218.",
            "But the exact bridge_output_schema object for the full-C_v1 inversion-sensitive pair12 to chart-sensitive transported section bridge still remains unexported.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether the frozen T218 future-only support interface, even if exported, would already suffice for inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_v1 or would still require one additional full_C_v1 transported-section lift.",
        "hard_limits": [
            "Does not claim that the exact bridge_output_schema is already exported.",
            "Does not claim that the T218 future-only support interface target already suffices for full C_v1.",
            "Does not claim that T176 is discharged.",
            "Does not claim that T177 is discharged.",
            "Does not claim that T185 is discharged.",
            "Does not claim that F647 is admissible S_sel_int.",
            "Does not claim that the rooted w_break sign lift is a strict physical orientation datum.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F946",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
