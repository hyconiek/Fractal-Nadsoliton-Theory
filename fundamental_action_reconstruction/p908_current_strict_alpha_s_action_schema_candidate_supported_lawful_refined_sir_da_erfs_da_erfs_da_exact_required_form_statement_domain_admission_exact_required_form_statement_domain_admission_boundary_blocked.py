#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P907 = GENERATED / "p907_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_rule_for_pair12_provider_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_interface_target_freeze_required.json"
IN_F907 = GENERATED / "f907_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_adapter_or_carrier_rule_target_packet.json"
IN_F896 = GENERATED / "f896_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"

OUT = GENERATED / "p908_current_strict_alpha_s_action_schema_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_boundary_blocked.json"
OUT_SUMMARY = GENERATED / "p908_current_strict_alpha_s_action_schema_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_boundary_blocked_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P907, IN_F907, IN_F896, IN_P764, IN_P788]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P908",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p907 = load_json(IN_P907)
    f907 = load_json(IN_F907)
    f896 = load_json(IN_F896)
    p764 = load_json(IN_P764)
    p788 = load_json(IN_P788)

    p907_theorem = p907.get("theorem_result") or {}
    f907_target = f907.get("target_object") or {}
    f907_refs = f907.get("target_refs") or {}
    f896_target = f896.get("target_object") or {}
    f896_refs = f896.get("target_refs") or {}
    p764_theorem = p764.get("theorem_result") or {}

    required_fields = [
        item.get("name")
        for item in (f907_target.get("required_fields") or [])
        if isinstance(item, dict)
    ]

    action_schema_checks = [
        {
            "id": "f907_already_names_action_schema_clause",
            "pass": "adapter_or_carrier_identification_action_schema" in required_fields,
            "details": "F907 already isolates adapter_or_carrier_identification_action_schema as an exact clause of the missing deeper lawful refined rule target.",
        },
        {
            "id": "p764_exports_exact_own_lane_typed_interface_target",
            "pass": (
                p764.get("status")
                == "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED"
                and p764_theorem.get("t218_target_exported_on_current_repo_state") is True
                and p764_theorem.get("current_t218_target_freezes_exact_t216_immediate_missing_interface") is True
                and p764_theorem.get("current_t218_target_is_future_route_only") is True
            ),
            "details": "P764 still exports one exact own-lane typed interface target, so the action-schema side keeps a real structural carrier candidate even though it is not an alpha_s export.",
        },
        {
            "id": "old_f896_preserves_action_schema_class_on_same_provider_lane",
            "pass": (
                f896.get("status")
                == "F896_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and any(
                    isinstance(item, dict) and item.get("name") == "adapter_or_carrier_identification_action_schema"
                    for item in (f896_target.get("required_fields") or [])
                )
            ),
            "details": "Old F896 preserves the same action-schema class on the same provider lane, but only for the shallower lawful refined rule target.",
        },
    ]

    boundary_checks = [
        {
            "id": "f907_already_names_lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_clause",
            "pass": "lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_ref" in required_fields,
            "details": "F907 already isolates the lawful refined deeper domain-admission / nonidentification boundary as an exact clause of the missing rule target.",
        },
        {
            "id": "p788_still_reports_no_generic_alpha_s_adapter",
            "pass": (
                p788.get("status")
                == "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE"
                and any(
                    isinstance(item, dict)
                    and item.get("id") == "no_exported_dimensionless_alpha_s_adapter_detected"
                    and item.get("pass") is True
                    for item in (p788.get("checks") or [])
                )
            ),
            "details": "P788 still blocks generic alpha_s adapter/export support on the current repo state.",
        },
        {
            "id": "p907_still_reports_no_exact_rule_for_new_deeper_domain_admission_lane",
            "pass": (
                p907.get("status")
                == "P907_CURRENT_STRICT_ALPHA_S_NO_SHIFT_INTERFACE_ADAPTER_OR_CARRIER_RULE_FOR_PAIR12_PROVIDER_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_INTERFACE_TARGET_FREEZE_REQUIRED"
                and p907_theorem.get("exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f906_target")
                is False
                and ((p907.get("exact_missing_rule_target_candidate") or {}).get("repo_scan_hits_for_exact_rule") or [])
                == []
            ),
            "details": "P907 still shows that no exact rule is exported for the deeper lawful refined domain-admission lane into the current alpha_s problem.",
        },
        {
            "id": "old_f896_is_not_reusable_as_lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_for_new_lane",
            "pass": (
                f896_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_v1"
                and f896_refs.get("shift_to_lawful_refined_deeper_exact_required_form_statement_interface_target_ref")
                != f907_refs.get("shift_to_lawful_refined_deeper_domain_admission_exact_required_form_statement_interface_target_ref")
                and f896_refs.get("downstream_lawful_refined_deeper_exact_required_form_statement_target_ref")
                != f907_refs.get("downstream_lawful_refined_deeper_domain_admission_exact_required_form_statement_target_ref")
            ),
            "details": "Old F896 remains shallower-lane specific, so it cannot silently discharge the deeper lawful refined domain-admission boundary for the new lane.",
        },
    ]

    failed_action = [item["id"] for item in action_schema_checks if not item["pass"]]
    failed_boundary = [item["id"] for item in boundary_checks if not item["pass"]]

    action_schema_clause_status = (
        "candidate_supported_not_yet_exported" if not failed_action else "requires_review"
    )
    lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_clause_status = (
        "blocked_nonexport" if not failed_boundary else "requires_review"
    )

    status = (
        "P908_CURRENT_STRICT_ALPHA_S_ACTION_SCHEMA_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_BOUNDARY_BLOCKED"
        if not failed_action and not failed_boundary
        else "P908_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P908",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p907_rule_absence_audit": rel(IN_P907),
            "f907_rule_target_packet": rel(IN_F907),
            "f896_old_same_provider_rule_target_packet": rel(IN_F896),
            "p764_own_lane_interface_target_probe": rel(IN_P764),
            "p788_alpha_s_adapter_probe": rel(IN_P788),
        },
        "clause_split_audit": {
            "action_schema_clause_status": action_schema_clause_status,
            "lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_clause_status": lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_clause_status,
            "sharp_blocker_clause": "lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary"
            if action_schema_clause_status == "candidate_supported_not_yet_exported"
            and lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_clause_status == "blocked_nonexport"
            else None,
        },
        "action_schema_checks": action_schema_checks,
        "lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_checks": boundary_checks,
        "failed_action_schema_checks": failed_action,
        "failed_lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_checks": failed_boundary,
        "current_honest_reading": [
            "The action-schema side now has only candidate support: the repo still exports one exact own-lane typed interface target on T213/T216, and old F896 preserves the same action-schema class on the same provider lane.",
            "But lawful refined deeper domain-admission exact-required-form-statement domain admission into the current alpha_s problem remains blocked: no generic alpha_s adapter is exported, no exact rule exists for the deeper new lane, and old F896 is not silently reusable as the needed lawful boundary.",
            "So the sharp blocker is now the missing lawful refined deeper domain-admission exact-required-form-statement domain-admission / nonidentification boundary for the frozen F907 rule target, not generic action-schema shape.",
        ],
        "recommended_next_packet": {
            "id": "F908_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined deeper domain-admission or nonidentification boundary object still missing before the new lane can move closer to an export-grade rule.",
            "minimum_fields": [
                "shift_to_lawful_refined_deeper_domain_admission_exact_required_form_statement_adapter_or_carrier_identification_rule_target_ref",
                "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref",
                "candidate_action_schema_support_ref",
                "downstream_lawful_refined_deeper_domain_admission_exact_required_form_statement_target_ref",
                "generic_alpha_s_adapter_absence_ref",
                "old_same_provider_lane_nontransfer_ref",
                "boundary_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P908",
        "status": status,
        "as_of": AS_OF,
        "action_schema_clause_status": action_schema_clause_status,
        "lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_clause_status": lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_boundary_clause_status,
        "sharp_blocker_clause": artifact["clause_split_audit"]["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
