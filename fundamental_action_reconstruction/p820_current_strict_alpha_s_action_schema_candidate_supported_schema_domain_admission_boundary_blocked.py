#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P819 = GENERATED / "p819_current_strict_alpha_s_no_shift_interface_adapter_or_carrier_identification_rule_for_pair12_source_side_branch_selection_provider_interface_target_exported_target_freeze_required.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F810 = GENERATED / "f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"

OUT = GENERATED / "p820_current_strict_alpha_s_action_schema_candidate_supported_schema_domain_admission_boundary_blocked.json"
OUT_SUMMARY = GENERATED / "p820_current_strict_alpha_s_action_schema_candidate_supported_schema_domain_admission_boundary_blocked_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P819, IN_F819, IN_F810, IN_P764, IN_P788]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P820",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p819 = load_json(IN_P819)
    f819 = load_json(IN_F819)
    f810 = load_json(IN_F810)
    p764 = load_json(IN_P764)
    p788 = load_json(IN_P788)

    f819_target = f819.get("target_object") or {}
    f810_target = f810.get("target_object") or {}
    p764_theorem = p764.get("theorem_result") or {}

    required_fields = [
        item.get("name")
        for item in (f819_target.get("required_fields") or [])
        if isinstance(item, dict)
    ]

    action_schema_checks = [
        {
            "id": "f819_already_names_action_schema_clause",
            "pass": "adapter_or_carrier_identification_action_schema" in required_fields,
            "details": "F819 already isolates adapter_or_carrier_identification_action_schema as an exact clause of the missing rule target.",
        },
        {
            "id": "p764_exports_exact_own_lane_typed_interface_target",
            "pass": (
                p764.get("status")
                == "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED"
                and p764_theorem.get("t218_target_exported_on_current_repo_state") is True
                and p764_theorem.get("current_t218_target_freezes_exact_t216_immediate_missing_interface")
                is True
                and p764_theorem.get("current_t218_target_is_future_route_only") is True
            ),
            "details": "P764 already exports one exact own-lane typed interface target, so the action-schema side has a real structural carrier candidate even though it is not an alpha_s export.",
        },
        {
            "id": "old_f810_preserves_action_schema_class_only_as_template",
            "pass": (
                f810.get("status")
                == "F810_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and any(
                    isinstance(item, dict) and item.get("name") == "adapter_action_schema"
                    for item in (f810_target.get("required_fields") or [])
                )
            ),
            "details": "Old F810 preserves the adapter-action-schema class as a lawful template category, but only on the old Shannon lane.",
        },
    ]

    schema_domain_boundary_checks = [
        {
            "id": "f819_already_names_schema_domain_boundary_clause",
            "pass": "schema_domain_admission_or_nonidentification_boundary_ref" in required_fields,
            "details": "F819 already isolates the schema-domain admission/nonidentification boundary as an exact clause of the missing rule target.",
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
            "id": "p819_still_reports_no_exact_rule_for_new_lane",
            "pass": (
                p819.get("status")
                == "P819_CURRENT_STRICT_ALPHA_S_NO_SHIFT_INTERFACE_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_FOR_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
                and (p819.get("theorem_result") or {}).get(
                    "exact_shift_interface_adapter_or_carrier_identification_rule_exported_for_f818_target"
                )
                is False
                and ((p819.get("exact_missing_rule_target_candidate") or {}).get("repo_scan_hits_for_exact_rule") or [])
                == []
            ),
            "details": "P819 still shows that no exact rule is exported for the new T213/T216 lane into the alpha_s schema problem.",
        },
        {
            "id": "old_f810_is_not_reusable_as_schema_domain_boundary_for_new_lane",
            "pass": (
                ((f810.get("target_refs") or {}).get("provider_shift_candidate_reference_lane_ref"))
                == "alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_v1"
            ),
            "details": "Old F810 remains lane-specific to strict-source Shannon, so it cannot discharge schema-domain admission for the new lane.",
        },
    ]

    failed_action = [item["id"] for item in action_schema_checks if not item["pass"]]
    failed_boundary = [item["id"] for item in schema_domain_boundary_checks if not item["pass"]]

    action_schema_clause_status = (
        "candidate_supported_not_yet_exported" if not failed_action else "requires_review"
    )
    schema_domain_boundary_clause_status = (
        "blocked_nonexport" if not failed_boundary else "requires_review"
    )

    status = (
        "P820_CURRENT_STRICT_ALPHA_S_ACTION_SCHEMA_CANDIDATE_SUPPORTED_SCHEMA_DOMAIN_ADMISSION_BOUNDARY_BLOCKED"
        if not failed_action and not failed_boundary
        else "P820_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P820",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p819_rule_absence_audit": rel(IN_P819),
            "f819_rule_target_packet": rel(IN_F819),
            "f810_old_shannon_rule_target_packet": rel(IN_F810),
            "p764_own_lane_interface_target_probe": rel(IN_P764),
            "p788_alpha_s_adapter_probe": rel(IN_P788),
        },
        "clause_split_audit": {
            "action_schema_clause_status": action_schema_clause_status,
            "schema_domain_admission_boundary_clause_status": schema_domain_boundary_clause_status,
            "sharp_blocker_clause": "schema_domain_admission_boundary"
            if action_schema_clause_status == "candidate_supported_not_yet_exported"
            and schema_domain_boundary_clause_status == "blocked_nonexport"
            else None,
        },
        "action_schema_checks": action_schema_checks,
        "schema_domain_boundary_checks": schema_domain_boundary_checks,
        "failed_action_schema_checks": failed_action,
        "failed_schema_domain_boundary_checks": failed_boundary,
        "current_honest_reading": [
            "The action-schema side now has candidate support: the repo already exports one exact own-lane typed interface target on T213/T216, and the old Shannon lane preserves adapter-action-schema as a lawful template class.",
            "But lawful schema-domain admission into the current alpha_s schema problem remains blocked: no generic alpha_s adapter is exported, no exact rule exists for the new lane, and old F810 is not silently reusable.",
            "So the sharp blocker is now the missing schema-domain admission / nonidentification boundary, not generic action-schema shape.",
        ],
        "recommended_next_packet": {
            "id": "F820_CURRENT_STRICT_ALPHA_S_SCHEMA_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET",
            "goal": "Freeze the exact schema-domain admission or nonidentification boundary object still missing before the new lane can move closer to an export-grade rule.",
            "minimum_fields": [
                "shift_to_schema_adapter_or_carrier_identification_rule_target_ref",
                "different_selection_provider_class_shift_candidate_reference_lane_ref",
                "candidate_action_schema_support_ref",
                "generic_alpha_s_adapter_absence_ref",
                "old_shannon_nontransfer_ref",
                "boundary_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P820",
        "status": status,
        "as_of": AS_OF,
        "action_schema_clause_status": action_schema_clause_status,
        "schema_domain_admission_boundary_clause_status": schema_domain_boundary_clause_status,
        "sharp_blocker_clause": artifact["clause_split_audit"]["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
