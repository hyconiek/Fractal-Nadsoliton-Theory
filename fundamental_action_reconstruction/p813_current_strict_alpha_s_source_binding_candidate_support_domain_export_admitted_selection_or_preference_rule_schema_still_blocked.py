#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F812 = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json"
IN_F811 = GENERATED / "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.json"
IN_P811 = GENERATED / "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json"
IN_P812 = GENERATED / "p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required.json"
IN_P755 = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe.json"
IN_P756 = GENERATED / "p756_current_strict_t210_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_audit_probe.json"

OUT = GENERATED / "p813_current_strict_alpha_s_source_binding_candidate_support_domain_export_admitted_selection_or_preference_rule_schema_still_blocked.json"
OUT_SUMMARY = GENERATED / "p813_current_strict_alpha_s_source_binding_candidate_support_domain_export_admitted_selection_or_preference_rule_schema_still_blocked_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F812, IN_F811, IN_P811, IN_P812, IN_P755, IN_P756]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P813",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f812 = load_json(IN_F812)
    f811 = load_json(IN_F811)
    p811 = load_json(IN_P811)
    p812 = load_json(IN_P812)
    p755 = load_json(IN_P755)
    p756 = load_json(IN_P756)

    f812_target = f812.get("target_object") or {}
    f812_refs = f812.get("target_refs") or {}
    f811_refs = f811.get("target_refs") or {}
    p811_support = p811.get("support_objects") or {}
    p811_theorem = p811.get("theorem_result") or {}
    p812_theorem = p812.get("theorem_result") or {}
    p755_theorem = p755.get("theorem_result") or {}
    p756_theorem = p756.get("theorem_result") or {}

    domain_members = [
        {
            "member_id": "current_source_candidate_support",
            "object_ref": p811_support.get("current_source_candidate_support_ref"),
            "member_grade": "current_candidate_only",
            "support_status": "exported_current_support",
        },
        {
            "member_id": "lawful_future_entry_object_support",
            "object_ref": p811_support.get("lawful_future_entry_object_support_ref"),
            "member_grade": "future_only_entry_target",
            "support_status": "exported_future_only_support",
        },
    ]

    checks = [
        {
            "id": "f812_requires_candidate_domain_and_rule_schema",
            "pass": (
                f812.get("status")
                == "F812_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and f812_target.get("object_id")
                == "alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "candidate_source_support_domain_ref"
                    for item in (f812_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "selection_or_preference_rule_schema"
                    for item in (f812_target.get("required_fields") or [])
                )
            ),
            "details": "F812 already freezes the exact missing selection/preference-rule target and requires both candidate_source_support_domain_ref and selection_or_preference_rule_schema.",
        },
        {
            "id": "f811_and_p811_already_make_the_support_pair_explicit",
            "pass": (
                f811.get("status")
                == "F811_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_TARGET_PACKET_NO_FALSE_PASS"
                and p811.get("status")
                == "P811_CURRENT_STRICT_ALPHA_S_SOURCE_SUPPORT_PRESENT_EXACT_SOURCE_BINDING_UNEXPORTED_ADAPTER_ACTION_SCHEMA_BLOCKED_SOURCE_BINDING_TARGET_FREEZE_REQUIRED"
                and p811_theorem.get("current_source_candidate_support_present") is True
                and p811_theorem.get("lawful_future_entry_object_support_present") is True
                and f811_refs.get("current_source_candidate_support_ref")
                == p811_support.get("current_source_candidate_support_ref")
                and f811_refs.get("lawful_future_entry_object_support_ref")
                == p811_support.get("lawful_future_entry_object_support_ref")
            ),
            "details": "F811/P811 already expose the exact two support objects needed for a finite candidate support domain.",
        },
        {
            "id": "future_entry_support_remains_future_only_not_actual_entry",
            "pass": (
                p755.get("status") == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_EXPORTED"
                and p755_theorem.get("current_t209_target_is_future_only") is True
                and p756.get("status")
                == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
                and p756_theorem.get("current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object")
                is False
            ),
            "details": "The future entry support remains future-only and explicitly below actual realization.",
        },
        {
            "id": "finite_two_support_domain_can_be_exported_without_selecting_one_source",
            "pass": (
                len(domain_members) == 2
                and all(member["object_ref"] for member in domain_members)
                and p812_theorem.get("exact_source_binding_selection_or_preference_rule_exported_for_f811_target")
                is False
            ),
            "details": "The current support pair is explicit enough to be exported as one finite support domain without thereby selecting one source.",
        },
        {
            "id": "selection_or_preference_rule_schema_remains_blocked",
            "pass": (
                p812.get("status")
                == "P812_CURRENT_STRICT_ALPHA_S_NO_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_FOR_F811_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
                and p812_theorem.get("exact_source_binding_selection_or_preference_rule_exported_for_f811_target")
                is False
            ),
            "details": "The selection/preference rule schema itself remains blocked and unexported on current repo state.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    clause_split = {
        "candidate_source_support_domain_clause_status": "export_admitted" if all_pass else "review_required",
        "selection_or_preference_rule_schema_clause_status": "blocked_nonexport" if all_pass else "review_required",
        "sharp_blocker_clause": "selection_or_preference_rule_schema" if all_pass else "review_required",
    }

    theorem_result = {
        "candidate_source_support_domain_exportable_now": checks[0]["pass"] and checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"],
        "selection_or_preference_rule_schema_exported_now": False if all_pass else None,
        "next_honest_move_is_export_support_domain_leave_rule_schema_blocked": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P813_CURRENT_STRICT_ALPHA_S_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED"
        if all_pass
        else "P813_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P813",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f812_selection_or_preference_rule_target_packet": rel(IN_F812),
            "f811_source_binding_target_packet": rel(IN_F811),
            "p811_source_binding_probe": rel(IN_P811),
            "p812_selection_or_preference_rule_probe": rel(IN_P812),
            "p755_future_entry_target_probe": rel(IN_P755),
            "p756_future_entry_actual_realization_nonexport_audit": rel(IN_P756),
        },
        "support_domain_candidate": {
            "domain_id": "alpha_s_strict_source_shannon_source_binding_candidate_support_domain_v1",
            "member_count": len(domain_members),
            "members": domain_members,
        },
        "clause_split": clause_split,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact two support objects needed for source binding are already explicit on current repo state.",
            "That pair can be exported as a finite candidate support domain without yet choosing one source.",
            "What remains missing is only the selection/preference rule schema that would act on that domain.",
        ],
        "recommended_next_packet": {
            "id": "F813_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_PACKET",
            "goal": "Export the admitted finite candidate support domain explicitly while leaving the selection/preference rule schema unresolved.",
            "export_object_id": "alpha_s_strict_source_shannon_source_binding_candidate_support_domain_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P813",
        "status": status,
        "as_of": AS_OF,
        "candidate_source_support_domain_exportable_now": theorem_result["candidate_source_support_domain_exportable_now"],
        "selection_or_preference_rule_schema_exported_now": theorem_result["selection_or_preference_rule_schema_exported_now"],
        "next_honest_move_is_export_support_domain_leave_rule_schema_blocked": theorem_result[
            "next_honest_move_is_export_support_domain_leave_rule_schema_blocked"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
