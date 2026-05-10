#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
IN_F813 = GENERATED / "f813_current_strict_alpha_s_strict_source_shannon_source_binding_candidate_support_domain_packet.json"
IN_F811 = GENERATED / "f811_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_target_packet.json"
IN_P811 = GENERATED / "p811_current_strict_alpha_s_source_support_present_exact_source_binding_unexported_adapter_action_schema_blocked_source_binding_target_freeze_required.json"
IN_P755 = GENERATED / "p755_current_strict_t209_t26_component2_minimal_designated_pair12_noncyclic_entry_object_target_probe.json"
IN_P756 = GENERATED / "p756_current_strict_t210_t26_component2_minimal_designated_pair12_noncyclic_entry_object_actual_realization_nonexport_audit_probe.json"

OUT = GENERATED / "p815_current_strict_alpha_s_domain_member_grade_handling_clause_export_admitted_selection_or_preference_rule_schema_still_blocked.json"
OUT_SUMMARY = GENERATED / "p815_current_strict_alpha_s_domain_member_grade_handling_clause_export_admitted_selection_or_preference_rule_schema_still_blocked_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F814, IN_F813, IN_F811, IN_P811, IN_P755, IN_P756]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P815",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f814 = load_json(IN_F814)
    f813 = load_json(IN_F813)
    f811 = load_json(IN_F811)
    p811 = load_json(IN_P811)
    p755 = load_json(IN_P755)
    p756 = load_json(IN_P756)

    f814_target = f814.get("target_object") or {}
    f813_export = f813.get("exported_object") or {}
    f811_refs = f811.get("target_refs") or {}
    p811_support = p811.get("support_objects") or {}
    p755_theorem = p755.get("theorem_result") or {}
    p756_theorem = p756.get("theorem_result") or {}

    clause_fields = {
        "candidate_support_member_id": "current_source_candidate_support",
        "candidate_support_ref": f811_refs.get("current_source_candidate_support_ref"),
        "candidate_support_grade": "current_candidate_only",
        "future_entry_support_member_id": "lawful_future_entry_object_support",
        "future_entry_support_ref": f811_refs.get("lawful_future_entry_object_support_ref"),
        "future_entry_support_grade": "future_only_entry_target",
        "no_silent_grade_promotion": True,
        "candidate_support_stays_below_actual_source_binding": True,
        "future_entry_support_stays_below_actual_entry_and_actual_source_binding": True,
    }

    checks = [
        {
            "id": "f814_requires_grade_clause_and_selection_schema",
            "pass": (
                f814.get("status")
                == "F814_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f814_target.get("object_id")
                == "alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "domain_member_grade_handling_clause"
                    for item in (f814_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "selection_or_preference_rule_schema"
                    for item in (f814_target.get("required_fields") or [])
                )
            ),
            "details": "F814 already freezes the exact schema-only target and requires both domain_member_grade_handling_clause and selection_or_preference_rule_schema.",
        },
        {
            "id": "f813_exports_exact_domain_with_explicit_member_grades",
            "pass": (
                f813.get("status")
                == "F813_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_PACKET_NO_FALSE_PASS"
                and f813_export.get("object_id")
                == "alpha_s_strict_source_shannon_source_binding_candidate_support_domain_v1"
                and f813_export.get("member_count") == 2
                and (f813_export.get("members") or [])[0].get("member_grade") == "current_candidate_only"
                and (f813_export.get("members") or [])[1].get("member_grade") == "future_only_entry_target"
            ),
            "details": "F813 already exports the exact finite support domain and keeps both member grades explicit.",
        },
        {
            "id": "f321_side_support_still_candidate_only",
            "pass": (
                p811.get("status")
                == "P811_CURRENT_STRICT_ALPHA_S_SOURCE_SUPPORT_PRESENT_EXACT_SOURCE_BINDING_UNEXPORTED_ADAPTER_ACTION_SCHEMA_BLOCKED_SOURCE_BINDING_TARGET_FREEZE_REQUIRED"
                and p811_support.get("current_source_candidate_support_ref")
                == "BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_population_refinement_candidate_v1"
                and f811_refs.get("current_source_candidate_support_ref")
                == p811_support.get("current_source_candidate_support_ref")
            ),
            "details": "The current strict-source Shannon support remains the F321 candidate-only object and is not promoted into actual source binding.",
        },
        {
            "id": "future_entry_support_still_future_only_and_below_actual_realization",
            "pass": (
                p755.get("status") == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_TARGET_EXPORTED"
                and p755_theorem.get("current_t209_target_is_future_only") is True
                and p756.get("status")
                == "PASS_STRICT_T26_COMPONENT2_MINIMAL_DESIGNATED_PAIR12_NONCYCLIC_ENTRY_OBJECT_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
                and p756_theorem.get("current_repo_still_does_not_export_actual_realization_of_t209_target") is True
                and p756_theorem.get("current_repo_has_exported_actual_t26_component2_minimal_designated_pair12_noncyclic_entry_object")
                is False
            ),
            "details": "The lawful future entry support remains future-only and below actual realization on current repo state.",
        },
        {
            "id": "grade_handling_clause_can_be_exported_without_selecting_one_source",
            "pass": (
                clause_fields["candidate_support_ref"]
                == "BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_population_refinement_candidate_v1"
                and clause_fields["future_entry_support_ref"]
                == "StrictT26Component2MinimalDesignatedPair12NoncyclicEntryObjectTarget_strict_v1"
                and clause_fields["no_silent_grade_promotion"] is True
            ),
            "details": "The current exports are explicit enough to support one exact grade-handling clause without thereby selecting one source or exporting a source-binding relation.",
        },
        {
            "id": "selection_or_preference_rule_schema_remains_blocked",
            "pass": (
                f813_export.get("unresolved_selection_or_preference_rule_schema_ref")
                == "alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1"
            ),
            "details": "The selection/preference rule schema itself remains blocked and unexported on current repo state.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "domain_member_grade_handling_clause_exportable_now": (
            checks[0]["pass"] and checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"] and checks[4]["pass"]
        ),
        "selection_or_preference_rule_schema_exported_now": False if all_pass else None,
        "next_honest_move_is_export_grade_clause_leave_selection_schema_blocked": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P815_CURRENT_STRICT_ALPHA_S_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED"
        if all_pass
        else "P815_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P815",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f814_schema_only_target_packet": rel(IN_F814),
            "f813_candidate_support_domain_packet": rel(IN_F813),
            "f811_source_binding_target_packet": rel(IN_F811),
            "p811_source_support_probe": rel(IN_P811),
            "p755_future_entry_target_probe": rel(IN_P755),
            "p756_future_entry_actual_realization_nonexport_audit": rel(IN_P756),
        },
        "grade_clause_candidate": {
            "clause_id": "alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_v1",
            "candidate_source_support_domain_ref": f813_export.get("object_id"),
            "fields": clause_fields,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact support domain already keeps both member grades explicit.",
            "Those grade facts are now strong enough to export one exact no-grade-promotion clause on the F813 domain.",
            "What remains missing is still only the selection/preference rule schema acting on that already grade-disciplined domain.",
        ],
        "recommended_next_packet": {
            "id": "F815_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_PACKET",
            "goal": "Export the exact grade-handling clause while leaving the source-binding selection/preference schema unresolved.",
            "export_object_id": "alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P815",
        "status": status,
        "as_of": AS_OF,
        "domain_member_grade_handling_clause_exportable_now": theorem_result[
            "domain_member_grade_handling_clause_exportable_now"
        ],
        "selection_or_preference_rule_schema_exported_now": theorem_result[
            "selection_or_preference_rule_schema_exported_now"
        ],
        "next_honest_move_is_export_grade_clause_leave_selection_schema_blocked": theorem_result[
            "next_honest_move_is_export_grade_clause_leave_selection_schema_blocked"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
