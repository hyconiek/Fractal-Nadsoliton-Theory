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
IN_F815 = GENERATED / "f815_current_strict_alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_packet.json"
IN_P814 = GENERATED / "p814_current_strict_alpha_s_no_selection_or_preference_rule_schema_for_f812_target_on_exported_f813_support_domain_narrow_schema_target_freeze_required.json"
IN_P815 = GENERATED / "p815_current_strict_alpha_s_domain_member_grade_handling_clause_export_admitted_selection_or_preference_rule_schema_still_blocked.json"
IN_P791 = GENERATED / "p791_current_strict_alpha_s_selection_principle_reuse_audit_probe.json"
IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"

OUT = GENERATED / "p816_current_strict_alpha_s_selection_or_preference_rule_schema_same_lane_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p816_current_strict_alpha_s_selection_or_preference_rule_schema_same_lane_exhaustion_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F814, IN_F815, IN_P814, IN_P815, IN_P791, IN_P792]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P816",
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
    f815 = load_json(IN_F815)
    p814 = load_json(IN_P814)
    p815 = load_json(IN_P815)
    p791 = load_json(IN_P791)
    p792 = load_json(IN_P792)

    f814_target = f814.get("target_object") or {}
    f815_export = f815.get("exported_object") or {}
    p814_theorem = p814.get("theorem_result") or {}
    p815_theorem = p815.get("theorem_result") or {}
    p792_rule = p792.get("probe_local_order_rule_tuple") or {}

    checks = [
        {
            "id": "f814_schema_only_target_already_frozen",
            "pass": (
                f814.get("status")
                == "F814_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f814_target.get("object_id")
                == "alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_v1"
            ),
            "details": "F814 already freezes the exact schema-only target on the F813 source-binding domain.",
        },
        {
            "id": "f815_grade_clause_already_exported",
            "pass": (
                f815.get("status")
                == "F815_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_PACKET_NO_FALSE_PASS"
                and f815_export.get("object_id")
                == "alpha_s_strict_source_shannon_source_binding_domain_member_grade_handling_clause_v1"
            ),
            "details": "F815 already exports the exact grade-handling clause below the missing schema blocker.",
        },
        {
            "id": "no_current_export_supplies_selection_or_preference_rule_schema",
            "pass": (
                p814.get("status")
                == "P814_CURRENT_STRICT_ALPHA_S_NO_SELECTION_OR_PREFERENCE_RULE_SCHEMA_FOR_F812_TARGET_ON_EXPORTED_F813_SUPPORT_DOMAIN_NARROW_SCHEMA_TARGET_FREEZE_REQUIRED"
                and p814_theorem.get("exact_selection_or_preference_rule_schema_exported_on_f813_domain") is False
                and p815.get("status")
                == "P815_CURRENT_STRICT_ALPHA_S_DOMAIN_MEMBER_GRADE_HANDLING_CLAUSE_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED"
                and p815_theorem.get("selection_or_preference_rule_schema_exported_now") is False
            ),
            "details": "No current export supplies the missing selection_or_preference_rule_schema on the already grade-disciplined F813 domain.",
        },
        {
            "id": "foreign_domain_selection_template_reuse_still_blocked",
            "pass": (
                p791.get("status")
                == "P791_CURRENT_STRICT_ALPHA_S_SELECTION_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
                and p791.get("reusable_selection_principle_found") is False
            ),
            "details": "The real strict selection template family already known in the repo remains foreign-domain only and still cannot be silently reused here.",
        },
        {
            "id": "probe_local_family_order_rule_reuse_still_blocked",
            "pass": (
                p792.get("status")
                == "P792_CURRENT_STRICT_ALPHA_S_PROBE_LOCAL_ORDER_RULE_UNIQUE_WINNER_NONEXPORT"
                and p792_rule.get("export_status") == "probe_local_only"
                and p814_theorem.get("probe_local_family_order_rule_reusable_for_source_binding_schema") is False
            ),
            "details": "The only close alpha_s-side order rule remains probe-local on the family domain and is still not reusable as the missing source-binding schema.",
        },
        {
            "id": "no_residual_passive_same_lane_loophole_below_schema_blocker",
            "pass": (
                f815_export.get("scope") == "strict_source_side_domain_grade_discipline_only"
                and p815_theorem.get("domain_member_grade_handling_clause_exportable_now") is True
                and p815_theorem.get("next_honest_move_is_export_grade_clause_leave_selection_schema_blocked") is True
            ),
            "details": "After domain export and grade-clause export, no residual passive same-lane loophole remains below the current schema blocker on this lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "selection_or_preference_rule_schema_exported_on_current_repo_state": False if all_pass else None,
        "same_lane_passive_groundwork_exhausted": True if all_pass else None,
        "next_honest_move_requires_continuation_boundary_export": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P816_CURRENT_STRICT_ALPHA_S_SELECTION_OR_PREFERENCE_RULE_SCHEMA_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        if all_pass
        else "P816_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P816",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f814_schema_only_target_packet": rel(IN_F814),
            "f815_grade_clause_packet": rel(IN_F815),
            "p814_schema_absence_probe": rel(IN_P814),
            "p815_grade_clause_probe": rel(IN_P815),
            "p791_foreign_template_reuse_audit": rel(IN_P791),
            "p792_probe_local_order_rule_probe": rel(IN_P792),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact schema-only target is already frozen and the exact grade clause is already exported.",
            "No current export supplies the missing selection/preference rule schema, and the obvious foreign-domain and probe-local reuse routes remain blocked.",
            "So no residual passive same-lane loophole remains below the current schema blocker on this lane.",
        ],
        "recommended_next_packet": {
            "id": "F816_CURRENT_STRICT_ALPHA_S_SELECTION_OR_PREFERENCE_RULE_SCHEMA_CONTINUATION_BOUNDARY_PACKET",
            "goal": "Export the continuation boundary after the current same-lane passive groundwork is exhausted but the missing schema still remains.",
            "export_object_id": "alpha_s_selection_or_preference_rule_schema_continuation_boundary_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P816",
        "status": status,
        "as_of": AS_OF,
        "selection_or_preference_rule_schema_exported_on_current_repo_state": theorem_result[
            "selection_or_preference_rule_schema_exported_on_current_repo_state"
        ],
        "same_lane_passive_groundwork_exhausted": theorem_result["same_lane_passive_groundwork_exhausted"],
        "next_honest_move_requires_continuation_boundary_export": theorem_result[
            "next_honest_move_requires_continuation_boundary_export"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
