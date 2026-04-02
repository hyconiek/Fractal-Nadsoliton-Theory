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
IN_F813 = GENERATED / "f813_current_strict_alpha_s_strict_source_shannon_source_binding_candidate_support_domain_packet.json"
IN_P812 = GENERATED / "p812_current_strict_alpha_s_no_source_binding_selection_or_preference_rule_for_f811_target_exported_target_freeze_required.json"
IN_P813 = GENERATED / "p813_current_strict_alpha_s_source_binding_candidate_support_domain_export_admitted_selection_or_preference_rule_schema_still_blocked.json"
IN_P791 = GENERATED / "p791_current_strict_alpha_s_selection_principle_reuse_audit_probe.json"
IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"

OUT = GENERATED / "p814_current_strict_alpha_s_no_selection_or_preference_rule_schema_for_f812_target_on_exported_f813_support_domain_narrow_schema_target_freeze_required.json"
OUT_SUMMARY = GENERATED / "p814_current_strict_alpha_s_no_selection_or_preference_rule_schema_for_f812_target_on_exported_f813_support_domain_narrow_schema_target_freeze_required_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F812, IN_F813, IN_P812, IN_P813, IN_P791, IN_P792]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P814",
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
    f813 = load_json(IN_F813)
    p812 = load_json(IN_P812)
    p813 = load_json(IN_P813)
    p791 = load_json(IN_P791)
    p792 = load_json(IN_P792)

    f812_target = f812.get("target_object") or {}
    f813_export = f813.get("exported_object") or {}
    p812_theorem = p812.get("theorem_result") or {}
    p813_theorem = p813.get("theorem_result") or {}
    p792_rule = p792.get("probe_local_order_rule_tuple") or {}

    checks = [
        {
            "id": "f812_target_still_requires_selection_or_preference_rule_schema",
            "pass": (
                f812.get("status")
                == "F812_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and f812_target.get("object_id")
                == "alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "selection_or_preference_rule_schema"
                    for item in (f812_target.get("required_fields") or [])
                )
            ),
            "details": "F812 still freezes the exact source-binding selection/preference-rule target and still requires selection_or_preference_rule_schema.",
        },
        {
            "id": "f813_now_exports_exact_candidate_support_domain_for_f812_route",
            "pass": (
                f813.get("status")
                == "F813_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_PACKET_NO_FALSE_PASS"
                and f813_export.get("object_id")
                == "alpha_s_strict_source_shannon_source_binding_candidate_support_domain_v1"
                and f813_export.get("source_binding_selection_or_preference_rule_target_ref")
                == f812_target.get("object_id")
                and f813_export.get("support_domain_kind") == "finite_two_support_domain"
            ),
            "details": "F813 now exports the exact finite candidate support domain on which the frozen F812 rule target would act.",
        },
        {
            "id": "no_current_export_fills_selection_or_preference_rule_schema_on_exported_domain",
            "pass": (
                p812.get("status")
                == "P812_CURRENT_STRICT_ALPHA_S_NO_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_FOR_F811_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
                and p812_theorem.get("exact_source_binding_selection_or_preference_rule_exported_for_f811_target")
                is False
                and p813.get("status")
                == "P813_CURRENT_STRICT_ALPHA_S_SOURCE_BINDING_CANDIDATE_SUPPORT_DOMAIN_EXPORT_ADMITTED_SELECTION_OR_PREFERENCE_RULE_SCHEMA_STILL_BLOCKED"
                and p813_theorem.get("candidate_source_support_domain_exportable_now") is True
                and p813_theorem.get("selection_or_preference_rule_schema_exported_now") is False
                and f813_export.get("unresolved_selection_or_preference_rule_schema_ref")
                == f812_target.get("object_id")
            ),
            "details": "The support-domain slot is now filled by F813, but the selection/preference rule schema slot remains unfilled on that exact domain.",
        },
        {
            "id": "foreign_domain_selection_template_reuse_still_blocked",
            "pass": (
                p791.get("status")
                == "P791_CURRENT_STRICT_ALPHA_S_SELECTION_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
                and p791.get("reusable_selection_principle_found") is False
            ),
            "details": "The real strict selection template family already known in the repo remains foreign-domain only and cannot be silently reused here.",
        },
        {
            "id": "probe_local_family_order_rule_still_domain_mismatched_and_nonexport",
            "pass": (
                p792.get("status")
                == "P792_CURRENT_STRICT_ALPHA_S_PROBE_LOCAL_ORDER_RULE_UNIQUE_WINNER_NONEXPORT"
                and p792_rule.get("export_status") == "probe_local_only"
                and p792.get("unique_winner_exists") is True
            ),
            "details": "The only close alpha_s-side order rule remains probe-local on the family domain and is not an exported source-binding schema for the F813 support domain.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "candidate_source_support_domain_now_exactly_exported": checks[1]["pass"],
        "exact_selection_or_preference_rule_schema_exported_on_f813_domain": False if all_pass else None,
        "foreign_domain_selection_template_reusable": False if all_pass else None,
        "probe_local_family_order_rule_reusable_for_source_binding_schema": False if all_pass else None,
        "next_honest_move_requires_freezing_exact_schema_only_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P814_CURRENT_STRICT_ALPHA_S_NO_SELECTION_OR_PREFERENCE_RULE_SCHEMA_FOR_F812_TARGET_ON_EXPORTED_F813_SUPPORT_DOMAIN_NARROW_SCHEMA_TARGET_FREEZE_REQUIRED"
        if all_pass
        else "P814_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P814",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f812_selection_or_preference_rule_target_packet": rel(IN_F812),
            "f813_candidate_support_domain_packet": rel(IN_F813),
            "p812_rule_absence_probe": rel(IN_P812),
            "p813_support_domain_probe": rel(IN_P813),
            "p791_selection_template_reuse_audit": rel(IN_P791),
            "p792_family_order_rule_probe": rel(IN_P792),
        },
        "exact_missing_schema_target_candidate": {
            "candidate_id": "alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_missing_v1",
            "source_binding_selection_or_preference_rule_target_ref": f812_target.get("object_id"),
            "candidate_source_support_domain_ref": f813_export.get("object_id"),
            "repo_state_schema_status": "absent_on_exported_support_domain" if all_pass else "review_required",
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact candidate support domain is now exported by F813 for the frozen F812 route.",
            "No current export fills the remaining selection_or_preference_rule_schema slot on that exact domain.",
            "So the blocker is now schema-only: the next honest move is to freeze the exact missing schema target, not to pretend that an existing selection template can be reused here.",
        ],
        "recommended_next_packet": {
            "id": "F814_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact schema-only target that remains missing now that the F813 support domain is explicit.",
            "export_object_id": "alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P814",
        "status": status,
        "as_of": AS_OF,
        "candidate_source_support_domain_now_exactly_exported": theorem_result[
            "candidate_source_support_domain_now_exactly_exported"
        ],
        "exact_selection_or_preference_rule_schema_exported_on_f813_domain": theorem_result[
            "exact_selection_or_preference_rule_schema_exported_on_f813_domain"
        ],
        "next_honest_move_requires_freezing_exact_schema_only_target": theorem_result[
            "next_honest_move_requires_freezing_exact_schema_only_target"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
