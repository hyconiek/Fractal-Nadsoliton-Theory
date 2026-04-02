#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P845 = GENERATED / "p845_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_F845 = GENERATED / "f845_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F844 = GENERATED / "f844_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F843 = GENERATED / "f843_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_target_packet.json"
IN_F842 = GENERATED / "f842_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F841 = GENERATED / "f841_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F840 = GENERATED / "f840_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F836 = GENERATED / "f836_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p846_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
OUT_SUMMARY = GENERATED / "p846_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def has_required_field(obj: dict[str, Any], name: str) -> bool:
    return any(
        isinstance(item, dict) and item.get("name") == name
        for item in (obj.get("required_fields") or [])
    )


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P845, IN_F845, IN_F844, IN_F843, IN_F842, IN_F841, IN_F840, IN_F836]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P846",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p845 = load_json(IN_P845)
    f845 = load_json(IN_F845)
    f844 = load_json(IN_F844)
    f843 = load_json(IN_F843)
    f842 = load_json(IN_F842)
    f841 = load_json(IN_F841)
    f840 = load_json(IN_F840)
    f836 = load_json(IN_F836)

    p845_theorem = p845.get("theorem_result") or {}
    f845_target = f845.get("target_object") or {}
    f845_refs = f845.get("target_refs") or {}
    f844_target = f844.get("target_object") or {}
    f843_target = f843.get("target_object") or {}
    f842_target = f842.get("target_object") or {}
    f841_target = f841.get("target_object") or {}
    f840_target = f840.get("target_object") or {}
    f836_target = f836.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_exact_required_form_statement_domain_admission_refined_exact_statement_required_form"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p845_", "f845_", "p846_", "f846_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f845_already_names_exact_missing_form_field",
            "pass": (
                f845.get("status")
                == "F845_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f845_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_output_schema_statement_target_v1"
                and has_required_field(f845_target, "exact_statement_required_form_ref")
            ),
            "details": "F845 already isolates exact_statement_required_form_ref as one exact missing field of the refined lawful exact-statement target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_form_class_support",
            "pass": (
                has_required_field(f844_target, "exact_output_schema_statement")
                and has_required_field(f843_target, "lawful_exact_required_form_statement_domain_admission_output_schema")
                and has_required_field(f842_target, "boundary_output_schema")
                and has_required_field(f841_target, "selected_interface_output_schema")
                and has_required_field(f840_target, "exact_interface_output_schema")
                and f836_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve refined output or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact required-form export for the refined new lane.",
        },
        {
            "id": "p845_already_keeps_refined_exact_statement_unexported",
            "pass": (
                p845.get("status")
                == "P845_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p845_theorem.get("exact_output_schema_statement_exported_now") is False
            ),
            "details": "P845 already keeps the exact refined new-lane statement object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_statement_required_form_ref outside the new frozen F845 lineage itself.",
        },
        {
            "id": "neighboring_form_support_remains_nonidentical_to_refined_new_lane_form",
            "pass": (
                f845_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_output_schema_statement_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_target_v1"
                in (f845_refs.get("statement_class_candidate_support_refs") or [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f845_refs.get("statement_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f845_refs.get("neighboring_statement_slot_refs") or [])
            ),
            "details": "F845 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the refined new lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "statement_form_class_candidate_supported_now": (
            checks[1]["pass"] and checks[4]["pass"]
        ),
        "exact_statement_required_form_exported_now": False if all_pass else None,
        "sharp_blocker_field": "exact_statement_required_form_ref" if all_pass else None,
        "next_honest_move_is_freeze_exact_required_form_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P846_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        if all_pass
        else "P846_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P846",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p845_exact_statement_probe": rel(IN_P845),
            "f845_refined_exact_statement_target_packet": rel(IN_F845),
            "f844_refined_output_schema_target_packet": rel(IN_F844),
            "f843_refined_lawful_admission_target_packet": rel(IN_F843),
            "f842_boundary_target_packet": rel(IN_F842),
            "f841_rule_target_packet": rel(IN_F841),
            "f840_interface_target_packet": rel(IN_F840),
            "f836_neighboring_exact_required_form_target_packet": rel(IN_F836),
        },
        "repo_scan_object_hits_for_exact_statement_required_form_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "statement_form_class_support_stack": {
            "candidate_support_refs": [
                "exact_output_schema_statement",
                "lawful_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring form-class supports only; none exports the exact refined statement-required form needed by F845.",
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the refined new lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F845 for the refined lawful T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form object itself.",
        ],
        "recommended_next_packet": {
            "id": "F846_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET",
            "goal": "Freeze the exact refined statement-required form object still missing after statement-form class support is present only at candidate grade.",
            "minimum_fields": [
                "exact_output_schema_statement_target_ref",
                "statement_form_class_candidate_support_refs",
                "neighboring_form_slot_refs",
                "exact_required_form_statement_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P846",
        "status": status,
        "as_of": AS_OF,
        "statement_form_class_candidate_supported_now": theorem_result[
            "statement_form_class_candidate_supported_now"
        ],
        "exact_statement_required_form_exported_now": theorem_result[
            "exact_statement_required_form_exported_now"
        ],
        "sharp_blocker_field": theorem_result["sharp_blocker_field"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
