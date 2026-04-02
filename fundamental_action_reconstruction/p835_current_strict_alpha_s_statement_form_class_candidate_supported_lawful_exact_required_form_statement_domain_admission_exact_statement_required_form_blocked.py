#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P834 = GENERATED / "p834_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_F834 = GENERATED / "f834_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F833 = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F832 = GENERATED / "f832_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_target_packet.json"
IN_F831 = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F829 = GENERATED / "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"

OUT = GENERATED / "p835_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
OUT_SUMMARY = GENERATED / "p835_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked_summary.json"


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

    prereq = [
        IN_P834,
        IN_F834,
        IN_F833,
        IN_F832,
        IN_F831,
        IN_F830,
        IN_F829,
        IN_F825,
        IN_F822,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P835",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p834 = load_json(IN_P834)
    f834 = load_json(IN_F834)
    f833 = load_json(IN_F833)
    f832 = load_json(IN_F832)
    f831 = load_json(IN_F831)
    f830 = load_json(IN_F830)
    f829 = load_json(IN_F829)
    f825 = load_json(IN_F825)
    f822 = load_json(IN_F822)

    p834_theorem = p834.get("theorem_result") or {}
    f834_target = f834.get("target_object") or {}
    f833_target = f833.get("target_object") or {}
    f833_refs = f833.get("target_refs") or {}
    f832_target = f832.get("target_object") or {}
    f831_target = f831.get("target_object") or {}
    f830_target = f830.get("target_object") or {}
    f829_target = f829.get("target_object") or {}
    f825_target = f825.get("target_object") or {}
    f822_target = f822.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form"
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p834_", "f834_", "p835_", "f835_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f834_already_names_exact_missing_form_field",
            "pass": (
                f834.get("status")
                == "F834_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f834_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and has_required_field(f834_target, "exact_statement_required_form_ref")
            ),
            "details": "F834 already isolates exact_statement_required_form_ref as one exact missing field of the statement target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_form_class_support",
            "pass": (
                has_required_field(f833_target, "exact_output_schema_statement")
                and has_required_field(f832_target, "lawful_exact_required_form_statement_domain_admission_output_schema")
                and has_required_field(f831_target, "boundary_output_schema")
                and has_required_field(f830_target, "selected_interface_output_schema")
                and has_required_field(f829_target, "exact_interface_output_schema")
                and f825_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                and has_required_field(f822_target, "exact_output_schema_statement")
            ),
            "details": "Neighboring lanes preserve output/form structure only through target fields and neighboring required-form objects, not through an exact required-form export for the new lane.",
        },
        {
            "id": "p834_already_keeps_exact_statement_unexported",
            "pass": (
                p834.get("status")
                == "P834_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p834_theorem.get("exact_output_schema_statement_exported_now") is False
            ),
            "details": "P834 already keeps the exact new-lane statement object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_exact_required_form_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_statement_required_form_ref outside the new frozen F834 lineage itself.",
        },
        {
            "id": "neighboring_form_support_remains_nonidentical_to_new_lane_form",
            "pass": (
                f834_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1"
                in (f834.get("target_refs") or {}).get("statement_class_candidate_support_refs", [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                in (f834.get("target_refs") or {}).get("statement_class_candidate_support_refs", [])
                and "exact_output_schema_statement"
                in (f834.get("target_refs") or {}).get("neighboring_statement_slot_refs", [])
            ),
            "details": "F834 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the new lane.",
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
        "P835_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        if all_pass
        else "P835_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P835",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p834_exact_statement_probe": rel(IN_P834),
            "f834_exact_statement_target_packet": rel(IN_F834),
            "f833_output_schema_target_packet": rel(IN_F833),
            "f832_lawful_admission_target_packet": rel(IN_F832),
            "f831_boundary_target_packet": rel(IN_F831),
            "f830_rule_target_packet": rel(IN_F830),
            "f829_interface_target_packet": rel(IN_F829),
            "f825_neighboring_exact_required_form_target_packet": rel(IN_F825),
            "f822_neighboring_schema_output_target_packet": rel(IN_F822),
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
                "lawful_schema_domain_admission_exact_required_form_statement_target",
            ],
            "support_grade": "candidate_supported_not_yet-exported",
            "scope_limit": "These are neighboring form-class supports only; none exports the exact statement-required form needed by F834.",
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the new lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F834 for the new T213/T216 -> alpha_s exact-required-form-statement route.",
            "So the sharp blocker is now the still-missing exact required-form object itself.",
        ],
        "recommended_next_packet": {
            "id": "F835_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET",
            "goal": "Freeze the exact statement-required form object still missing after statement-form class support is present only at candidate grade.",
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
        "stage": "P835",
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
