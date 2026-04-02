#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P879 = GENERATED / "p879_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F879 = GENERATED / "f879_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F878 = GENERATED / "f878_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F877 = GENERATED / "f877_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F876 = GENERATED / "f876_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F875 = GENERATED / "f875_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F874 = GENERATED / "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F873 = GENERATED / "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p880_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p880_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked_summary.json"


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

    prereq = [IN_P879, IN_F879, IN_F878, IN_F877, IN_F876, IN_F875, IN_F874, IN_F873, IN_F869]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P880",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p879 = load_json(IN_P879)
    f879 = load_json(IN_F879)
    f878 = load_json(IN_F878)
    f877 = load_json(IN_F877)
    f876 = load_json(IN_F876)
    f875 = load_json(IN_F875)
    f874 = load_json(IN_F874)
    f873 = load_json(IN_F873)
    f869 = load_json(IN_F869)

    p879_theorem = p879.get("theorem_result") or {}
    f879_target = f879.get("target_object") or {}
    f879_refs = f879.get("target_refs") or {}
    f878_target = f878.get("target_object") or {}
    f877_target = f877.get("target_object") or {}
    f876_target = f876.get("target_object") or {}
    f875_target = f875.get("target_object") or {}
    f874_target = f874.get("target_object") or {}
    f873_target = f873.get("target_object") or {}
    f869_target = f869.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p879_", "f879_", "p880_", "f880_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f879_already_names_exact_missing_required_form_statement_field",
            "pass": (
                f879.get("status")
                == "F879_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f879_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f879_target, "exact_required_form_statement_ref")
            ),
            "details": "F879 already isolates exact_required_form_statement_ref as one exact missing field of the lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission exact-statement-required-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_required_form_statement_class_support",
            "pass": (
                has_required_field(f878_target, "exact_statement_required_form_ref")
                and has_required_field(f877_target, "exact_output_schema_statement")
                and has_required_field(
                    f876_target,
                    "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                )
                and has_required_field(f875_target, "boundary_output_schema")
                and has_required_field(f874_target, "selected_interface_output_schema")
                and has_required_field(f873_target, "exact_interface_output_schema")
                and f869_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve lawful refined statement or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact required-form-statement export for the new shift-interface-rule domain-admission lane.",
        },
        {
            "id": "p879_already_keeps_lawful_refined_exact_statement_required_form_unexported",
            "pass": (
                p879.get("status")
                == "P879_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p879_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P879 already keeps the exact lawful refined new-lane statement-required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_required_form_statement_ref outside the new frozen F879 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_required_form_statement",
            "pass": (
                f879_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                == f879_refs.get("exact_output_schema_statement_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_v1"
                in (f879_refs.get("statement_form_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f879_refs.get("neighboring_form_slot_refs") or [])
                and "exact_required_form_statement"
                in (f879_refs.get("neighboring_form_slot_refs") or [])
            ),
            "details": "F879 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the lawful refined shift-interface-rule domain-admission lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "required_form_statement_class_candidate_supported_now": (
            checks[1]["pass"] and checks[4]["pass"]
        ),
        "exact_required_form_statement_exported_now": False if all_pass else None,
        "sharp_blocker_field": "exact_required_form_statement_ref" if all_pass else None,
        "next_honest_move_is_freeze_exact_required_form_statement_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P880_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P880_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P880",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p879_required_form_probe": rel(IN_P879),
            "f879_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet": rel(IN_F879),
            "f878_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F878),
            "f877_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet": rel(IN_F877),
            "f876_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet": rel(IN_F876),
            "f875_boundary_target_packet": rel(IN_F875),
            "f874_rule_target_packet": rel(IN_F874),
            "f873_interface_target_packet": rel(IN_F873),
            "f869_neighboring_exact_required_form_statement_target_packet": rel(IN_F869),
        },
        "repo_scan_object_hits_for_exact_required_form_statement_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "required_form_statement_class_support_stack": {
            "candidate_support_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission required-form statement needed by F879.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the lawful refined shift-interface-rule domain-admission lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F879 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F880_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission required-form statement object still missing after required-form-statement class support is present only at candidate grade.",
            "minimum_fields": [
                "exact_statement_required_form_target_ref",
                "required_form_statement_class_candidate_support_refs",
                "neighboring_statement_or_form_slot_refs",
                "exact_required_form_statement_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P880",
        "status": status,
        "as_of": AS_OF,
        "required_form_statement_class_candidate_supported_now": theorem_result[
            "required_form_statement_class_candidate_supported_now"
        ],
        "exact_required_form_statement_exported_now": theorem_result[
            "exact_required_form_statement_exported_now"
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
