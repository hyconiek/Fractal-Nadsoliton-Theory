#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P878 = GENERATED / "p878_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_F878 = GENERATED / "f878_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F877 = GENERATED / "f877_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F876 = GENERATED / "f876_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F875 = GENERATED / "f875_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F874 = GENERATED / "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F873 = GENERATED / "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p879_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
OUT_SUMMARY = GENERATED / "p879_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked_summary.json"


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

    prereq = [IN_P878, IN_F878, IN_F877, IN_F876, IN_F875, IN_F874, IN_F873, IN_F869]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P879",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p878 = load_json(IN_P878)
    f878 = load_json(IN_F878)
    f877 = load_json(IN_F877)
    f876 = load_json(IN_F876)
    f875 = load_json(IN_F875)
    f874 = load_json(IN_F874)
    f873 = load_json(IN_F873)
    f869 = load_json(IN_F869)

    p878_theorem = p878.get("theorem_result") or {}
    f878_target = f878.get("target_object") or {}
    f878_refs = f878.get("target_refs") or {}
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
        "domain_admission_exact_statement_required_form"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p878_", "f878_", "p879_", "f879_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f878_already_names_exact_missing_form_field",
            "pass": (
                f878.get("status")
                == "F878_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f878_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and has_required_field(f878_target, "exact_statement_required_form_ref")
            ),
            "details": "F878 already isolates exact_statement_required_form_ref as one exact missing field of the lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission exact-statement target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_form_class_support",
            "pass": (
                has_required_field(f877_target, "exact_output_schema_statement")
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
            "details": "Neighboring lanes preserve lawful refined output or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact statement-required-form export for the new shift-interface-rule lane.",
        },
        {
            "id": "p878_already_keeps_lawful_refined_exact_statement_unexported",
            "pass": (
                p878.get("status")
                == "P878_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p878_theorem.get("exact_output_schema_statement_exported_now") is False
            ),
            "details": "P878 already keeps the exact lawful refined new-lane statement object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_statement_required_form_ref outside the new frozen F878 lineage itself.",
        },
        {
            "id": "neighboring_form_support_remains_nonidentical_to_new_lane_form",
            "pass": (
                f878_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                == f878_refs.get("lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_v1"
                in (f878_refs.get("statement_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f878_refs.get("neighboring_statement_slot_refs") or [])
                and "exact_required_form_statement"
                in (f878_refs.get("neighboring_statement_slot_refs") or [])
            ),
            "details": "F878 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the lawful refined shift-interface-rule lane.",
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
        "P879_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        if all_pass
        else "P879_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P879",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p878_exact_statement_probe": rel(IN_P878),
            "f878_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F878),
            "f877_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet": rel(IN_F877),
            "f876_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet": rel(IN_F876),
            "f875_boundary_target_packet": rel(IN_F875),
            "f874_rule_target_packet": rel(IN_F874),
            "f873_interface_target_packet": rel(IN_F873),
            "f869_neighboring_exact_required_form_statement_target_packet": rel(IN_F869),
        },
        "repo_scan_object_hits_for_exact_statement_required_form_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "statement_form_class_support_stack": {
            "candidate_support_refs": [
                "exact_output_schema_statement",
                "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring form-class supports only; none exports the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission statement-required form needed by F878.",
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the lawful refined shift-interface-rule lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F878 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form object itself.",
        ],
        "recommended_next_packet": {
            "id": "F879_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission statement-required form object still missing after statement-form class support is present only at candidate grade.",
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
        "stage": "P879",
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
