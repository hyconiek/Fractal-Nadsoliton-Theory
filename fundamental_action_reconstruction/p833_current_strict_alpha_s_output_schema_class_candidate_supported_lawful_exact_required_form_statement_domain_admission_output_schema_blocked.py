#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P832 = GENERATED / "p832_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_exact_required_form_statement_domain_admission_blocked.json"
IN_F832 = GENERATED / "f832_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_target_packet.json"
IN_F831 = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F829 = GENERATED / "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p833_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p833_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_output_schema_blocked_summary.json"


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

    prereq = [IN_P832, IN_F832, IN_F831, IN_F830, IN_F829, IN_F822, IN_F825]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P833",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p832 = load_json(IN_P832)
    f832 = load_json(IN_F832)
    f831 = load_json(IN_F831)
    f830 = load_json(IN_F830)
    f829 = load_json(IN_F829)
    f822 = load_json(IN_F822)
    f825 = load_json(IN_F825)

    p832_theorem = p832.get("theorem_result") or {}
    f832_target = f832.get("target_object") or {}
    f831_target = f831.get("target_object") or {}
    f830_target = f830.get("target_object") or {}
    f829_target = f829.get("target_object") or {}
    f822_target = f822.get("target_object") or {}
    f825_target = f825.get("target_object") or {}
    f825_refs = f825.get("target_refs") or {}

    checks = [
        {
            "id": "f832_already_names_exact_missing_output_schema_field",
            "pass": (
                f832.get("status")
                == "F832_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
                and f832_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_target_v1"
                and has_required_field(
                    f832_target,
                    "lawful_exact_required_form_statement_domain_admission_output_schema",
                )
            ),
            "details": "F832 already isolates lawful_exact_required_form_statement_domain_admission_output_schema as one exact missing field of the lawful-admission target.",
        },
        {
            "id": "f831_preserves_combined_boundary_output_schema_class",
            "pass": (
                f831.get("status")
                == "F831_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f831_target, "boundary_output_schema")
            ),
            "details": "F831 preserves one combined boundary-output-schema class at the earlier admission/nonidentification layer.",
        },
        {
            "id": "f830_and_f829_preserve_upstream_output_schema_classes",
            "pass": (
                f830.get("status")
                == "F830_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f830_target, "selected_interface_output_schema")
                and f829.get("status")
                == "F829_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f829_target, "exact_interface_output_schema")
            ),
            "details": "F830 and F829 preserve output-schema classes on the immediate upstream rule/interface stack without exporting the new lawful-admission output schema.",
        },
        {
            "id": "f822_and_f825_preserve_neighboring_output_schema_or_statement_classes",
            "pass": (
                f822.get("status")
                == "F822_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f822_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1"
                and has_required_field(f822_target, "exact_output_schema_statement")
                and f825.get("status")
                == "F825_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f825_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                and "lawful_schema_domain_admission_output_schema"
                in (f825_refs.get("neighboring_statement_or_form_slot_refs") or [])
                and "exact_output_schema_statement"
                in (f825_refs.get("neighboring_statement_or_form_slot_refs") or [])
            ),
            "details": "F822 and F825 preserve neighboring output-schema or output-statement classes near the same downstream exact-required-form-statement lane, but not the new lawful exact-required-form-statement domain-admission output schema.",
        },
        {
            "id": "p832_already_keeps_lawful_exact_required_form_statement_domain_admission_unexported",
            "pass": (
                p832.get("status")
                == "P832_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_BLOCKED"
                and p832_theorem.get("lawful_exact_required_form_statement_domain_admission_exported_now")
                is False
            ),
            "details": "P832 already keeps the lawful exact-required-form-statement domain-admission object itself unexported on the current repo state.",
        },
        {
            "id": "no_exact_new_lane_output_schema_export_present",
            "pass": (
                f832_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_target_v1"
                and f822_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_v1"
                and f825_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_v1"
            ),
            "details": "No current export names one exact lawful_exact_required_form_statement_domain_admission_output_schema object for the new T213/T216 -> alpha_s exact-required-form-statement lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "output_schema_class_candidate_supported_now": (
            checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"]
        ),
        "lawful_exact_required_form_statement_domain_admission_output_schema_exported_now": False
        if all_pass
        else None,
        "sharp_blocker_field": (
            "lawful_exact_required_form_statement_domain_admission_output_schema"
            if all_pass
            else None
        ),
        "next_honest_move_is_freeze_exact_output_schema_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P833_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P833_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P833",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p832_lawful_exact_required_form_statement_domain_admission_probe": rel(IN_P832),
            "f832_lawful_exact_required_form_statement_domain_admission_target_packet": rel(IN_F832),
            "f831_combined_boundary_target_packet": rel(IN_F831),
            "f830_upstream_rule_target_packet": rel(IN_F830),
            "f829_upstream_interface_target_packet": rel(IN_F829),
            "f822_neighboring_schema_output_target_packet": rel(IN_F822),
            "f825_downstream_exact_required_form_statement_target_packet": rel(IN_F825),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "candidate_support_refs": [
                "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_target_v1",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These objects preserve output-schema class or neighboring output-statement class only; none exports the exact lawful exact-required-form-statement domain-admission output schema required by F832.",
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class around the new lane at upstream, boundary, and neighboring downstream levels.",
            "But no current export supplies the exact lawful_exact_required_form_statement_domain_admission_output_schema required by F832 for the new T213/T216 -> alpha_s exact-required-form-statement route.",
            "So the sharp blocker is now the still-missing exact output-schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F833_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact lawful exact-required-form-statement domain-admission output-schema object still missing after output-schema class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_exact_required_form_statement_domain_admission_target_ref",
                "output_schema_class_candidate_support_refs",
                "upstream_rule_or_interface_output_refs",
                "neighboring_output_schema_or_statement_refs",
                "exact_output_schema_statement",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P833",
        "status": status,
        "as_of": AS_OF,
        "output_schema_class_candidate_supported_now": theorem_result[
            "output_schema_class_candidate_supported_now"
        ],
        "lawful_exact_required_form_statement_domain_admission_output_schema_exported_now": theorem_result[
            "lawful_exact_required_form_statement_domain_admission_output_schema_exported_now"
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
