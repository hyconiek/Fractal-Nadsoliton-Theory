#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P909 = GENERATED / "p909_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_blocked.json"
IN_F909 = GENERATED / "f909_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F908 = GENERATED / "f908_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F907 = GENERATED / "f907_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_adapter_or_carrier_rule_target_packet.json"
IN_F906 = GENERATED / "f906_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_interface_target_packet.json"
IN_F899 = GENERATED / "f899_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F902 = GENERATED / "f902_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p910_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p910_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_blocked_summary.json"


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

    prereq = [IN_P909, IN_F909, IN_F908, IN_F907, IN_F906, IN_F899, IN_F902]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P910",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p909 = load_json(IN_P909)
    f909 = load_json(IN_F909)
    f908 = load_json(IN_F908)
    f907 = load_json(IN_F907)
    f906 = load_json(IN_F906)
    f899 = load_json(IN_F899)
    f902 = load_json(IN_F902)

    p909_theorem = p909.get("theorem_result") or {}
    f909_target = f909.get("target_object") or {}
    f908_target = f908.get("target_object") or {}
    f907_target = f907.get("target_object") or {}
    f906_target = f906.get("target_object") or {}
    f899_target = f899.get("target_object") or {}
    f902_target = f902.get("target_object") or {}
    f902_refs = f902.get("target_refs") or {}

    checks = [
        {
            "id": "f909_already_names_exact_missing_output_schema_field",
            "pass": (
                f909.get("status")
                == "F909_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
                and f909_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_target_v1"
                and has_required_field(
                    f909_target,
                    "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                )
            ),
            "details": "F909 already isolates the deeper lawful domain-admission output schema as one exact missing field of the lawful refined admission target.",
        },
        {
            "id": "f908_preserves_combined_boundary_output_schema_class",
            "pass": (
                f908.get("status")
                == "F908_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f908_target, "boundary_output_schema")
            ),
            "details": "F908 preserves one combined lawful refined boundary-output-schema class at the current deeper admission/nonidentification layer.",
        },
        {
            "id": "f907_and_f906_preserve_upstream_output_schema_classes",
            "pass": (
                f907.get("status")
                == "F907_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_ADAPTER_OR_CARRIER_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f907_target, "selected_interface_output_schema")
                and f906.get("status")
                == "F906_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f906_target, "exact_interface_output_schema")
            ),
            "details": "F907 and F906 preserve output-schema classes on the immediate upstream lawful refined rule/interface stack without exporting the new deeper lawful refined output schema.",
        },
        {
            "id": "f899_and_f902_preserve_neighboring_output_schema_or_statement_classes",
            "pass": (
                f899.get("status")
                == "F899_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f899_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and has_required_field(f899_target, "exact_output_schema_statement")
                and f902.get("status")
                == "F902_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f902_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                and "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema"
                in (f902_refs.get("neighboring_statement_or_form_slot_refs") or [])
                and "exact_output_schema_statement"
                in (f902_refs.get("neighboring_statement_or_form_slot_refs") or [])
            ),
            "details": "F899 and F902 preserve neighboring output-schema or output-statement classes near the same lawful refined lane family, but not the new deeper lawful refined output schema.",
        },
        {
            "id": "p909_already_keeps_lawful_refined_deeper_domain_admission_unexported",
            "pass": (
                p909.get("status")
                == "P909_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_BLOCKED"
                and p909_theorem.get("lawful_refined_deeper_domain_admission_exact_required_form_statement_domain_admission_exported_now")
                is False
            ),
            "details": "P909 already keeps the deeper lawful refined domain-admission object itself unexported on the current repo state.",
        },
        {
            "id": "no_exact_new_lawful_refined_deeper_domain_admission_output_schema_export_present",
            "pass": (
                f909_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_target_v1"
                and f899_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and f902_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
            ),
            "details": "No current export names one exact lawful refined deeper domain-admission output-schema object for the new lawful T213/T216 -> alpha_s lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "output_schema_class_candidate_supported_now": (
            checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"]
        ),
        "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema_exported_now": False
        if all_pass
        else None,
        "sharp_blocker_field": (
            "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema"
            if all_pass
            else None
        ),
        "next_honest_move_is_freeze_exact_output_schema_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P910_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P910_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P910",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p909_lawful_refined_deeper_domain_admission_probe": rel(IN_P909),
            "f909_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F909),
            "f908_combined_boundary_target_packet": rel(IN_F908),
            "f907_upstream_rule_target_packet": rel(IN_F907),
            "f906_upstream_interface_target_packet": rel(IN_F906),
            "f899_neighboring_output_target_packet": rel(IN_F899),
            "f902_neighboring_exact_required_form_statement_target_packet": rel(IN_F902),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "candidate_support_refs": [
                f908_target.get("object_id"),
                f907_target.get("object_id"),
                f906_target.get("object_id"),
                f899_target.get("object_id"),
                f902_target.get("object_id"),
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These objects preserve output-schema class or neighboring output-statement class only; none exports the exact lawful refined deeper domain-admission output schema required by F909.",
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class around the new deeper lawful refined domain-admission lane at upstream, boundary, and neighboring downstream levels.",
            "But no current export supplies the exact lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema required by F909 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact output-schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F910_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined deeper domain-admission output-schema object still missing after output-schema class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_refined_deeper_exact_required_form_statement_domain_admission_target_ref",
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
        "stage": "P910",
        "status": status,
        "as_of": AS_OF,
        "output_schema_class_candidate_supported_now": theorem_result[
            "output_schema_class_candidate_supported_now"
        ],
        "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema_exported_now": theorem_result[
            "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema_exported_now"
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
