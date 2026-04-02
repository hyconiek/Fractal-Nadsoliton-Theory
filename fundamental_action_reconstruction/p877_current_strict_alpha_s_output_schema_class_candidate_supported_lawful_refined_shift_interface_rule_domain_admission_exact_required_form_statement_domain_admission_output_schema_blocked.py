#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P876 = GENERATED / "p876_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_blocked.json"
IN_F876 = GENERATED / "f876_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F875 = GENERATED / "f875_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F874 = GENERATED / "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F873 = GENERATED / "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json"
IN_F866 = GENERATED / "f866_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p877_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p877_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_blocked_summary.json"


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

    prereq = [IN_P876, IN_F876, IN_F875, IN_F874, IN_F873, IN_F866, IN_F869]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P877",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p876 = load_json(IN_P876)
    f876 = load_json(IN_F876)
    f875 = load_json(IN_F875)
    f874 = load_json(IN_F874)
    f873 = load_json(IN_F873)
    f866 = load_json(IN_F866)
    f869 = load_json(IN_F869)

    p876_theorem = p876.get("theorem_result") or {}
    f876_target = f876.get("target_object") or {}
    f875_target = f875.get("target_object") or {}
    f874_target = f874.get("target_object") or {}
    f873_target = f873.get("target_object") or {}
    f866_target = f866.get("target_object") or {}
    f869_target = f869.get("target_object") or {}
    f869_refs = f869.get("target_refs") or {}

    checks = [
        {
            "id": "f876_already_names_exact_missing_output_schema_field",
            "pass": (
                f876.get("status")
                == "F876_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
                and f876_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_v1"
                and has_required_field(
                    f876_target,
                    "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                )
            ),
            "details": "F876 already isolates lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema as one exact missing field of the lawful refined admission target.",
        },
        {
            "id": "f875_preserves_combined_boundary_output_schema_class",
            "pass": (
                f875.get("status")
                == "F875_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f875_target, "boundary_output_schema")
            ),
            "details": "F875 preserves one combined lawful refined boundary-output-schema class at the earlier admission/nonidentification layer.",
        },
        {
            "id": "f874_and_f873_preserve_upstream_output_schema_classes",
            "pass": (
                f874.get("status")
                == "F874_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f874_target, "selected_interface_output_schema")
                and f873.get("status")
                == "F873_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f873_target, "exact_interface_output_schema")
            ),
            "details": "F874 and F873 preserve output-schema classes on the immediate upstream lawful refined rule/interface stack without exporting the new lawful refined output schema.",
        },
        {
            "id": "f866_and_f869_preserve_neighboring_output_schema_or_statement_classes",
            "pass": (
                f866.get("status")
                == "F866_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f866_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1"
                and has_required_field(f866_target, "exact_output_schema_statement")
                and f869.get("status")
                == "F869_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f869_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_v1"
                and "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema"
                in (f869_refs.get("neighboring_statement_or_form_slot_refs") or [])
                and "exact_output_schema_statement"
                in (f869_refs.get("neighboring_statement_or_form_slot_refs") or [])
            ),
            "details": "F866 and F869 preserve neighboring output-schema or output-statement classes near the same lawful refined lane family, but not the new shift-interface-rule domain-admission exact-required-form-statement domain-admission output schema.",
        },
        {
            "id": "p876_already_keeps_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_unexported",
            "pass": (
                p876.get("status")
                == "P876_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_BLOCKED"
                and p876_theorem.get(
                    "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exported_now"
                )
                is False
            ),
            "details": "P876 already keeps the lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission object itself unexported on the current repo state.",
        },
        {
            "id": "no_exact_new_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_export_present",
            "pass": (
                f876_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_v1"
                and f866_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and f869_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
            ),
            "details": "No current export names one exact lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema object for the new lawful T213/T216 -> alpha_s lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "output_schema_class_candidate_supported_now": (
            checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"]
        ),
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_exported_now": False
        if all_pass
        else None,
        "sharp_blocker_field": (
            "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema"
            if all_pass
            else None
        ),
        "next_honest_move_is_freeze_exact_output_schema_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P877_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P877_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P877",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p876_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_probe": rel(IN_P876),
            "f876_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet": rel(IN_F876),
            "f875_combined_boundary_target_packet": rel(IN_F875),
            "f874_upstream_rule_target_packet": rel(IN_F874),
            "f873_upstream_interface_target_packet": rel(IN_F873),
            "f866_neighboring_output_target_packet": rel(IN_F866),
            "f869_neighboring_exact_required_form_statement_target_packet": rel(IN_F869),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "candidate_support_refs": [
                f875_target.get("object_id"),
                f874_target.get("object_id"),
                f873_target.get("object_id"),
                f866_target.get("object_id"),
                f869_target.get("object_id"),
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These objects preserve output-schema class or neighboring output-statement class only; none exports the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission output schema required by F876.",
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class around the new lawful refined shift-interface-rule domain-admission exact-required-form-statement lane at upstream, boundary, and neighboring downstream levels.",
            "But no current export supplies the exact lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema required by F876 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact output-schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F877_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission output-schema object still missing after output-schema class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_ref",
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
        "stage": "P877",
        "status": status,
        "as_of": AS_OF,
        "output_schema_class_candidate_supported_now": theorem_result[
            "output_schema_class_candidate_supported_now"
        ],
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_exported_now": theorem_result[
            "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_exported_now"
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
