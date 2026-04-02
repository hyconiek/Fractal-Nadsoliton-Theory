#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P865 = GENERATED / "p865_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_blocked.json"
IN_F865 = GENERATED / "f865_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_packet.json"
IN_F864 = GENERATED / "f864_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F863 = GENERATED / "f863_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F862 = GENERATED / "f862_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_interface_target_packet.json"
IN_F855 = GENERATED / "f855_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F858 = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p866_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p866_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_blocked_summary.json"


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

    prereq = [IN_P865, IN_F865, IN_F864, IN_F863, IN_F862, IN_F855, IN_F858]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P866",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p865 = load_json(IN_P865)
    f865 = load_json(IN_F865)
    f864 = load_json(IN_F864)
    f863 = load_json(IN_F863)
    f862 = load_json(IN_F862)
    f855 = load_json(IN_F855)
    f858 = load_json(IN_F858)

    p865_theorem = p865.get("theorem_result") or {}
    f865_target = f865.get("target_object") or {}
    f864_target = f864.get("target_object") or {}
    f863_target = f863.get("target_object") or {}
    f862_target = f862.get("target_object") or {}
    f855_target = f855.get("target_object") or {}
    f858_target = f858.get("target_object") or {}
    f858_refs = f858.get("target_refs") or {}

    checks = [
        {
            "id": "f865_already_names_exact_missing_output_schema_field",
            "pass": (
                f865.get("status")
                == "F865_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
                and f865_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_v1"
                and has_required_field(
                    f865_target,
                    "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema",
                )
            ),
            "details": "F865 already isolates lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema as one exact missing field of the lawful refined admission target.",
        },
        {
            "id": "f864_preserves_combined_boundary_output_schema_class",
            "pass": (
                f864.get("status")
                == "F864_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f864_target, "boundary_output_schema")
            ),
            "details": "F864 preserves one combined lawful refined boundary-output-schema class at the earlier admission/nonidentification layer.",
        },
        {
            "id": "f863_and_f862_preserve_upstream_output_schema_classes",
            "pass": (
                f863.get("status")
                == "F863_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f863_target, "selected_interface_output_schema")
                and f862.get("status")
                == "F862_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f862_target, "exact_interface_output_schema")
            ),
            "details": "F863 and F862 preserve output-schema classes on the immediate upstream lawful refined rule/interface stack without exporting the new lawful refined output schema.",
        },
        {
            "id": "f855_and_f858_preserve_neighboring_output_schema_or_statement_classes",
            "pass": (
                f855.get("status")
                == "F855_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f855_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and has_required_field(f855_target, "exact_output_schema_statement")
                and f858.get("status")
                == "F858_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f858_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                and "lawful_refined_exact_required_form_statement_domain_admission_output_schema"
                in (f858_refs.get("neighboring_statement_or_form_slot_refs") or [])
                and "exact_output_schema_statement"
                in (f858_refs.get("neighboring_statement_or_form_slot_refs") or [])
            ),
            "details": "F855 and F858 preserve neighboring output-schema or output-statement classes near the same lawful refined lane family, but not the new shift-interface-rule domain-admission output schema.",
        },
        {
            "id": "p865_already_keeps_lawful_refined_shift_interface_rule_domain_admission_unexported",
            "pass": (
                p865.get("status")
                == "P865_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_BLOCKED"
                and p865_theorem.get(
                    "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exported_now"
                )
                is False
            ),
            "details": "P865 already keeps the lawful refined shift-interface-rule domain-admission object itself unexported on the current repo state.",
        },
        {
            "id": "no_exact_new_lawful_refined_shift_interface_rule_lane_output_schema_export_present",
            "pass": (
                f865_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_v1"
                and f855_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1"
                and f858_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1"
            ),
            "details": "No current export names one exact lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema object for the new lawful T213/T216 -> alpha_s lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "output_schema_class_candidate_supported_now": (
            checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"]
        ),
        "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_exported_now": False
        if all_pass
        else None,
        "sharp_blocker_field": (
            "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema"
            if all_pass
            else None
        ),
        "next_honest_move_is_freeze_exact_output_schema_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P866_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P866_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P866",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p865_lawful_refined_shift_interface_rule_domain_admission_probe": rel(IN_P865),
            "f865_lawful_refined_shift_interface_rule_domain_admission_target_packet": rel(IN_F865),
            "f864_combined_boundary_target_packet": rel(IN_F864),
            "f863_upstream_rule_target_packet": rel(IN_F863),
            "f862_upstream_interface_target_packet": rel(IN_F862),
            "f855_neighboring_output_target_packet": rel(IN_F855),
            "f858_neighboring_exact_required_form_statement_target_packet": rel(IN_F858),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "candidate_support_refs": [
                f864_target.get("object_id"),
                f863_target.get("object_id"),
                f862_target.get("object_id"),
                f855_target.get("object_id"),
                f858_target.get("object_id"),
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These objects preserve output-schema class or neighboring output-statement class only; none exports the exact lawful refined shift-interface-rule domain-admission output schema required by F865.",
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class around the new lawful refined shift-interface-rule lane at upstream, boundary, and neighboring downstream levels.",
            "But no current export supplies the exact lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema required by F865 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact output-schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F866_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined exact-required-form-statement shift-interface-rule domain-admission output-schema object still missing after output-schema class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_ref",
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
        "stage": "P866",
        "status": status,
        "as_of": AS_OF,
        "output_schema_class_candidate_supported_now": theorem_result[
            "output_schema_class_candidate_supported_now"
        ],
        "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_exported_now": theorem_result[
            "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_exported_now"
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
