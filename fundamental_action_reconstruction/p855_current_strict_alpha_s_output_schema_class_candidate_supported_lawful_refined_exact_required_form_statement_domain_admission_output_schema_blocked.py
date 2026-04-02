#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P854 = GENERATED / "p854_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_blocked.json"
IN_F854 = GENERATED / "f854_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_target_packet.json"
IN_F853 = GENERATED / "f853_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F852 = GENERATED / "f852_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F851 = GENERATED / "f851_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_interface_target_packet.json"
IN_F844 = GENERATED / "f844_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F847 = GENERATED / "f847_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p855_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p855_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_output_schema_blocked_summary.json"


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

    prereq = [IN_P854, IN_F854, IN_F853, IN_F852, IN_F851, IN_F844, IN_F847]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P855",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p854 = load_json(IN_P854)
    f854 = load_json(IN_F854)
    f853 = load_json(IN_F853)
    f852 = load_json(IN_F852)
    f851 = load_json(IN_F851)
    f844 = load_json(IN_F844)
    f847 = load_json(IN_F847)

    p854_theorem = p854.get("theorem_result") or {}
    f854_target = f854.get("target_object") or {}
    f853_target = f853.get("target_object") or {}
    f852_target = f852.get("target_object") or {}
    f851_target = f851.get("target_object") or {}
    f844_target = f844.get("target_object") or {}
    f847_target = f847.get("target_object") or {}
    f847_refs = f847.get("target_refs") or {}

    checks = [
        {
            "id": "f854_already_names_exact_missing_output_schema_field",
            "pass": (
                f854.get("status")
                == "F854_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
                and f854_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_target_v1"
                and has_required_field(
                    f854_target,
                    "lawful_refined_exact_required_form_statement_domain_admission_output_schema",
                )
            ),
            "details": "F854 already isolates lawful_refined_exact_required_form_statement_domain_admission_output_schema as one exact missing field of the lawful refined admission target.",
        },
        {
            "id": "f853_preserves_combined_boundary_output_schema_class",
            "pass": (
                f853.get("status")
                == "F853_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f853_target, "boundary_output_schema")
            ),
            "details": "F853 preserves one combined lawful refined boundary-output-schema class at the earlier admission/nonidentification layer.",
        },
        {
            "id": "f852_and_f851_preserve_upstream_output_schema_classes",
            "pass": (
                f852.get("status")
                == "F852_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_REFINED_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f852_target, "selected_interface_output_schema")
                and f851.get("status")
                == "F851_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_REFINED_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f851_target, "exact_interface_output_schema")
            ),
            "details": "F852 and F851 preserve output-schema classes on the immediate upstream refined rule/interface stack without exporting the new lawful refined output schema.",
        },
        {
            "id": "f844_and_f847_preserve_neighboring_output_schema_or_statement_classes",
            "pass": (
                f844.get("status")
                == "F844_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f844_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_output_schema_target_v1"
                and has_required_field(f844_target, "exact_output_schema_statement")
                and f847.get("status")
                == "F847_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f847_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_target_v1"
                and "lawful_exact_required_form_statement_domain_admission_output_schema"
                in (f847_refs.get("neighboring_statement_or_form_slot_refs") or [])
                and "exact_output_schema_statement"
                in (f847_refs.get("neighboring_statement_or_form_slot_refs") or [])
            ),
            "details": "F844 and F847 preserve neighboring output-schema or output-statement classes near the same refined lawful lane, but not the new lawful refined admission output schema.",
        },
        {
            "id": "p854_already_keeps_lawful_refined_exact_required_form_statement_domain_admission_unexported",
            "pass": (
                p854.get("status")
                == "P854_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_BLOCKED"
                and p854_theorem.get("lawful_refined_exact_required_form_statement_domain_admission_exported_now")
                is False
            ),
            "details": "P854 already keeps the lawful refined exact-required-form-statement domain-admission object itself unexported on the current repo state.",
        },
        {
            "id": "no_exact_new_refined_lawful_lane_output_schema_export_present",
            "pass": (
                f854_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_target_v1"
                and f844_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and f847_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_v1"
            ),
            "details": "No current export names one exact lawful_refined_exact_required_form_statement_domain_admission_output_schema object for the new lawful T213/T216 -> alpha_s lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "output_schema_class_candidate_supported_now": (
            checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"]
        ),
        "lawful_refined_exact_required_form_statement_domain_admission_output_schema_exported_now": False
        if all_pass
        else None,
        "sharp_blocker_field": (
            "lawful_refined_exact_required_form_statement_domain_admission_output_schema"
            if all_pass
            else None
        ),
        "next_honest_move_is_freeze_exact_output_schema_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P855_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P855_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P855",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p854_lawful_refined_exact_required_form_statement_domain_admission_probe": rel(IN_P854),
            "f854_lawful_refined_exact_required_form_statement_domain_admission_target_packet": rel(IN_F854),
            "f853_combined_boundary_target_packet": rel(IN_F853),
            "f852_upstream_rule_target_packet": rel(IN_F852),
            "f851_upstream_interface_target_packet": rel(IN_F851),
            "f844_old_refined_lawful_output_target_packet": rel(IN_F844),
            "f847_downstream_refined_exact_required_form_statement_target_packet": rel(IN_F847),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "candidate_support_refs": [
                f853_target.get("object_id"),
                f852_target.get("object_id"),
                f851_target.get("object_id"),
                f844_target.get("object_id"),
                f847_target.get("object_id"),
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These objects preserve output-schema class or neighboring output-statement class only; none exports the exact lawful refined exact-required-form-statement domain-admission output schema required by F854.",
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class around the new lawful refined lane at upstream, boundary, and neighboring downstream levels.",
            "But no current export supplies the exact lawful_refined_exact_required_form_statement_domain_admission_output_schema required by F854 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact output-schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F855_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined exact-required-form-statement domain-admission output-schema object still missing after output-schema class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_refined_exact_required_form_statement_domain_admission_target_ref",
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
        "stage": "P855",
        "status": status,
        "as_of": AS_OF,
        "output_schema_class_candidate_supported_now": theorem_result[
            "output_schema_class_candidate_supported_now"
        ],
        "lawful_refined_exact_required_form_statement_domain_admission_output_schema_exported_now": theorem_result[
            "lawful_refined_exact_required_form_statement_domain_admission_output_schema_exported_now"
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
