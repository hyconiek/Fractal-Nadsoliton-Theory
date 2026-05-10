#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P843 = GENERATED / "p843_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_exact_required_form_statement_domain_admission_blocked.json"
IN_F843 = GENERATED / "f843_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_target_packet.json"
IN_F842 = GENERATED / "f842_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F841 = GENERATED / "f841_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F840 = GENERATED / "f840_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F833 = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F836 = GENERATED / "f836_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p844_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p844_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_output_schema_blocked_summary.json"


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

    prereq = [IN_P843, IN_F843, IN_F842, IN_F841, IN_F840, IN_F833, IN_F836]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P844",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p843 = load_json(IN_P843)
    f843 = load_json(IN_F843)
    f842 = load_json(IN_F842)
    f841 = load_json(IN_F841)
    f840 = load_json(IN_F840)
    f833 = load_json(IN_F833)
    f836 = load_json(IN_F836)

    p843_theorem = p843.get("theorem_result") or {}
    f843_target = f843.get("target_object") or {}
    f842_target = f842.get("target_object") or {}
    f841_target = f841.get("target_object") or {}
    f840_target = f840.get("target_object") or {}
    f833_target = f833.get("target_object") or {}
    f836_target = f836.get("target_object") or {}
    f836_refs = f836.get("target_refs") or {}

    checks = [
        {
            "id": "f843_already_names_exact_missing_output_schema_field",
            "pass": (
                f843.get("status")
                == "F843_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
                and f843_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_target_v1"
                and has_required_field(
                    f843_target,
                    "lawful_exact_required_form_statement_domain_admission_output_schema",
                )
            ),
            "details": "F843 already isolates lawful_exact_required_form_statement_domain_admission_output_schema as one exact missing field of the lawful-admission target.",
        },
        {
            "id": "f842_preserves_combined_boundary_output_schema_class",
            "pass": (
                f842.get("status")
                == "F842_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f842_target, "boundary_output_schema")
            ),
            "details": "F842 preserves one combined lawful boundary-output-schema class at the earlier admission/nonidentification layer.",
        },
        {
            "id": "f841_and_f840_preserve_upstream_output_schema_classes",
            "pass": (
                f841.get("status")
                == "F841_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f841_target, "selected_interface_output_schema")
                and f840.get("status")
                == "F840_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f840_target, "exact_interface_output_schema")
            ),
            "details": "F841 and F840 preserve output-schema classes on the immediate upstream rule/interface stack without exporting the new lawful-admission output schema.",
        },
        {
            "id": "f833_and_f836_preserve_neighboring_output_schema_or_statement_classes",
            "pass": (
                f833.get("status")
                == "F833_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f833_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and has_required_field(f833_target, "exact_output_schema_statement")
                and f836.get("status")
                == "F836_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f836_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                and "lawful_exact_required_form_statement_domain_admission_output_schema"
                in (f836_refs.get("neighboring_statement_or_form_slot_refs") or [])
                and "exact_output_schema_statement"
                in (f836_refs.get("neighboring_statement_or_form_slot_refs") or [])
            ),
            "details": "F833 and F836 preserve neighboring output-schema or output-statement classes near the same downstream lawful lane, but not the new refined lawful-admission output schema.",
        },
        {
            "id": "p843_already_keeps_lawful_exact_required_form_statement_domain_admission_unexported",
            "pass": (
                p843.get("status")
                == "P843_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_BLOCKED"
                and p843_theorem.get("lawful_exact_required_form_statement_domain_admission_exported_now")
                is False
            ),
            "details": "P843 already keeps the lawful exact-required-form-statement domain-admission object itself unexported on the current repo state.",
        },
        {
            "id": "no_exact_new_lawful_lane_output_schema_export_present",
            "pass": (
                f843_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_target_v1"
                and f833_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_output_schema_v1"
                and f836_target.get("object_id")
                != "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_output_schema_v1"
            ),
            "details": "No current export names one exact lawful_exact_required_form_statement_domain_admission_output_schema object for the new lawful T213/T216 -> alpha_s lane.",
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
        "P844_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P844_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P844",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p843_lawful_exact_required_form_statement_domain_admission_probe": rel(IN_P843),
            "f843_lawful_exact_required_form_statement_domain_admission_target_packet": rel(IN_F843),
            "f842_combined_boundary_target_packet": rel(IN_F842),
            "f841_upstream_rule_target_packet": rel(IN_F841),
            "f840_upstream_interface_target_packet": rel(IN_F840),
            "f833_old_lawful_output_target_packet": rel(IN_F833),
            "f836_downstream_exact_required_form_statement_target_packet": rel(IN_F836),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "candidate_support_refs": [
                f842_target.get("object_id"),
                f841_target.get("object_id"),
                f840_target.get("object_id"),
                f833_target.get("object_id"),
                f836_target.get("object_id"),
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These objects preserve output-schema class or neighboring output-statement class only; none exports the exact lawful exact-required-form-statement domain-admission output schema required by F843.",
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class around the new lawful lane at upstream, boundary, and neighboring downstream levels.",
            "But no current export supplies the exact lawful_exact_required_form_statement_domain_admission_output_schema required by F843 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact output-schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F844_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET",
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
        "stage": "P844",
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
