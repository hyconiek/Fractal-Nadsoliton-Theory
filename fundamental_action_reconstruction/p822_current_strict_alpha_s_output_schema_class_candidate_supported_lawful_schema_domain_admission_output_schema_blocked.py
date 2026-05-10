#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P821 = GENERATED / "p821_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_schema_domain_admission_blocked.json"
IN_F821 = GENERATED / "f821_current_strict_alpha_s_lawful_schema_domain_admission_target_packet.json"
IN_F820 = GENERATED / "f820_current_strict_alpha_s_schema_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F818 = GENERATED / "f818_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_packet.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"

OUT = GENERATED / "p822_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_schema_domain_admission_output_schema_blocked.json"
OUT_SUMMARY = GENERATED / "p822_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_schema_domain_admission_output_schema_blocked_summary.json"


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

    prereq = [IN_P821, IN_F821, IN_F820, IN_F819, IN_F818, IN_F814]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P822",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p821 = load_json(IN_P821)
    f821 = load_json(IN_F821)
    f820 = load_json(IN_F820)
    f819 = load_json(IN_F819)
    f818 = load_json(IN_F818)
    f814 = load_json(IN_F814)

    p821_theorem = p821.get("theorem_result") or {}
    f821_target = f821.get("target_object") or {}
    f820_target = f820.get("target_object") or {}
    f819_target = f819.get("target_object") or {}
    f818_target = f818.get("target_object") or {}
    f814_target = f814.get("target_object") or {}

    checks = [
        {
            "id": "f821_already_names_exact_missing_output_schema_field",
            "pass": (
                f821.get("status")
                == "F821_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
                and f821_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_target_v1"
                and has_required_field(f821_target, "lawful_schema_domain_admission_output_schema")
            ),
            "details": "F821 already isolates lawful_schema_domain_admission_output_schema as one exact missing field of the lawful-admission target.",
        },
        {
            "id": "f820_preserves_combined_boundary_output_schema_class",
            "pass": (
                f820.get("status")
                == "F820_EXECUTED_CURRENT_STRICT_ALPHA_S_SCHEMA_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f820_target, "boundary_output_schema")
            ),
            "details": "F820 preserves one combined boundary-output-schema class at the earlier admission/nonidentification layer.",
        },
        {
            "id": "f819_and_f818_preserve_upstream_output_schema_classes",
            "pass": (
                f819.get("status")
                == "F819_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_SCHEMA_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f819_target, "selected_interface_output_schema")
                and f818.get("status")
                == "F818_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_SCHEMA_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f818_target, "exact_interface_output_schema")
            ),
            "details": "F819 and F818 preserve output-schema classes on the immediate upstream rule/interface stack without exporting the new lawful-admission output schema.",
        },
        {
            "id": "f814_preserves_downstream_output_schema_class",
            "pass": (
                f814.get("status")
                == "F814_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SOURCE_BINDING_SELECTION_OR_PREFERENCE_RULE_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and has_required_field(f814_target, "selected_source_binding_output_schema")
            ),
            "details": "F814 preserves output-schema class on the downstream schema lane, again without discharging the new lawful-admission output schema.",
        },
        {
            "id": "p821_already_keeps_lawful_admission_unexported",
            "pass": (
                p821.get("status")
                == "P821_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_SCHEMA_DOMAIN_ADMISSION_BLOCKED"
                and p821_theorem.get("lawful_schema_domain_admission_exported_now") is False
            ),
            "details": "P821 already keeps the lawful-admission object itself unexported on the current repo state.",
        },
        {
            "id": "no_exact_new_lane_output_schema_export_present",
            "pass": (
                f821_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_target_v1"
                and not any(
                    target.get("object_id")
                    == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_v1"
                    for target in [f820_target, f819_target, f818_target, f814_target]
                )
            ),
            "details": "No current export names one exact lawful_schema_domain_admission_output_schema object for the new T213/T216 -> alpha_s schema lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "output_schema_class_candidate_supported_now": (
            checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"]
        ),
        "lawful_schema_domain_admission_output_schema_exported_now": False if all_pass else None,
        "sharp_blocker_field": "lawful_schema_domain_admission_output_schema" if all_pass else None,
        "next_honest_move_is_freeze_exact_output_schema_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P822_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_SCHEMA_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        if all_pass
        else "P822_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P822",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p821_lawful_admission_probe": rel(IN_P821),
            "f821_lawful_admission_target_packet": rel(IN_F821),
            "f820_combined_boundary_target_packet": rel(IN_F820),
            "f819_upstream_rule_target_packet": rel(IN_F819),
            "f818_upstream_interface_target_packet": rel(IN_F818),
            "f814_downstream_schema_target_packet": rel(IN_F814),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "output_schema_class_support_stack": {
            "candidate_support_refs": [
                "alpha_s_pair12_source_side_branch_selection_provider_schema_domain_admission_or_nonidentification_boundary_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_v1",
                "alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_v1",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These objects preserve output-schema class only; none exports the exact lawful schema-domain admission output schema required by F821.",
        },
        "current_honest_reading": [
            "The repo already preserves output-schema class on nearby upstream and downstream targets around the new lane.",
            "But no current export supplies the exact lawful_schema_domain_admission_output_schema required by F821 for the new T213/T216 -> alpha_s schema route.",
            "So the sharp blocker is now the still-missing exact output-schema object itself.",
        ],
        "recommended_next_packet": {
            "id": "F822_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET",
            "goal": "Freeze the exact lawful schema-domain admission output-schema object still missing after output-schema class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_schema_domain_admission_target_ref",
                "output_schema_class_candidate_support_refs",
                "upstream_rule_or_interface_output_refs",
                "downstream_schema_output_ref",
                "exact_output_schema_statement",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P822",
        "status": status,
        "as_of": AS_OF,
        "output_schema_class_candidate_supported_now": theorem_result[
            "output_schema_class_candidate_supported_now"
        ],
        "lawful_schema_domain_admission_output_schema_exported_now": theorem_result[
            "lawful_schema_domain_admission_output_schema_exported_now"
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
