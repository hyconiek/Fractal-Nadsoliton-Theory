#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P835 = GENERATED / "p835_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F835 = GENERATED / "f835_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F834 = GENERATED / "f834_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F833 = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F832 = GENERATED / "f832_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_target_packet.json"
IN_F831 = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F829 = GENERATED / "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"

OUT = GENERATED / "p836_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p836_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked_summary.json"


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
        IN_P835,
        IN_F835,
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
            "stage": "P836",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p835 = load_json(IN_P835)
    f835 = load_json(IN_F835)
    f834 = load_json(IN_F834)
    f833 = load_json(IN_F833)
    f832 = load_json(IN_F832)
    f831 = load_json(IN_F831)
    f830 = load_json(IN_F830)
    f829 = load_json(IN_F829)
    f825 = load_json(IN_F825)
    f822 = load_json(IN_F822)

    p835_theorem = p835.get("theorem_result") or {}
    f835_target = f835.get("target_object") or {}
    f834_target = f834.get("target_object") or {}
    f833_target = f833.get("target_object") or {}
    f832_target = f832.get("target_object") or {}
    f831_target = f831.get("target_object") or {}
    f830_target = f830.get("target_object") or {}
    f829_target = f829.get("target_object") or {}
    f825_target = f825.get("target_object") or {}
    f822_target = f822.get("target_object") or {}

    support_refs = (f835.get("target_refs") or {}).get("statement_form_class_candidate_support_refs") or []
    slot_refs = (f835.get("target_refs") or {}).get("neighboring_form_slot_refs") or []

    token_hits: list[str] = []
    lane_token = "lawful_exact_required_form_statement_domain_admission_exact_required_form_statement"
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p835_", "f835_", "p836_", "f836_")):
            continue
        text = path.read_text(encoding="utf-8")
        if lane_token in text:
            token_hits.append(rel(path))

    checks = [
        {
            "id": "f835_already_names_exact_missing_statement_field",
            "pass": (
                f835.get("status")
                == "F835_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f835_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f835_target, "exact_required_form_statement_ref")
            ),
            "details": "F835 already isolates exact_required_form_statement_ref as one exact missing field of the current required-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_statement_or_form_class_support",
            "pass": (
                f834_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and f833_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and f832_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_target_v1"
                and f831_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_v1"
                and f830_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_v1"
                and f829_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_v1"
                and f825_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                and f822_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1"
                and len(support_refs) >= 7
            ),
            "details": "Neighboring lanes preserve only statement or form-class support through current-lane neighboring targets and old-lane neighboring targets, not through an exact required-form statement export for the new lane.",
        },
        {
            "id": "p835_already_keeps_exact_statement_required_form_unexported",
            "pass": (
                p835.get("status")
                == "P835_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p835_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P835 already keeps the exact statement-required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_current_lane_exact_required_form_statement_export",
            "pass": token_hits == [],
            "details": "Repo scan finds no generated artifact exporting a lawful exact-required-form-statement domain-admission exact required-form statement outside the new frozen F835 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_statement",
            "pass": (
                len(support_refs) >= 7
                and len(slot_refs) >= 6
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_required_form_statement_target_v1"
                in support_refs
                and "exact_output_schema_statement" in slot_refs
                and "lawful_exact_required_form_statement_domain_admission_output_schema" in slot_refs
            ),
            "details": "F835 already records neighboring statement/form supports only as nonidentical candidate support, not as silent discharge of the new lane.",
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
        "P836_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P836_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P836",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p835_required_form_probe": rel(IN_P835),
            "f835_required_form_target_packet": rel(IN_F835),
            "f834_exact_statement_target_packet": rel(IN_F834),
            "f833_output_schema_target_packet": rel(IN_F833),
            "f832_lawful_admission_target_packet": rel(IN_F832),
            "f831_boundary_target_packet": rel(IN_F831),
            "f830_rule_target_packet": rel(IN_F830),
            "f829_interface_target_packet": rel(IN_F829),
            "f825_neighboring_exact_required_form_target_packet": rel(IN_F825),
            "f822_neighboring_schema_output_target_packet": rel(IN_F822),
        },
        "repo_scan_token_hits_for_current_lane_exact_required_form_statement_ref": token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "required_form_statement_class_support_stack": {
            "candidate_support_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "lawful_schema_domain_admission_exact_required_form_statement_target",
                "lawful_schema_domain_admission_output_schema",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact required-form statement needed by F835.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the new lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F835 for the new T213/T216 -> alpha_s exact-required-form-statement route.",
            "So the sharp blocker is now the still-missing exact required-form statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F836_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact required-form statement object still missing after required-form-statement class support is present only at candidate grade on the current lawful exact-required-form-statement domain-admission lane.",
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
        "stage": "P836",
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
