#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P833 = GENERATED / "p833_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_output_schema_blocked.json"
IN_F833 = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F831 = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F829 = GENERATED / "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p834_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p834_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked_summary.json"


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

    prereq = [IN_P833, IN_F833, IN_F831, IN_F830, IN_F829, IN_F822, IN_F825]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P834",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p833 = load_json(IN_P833)
    f833 = load_json(IN_F833)
    f831 = load_json(IN_F831)
    f830 = load_json(IN_F830)
    f829 = load_json(IN_F829)
    f822 = load_json(IN_F822)
    f825 = load_json(IN_F825)

    p833_theorem = p833.get("theorem_result") or {}
    f833_target = f833.get("target_object") or {}
    f831_target = f831.get("target_object") or {}
    f830_target = f830.get("target_object") or {}
    f829_target = f829.get("target_object") or {}
    f822_target = f822.get("target_object") or {}
    f825_target = f825.get("target_object") or {}
    f825_refs = f825.get("target_refs") or {}

    object_hits: list[str] = []
    exact_object_token = "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement"
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p833_", "f833_", "p834_", "f834_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f833_already_names_exact_missing_statement_field",
            "pass": (
                f833.get("status")
                == "F833_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f833_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and has_required_field(f833_target, "exact_output_schema_statement")
            ),
            "details": "F833 already isolates exact_output_schema_statement as one exact missing field of the output-schema target.",
        },
        {
            "id": "nearby_targets_preserve_statement_slots",
            "pass": (
                has_required_field(f831_target, "boundary_output_schema")
                and has_required_field(f830_target, "selected_interface_output_schema")
                and has_required_field(f829_target, "exact_interface_output_schema")
                and has_required_field(f822_target, "exact_output_schema_statement")
                and "exact_output_schema_statement"
                in (f825_refs.get("neighboring_statement_or_form_slot_refs") or [])
            ),
            "details": "Nearby upstream and downstream targets do preserve output-schema statement slots on neighboring lanes.",
        },
        {
            "id": "p833_already_keeps_exact_output_schema_unexported",
            "pass": (
                p833.get("status")
                == "P833_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
                and p833_theorem.get("lawful_exact_required_form_statement_domain_admission_output_schema_exported_now")
                is False
            ),
            "details": "P833 already keeps the exact new-lane output-schema object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_exact_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting the exact new-lane output-schema statement outside the new frozen F833 target lineage itself.",
        },
        {
            "id": "neighboring_statement_slots_remain_nonidentical_to_new_lane_statement",
            "pass": (
                f833_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and "boundary_output_schema" in (f833.get("target_refs") or {}).get("upstream_rule_or_interface_output_refs", [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1"
                in (f833.get("target_refs") or {}).get("neighboring_output_schema_or_statement_refs", [])
                and "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_target_v1"
                in (f833.get("target_refs") or {}).get("neighboring_output_schema_or_statement_refs", [])
            ),
            "details": "F833 already records neighboring statement slots only as nonidentical support refs, not as silent discharge of the new lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "output_schema_statement_class_candidate_supported_now": (
            checks[1]["pass"] and checks[4]["pass"]
        ),
        "exact_output_schema_statement_exported_now": False if all_pass else None,
        "sharp_blocker_field": "exact_output_schema_statement" if all_pass else None,
        "next_honest_move_is_freeze_exact_statement_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P834_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        if all_pass
        else "P834_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P834",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p833_output_schema_probe": rel(IN_P833),
            "f833_output_schema_target_packet": rel(IN_F833),
            "f831_boundary_target_packet": rel(IN_F831),
            "f830_rule_target_packet": rel(IN_F830),
            "f829_interface_target_packet": rel(IN_F829),
            "f822_neighboring_schema_output_target_packet": rel(IN_F822),
            "f825_downstream_exact_required_form_statement_target_packet": rel(IN_F825),
        },
        "repo_scan_object_hits_for_exact_output_schema_statement": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "statement_class_support_stack": {
            "candidate_support_refs": [
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_output_schema_statement",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement slots only; none exports the exact output-schema statement required by F833.",
        },
        "current_honest_reading": [
            "The repo preserves statement-level output slots around the new lane, but only as neighboring target fields or neighboring target refs.",
            "No current export yet names the exact output-schema statement required by F833 for the new T213/T216 -> alpha_s exact-required-form-statement route.",
            "So the sharp blocker is now the still-missing exact statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F834_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact output-schema statement object still missing after statement-class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_exact_required_form_statement_domain_admission_output_schema_target_ref",
                "statement_class_candidate_support_refs",
                "neighboring_statement_slot_refs",
                "exact_statement_required_form_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P834",
        "status": status,
        "as_of": AS_OF,
        "output_schema_statement_class_candidate_supported_now": theorem_result[
            "output_schema_statement_class_candidate_supported_now"
        ],
        "exact_output_schema_statement_exported_now": theorem_result[
            "exact_output_schema_statement_exported_now"
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
