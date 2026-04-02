#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P912 = GENERATED / "p912_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F912 = GENERATED / "f912_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F911 = GENERATED / "f911_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F910 = GENERATED / "f910_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F909 = GENERATED / "f909_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F908 = GENERATED / "f908_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F907 = GENERATED / "f907_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_adapter_or_carrier_rule_target_packet.json"
IN_F906 = GENERATED / "f906_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_interface_target_packet.json"
IN_F902 = GENERATED / "f902_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p913_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
# Technical shortening only: the full summary basename would exceed filesystem limits.
OUT_SUMMARY = GENERATED / "p913_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_blocked_summary.json"


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

    prereq = [IN_P912, IN_F912, IN_F911, IN_F910, IN_F909, IN_F908, IN_F907, IN_F906, IN_F902]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P913",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p912 = load_json(IN_P912)
    f912 = load_json(IN_F912)
    f911 = load_json(IN_F911)
    f910 = load_json(IN_F910)
    f909 = load_json(IN_F909)
    f908 = load_json(IN_F908)
    f907 = load_json(IN_F907)
    f906 = load_json(IN_F906)
    f902 = load_json(IN_F902)

    p912_theorem = p912.get("theorem_result") or {}
    f912_target = f912.get("target_object") or {}
    f912_refs = f912.get("target_refs") or {}
    f911_target = f911.get("target_object") or {}
    f910_target = f910.get("target_object") or {}
    f909_target = f909.get("target_object") or {}
    f908_target = f908.get("target_object") or {}
    f907_target = f907.get("target_object") or {}
    f906_target = f906.get("target_object") or {}
    f902_target = f902.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p912_", "f912_", "p913_", "f913_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f912_already_names_exact_missing_required_form_statement_field",
            "pass": (
                f912.get("status")
                == "F912_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f912_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f912_target, "exact_required_form_statement_ref")
            ),
            "details": "F912 already isolates exact_required_form_statement_ref as one exact missing field of the deeper lawful refined domain-admission exact-statement-required-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_required_form_statement_class_support",
            "pass": (
                has_required_field(f911_target, "exact_statement_required_form_ref")
                and has_required_field(f910_target, "exact_output_schema_statement")
                and has_required_field(
                    f909_target,
                    "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                )
                and has_required_field(f908_target, "boundary_output_schema")
                and has_required_field(f907_target, "selected_interface_output_schema")
                and has_required_field(f906_target, "exact_interface_output_schema")
                and f902_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve lawful refined statement or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact required-form-statement export for the new deeper lane.",
        },
        {
            "id": "p912_already_keeps_lawful_refined_deeper_exact_statement_required_form_unexported",
            "pass": (
                p912.get("status")
                == "P912_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p912_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P912 already keeps the exact deeper lawful refined new-lane statement-required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_required_form_statement_ref outside the new frozen F912 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_required_form_statement",
            "pass": (
                f912_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                == f912_refs.get("exact_output_schema_statement_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f912_refs.get("statement_form_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f912_refs.get("neighboring_form_slot_refs") or [])
                and "exact_required_form_statement"
                in (f912_refs.get("neighboring_form_slot_refs") or [])
            ),
            "details": "F912 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the deeper lawful refined lane.",
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
        "P913_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P913_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P913",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p912_required_form_probe": rel(IN_P912),
            "f912_lawful_refined_deeper_domain_admission_exact_statement_required_form_target_packet": rel(IN_F912),
            "f911_lawful_refined_deeper_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F911),
            "f910_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F910),
            "f909_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F909),
            "f908_boundary_target_packet": rel(IN_F908),
            "f907_rule_target_packet": rel(IN_F907),
            "f906_interface_target_packet": rel(IN_F906),
            "f902_neighboring_exact_required_form_statement_target_packet": rel(IN_F902),
        },
        "repo_scan_object_hits_for_exact_required_form_statement_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "required_form_statement_class_support_stack": {
            "candidate_support_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact deeper lawful refined domain-admission required-form statement needed by F912.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the deeper lawful refined domain-admission lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F912 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form-statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F913_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact deeper lawful refined domain-admission required-form statement object still missing after required-form-statement class support is present only at candidate grade.",
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
        "stage": "P913",
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
