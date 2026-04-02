#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P911 = GENERATED / "p911_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_F911 = GENERATED / "f911_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F910 = GENERATED / "f910_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F909 = GENERATED / "f909_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F908 = GENERATED / "f908_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F907 = GENERATED / "f907_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_adapter_or_carrier_rule_target_packet.json"
IN_F906 = GENERATED / "f906_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_interface_target_packet.json"
IN_F902 = GENERATED / "f902_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p912_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
OUT_SUMMARY = GENERATED / "p912_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked_summary.json"


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

    prereq = [IN_P911, IN_F911, IN_F910, IN_F909, IN_F908, IN_F907, IN_F906, IN_F902]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P912",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p911 = load_json(IN_P911)
    f911 = load_json(IN_F911)
    f910 = load_json(IN_F910)
    f909 = load_json(IN_F909)
    f908 = load_json(IN_F908)
    f907 = load_json(IN_F907)
    f906 = load_json(IN_F906)
    f902 = load_json(IN_F902)

    p911_theorem = p911.get("theorem_result") or {}
    f911_target = f911.get("target_object") or {}
    f911_refs = f911.get("target_refs") or {}
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
        "domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p911_", "f911_", "p912_", "f912_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f911_already_names_exact_missing_form_field",
            "pass": (
                f911.get("status")
                == "F911_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f911_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and has_required_field(f911_target, "exact_statement_required_form_ref")
            ),
            "details": "F911 already isolates exact_statement_required_form_ref as one exact missing field of the deeper lawful refined domain-admission exact-output-schema-statement target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_form_class_support",
            "pass": (
                has_required_field(f910_target, "exact_output_schema_statement")
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
            "details": "Neighboring lanes preserve lawful refined output or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact statement-required-form export for the new deeper lane.",
        },
        {
            "id": "p911_already_keeps_lawful_refined_deeper_exact_statement_unexported",
            "pass": (
                p911.get("status")
                == "P911_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p911_theorem.get("exact_output_schema_statement_exported_now") is False
            ),
            "details": "P911 already keeps the exact deeper lawful refined new-lane statement object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_statement_required_form_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_statement_required_form_ref outside the new frozen F911 lineage itself.",
        },
        {
            "id": "neighboring_form_support_remains_nonidentical_to_new_lane_form",
            "pass": (
                f911_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                == f911_refs.get("lawful_refined_deeper_domain_admission_output_schema_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f911_refs.get("statement_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f911_refs.get("neighboring_statement_slot_refs") or [])
                and "exact_required_form_statement"
                in (f911_refs.get("neighboring_statement_slot_refs") or [])
            ),
            "details": "F911 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the deeper lawful refined lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "statement_form_class_candidate_supported_now": (
            checks[1]["pass"] and checks[4]["pass"]
        ),
        "exact_statement_required_form_exported_now": False if all_pass else None,
        "sharp_blocker_field": "exact_statement_required_form_ref" if all_pass else None,
        "next_honest_move_is_freeze_exact_required_form_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P912_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        if all_pass
        else "P912_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P912",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p911_exact_statement_probe": rel(IN_P911),
            "f911_lawful_refined_deeper_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F911),
            "f910_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F910),
            "f909_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F909),
            "f908_boundary_target_packet": rel(IN_F908),
            "f907_rule_target_packet": rel(IN_F907),
            "f906_interface_target_packet": rel(IN_F906),
            "f902_neighboring_exact_required_form_statement_target_packet": rel(IN_F902),
        },
        "repo_scan_object_hits_for_exact_statement_required_form_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "statement_form_class_support_stack": {
            "candidate_support_refs": [
                "exact_output_schema_statement",
                "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring form-class supports only; none exports the exact deeper lawful refined domain-admission statement-required form needed by F911.",
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the deeper lawful refined domain-admission lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F911 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form object itself.",
        ],
        "recommended_next_packet": {
            "id": "F912_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET",
            "goal": "Freeze the exact deeper lawful refined domain-admission statement-required form object still missing after statement-form class support is present only at candidate grade.",
            "minimum_fields": [
                "exact_output_schema_statement_target_ref",
                "statement_form_class_candidate_support_refs",
                "neighboring_form_slot_refs",
                "exact_required_form_statement_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P912",
        "status": status,
        "as_of": AS_OF,
        "statement_form_class_candidate_supported_now": theorem_result[
            "statement_form_class_candidate_supported_now"
        ],
        "exact_statement_required_form_exported_now": theorem_result[
            "exact_statement_required_form_exported_now"
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
