#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P856 = GENERATED / "p856_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_F856 = GENERATED / "f856_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F855 = GENERATED / "f855_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F854 = GENERATED / "f854_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_target_packet.json"
IN_F853 = GENERATED / "f853_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F852 = GENERATED / "f852_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F851 = GENERATED / "f851_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_interface_target_packet.json"
IN_F847 = GENERATED / "f847_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p857_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
OUT_SUMMARY = GENERATED / "p857_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked_summary.json"


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

    prereq = [IN_P856, IN_F856, IN_F855, IN_F854, IN_F853, IN_F852, IN_F851, IN_F847]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P857",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p856 = load_json(IN_P856)
    f856 = load_json(IN_F856)
    f855 = load_json(IN_F855)
    f854 = load_json(IN_F854)
    f853 = load_json(IN_F853)
    f852 = load_json(IN_F852)
    f851 = load_json(IN_F851)
    f847 = load_json(IN_F847)

    p856_theorem = p856.get("theorem_result") or {}
    f856_target = f856.get("target_object") or {}
    f856_refs = f856.get("target_refs") or {}
    f855_target = f855.get("target_object") or {}
    f854_target = f854.get("target_object") or {}
    f853_target = f853.get("target_object") or {}
    f852_target = f852.get("target_object") or {}
    f851_target = f851.get("target_object") or {}
    f847_target = f847.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p856_", "f856_", "p857_", "f857_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f856_already_names_exact_missing_form_field",
            "pass": (
                f856.get("status")
                == "F856_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f856_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and has_required_field(f856_target, "exact_statement_required_form_ref")
            ),
            "details": "F856 already isolates exact_statement_required_form_ref as one exact missing field of the refined lawful exact-statement target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_form_class_support",
            "pass": (
                has_required_field(f855_target, "exact_output_schema_statement")
                and has_required_field(
                    f854_target,
                    "lawful_refined_exact_required_form_statement_domain_admission_output_schema",
                )
                and has_required_field(f853_target, "boundary_output_schema")
                and has_required_field(f852_target, "selected_interface_output_schema")
                and has_required_field(f851_target, "exact_interface_output_schema")
                and f847_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve refined output or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact statement-required-form export for the refined new lane.",
        },
        {
            "id": "p856_already_keeps_refined_exact_statement_unexported",
            "pass": (
                p856.get("status")
                == "P856_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p856_theorem.get("exact_output_schema_statement_exported_now") is False
            ),
            "details": "P856 already keeps the exact refined new-lane statement object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_statement_required_form_ref outside the new frozen F856 lineage itself.",
        },
        {
            "id": "neighboring_form_support_remains_nonidentical_to_refined_new_lane_form",
            "pass": (
                f856_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_v1"
                == f856_refs.get("lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_target_v1"
                in (f856_refs.get("statement_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f856_refs.get("neighboring_statement_slot_refs") or [])
                and "exact_required_form_statement"
                in (f856_refs.get("neighboring_statement_slot_refs") or [])
            ),
            "details": "F856 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the refined new lane.",
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
        "P857_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        if all_pass
        else "P857_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P857",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p856_exact_statement_probe": rel(IN_P856),
            "f856_refined_exact_statement_target_packet": rel(IN_F856),
            "f855_refined_output_schema_target_packet": rel(IN_F855),
            "f854_refined_lawful_admission_target_packet": rel(IN_F854),
            "f853_boundary_target_packet": rel(IN_F853),
            "f852_rule_target_packet": rel(IN_F852),
            "f851_interface_target_packet": rel(IN_F851),
            "f847_neighboring_exact_required_form_statement_target_packet": rel(IN_F847),
        },
        "repo_scan_object_hits_for_exact_statement_required_form_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "statement_form_class_support_stack": {
            "candidate_support_refs": [
                "exact_output_schema_statement",
                "lawful_refined_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring form-class supports only; none exports the exact refined statement-required form needed by F856.",
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the refined new lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F856 for the refined lawful T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form object itself.",
        ],
        "recommended_next_packet": {
            "id": "F857_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET",
            "goal": "Freeze the exact refined statement-required form object still missing after statement-form class support is present only at candidate grade.",
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
        "stage": "P857",
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
