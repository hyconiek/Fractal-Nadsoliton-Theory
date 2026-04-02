#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P868 = GENERATED / "p868_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_blocked.json"
IN_F868 = GENERATED / "f868_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_packet.json"
IN_F867 = GENERATED / "f867_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F866 = GENERATED / "f866_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_packet.json"
IN_F865 = GENERATED / "f865_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_packet.json"
IN_F864 = GENERATED / "f864_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F863 = GENERATED / "f863_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F862 = GENERATED / "f862_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_interface_target_packet.json"
IN_F858 = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p869_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p869_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_blocked_summary.json"


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

    prereq = [IN_P868, IN_F868, IN_F867, IN_F866, IN_F865, IN_F864, IN_F863, IN_F862, IN_F858]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P869",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p868 = load_json(IN_P868)
    f868 = load_json(IN_F868)
    f867 = load_json(IN_F867)
    f866 = load_json(IN_F866)
    f865 = load_json(IN_F865)
    f864 = load_json(IN_F864)
    f863 = load_json(IN_F863)
    f862 = load_json(IN_F862)
    f858 = load_json(IN_F858)

    p868_theorem = p868.get("theorem_result") or {}
    f868_target = f868.get("target_object") or {}
    f868_refs = f868.get("target_refs") or {}
    f867_target = f867.get("target_object") or {}
    f866_target = f866.get("target_object") or {}
    f865_target = f865.get("target_object") or {}
    f864_target = f864.get("target_object") or {}
    f863_target = f863.get("target_object") or {}
    f862_target = f862.get("target_object") or {}
    f858_target = f858.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p868_", "f868_", "p869_", "f869_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f868_already_names_exact_missing_required_form_statement_field",
            "pass": (
                f868.get("status")
                == "F868_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f868_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f868_target, "exact_required_form_statement_ref")
            ),
            "details": "F868 already isolates exact_required_form_statement_ref as one exact missing field of the lawful refined shift-interface-rule exact-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_required_form_statement_class_support",
            "pass": (
                has_required_field(f867_target, "exact_statement_required_form_ref")
                and has_required_field(
                    f866_target,
                    "exact_output_schema_statement",
                )
                and has_required_field(
                    f865_target,
                    "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema",
                )
                and has_required_field(f864_target, "boundary_output_schema")
                and has_required_field(f863_target, "selected_interface_output_schema")
                and has_required_field(f862_target, "exact_interface_output_schema")
                and f858_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve lawful refined statement or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact required-form-statement export for the new shift-interface-rule lane.",
        },
        {
            "id": "p868_already_keeps_lawful_refined_exact_statement_required_form_unexported",
            "pass": (
                p868.get("status")
                == "P868_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p868_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P868 already keeps the exact lawful refined new-lane statement-required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_required_form_statement_ref outside the new frozen F868 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_required_form_statement",
            "pass": (
                f868_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1"
                == f868_refs.get("exact_output_schema_statement_target_ref") or True
            ),
            "details": "Neighboring statement/form supports remain only candidate support and do not silently discharge the new lawful refined shift-interface-rule lane.",
        },
    ]

    # Strengthen the final check with explicit refs from F868 after the basic object-id check above.
    checks[-1]["pass"] = (
        f868_target.get("object_id")
        == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_v1"
        and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_target_v1"
        == f868_refs.get("exact_output_schema_statement_target_ref")
        and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
        in (f868_refs.get("statement_form_class_candidate_support_refs") or [])
        and "exact_output_schema_statement" in (f868_refs.get("neighboring_form_slot_refs") or [])
        and "exact_required_form_statement" in (f868_refs.get("neighboring_form_slot_refs") or [])
    )

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
        "P869_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P869_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P869",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p868_required_form_probe": rel(IN_P868),
            "f868_lawful_refined_shift_interface_rule_required_form_target_packet": rel(IN_F868),
            "f867_lawful_refined_shift_interface_rule_exact_statement_target_packet": rel(IN_F867),
            "f866_lawful_refined_shift_interface_rule_output_schema_target_packet": rel(IN_F866),
            "f865_lawful_refined_shift_interface_rule_admission_target_packet": rel(IN_F865),
            "f864_boundary_target_packet": rel(IN_F864),
            "f863_rule_target_packet": rel(IN_F863),
            "f862_interface_target_packet": rel(IN_F862),
            "f858_neighboring_exact_required_form_statement_target_packet": rel(IN_F858),
        },
        "repo_scan_object_hits_for_exact_required_form_statement_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "required_form_statement_class_support_stack": {
            "candidate_support_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact lawful refined shift-interface-rule required-form statement needed by F868.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the lawful refined shift-interface-rule lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F868 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F869_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined shift-interface-rule required-form statement object still missing after required-form-statement class support is present only at candidate grade.",
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
        "stage": "P869",
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
