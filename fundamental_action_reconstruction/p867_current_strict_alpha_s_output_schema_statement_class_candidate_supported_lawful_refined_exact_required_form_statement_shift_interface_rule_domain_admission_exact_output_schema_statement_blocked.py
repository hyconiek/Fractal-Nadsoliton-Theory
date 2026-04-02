#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P866 = GENERATED / "p866_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_blocked.json"
IN_F866 = GENERATED / "f866_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_packet.json"
IN_F864 = GENERATED / "f864_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F863 = GENERATED / "f863_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F862 = GENERATED / "f862_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_interface_target_packet.json"
IN_F855 = GENERATED / "f855_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F858 = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p867_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p867_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_blocked_summary.json"


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

    prereq = [IN_P866, IN_F866, IN_F864, IN_F863, IN_F862, IN_F855, IN_F858]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P867",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p866 = load_json(IN_P866)
    f866 = load_json(IN_F866)
    f864 = load_json(IN_F864)
    f863 = load_json(IN_F863)
    f862 = load_json(IN_F862)
    f855 = load_json(IN_F855)
    f858 = load_json(IN_F858)

    p866_theorem = p866.get("theorem_result") or {}
    f866_target = f866.get("target_object") or {}
    f864_target = f864.get("target_object") or {}
    f863_target = f863.get("target_object") or {}
    f862_target = f862.get("target_object") or {}
    f855_target = f855.get("target_object") or {}
    f858_target = f858.get("target_object") or {}
    f866_refs = f866.get("target_refs") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p866_", "f866_", "p867_", "f867_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f866_already_names_exact_missing_statement_field",
            "pass": (
                f866.get("status")
                == "F866_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f866_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1"
                and has_required_field(f866_target, "exact_output_schema_statement")
            ),
            "details": "F866 already isolates exact_output_schema_statement as one exact missing field of the lawful refined shift-interface-rule output-schema target.",
        },
        {
            "id": "nearby_targets_preserve_statement_slots",
            "pass": (
                has_required_field(f864_target, "boundary_output_schema")
                and has_required_field(f863_target, "selected_interface_output_schema")
                and has_required_field(f862_target, "exact_interface_output_schema")
                and has_required_field(f855_target, "exact_output_schema_statement")
                and "exact_output_schema_statement"
                in (f866_refs.get("neighboring_output_schema_or_statement_refs") or [])
            ),
            "details": "Nearby upstream and downstream targets do preserve output-schema statement slots on neighboring lanes.",
        },
        {
            "id": "p866_already_keeps_lawful_refined_shift_interface_rule_output_schema_unexported",
            "pass": (
                p866.get("status")
                == "P866_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
                and p866_theorem.get(
                    "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_exported_now"
                )
                is False
            ),
            "details": "P866 already keeps the lawful refined shift-interface-rule output-schema object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_new_exact_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting the exact lawful refined shift-interface-rule output-schema statement outside the new frozen F866 target lineage itself.",
        },
        {
            "id": "neighboring_statement_slots_remain_nonidentical_to_new_lane_statement",
            "pass": (
                f866_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1"
                and "boundary_output_schema"
                in (f866_refs.get("upstream_rule_or_interface_output_refs") or [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_v1"
                in (f866_refs.get("neighboring_output_schema_or_statement_refs") or [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f866_refs.get("neighboring_output_schema_or_statement_refs") or [])
                and "exact_output_schema_statement"
                in (f866_refs.get("neighboring_output_schema_or_statement_refs") or [])
            ),
            "details": "F866 already records neighboring statement slots only as nonidentical support refs, not as silent discharge of the new lawful refined shift-interface-rule lane.",
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
        "P867_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        if all_pass
        else "P867_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P867",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p866_output_schema_probe": rel(IN_P866),
            "f866_lawful_refined_shift_interface_rule_output_schema_target_packet": rel(IN_F866),
            "f864_boundary_target_packet": rel(IN_F864),
            "f863_rule_target_packet": rel(IN_F863),
            "f862_interface_target_packet": rel(IN_F862),
            "f855_neighboring_output_target_packet": rel(IN_F855),
            "f858_neighboring_exact_required_form_statement_target_packet": rel(IN_F858),
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
            "scope_limit": "These are neighboring statement slots only; none exports the exact lawful refined shift-interface-rule output-schema statement required by F866.",
        },
        "current_honest_reading": [
            "The repo preserves statement-level output slots around the lawful refined shift-interface-rule lane, but only as neighboring target fields or neighboring target refs.",
            "No current export yet names the exact output-schema statement required by F866 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F867_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined shift-interface-rule output-schema statement object still missing after statement-class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_ref",
                "statement_class_candidate_support_refs",
                "neighboring_statement_slot_refs",
                "exact_statement_required_form_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P867",
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
