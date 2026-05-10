#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P877 = GENERATED / "p877_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_blocked.json"
IN_F877 = GENERATED / "f877_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F875 = GENERATED / "f875_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F874 = GENERATED / "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F873 = GENERATED / "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json"
IN_F866 = GENERATED / "f866_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p878_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p878_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked_summary.json"


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

    prereq = [IN_P877, IN_F877, IN_F875, IN_F874, IN_F873, IN_F866, IN_F869]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P878",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p877 = load_json(IN_P877)
    f877 = load_json(IN_F877)
    f875 = load_json(IN_F875)
    f874 = load_json(IN_F874)
    f873 = load_json(IN_F873)
    f866 = load_json(IN_F866)
    f869 = load_json(IN_F869)

    p877_theorem = p877.get("theorem_result") or {}
    f877_target = f877.get("target_object") or {}
    f875_target = f875.get("target_object") or {}
    f874_target = f874.get("target_object") or {}
    f873_target = f873.get("target_object") or {}
    f866_target = f866.get("target_object") or {}
    f869_target = f869.get("target_object") or {}
    f877_refs = f877.get("target_refs") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_output_schema_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p877_", "f877_", "p878_", "f878_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f877_already_names_exact_missing_statement_field",
            "pass": (
                f877.get("status")
                == "F877_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f877_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and has_required_field(f877_target, "exact_output_schema_statement")
            ),
            "details": "F877 already isolates exact_output_schema_statement as one exact missing field of the lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission output-schema target.",
        },
        {
            "id": "nearby_targets_preserve_statement_slots",
            "pass": (
                has_required_field(f875_target, "boundary_output_schema")
                and has_required_field(f874_target, "selected_interface_output_schema")
                and has_required_field(f873_target, "exact_interface_output_schema")
                and has_required_field(f866_target, "exact_output_schema_statement")
                and "exact_output_schema_statement"
                in (f877_refs.get("neighboring_output_schema_or_statement_refs") or [])
            ),
            "details": "Nearby upstream and downstream targets do preserve output-schema statement slots on neighboring lanes.",
        },
        {
            "id": "p877_already_keeps_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_unexported",
            "pass": (
                p877.get("status")
                == "P877_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
                and p877_theorem.get(
                    "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_exported_now"
                )
                is False
            ),
            "details": "P877 already keeps the lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission output-schema object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_new_exact_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission output-schema statement outside the new frozen F877 target lineage itself.",
        },
        {
            "id": "neighboring_statement_slots_remain_nonidentical_to_new_lane_statement",
            "pass": (
                f877_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and "boundary_output_schema"
                in (f877_refs.get("upstream_rule_or_interface_output_refs") or [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1"
                in (f877_refs.get("neighboring_output_schema_or_statement_refs") or [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_v1"
                in (f877_refs.get("neighboring_output_schema_or_statement_refs") or [])
                and "exact_output_schema_statement"
                in (f877_refs.get("neighboring_output_schema_or_statement_refs") or [])
            ),
            "details": "F877 already records neighboring statement slots only as nonidentical support refs, not as silent discharge of the new lawful refined shift-interface-rule lane.",
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
        "P878_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        if all_pass
        else "P878_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P878",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p877_output_schema_probe": rel(IN_P877),
            "f877_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet": rel(IN_F877),
            "f875_boundary_target_packet": rel(IN_F875),
            "f874_rule_target_packet": rel(IN_F874),
            "f873_interface_target_packet": rel(IN_F873),
            "f866_neighboring_output_target_packet": rel(IN_F866),
            "f869_neighboring_exact_required_form_statement_target_packet": rel(IN_F869),
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
            "scope_limit": "These are neighboring statement slots only; none exports the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission output-schema statement required by F877.",
        },
        "current_honest_reading": [
            "The repo preserves statement-level output slots around the lawful refined shift-interface-rule domain-admission exact-required-form-statement lane, but only as neighboring target fields or neighboring target refs.",
            "No current export yet names the exact output-schema statement required by F877 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F878_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission output-schema statement object still missing after statement-class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_ref",
                "statement_class_candidate_support_refs",
                "neighboring_statement_slot_refs",
                "exact_statement_required_form_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P878",
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
