#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P889 = GENERATED / "p889_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_output_schema_statement_blocked.json"
IN_F889 = GENERATED / "f889_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_output_schema_statement_target_packet.json"
IN_F888 = GENERATED / "f888_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F887 = GENERATED / "f887_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F886 = GENERATED / "f886_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F885 = GENERATED / "f885_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F884 = GENERATED / "f884_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_interface_target_packet.json"
IN_F880 = GENERATED / "f880_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p890_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_statement_required_form_blocked.json"
OUT_SUMMARY = GENERATED / "p890_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_statement_required_form_blocked_summary.json"


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

    prereq = [IN_P889, IN_F889, IN_F888, IN_F887, IN_F886, IN_F885, IN_F884, IN_F880]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P890",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p889 = load_json(IN_P889)
    f889 = load_json(IN_F889)
    f888 = load_json(IN_F888)
    f887 = load_json(IN_F887)
    f886 = load_json(IN_F886)
    f885 = load_json(IN_F885)
    f884 = load_json(IN_F884)
    f880 = load_json(IN_F880)

    p889_theorem = p889.get("theorem_result") or {}
    f889_target = f889.get("target_object") or {}
    f889_refs = f889.get("target_refs") or {}
    f888_target = f888.get("target_object") or {}
    f887_target = f887.get("target_object") or {}
    f886_target = f886.get("target_object") or {}
    f885_target = f885.get("target_object") or {}
    f884_target = f884.get("target_object") or {}
    f880_target = f880.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p889_", "f889_", "p890_", "f890_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f889_already_names_exact_missing_form_field",
            "pass": (
                f889.get("status")
                == "F889_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f889_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and has_required_field(f889_target, "exact_statement_required_form_ref")
            ),
            "details": "F889 already isolates exact_statement_required_form_ref as one exact missing field of the deeper lawful refined domain-admission exact-output-schema statement target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_form_class_support",
            "pass": (
                has_required_field(f888_target, "exact_output_schema_statement")
                and has_required_field(
                    f887_target,
                    "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                )
                and has_required_field(f886_target, "boundary_output_schema")
                and has_required_field(f885_target, "selected_interface_output_schema")
                and has_required_field(f884_target, "exact_interface_output_schema")
                and f880_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve lawful refined output or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact statement-required-form export for the new deeper lane.",
        },
        {
            "id": "p889_already_keeps_lawful_refined_deeper_exact_statement_unexported",
            "pass": (
                p889.get("status")
                == "P889_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
                and p889_theorem.get("exact_output_schema_statement_exported_now") is False
            ),
            "details": "P889 already keeps the exact deeper lawful refined new-lane statement object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_statement_required_form_ref outside the new frozen F889 lineage itself.",
        },
        {
            "id": "neighboring_form_support_remains_nonidentical_to_new_lane_form",
            "pass": (
                f889_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                == f889_refs.get("lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f889_refs.get("statement_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f889_refs.get("neighboring_statement_slot_refs") or [])
                and "exact_required_form_statement"
                in (f889_refs.get("neighboring_statement_slot_refs") or [])
            ),
            "details": "F889 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the deeper lawful refined lane.",
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
        "P890_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        if all_pass
        else "P890_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P890",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p889_exact_statement_probe": rel(IN_P889),
            "f889_lawful_refined_deeper_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F889),
            "f888_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F888),
            "f887_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F887),
            "f886_boundary_target_packet": rel(IN_F886),
            "f885_rule_target_packet": rel(IN_F885),
            "f884_interface_target_packet": rel(IN_F884),
            "f880_neighboring_exact_required_form_statement_target_packet": rel(IN_F880),
        },
        "repo_scan_object_hits_for_exact_statement_required_form_ref": object_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "statement_form_class_support_stack": {
            "candidate_support_refs": [
                "exact_output_schema_statement",
                "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring form-class supports only; none exports the exact deeper lawful refined domain-admission statement-required form needed by F889.",
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the deeper lawful refined domain-admission lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F889 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form object itself.",
        ],
        "recommended_next_packet": {
            "id": "F890_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET",
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
        "stage": "P890",
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
