#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P901 = GENERATED / "p901_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F901 = GENERATED / "f901_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F900 = GENERATED / "f900_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F899 = GENERATED / "f899_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F898 = GENERATED / "f898_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_target_packet.json"
IN_F897 = GENERATED / "f897_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F896 = GENERATED / "f896_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F895 = GENERATED / "f895_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_interface_target_packet.json"
IN_F891 = GENERATED / "f891_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p902_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p902_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked_summary.json"


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

    prereq = [IN_P901, IN_F901, IN_F900, IN_F899, IN_F898, IN_F897, IN_F896, IN_F895, IN_F891]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P902",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p901 = load_json(IN_P901)
    f901 = load_json(IN_F901)
    f900 = load_json(IN_F900)
    f899 = load_json(IN_F899)
    f898 = load_json(IN_F898)
    f897 = load_json(IN_F897)
    f896 = load_json(IN_F896)
    f895 = load_json(IN_F895)
    f891 = load_json(IN_F891)

    p901_theorem = p901.get("theorem_result") or {}
    f901_target = f901.get("target_object") or {}
    f901_refs = f901.get("target_refs") or {}
    f900_target = f900.get("target_object") or {}
    f899_target = f899.get("target_object") or {}
    f898_target = f898.get("target_object") or {}
    f897_target = f897.get("target_object") or {}
    f896_target = f896.get("target_object") or {}
    f895_target = f895.get("target_object") or {}
    f891_target = f891.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p901_", "f901_", "p902_", "f902_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f901_already_names_exact_missing_required_form_statement_field",
            "pass": (
                f901.get("status")
                == "F901_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f901_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f901_target, "exact_required_form_statement_ref")
            ),
            "details": "F901 already isolates exact_required_form_statement_ref as one exact missing field of the deeper lawful refined domain-admission exact-statement-required-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_required_form_statement_class_support",
            "pass": (
                has_required_field(f900_target, "exact_statement_required_form_ref")
                and has_required_field(f899_target, "exact_output_schema_statement")
                and has_required_field(
                    f898_target,
                    "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                )
                and has_required_field(f897_target, "boundary_output_schema")
                and has_required_field(f896_target, "selected_interface_output_schema")
                and has_required_field(f895_target, "exact_interface_output_schema")
                and f891_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve lawful refined statement or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact required-form-statement export for the new deeper lane.",
        },
        {
            "id": "p901_already_keeps_lawful_refined_deeper_exact_statement_required_form_unexported",
            "pass": (
                p901.get("status")
                == "P901_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p901_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P901 already keeps the exact deeper lawful refined new-lane statement-required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_required_form_statement_ref outside the new frozen F901 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_required_form_statement",
            "pass": (
                f901_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement.domain_admission_exact_output_schema_statement_target_v1".replace(".", "_")
                == f901_refs.get("exact_output_schema_statement_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f901_refs.get("statement_form_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f901_refs.get("neighboring_form_slot_refs") or [])
                and "exact_required_form_statement"
                in (f901_refs.get("neighboring_form_slot_refs") or [])
            ),
            "details": "F901 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the deeper lawful refined lane.",
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
        "P902_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P902_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P902",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p901_required_form_probe": rel(IN_P901),
            "f901_lawful_refined_deeper_domain_admission_exact_statement_required_form_target_packet": rel(IN_F901),
            "f900_lawful_refined_deeper_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F900),
            "f899_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F899),
            "f898_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F898),
            "f897_boundary_target_packet": rel(IN_F897),
            "f896_rule_target_packet": rel(IN_F896),
            "f895_interface_target_packet": rel(IN_F895),
            "f891_neighboring_exact_required_form_statement_target_packet": rel(IN_F891),
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
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact deeper lawful refined domain-admission required-form statement needed by F901.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the deeper lawful refined domain-admission lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F901 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form-statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F902_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
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
        "stage": "P902",
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
