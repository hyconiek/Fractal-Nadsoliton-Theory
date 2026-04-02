#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P923 = GENERATED / "p923_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_statement_required_form_blocked.json"
IN_F923 = GENERATED / "f923_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_statement_required_form_target_packet.json"
IN_F922 = GENERATED / "f922_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F921 = GENERATED / "f921_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet.json"
IN_F920 = GENERATED / "f920_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_target_packet.json"
IN_F919 = GENERATED / "f919_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F918 = GENERATED / "f918_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json"
IN_F917 = GENERATED / "f917_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.json"
IN_F913 = GENERATED / "f913_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p924_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p924_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_blocked_summary.json"


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

    prereq = [IN_P923, IN_F923, IN_F922, IN_F921, IN_F920, IN_F919, IN_F918, IN_F917, IN_F913]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P924",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p923 = load_json(IN_P923)
    f923 = load_json(IN_F923)
    f922 = load_json(IN_F922)
    f921 = load_json(IN_F921)
    f920 = load_json(IN_F920)
    f919 = load_json(IN_F919)
    f918 = load_json(IN_F918)
    f917 = load_json(IN_F917)
    f913 = load_json(IN_F913)

    p923_theorem = p923.get("theorem_result") or {}
    f923_target = f923.get("target_object") or {}
    f923_refs = f923.get("target_refs") or {}
    f922_target = f922.get("target_object") or {}
    f921_target = f921.get("target_object") or {}
    f920_target = f920.get("target_object") or {}
    f919_target = f919.get("target_object") or {}
    f918_target = f918.get("target_object") or {}
    f917_target = f917.get("target_object") or {}
    f913_target = f913.get("target_object") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p923_", "f923_", "p924_", "f924_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f923_already_names_exact_missing_required_form_statement_field",
            "pass": (
                f923.get("status")
                == "F923_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f923_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f923_target, "exact_required_form_statement_ref")
            ),
            "details": "F923 already isolates exact_required_form_statement_ref as one exact missing field of the deeper lawful refined domain-admission exact-statement-required-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_required_form_statement_class_support",
            "pass": (
                has_required_field(f922_target, "exact_statement_required_form_ref")
                and has_required_field(f921_target, "exact_output_schema_statement")
                and has_required_field(
                    f920_target,
                    "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                )
                and has_required_field(f919_target, "boundary_output_schema")
                and has_required_field(f918_target, "selected_interface_output_schema")
                and has_required_field(f917_target, "exact_interface_output_schema")
                and f913_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "Neighboring lanes preserve lawful refined statement or form-like structure only through target fields and neighboring exact-required-form objects, not through an exact required-form-statement export for the new deeper lane.",
        },
        {
            "id": "p923_already_keeps_lawful_refined_deeper_exact_statement_required_form_unexported",
            "pass": (
                p923.get("status")
                == "P923_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p923_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P923 already keeps the exact deeper lawful refined new-lane statement-required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_refined_exact_required_form_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_required_form_statement_ref outside the new frozen F923 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_required_form_statement",
            "pass": (
                f923_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                == f923_refs.get("exact_output_schema_statement_target_ref")
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f923_refs.get("statement_form_class_candidate_support_refs") or [])
                and "exact_output_schema_statement"
                in (f923_refs.get("neighboring_form_slot_refs") or [])
                and "exact_required_form_statement"
                in (f923_refs.get("neighboring_form_slot_refs") or [])
            ),
            "details": "F923 already records neighboring slots and required-form-like supports only as nonidentical candidate support, not as silent discharge of the deeper lawful refined lane.",
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
        "P924_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P924_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P924",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p923_required_form_probe": rel(IN_P923),
            "f923_lawful_refined_deeper_domain_admission_exact_statement_required_form_target_packet": rel(IN_F923),
            "f922_lawful_refined_deeper_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F922),
            "f921_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F921),
            "f920_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F920),
            "f919_boundary_target_packet": rel(IN_F919),
            "f918_rule_target_packet": rel(IN_F918),
            "f917_interface_target_packet": rel(IN_F917),
            "f913_neighboring_exact_required_form_statement_target_packet": rel(IN_F913),
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
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact deeper lawful refined domain-admission required-form statement needed by F923.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the deeper lawful refined domain-admission lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F923 for the lawful refined T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form-statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F924_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
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
        "stage": "P924",
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
