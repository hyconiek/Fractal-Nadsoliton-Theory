#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P921 = GENERATED / "p921_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_blocked.json"
IN_F921 = GENERATED / "f921_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet.json"
IN_F919 = GENERATED / "f919_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F918 = GENERATED / "f918_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json"
IN_F917 = GENERATED / "f917_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.json"
IN_F910 = GENERATED / "f910_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F913 = GENERATED / "f913_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p922_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_output_schema_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p922_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_output_schema_statement_blocked_summary.json"


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

    prereq = [IN_P921, IN_F921, IN_F919, IN_F918, IN_F917, IN_F910, IN_F913]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P922",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p921 = load_json(IN_P921)
    f921 = load_json(IN_F921)
    f919 = load_json(IN_F919)
    f918 = load_json(IN_F918)
    f917 = load_json(IN_F917)
    f910 = load_json(IN_F910)
    f913 = load_json(IN_F913)

    p921_theorem = p921.get("theorem_result") or {}
    f921_target = f921.get("target_object") or {}
    f919_target = f919.get("target_object") or {}
    f918_target = f918.get("target_object") or {}
    f917_target = f917.get("target_object") or {}
    f910_target = f910.get("target_object") or {}
    f913_target = f913.get("target_object") or {}
    f921_refs = f921.get("target_refs") or {}

    object_hits: list[str] = []
    exact_object_token = (
        "alpha_s_pair12_source_side_branch_selection_provider_"
        "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_"
        "domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement"
    )
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p921_", "f921_", "p922_", "f922_")):
            continue
        text = path.read_text(encoding="utf-8")
        if exact_object_token in text:
            object_hits.append(rel(path))

    checks = [
        {
            "id": "f921_already_names_exact_missing_statement_field",
            "pass": (
                f921.get("status")
                == "F921_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f921_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and has_required_field(f921_target, "exact_output_schema_statement")
            ),
            "details": "F921 already isolates exact_output_schema_statement as one exact missing field of the deeper lawful refined domain-admission output-schema target.",
        },
        {
            "id": "nearby_targets_preserve_statement_slots",
            "pass": (
                has_required_field(f919_target, "boundary_output_schema")
                and has_required_field(f918_target, "selected_interface_output_schema")
                and has_required_field(f917_target, "exact_interface_output_schema")
                and has_required_field(f910_target, "exact_output_schema_statement")
                and "exact_output_schema_statement"
                in (f921_refs.get("neighboring_output_schema_or_statement_refs") or [])
            ),
            "details": "Nearby upstream and downstream targets do preserve output-schema statement slots on neighboring lanes.",
        },
        {
            "id": "p921_already_keeps_lawful_refined_deeper_domain_admission_output_schema_unexported",
            "pass": (
                p921.get("status")
                == "P921_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
                and p921_theorem.get(
                    "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema_exported_now"
                )
                is False
            ),
            "details": "P921 already keeps the deeper lawful refined domain-admission output-schema object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_new_exact_statement_export",
            "pass": object_hits == [],
            "details": "Repo scan finds no generated artifact exporting the exact deeper lawful refined domain-admission output-schema statement outside the new frozen F921 target lineage itself.",
        },
        {
            "id": "neighboring_statement_slots_remain_nonidentical_to_new_lane_statement",
            "pass": (
                f921_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and "boundary_output_schema"
                in (f921_refs.get("upstream_rule_or_interface_output_refs") or [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1"
                in (f921_refs.get("neighboring_output_schema_or_statement_refs") or [])
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
                in (f921_refs.get("neighboring_output_schema_or_statement_refs") or [])
                and "exact_output_schema_statement"
                in (f921_refs.get("neighboring_output_schema_or_statement_refs") or [])
                and f913_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "F921 already records neighboring statement slots only as nonidentical support refs, not as silent discharge of the new deeper lawful refined lane.",
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
        "P922_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        if all_pass
        else "P922_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P922",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p921_output_schema_probe": rel(IN_P921),
            "f921_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F921),
            "f919_boundary_target_packet": rel(IN_F919),
            "f918_rule_target_packet": rel(IN_F918),
            "f917_interface_target_packet": rel(IN_F917),
            "f910_neighboring_output_target_packet": rel(IN_F910),
            "f913_neighboring_exact_required_form_statement_target_packet": rel(IN_F913),
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
            "scope_limit": "These are neighboring statement slots only; none exports the exact deeper lawful refined domain-admission output-schema statement required by F921.",
        },
        "current_honest_reading": [
            "The repo preserves statement-level output slots around the deeper lawful refined domain-admission lane, but only as neighboring target fields or neighboring target refs.",
            "No current export yet names the exact output-schema statement required by F921 for the new T213/T216 -> alpha_s lawful route.",
            "So the sharp blocker is now the still-missing exact statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F922_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact deeper lawful refined domain-admission output-schema statement object still missing after statement-class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_refined_deeper_domain_admission_output_schema_target_ref",
                "statement_class_candidate_support_refs",
                "neighboring_statement_slot_refs",
                "exact_statement_required_form_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P922",
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
