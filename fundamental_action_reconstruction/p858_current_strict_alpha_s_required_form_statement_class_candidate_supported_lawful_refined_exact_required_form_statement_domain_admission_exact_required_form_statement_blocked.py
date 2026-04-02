#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P857 = GENERATED / "p857_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F857 = GENERATED / "f857_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F856 = GENERATED / "f856_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F855 = GENERATED / "f855_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F854 = GENERATED / "f854_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_target_packet.json"
IN_F853 = GENERATED / "f853_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F852 = GENERATED / "f852_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F851 = GENERATED / "f851_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_interface_target_packet.json"
IN_F847 = GENERATED / "f847_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "p858_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p858_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked_summary.json"


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

    prereq = [IN_P857, IN_F857, IN_F856, IN_F855, IN_F854, IN_F853, IN_F852, IN_F851, IN_F847]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P858",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p857 = load_json(IN_P857)
    f857 = load_json(IN_F857)
    f856 = load_json(IN_F856)
    f855 = load_json(IN_F855)
    f854 = load_json(IN_F854)
    f853 = load_json(IN_F853)
    f852 = load_json(IN_F852)
    f851 = load_json(IN_F851)
    f847 = load_json(IN_F847)

    p857_theorem = p857.get("theorem_result") or {}
    f857_target = f857.get("target_object") or {}
    f856_target = f856.get("target_object") or {}
    f855_target = f855.get("target_object") or {}
    f854_target = f854.get("target_object") or {}
    f853_target = f853.get("target_object") or {}
    f852_target = f852.get("target_object") or {}
    f851_target = f851.get("target_object") or {}
    f847_target = f847.get("target_object") or {}

    support_refs = (f857.get("target_refs") or {}).get("statement_form_class_candidate_support_refs") or []
    slot_refs = (f857.get("target_refs") or {}).get("neighboring_form_slot_refs") or []

    token_hits: list[str] = []
    lane_token = "lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement"
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p857_", "f857_", "p858_", "f858_")):
            continue
        text = path.read_text(encoding="utf-8")
        if lane_token in text:
            token_hits.append(rel(path))

    checks = [
        {
            "id": "f857_already_names_exact_missing_statement_field",
            "pass": (
                f857.get("status")
                == "F857_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f857_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f857_target, "exact_required_form_statement_ref")
            ),
            "details": "F857 already isolates exact_required_form_statement_ref as one exact missing field of the refined current required-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_statement_or_form_class_support",
            "pass": (
                f856_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1"
                and f855_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_v1"
                and f854_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_target_v1"
                and f853_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_v1"
                and f852_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_shift_interface_adapter_or_carrier_identification_rule_target_v1"
                and f851_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_shift_interface_target_v1"
                and f847_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_target_v1"
                and len(support_refs) >= 6
            ),
            "details": "Neighboring lanes preserve only statement or form-class support through refined-lane neighboring targets and old-lane neighboring targets, not through an exact required-form statement export for the refined lane.",
        },
        {
            "id": "p857_already_keeps_exact_statement_required_form_unexported",
            "pass": (
                p857.get("status")
                == "P857_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p857_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P857 already keeps the exact refined statement-required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_current_lane_exact_required_form_statement_export",
            "pass": token_hits == [],
            "details": "Repo scan finds no generated artifact exporting a refined lawful exact-required-form-statement domain-admission exact required-form statement outside the new frozen F857 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_statement",
            "pass": (
                len(support_refs) >= 6
                and len(slot_refs) >= 6
                and "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_required_form_statement_target_v1"
                in support_refs
                and "exact_output_schema_statement" in slot_refs
                and "lawful_refined_exact_required_form_statement_domain_admission_output_schema" in slot_refs
                and "exact_required_form_statement" in slot_refs
            ),
            "details": "F857 already records neighboring statement/form supports only as nonidentical candidate support, not as silent discharge of the refined new lane.",
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
        "P858_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P858_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P858",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p857_required_form_probe": rel(IN_P857),
            "f857_required_form_target_packet": rel(IN_F857),
            "f856_refined_exact_statement_target_packet": rel(IN_F856),
            "f855_refined_output_schema_target_packet": rel(IN_F855),
            "f854_refined_lawful_admission_target_packet": rel(IN_F854),
            "f853_boundary_target_packet": rel(IN_F853),
            "f852_rule_target_packet": rel(IN_F852),
            "f851_interface_target_packet": rel(IN_F851),
            "f847_neighboring_exact_required_form_statement_target_packet": rel(IN_F847),
        },
        "repo_scan_token_hits_for_current_lane_exact_required_form_statement_ref": token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "required_form_statement_class_support_stack": {
            "candidate_support_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_refined_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact refined required-form statement needed by F857.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the refined new lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F857 for the refined lawful T213/T216 -> alpha_s route.",
            "So the sharp blocker is now the still-missing exact required-form statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F858_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact refined required-form statement object still missing after required-form-statement class support is present only at candidate grade on the refined lawful exact-required-form-statement domain-admission lane.",
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
        "stage": "P858",
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
