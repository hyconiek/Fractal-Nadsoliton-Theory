#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P824 = GENERATED / "p824_current_strict_alpha_s_statement_form_class_candidate_supported_exact_statement_required_form_blocked.json"
IN_F824 = GENERATED / "f824_current_strict_alpha_s_exact_statement_required_form_target_packet.json"
IN_F823 = GENERATED / "f823_current_strict_alpha_s_exact_output_schema_statement_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
IN_F821 = GENERATED / "f821_current_strict_alpha_s_lawful_schema_domain_admission_target_packet.json"
IN_F820 = GENERATED / "f820_current_strict_alpha_s_schema_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F818 = GENERATED / "f818_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_packet.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
IN_F812 = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json"
IN_F789 = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet.json"

OUT = GENERATED / "p825_current_strict_alpha_s_required_form_statement_class_candidate_supported_exact_required_form_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p825_current_strict_alpha_s_required_form_statement_class_candidate_supported_exact_required_form_statement_blocked_summary.json"


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

    prereq = [
        IN_P824,
        IN_F824,
        IN_F823,
        IN_F822,
        IN_F821,
        IN_F820,
        IN_F819,
        IN_F818,
        IN_F814,
        IN_F812,
        IN_F789,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P825",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p824 = load_json(IN_P824)
    f824 = load_json(IN_F824)
    f823 = load_json(IN_F823)
    f822 = load_json(IN_F822)
    f821 = load_json(IN_F821)
    f820 = load_json(IN_F820)
    f819 = load_json(IN_F819)
    f818 = load_json(IN_F818)
    f814 = load_json(IN_F814)
    f812 = load_json(IN_F812)
    f789 = load_json(IN_F789)

    p824_theorem = p824.get("theorem_result") or {}
    f824_target = f824.get("target_object") or {}
    f823_target = f823.get("target_object") or {}
    f822_target = f822.get("target_object") or {}
    f821_target = f821.get("target_object") or {}
    f820_target = f820.get("target_object") or {}
    f819_target = f819.get("target_object") or {}
    f818_target = f818.get("target_object") or {}
    f814_target = f814.get("target_object") or {}
    f812_target = f812.get("target_object") or {}
    f789_target = f789.get("target_interface") or {}

    token_hits: list[str] = []
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p824_", "f824_", "p825_", "f825_")):
            continue
        text = path.read_text(encoding="utf-8")
        if "exact_required_form_statement_ref" in text:
            token_hits.append(rel(path))

    checks = [
        {
            "id": "f824_already_names_exact_missing_statement_field",
            "pass": (
                f824.get("status")
                == "F824_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
                and f824_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_statement_required_form_target_v1"
                and has_required_field(f824_target, "exact_required_form_statement_ref")
            ),
            "details": "F824 already isolates exact_required_form_statement_ref as one exact missing field of the required-form target.",
        },
        {
            "id": "neighboring_lanes_preserve_only_statement_or_form_class_support",
            "pass": (
                has_required_field(f823_target, "exact_statement_required_form_ref")
                and has_required_field(f822_target, "exact_output_schema_statement")
                and has_required_field(f821_target, "lawful_schema_domain_admission_output_schema")
                and has_required_field(f820_target, "boundary_output_schema")
                and has_required_field(f819_target, "selected_interface_output_schema")
                and has_required_field(f818_target, "exact_interface_output_schema")
                and has_required_field(f814_target, "selected_source_binding_output_schema")
                and has_required_field(f812_target, "selected_source_binding_output_schema")
                and has_required_field(f789_target, "normalization_rule_ref")
            ),
            "details": "Neighboring lanes preserve only statement or form-class structure through target fields and generic interface-form scaffolding, not through an exact required-form statement export for the new lane.",
        },
        {
            "id": "p824_already_keeps_exact_required_form_unexported",
            "pass": (
                p824.get("status")
                == "P824_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
                and p824_theorem.get("exact_statement_required_form_exported_now") is False
            ),
            "details": "P824 already keeps the exact required-form object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_exact_required_form_statement_export",
            "pass": token_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_required_form_statement_ref outside the new frozen F824 lineage itself.",
        },
        {
            "id": "neighboring_support_remains_nonidentical_to_new_lane_statement",
            "pass": (
                f824_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_exact_statement_required_form_target_v1"
                and len((f824.get("target_refs") or {}).get("statement_form_class_candidate_support_refs") or []) >= 8
            ),
            "details": "F824 already records neighboring statement/form supports only as nonidentical candidate support, not as silent discharge of the new lane.",
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
        "P825_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        if all_pass
        else "P825_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P825",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p824_required_form_probe": rel(IN_P824),
            "f824_required_form_target_packet": rel(IN_F824),
            "f823_exact_statement_target_packet": rel(IN_F823),
            "f822_output_schema_target_packet": rel(IN_F822),
            "f821_lawful_admission_target_packet": rel(IN_F821),
            "f820_boundary_target_packet": rel(IN_F820),
            "f819_rule_target_packet": rel(IN_F819),
            "f818_interface_target_packet": rel(IN_F818),
            "f814_downstream_schema_target_packet": rel(IN_F814),
            "f812_downstream_rule_target_packet": rel(IN_F812),
            "f789_generic_interface_target_packet": rel(IN_F789),
        },
        "repo_scan_token_hits_for_exact_required_form_statement_ref": token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "required_form_statement_class_support_stack": {
            "candidate_support_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_schema_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "selected_source_binding_output_schema",
                "normalization_rule_ref",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement/form-class supports only; none exports the exact required-form statement needed by F824.",
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the new lane, but only through neighboring target fields and generic interface scaffolding.",
            "No current export yet names the exact required-form statement required by F824 for the new T213/T216 -> alpha_s schema route.",
            "So the sharp blocker is now the still-missing exact required-form statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F825_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact required-form statement object still missing after required-form statement-class support is present only at candidate grade.",
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
        "stage": "P825",
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
