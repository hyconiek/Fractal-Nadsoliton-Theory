#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P822 = GENERATED / "p822_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_schema_domain_admission_output_schema_blocked.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
IN_F820 = GENERATED / "f820_current_strict_alpha_s_schema_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"
IN_F818 = GENERATED / "f818_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_interface_target_packet.json"
IN_F814 = GENERATED / "f814_current_strict_alpha_s_strict_source_shannon_source_binding_selection_or_preference_rule_schema_target_packet.json"
IN_F812 = GENERATED / "f812_current_strict_alpha_s_strict_source_shannon_shift_to_domain_source_binding_selection_or_preference_rule_target_packet.json"

OUT = GENERATED / "p823_current_strict_alpha_s_output_schema_statement_class_candidate_supported_exact_output_schema_statement_blocked.json"
OUT_SUMMARY = GENERATED / "p823_current_strict_alpha_s_output_schema_statement_class_candidate_supported_exact_output_schema_statement_blocked_summary.json"


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

    prereq = [IN_P822, IN_F822, IN_F820, IN_F819, IN_F818, IN_F814, IN_F812]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P823",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p822 = load_json(IN_P822)
    f822 = load_json(IN_F822)
    f820 = load_json(IN_F820)
    f819 = load_json(IN_F819)
    f818 = load_json(IN_F818)
    f814 = load_json(IN_F814)
    f812 = load_json(IN_F812)

    p822_theorem = p822.get("theorem_result") or {}
    f822_target = f822.get("target_object") or {}
    f820_target = f820.get("target_object") or {}
    f819_target = f819.get("target_object") or {}
    f818_target = f818.get("target_object") or {}
    f814_target = f814.get("target_object") or {}
    f812_target = f812.get("target_object") or {}

    token_hits: list[str] = []
    for path in sorted(GENERATED.glob("*.json")):
        name = path.name
        if name.startswith(("p822_", "f822_", "p823_", "f823_")):
            continue
        text = path.read_text(encoding="utf-8")
        if "exact_output_schema_statement" in text:
            token_hits.append(rel(path))

    checks = [
        {
            "id": "f822_already_names_exact_missing_statement_field",
            "pass": (
                f822.get("status")
                == "F822_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
                and f822_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1"
                and has_required_field(f822_target, "exact_output_schema_statement")
            ),
            "details": "F822 already isolates exact_output_schema_statement as one exact missing field of the output-schema target.",
        },
        {
            "id": "neighboring_targets_preserve_statement_slots",
            "pass": (
                has_required_field(f820_target, "boundary_output_schema")
                and has_required_field(f819_target, "selected_interface_output_schema")
                and has_required_field(f818_target, "exact_interface_output_schema")
                and has_required_field(f814_target, "selected_source_binding_output_schema")
                and has_required_field(f812_target, "selected_source_binding_output_schema")
            ),
            "details": "Nearby upstream and downstream targets do preserve output-schema statement slots on neighboring lanes.",
        },
        {
            "id": "p822_already_keeps_exact_output_schema_unexported",
            "pass": (
                p822.get("status")
                == "P822_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_SCHEMA_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
                and p822_theorem.get("lawful_schema_domain_admission_output_schema_exported_now") is False
            ),
            "details": "P822 already keeps the exact new-lane output-schema object unexported on the current repo state.",
        },
        {
            "id": "repo_scan_finds_no_exact_statement_export",
            "pass": token_hits == [],
            "details": "Repo scan finds no generated artifact exporting exact_output_schema_statement outside the new frozen F822 target lineage itself.",
        },
        {
            "id": "neighboring_statement_slots_remain_nonidentical_to_new_lane_statement",
            "pass": (
                f822_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1"
                and "boundary_output_schema" in (f822.get("target_refs") or {}).get("upstream_rule_or_interface_output_refs", [])
                and (f822.get("target_refs") or {}).get("downstream_schema_output_ref")
                == "selected_source_binding_output_schema"
            ),
            "details": "F822 already records neighboring statement slots only as nonidentical support refs, not as silent discharge of the new lane.",
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
        "P823_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        if all_pass
        else "P823_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P823",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p822_output_schema_probe": rel(IN_P822),
            "f822_output_schema_target_packet": rel(IN_F822),
            "f820_boundary_target_packet": rel(IN_F820),
            "f819_rule_target_packet": rel(IN_F819),
            "f818_interface_target_packet": rel(IN_F818),
            "f814_downstream_schema_target_packet": rel(IN_F814),
            "f812_downstream_rule_target_packet": rel(IN_F812),
        },
        "repo_scan_token_hits_for_exact_output_schema_statement": token_hits,
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "statement_class_support_stack": {
            "candidate_support_refs": [
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "selected_source_binding_output_schema",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "These are neighboring statement slots only; none exports the exact output-schema statement required by F822.",
        },
        "current_honest_reading": [
            "The repo preserves statement-level output slots around the new lane, but only as neighboring target fields.",
            "No current export yet names the exact output-schema statement required by F822 for the new T213/T216 -> alpha_s schema route.",
            "So the sharp blocker is now the still-missing exact statement object itself.",
        ],
        "recommended_next_packet": {
            "id": "F823_CURRENT_STRICT_ALPHA_S_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET",
            "goal": "Freeze the exact output-schema statement object still missing after statement-class support is present only at candidate grade.",
            "minimum_fields": [
                "lawful_schema_domain_admission_output_schema_target_ref",
                "statement_class_candidate_support_refs",
                "neighboring_statement_slot_refs",
                "exact_statement_required_form_ref",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P823",
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
