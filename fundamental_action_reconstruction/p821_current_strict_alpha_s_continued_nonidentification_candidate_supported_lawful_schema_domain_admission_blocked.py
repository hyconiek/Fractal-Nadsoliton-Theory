#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P820 = GENERATED / "p820_current_strict_alpha_s_action_schema_candidate_supported_schema_domain_admission_boundary_blocked.json"
IN_F820 = GENERATED / "f820_current_strict_alpha_s_schema_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_P808 = GENERATED / "p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_N50 = GENERATED / "n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
IN_N117 = GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"

OUT = GENERATED / "p821_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_schema_domain_admission_blocked.json"
OUT_SUMMARY = GENERATED / "p821_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_schema_domain_admission_blocked_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P820, IN_F820, IN_P808, IN_P788, IN_N50, IN_N117]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P821",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p820 = load_json(IN_P820)
    f820 = load_json(IN_F820)
    p808 = load_json(IN_P808)
    p788 = load_json(IN_P788)
    n50 = load_json(IN_N50)
    n117 = load_json(IN_N117)

    p820_split = p820.get("clause_split_audit") or {}
    f820_target = f820.get("target_object") or {}
    p808_theorem = p808.get("theorem_result") or {}
    n50_theorem = n50.get("theorem_result") or {}
    n117_theorem = n117.get("theorem_result") or {}

    checks = [
        {
            "id": "f820_freezes_combined_boundary_target",
            "pass": (
                f820.get("status")
                == "F820_EXECUTED_CURRENT_STRICT_ALPHA_S_SCHEMA_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
                and f820_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_schema_domain_admission_or_nonidentification_boundary_target_v1"
                and any(
                    isinstance(item, dict)
                    and item.get("name") == "old_shannon_nontransfer_ref"
                    for item in (f820_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict)
                    and item.get("name") == "boundary_output_schema"
                    for item in (f820_target.get("required_fields") or [])
                )
            ),
            "details": "F820 already freezes one exact combined boundary target with explicit room for either lawful admission or continued nonidentification output, but does not discharge either branch.",
        },
        {
            "id": "p820_keeps_lawful_schema_domain_admission_blocked",
            "pass": (
                p820.get("status")
                == "P820_CURRENT_STRICT_ALPHA_S_ACTION_SCHEMA_CANDIDATE_SUPPORTED_SCHEMA_DOMAIN_ADMISSION_BOUNDARY_BLOCKED"
                and p820_split.get("schema_domain_admission_boundary_clause_status") == "blocked_nonexport"
                and p820_split.get("sharp_blocker_clause") == "schema_domain_admission_boundary"
            ),
            "details": "P820 already keeps lawful schema-domain admission blocked and localizes the sharp blocker inside the larger F820 boundary target.",
        },
        {
            "id": "p808_supports_blocked_admission_and_foreign_domain_fence",
            "pass": (
                p808.get("status")
                == "P808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_ALPHA_S_DOMAIN_INTERFACE_BLOCKED"
                and p808_theorem.get("strict_source_shannon_provider_shift_candidate_reference_lane_admitted")
                is True
                and p808_theorem.get("alpha_s_domain_interface_exported") is False
                and any(
                    isinstance(item, dict)
                    and item.get("id") == "alpha_s_domain_admission_still_blocked"
                    and item.get("pass") is True
                    for item in (p808.get("checks") or [])
                )
            ),
            "details": "P808 provides real blocked-admission support with explicit foreign-domain fencing, which is compatible with candidate-only continued nonidentification support but not with lawful admission export.",
        },
        {
            "id": "p788_still_reports_no_generic_alpha_s_adapter",
            "pass": (
                p788.get("status")
                == "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE"
                and any(
                    isinstance(item, dict)
                    and item.get("id") == "no_exported_dimensionless_alpha_s_adapter_detected"
                    and item.get("pass") is True
                    for item in (p788.get("checks") or [])
                )
            ),
            "details": "P788 still blocks lawful schema-domain admission by keeping generic alpha_s adapter export absent on current repo state.",
        },
        {
            "id": "n50_is_exact_nonidentification_but_scope_mismatched_for_new_lane",
            "pass": (
                n50.get("status")
                == "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS"
                and n50.get("scope") == "current_legacy_to_strict_kernel_bridge_question_only"
                and n50_theorem.get("rigorous_nonidentification_on_current_repo_state") is True
            ),
            "details": "N50 is a real exact nonidentification theorem, but only on the legacy-vs-strict kernel bridge scope, so it can support nonidentification analogy only here.",
        },
        {
            "id": "n117_is_exact_nontransfer_but_scope_mismatched_for_new_lane",
            "pass": (
                n117.get("status")
                == "N117_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_PACKAGE_NONTRANSFER_THEOREM_NO_FALSE_PASS"
                and n117.get("scope") == "current_legacy_to_strict_package_question_only"
                and n117_theorem.get("legacy_to_strict_package_nontransfer_on_current_repo_state") is True
            ),
            "details": "N117 is a real exact nontransfer theorem, but only on the legacy-to-strict package scope, so it remains analogy-only for the new T213/T216 lane.",
        },
        {
            "id": "no_exact_new_lane_boundary_export_present",
            "pass": (
                p820_split.get("schema_domain_admission_boundary_clause_status") == "blocked_nonexport"
                and p808_theorem.get("alpha_s_domain_interface_exported") is False
                and n50.get("scope") != "current_t213_t216_to_alpha_s_schema_lane"
                and n117.get("scope") != "current_t213_t216_to_alpha_s_schema_lane"
            ),
            "details": "No current export names an exact lawful-admission object or an exact continued-nonidentification boundary object for the new T213/T216 -> alpha_s schema lane.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "lawful_schema_domain_admission_exported_now": False if all_pass else None,
        "continued_nonidentification_boundary_exported_now": False if all_pass else None,
        "continued_nonidentification_boundary_candidate_supported_now": (
            checks[0]["pass"] and checks[2]["pass"] and checks[4]["pass"] and checks[5]["pass"]
        ),
        "sharp_blocker_clause": (
            "lawful_schema_domain_admission"
            if all_pass
            else None
        ),
        "next_honest_move_is_freeze_lawful_schema_domain_admission_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P821_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_SCHEMA_DOMAIN_ADMISSION_BLOCKED"
        if all_pass
        else "P821_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P821",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p820_clause_split_audit": rel(IN_P820),
            "f820_boundary_target_packet": rel(IN_F820),
            "p808_blocked_admission_reference_probe": rel(IN_P808),
            "p788_generic_alpha_s_adapter_probe": rel(IN_P788),
            "n50_kernel_nonidentification_theorem": rel(IN_N50),
            "n117_package_nontransfer_theorem": rel(IN_N117),
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "candidate_support_stack": {
            "continued_nonidentification_candidate_support_refs": [
                "P808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_ALPHA_S_DOMAIN_INTERFACE_BLOCKED",
                "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS",
                "N117_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_PACKAGE_NONTRANSFER_THEOREM_NO_FALSE_PASS",
            ],
            "support_grade": "candidate_supported_not_yet_exported",
            "scope_limit": "N50 and N117 are analogy-only for this lane; they do not discharge exact boundary export on the new T213/T216 route.",
        },
        "current_honest_reading": [
            "No current export yet supplies lawful schema-domain admission for the new T213/T216 -> alpha_s schema lane.",
            "Continued nonidentification does have candidate-only support through blocked-admission evidence plus nontransfer analogy, but it is still not exported as an exact boundary object for this lane.",
            "So the sharp blocker is now the still-missing lawful schema-domain admission object, not the weaker continued-nonidentification branch.",
        ],
        "recommended_next_packet": {
            "id": "F821_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_TARGET_PACKET",
            "goal": "Freeze the exact lawful schema-domain admission object that remains missing after continued nonidentification is supported only at candidate grade.",
            "minimum_fields": [
                "schema_domain_or_nonidentification_boundary_target_ref",
                "continued_nonidentification_candidate_support_refs",
                "generic_alpha_s_adapter_absence_ref",
                "scope_mismatch_nonidentification_or_nontransfer_refs",
                "lawful_schema_domain_admission_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P821",
        "status": status,
        "as_of": AS_OF,
        "lawful_schema_domain_admission_exported_now": theorem_result[
            "lawful_schema_domain_admission_exported_now"
        ],
        "continued_nonidentification_boundary_exported_now": theorem_result[
            "continued_nonidentification_boundary_exported_now"
        ],
        "continued_nonidentification_boundary_candidate_supported_now": theorem_result[
            "continued_nonidentification_boundary_candidate_supported_now"
        ],
        "sharp_blocker_clause": theorem_result["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
