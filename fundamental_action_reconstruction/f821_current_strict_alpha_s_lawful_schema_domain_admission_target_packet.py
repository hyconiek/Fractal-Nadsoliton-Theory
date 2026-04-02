#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P821 = GENERATED / "p821_current_strict_alpha_s_continued_nonidentification_candidate_supported_lawful_schema_domain_admission_blocked.json"
IN_F820 = GENERATED / "f820_current_strict_alpha_s_schema_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_N50 = GENERATED / "n50_current_legacy_ontological_kernel_to_strict_gate_kernel_nonidentification_theorem_summary.json"
IN_N117 = GENERATED / "n117_current_legacy_to_strict_kernel_bridge_package_nontransfer_theorem_summary.json"

OUT = GENERATED / "f821_current_strict_alpha_s_lawful_schema_domain_admission_target_packet.json"
OUT_SUMMARY = GENERATED / "f821_current_strict_alpha_s_lawful_schema_domain_admission_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P821, IN_F820, IN_P788, IN_N50, IN_N117]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F821",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p821 = load_json(IN_P821)
    f820 = load_json(IN_F820)
    p788 = load_json(IN_P788)
    n50 = load_json(IN_N50)
    n117 = load_json(IN_N117)

    p821_theorem = p821.get("theorem_result") or {}
    f820_target = f820.get("target_object") or {}

    if (
        p821.get("status")
        == "P821_CURRENT_STRICT_ALPHA_S_CONTINUED_NONIDENTIFICATION_CANDIDATE_SUPPORTED_LAWFUL_SCHEMA_DOMAIN_ADMISSION_BLOCKED"
        and p821_theorem.get("lawful_schema_domain_admission_exported_now") is False
        and p821_theorem.get("continued_nonidentification_boundary_exported_now") is False
        and p821_theorem.get("continued_nonidentification_boundary_candidate_supported_now") is True
        and p821_theorem.get("sharp_blocker_clause") == "lawful_schema_domain_admission"
        and p821_theorem.get("next_honest_move_is_freeze_lawful_schema_domain_admission_target") is True
        and f820.get("status")
        == "F820_EXECUTED_CURRENT_STRICT_ALPHA_S_SCHEMA_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
        and p788.get("status")
        == "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE"
    ):
        status = "F821_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_SCHEMA_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F821_REQUIRES_REVIEW"

    artifact = {
        "stage": "F821",
        "packet_name": "CurrentStrictAlphaSLawfulSchemaDomainAdmissionTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p821_clause_split_audit": rel(IN_P821),
            "f820_combined_boundary_target_packet": rel(IN_F820),
            "p788_generic_alpha_s_adapter_probe": rel(IN_P788),
            "n50_kernel_nonidentification_theorem": rel(IN_N50),
            "n117_package_nontransfer_theorem": rel(IN_N117),
        },
        "why_this_packet_exists": [
            "F820 freezes the larger boundary target containing both lawful-admission and continued-nonidentification branches.",
            "P821 shows that continued nonidentification has only candidate-supported-not-yet-exported status, while lawful schema-domain admission remains the sharp blocker.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_target_v1",
            "goal": "Freeze the exact lawful schema-domain admission object still missing after the weaker continued-nonidentification branch remains only candidate-supported on current repo state.",
            "required_fields": [
                {
                    "name": "schema_domain_or_nonidentification_boundary_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F820 combined boundary target and not silently replace the problem.",
                },
                {
                    "name": "continued_nonidentification_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade support for the weaker branch and must not promote it into exact export.",
                },
                {
                    "name": "generic_alpha_s_adapter_absence_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that no generic alpha_s adapter is exported on current repo state.",
                },
                {
                    "name": "scope_mismatch_nonidentification_or_nontransfer_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit that N50/N117 remain scope-mismatched analogy-only support for this new lane.",
                },
                {
                    "name": "lawful_schema_domain_admission_output_schema",
                    "required": True,
                    "hard_limit": "Must state what lawful admission would output for the F819/F820 lane without claiming that such output already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, continued-nonidentification export, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "schema_domain_or_nonidentification_boundary_target_ref": f820_target.get("object_id"),
            "continued_nonidentification_candidate_support_refs": [
                "P808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_ALPHA_S_DOMAIN_INTERFACE_BLOCKED",
                "N50_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_NONIDENTIFICATION_THEOREM_NO_FALSE_PASS",
                "N117_DISCHARGED_CURRENT_LEGACY_TO_STRICT_KERNEL_BRIDGE_PACKAGE_NONTRANSFER_THEOREM_NO_FALSE_PASS",
            ],
            "generic_alpha_s_adapter_absence_ref": "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE",
            "scope_mismatch_nonidentification_or_nontransfer_refs": [
                "current_legacy_to_strict_kernel_bridge_question_only",
                "current_legacy_to_strict_package_question_only",
            ],
        },
        "current_honest_reading": [
            "Continued nonidentification now has only candidate-grade support for the new T213/T216 lane, not exact boundary export.",
            "Lawful schema-domain admission remains the exact missing object that would have to be exported before the new lane could lawfully enter the current alpha_s schema problem.",
            "F821 freezes that sharper blocker without pretending that admission or boundary export already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply lawful_schema_domain_admission_output_schema for alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that lawful schema-domain admission already exists.",
            "Does not claim that continued nonidentification is already exported as an exact boundary object for the new lane.",
            "Does not claim that N50 or N117 discharge the new T213/T216 -> alpha_s schema lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s schema domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F821",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
