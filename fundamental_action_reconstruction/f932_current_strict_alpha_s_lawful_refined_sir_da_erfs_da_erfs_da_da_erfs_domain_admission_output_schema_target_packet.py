#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P932 = GENERATED / "p932_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_blocked.json"
IN_F931 = GENERATED / "f931_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_target_packet.json"
IN_F930 = GENERATED / "f930_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F929 = GENERATED / "f929_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json"
IN_F928 = GENERATED / "f928_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.json"
IN_F921 = GENERATED / "f921_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet.json"
IN_F924 = GENERATED / "f924_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f932_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet.json"
OUT_SUMMARY = GENERATED / "f932_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P932, IN_F931, IN_F930, IN_F929, IN_F928, IN_F921, IN_F924]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F932",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p932 = load_json(IN_P932)
    f931 = load_json(IN_F931)
    f930 = load_json(IN_F930)
    f929 = load_json(IN_F929)
    f928 = load_json(IN_F928)
    f921 = load_json(IN_F921)
    f924 = load_json(IN_F924)

    p932_theorem = p932.get("theorem_result") or {}
    f931_target = f931.get("target_object") or {}
    f930_target = f930.get("target_object") or {}
    f929_target = f929.get("target_object") or {}
    f928_target = f928.get("target_object") or {}
    f921_target = f921.get("target_object") or {}
    f924_target = f924.get("target_object") or {}

    if (
        p932.get("status")
        == "P932_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        and p932_theorem.get("output_schema_class_candidate_supported_now") is True
        and p932_theorem.get(
            "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema_exported_now"
        )
        is False
        and p932_theorem.get("sharp_blocker_field")
        == "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema"
        and p932_theorem.get("next_honest_move_is_freeze_exact_output_schema_target") is True
        and f931.get("status")
        == "F931_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F932_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F932_REQUIRES_REVIEW"

    artifact = {
        "stage": "F932",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedSirDaErfsDaErfsDaDaErfsDomainAdmissionOutputSchemaTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p932_output_schema_probe": rel(IN_P932),
            "f931_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F931),
            "f930_combined_boundary_target_packet": rel(IN_F930),
            "f929_upstream_rule_target_packet": rel(IN_F929),
            "f928_upstream_interface_target_packet": rel(IN_F928),
            "f921_neighboring_output_target_packet": rel(IN_F921),
            "f924_neighboring_exact_required_form_statement_target_packet": rel(IN_F924),
        },
        "why_this_packet_exists": [
            "F931 already freezes the deeper lawful refined domain-admission object and names one exact missing output-schema field.",
            "P932 shows that output-schema class support exists on nearby lanes, but the exact output schema for the new deeper lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1",
            "goal": "Freeze the exact lawful refined deeper domain-admission output-schema object still missing for the new lawful T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "lawful_refined_deeper_exact_required_form_statement_domain_admission_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F931 lawful refined deeper admission target and not silently replace the problem.",
                },
                {
                    "name": "output_schema_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade output-schema class support from neighboring targets and must not promote it into exact discharge.",
                },
                {
                    "name": "upstream_rule_or_interface_output_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit the upstream output-schema classes that remain nonidentical to the new deeper lane output schema.",
                },
                {
                    "name": "neighboring_output_schema_or_statement_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit neighboring output-schema or output-statement classes without silently reusing them as the new deeper lawful refined admission output schema.",
                },
                {
                    "name": "exact_output_schema_statement",
                    "required": True,
                    "hard_limit": "Must state what exact lawful refined deeper admission output would be required for the new lane without claiming that such output already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring output-schema classes, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "lawful_refined_deeper_exact_required_form_statement_domain_admission_target_ref": f931_target.get("object_id"),
            "output_schema_class_candidate_support_refs": [
                f930_target.get("object_id"),
                f929_target.get("object_id"),
                f928_target.get("object_id"),
                f921_target.get("object_id"),
                f924_target.get("object_id"),
            ],
            "upstream_rule_or_interface_output_refs": [
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
            ],
            "neighboring_output_schema_or_statement_refs": [
                f921_target.get("object_id"),
                f924_target.get("object_id"),
                "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                "exact_output_schema_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves output-schema class around the new deeper lawful refined domain-admission lane, but only at upstream, boundary, and neighboring target level.",
            "No current export yet names the exact lawful refined deeper domain-admission output schema required by F931.",
            "F932 freezes that exact missing output-schema object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_output_schema_statement for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact lawful refined deeper domain-admission output schema already exists.",
            "Does not claim that any upstream or neighboring output-schema class silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F932",
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
