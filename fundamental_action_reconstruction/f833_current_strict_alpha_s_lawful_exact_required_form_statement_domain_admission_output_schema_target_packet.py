#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P833 = GENERATED / "p833_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_output_schema_blocked.json"
IN_F832 = GENERATED / "f832_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_target_packet.json"
IN_F831 = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F829 = GENERATED / "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
OUT_SUMMARY = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P833, IN_F832, IN_F831, IN_F830, IN_F829, IN_F822, IN_F825]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F833",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p833 = load_json(IN_P833)
    f832 = load_json(IN_F832)
    f831 = load_json(IN_F831)
    f830 = load_json(IN_F830)
    f829 = load_json(IN_F829)
    f822 = load_json(IN_F822)
    f825 = load_json(IN_F825)

    p833_theorem = p833.get("theorem_result") or {}
    f832_target = f832.get("target_object") or {}
    f831_target = f831.get("target_object") or {}
    f830_target = f830.get("target_object") or {}
    f829_target = f829.get("target_object") or {}
    f822_target = f822.get("target_object") or {}
    f825_target = f825.get("target_object") or {}

    if (
        p833.get("status")
        == "P833_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        and p833_theorem.get("output_schema_class_candidate_supported_now") is True
        and p833_theorem.get("lawful_exact_required_form_statement_domain_admission_output_schema_exported_now")
        is False
        and p833_theorem.get("sharp_blocker_field")
        == "lawful_exact_required_form_statement_domain_admission_output_schema"
        and p833_theorem.get("next_honest_move_is_freeze_exact_output_schema_target") is True
        and f832.get("status")
        == "F832_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F833_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F833_REQUIRES_REVIEW"

    artifact = {
        "stage": "F833",
        "packet_name": "CurrentStrictAlphaSLawfulExactRequiredFormStatementDomainAdmissionOutputSchemaTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p833_output_schema_probe": rel(IN_P833),
            "f832_lawful_exact_required_form_statement_domain_admission_target_packet": rel(IN_F832),
            "f831_combined_boundary_target_packet": rel(IN_F831),
            "f830_upstream_rule_target_packet": rel(IN_F830),
            "f829_upstream_interface_target_packet": rel(IN_F829),
            "f822_neighboring_schema_output_target_packet": rel(IN_F822),
            "f825_downstream_exact_required_form_statement_target_packet": rel(IN_F825),
        },
        "why_this_packet_exists": [
            "F832 already freezes the lawful exact-required-form-statement domain-admission object and names one exact missing output field.",
            "P833 shows that output-schema class support exists on nearby lanes, but the exact lawful exact-required-form-statement domain-admission output schema for the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_target_v1",
            "goal": "Freeze the exact lawful exact-required-form-statement domain-admission output-schema object still missing for the new T213/T216 -> alpha_s exact-required-form-statement lane.",
            "required_fields": [
                {
                    "name": "lawful_exact_required_form_statement_domain_admission_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F832 lawful-admission target and not silently replace the problem.",
                },
                {
                    "name": "output_schema_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade output-schema class support from neighboring targets and must not promote it into exact discharge.",
                },
                {
                    "name": "upstream_rule_or_interface_output_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit the upstream output-schema classes that remain nonidentical to the new lane output schema.",
                },
                {
                    "name": "neighboring_output_schema_or_statement_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit neighboring output-schema or output-statement classes without silently reusing them as the new lawful-admission output schema.",
                },
                {
                    "name": "exact_output_schema_statement",
                    "required": True,
                    "hard_limit": "Must state what exact lawful-admission output would be required for the new lane without claiming that such output already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring output-schema classes, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "lawful_exact_required_form_statement_domain_admission_target_ref": f832_target.get("object_id"),
            "output_schema_class_candidate_support_refs": [
                f831_target.get("object_id"),
                f830_target.get("object_id"),
                f829_target.get("object_id"),
                f822_target.get("object_id"),
                f825_target.get("object_id"),
            ],
            "upstream_rule_or_interface_output_refs": [
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
            ],
            "neighboring_output_schema_or_statement_refs": [
                "alpha_s_pair12_source_side_branch_selection_provider_lawful_schema_domain_admission_output_schema_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_target_v1",
                "lawful_schema_domain_admission_output_schema",
                "exact_output_schema_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves output-schema class around the new lane, but only at upstream, boundary, and neighboring target level.",
            "No current export yet names the exact lawful exact-required-form-statement domain-admission output schema required by F832.",
            "F833 freezes that exact missing output-schema object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_output_schema_statement for alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_output_schema_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact lawful exact-required-form-statement domain-admission output schema already exists.",
            "Does not claim that any upstream or neighboring output-schema class silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F833",
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
