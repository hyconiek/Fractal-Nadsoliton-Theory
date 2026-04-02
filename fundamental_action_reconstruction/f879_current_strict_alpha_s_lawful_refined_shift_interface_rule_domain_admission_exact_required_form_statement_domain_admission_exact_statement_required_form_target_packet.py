#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P879 = GENERATED / "p879_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F878 = GENERATED / "f878_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F877 = GENERATED / "f877_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F876 = GENERATED / "f876_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F875 = GENERATED / "f875_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F874 = GENERATED / "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F873 = GENERATED / "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f879_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
OUT_SUMMARY = GENERATED / "f879_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P879, IN_F878, IN_F877, IN_F876, IN_F875, IN_F874, IN_F873, IN_F869]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F879",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p879 = load_json(IN_P879)
    f878 = load_json(IN_F878)
    f877 = load_json(IN_F877)
    f876 = load_json(IN_F876)
    f875 = load_json(IN_F875)
    f874 = load_json(IN_F874)
    f873 = load_json(IN_F873)
    f869 = load_json(IN_F869)

    p879_theorem = p879.get("theorem_result") or {}
    f878_target = f878.get("target_object") or {}
    f877_target = f877.get("target_object") or {}
    f876_target = f876.get("target_object") or {}
    f875_target = f875.get("target_object") or {}
    f874_target = f874.get("target_object") or {}
    f873_target = f873.get("target_object") or {}
    f869_target = f869.get("target_object") or {}

    if (
        p879.get("status")
        == "P879_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        and p879_theorem.get("statement_form_class_candidate_supported_now") is True
        and p879_theorem.get("exact_statement_required_form_exported_now") is False
        and p879_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
        and p879_theorem.get("next_honest_move_is_freeze_exact_required_form_target") is True
        and f878.get("status")
        == "F878_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F879_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F879_REQUIRES_REVIEW"

    artifact = {
        "stage": "F879",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedShiftInterfaceRuleDomainAdmissionExactRequiredFormStatementDomainAdmissionExactStatementRequiredFormTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p879_required_form_probe": rel(IN_P879),
            "f878_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F878),
            "f877_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet": rel(IN_F877),
            "f876_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet": rel(IN_F876),
            "f875_boundary_target_packet": rel(IN_F875),
            "f874_rule_target_packet": rel(IN_F874),
            "f873_interface_target_packet": rel(IN_F873),
            "f869_neighboring_exact_required_form_statement_target_packet": rel(IN_F869),
        },
        "why_this_packet_exists": [
            "F878 already freezes the lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission exact-output-schema-statement object and names one exact missing form field.",
            "P879 shows that neighboring statement slots and neighboring required-form supports exist, but the exact required form needed by the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1",
            "goal": "Freeze the exact statement-required form object still missing for the lawful refined shift-interface-rule T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "exact_output_schema_statement_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F878 statement target and not silently replace the problem.",
                },
                {
                    "name": "statement_form_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring form-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_form_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring form slots remain nonidentical to the new-lane form.",
                },
                {
                    "name": "exact_required_form_statement_ref",
                    "required": True,
                    "hard_limit": "Must state what exact required form is needed for the new lane without claiming that such form already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring form slots, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "exact_output_schema_statement_target_ref": f878_target.get("object_id"),
            "statement_form_class_candidate_support_refs": [
                f877_target.get("object_id"),
                f876_target.get("object_id"),
                f875_target.get("object_id"),
                f874_target.get("object_id"),
                f873_target.get("object_id"),
                f869_target.get("object_id"),
            ],
            "neighboring_form_slot_refs": [
                "exact_output_schema_statement",
                "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the lawful refined shift-interface-rule lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F878.",
            "F879 freezes that exact missing required-form object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_required_form_statement_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact statement-required form already exists.",
            "Does not claim that any neighboring statement slot or neighboring form support silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F879",
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
