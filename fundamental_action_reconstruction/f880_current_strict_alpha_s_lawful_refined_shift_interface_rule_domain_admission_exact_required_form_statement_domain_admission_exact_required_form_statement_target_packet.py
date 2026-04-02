#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P880 = GENERATED / "p880_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
IN_F879 = GENERATED / "f879_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F878 = GENERATED / "f878_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F877 = GENERATED / "f877_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F876 = GENERATED / "f876_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet.json"
IN_F875 = GENERATED / "f875_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F874 = GENERATED / "f874_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F873 = GENERATED / "f873_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_interface_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f880_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f880_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P880, IN_F879, IN_F878, IN_F877, IN_F876, IN_F875, IN_F874, IN_F873, IN_F869]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F880",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p880 = load_json(IN_P880)
    f879 = load_json(IN_F879)
    f878 = load_json(IN_F878)
    f877 = load_json(IN_F877)
    f876 = load_json(IN_F876)
    f875 = load_json(IN_F875)
    f874 = load_json(IN_F874)
    f873 = load_json(IN_F873)
    f869 = load_json(IN_F869)

    p880_theorem = p880.get("theorem_result") or {}
    f879_target = f879.get("target_object") or {}
    f878_target = f878.get("target_object") or {}
    f877_target = f877.get("target_object") or {}
    f876_target = f876.get("target_object") or {}
    f875_target = f875.get("target_object") or {}
    f874_target = f874.get("target_object") or {}
    f873_target = f873.get("target_object") or {}
    f869_target = f869.get("target_object") or {}

    if (
        p880.get("status")
        == "P880_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        and p880_theorem.get("required_form_statement_class_candidate_supported_now") is True
        and p880_theorem.get("exact_required_form_statement_exported_now") is False
        and p880_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
        and p880_theorem.get("next_honest_move_is_freeze_exact_required_form_statement_target") is True
        and f879.get("status")
        == "F879_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F880_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F880_REQUIRES_REVIEW"

    artifact = {
        "stage": "F880",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedShiftInterfaceRuleDomainAdmissionExactRequiredFormStatementDomainAdmissionExactRequiredFormStatementTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p880_required_form_statement_probe": rel(IN_P880),
            "f879_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet": rel(IN_F879),
            "f878_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F878),
            "f877_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet": rel(IN_F877),
            "f876_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_target_packet": rel(IN_F876),
            "f875_boundary_target_packet": rel(IN_F875),
            "f874_rule_target_packet": rel(IN_F874),
            "f873_interface_target_packet": rel(IN_F873),
            "f869_neighboring_exact_required_form_statement_target_packet": rel(IN_F869),
        },
        "why_this_packet_exists": [
            "F879 already freezes the lawful refined shift-interface-rule domain-admission exact-required-form-statement domain-admission exact-statement-required-form object and names one exact missing required-form-statement field.",
            "P880 shows that neighboring statement slots and neighboring required-form supports exist, but the exact required-form statement needed by the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1",
            "goal": "Freeze the exact required-form-statement object still missing for the lawful refined shift-interface-rule T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "exact_statement_required_form_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F879 required-form target and not silently replace the problem.",
                },
                {
                    "name": "required_form_statement_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring statement/form-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_statement_or_form_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring statement or form slots remain nonidentical to the new-lane statement.",
                },
                {
                    "name": "exact_required_form_statement_ref",
                    "required": True,
                    "hard_limit": "Must state what exact required-form statement is needed for the new lane without claiming that the statement already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring statement or form slots, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "exact_statement_required_form_target_ref": f879_target.get("object_id"),
            "required_form_statement_class_candidate_support_refs": [
                f878_target.get("object_id"),
                f877_target.get("object_id"),
                f876_target.get("object_id"),
                f875_target.get("object_id"),
                f874_target.get("object_id"),
                f873_target.get("object_id"),
                f869_target.get("object_id"),
            ],
            "neighboring_statement_or_form_slot_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the lawful refined shift-interface-rule domain-admission lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F879.",
            "F880 freezes that exact missing required-form-statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow same-lane exhaustion-boundary audit testing whether any further passive split remains below alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact required-form statement already exists.",
            "Does not claim that any neighboring statement or form slot silently discharges the lawful refined new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F880",
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
