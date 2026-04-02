#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P889 = GENERATED / "p889_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_output_schema_statement_blocked.json"
IN_F888 = GENERATED / "f888_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F886 = GENERATED / "f886_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F885 = GENERATED / "f885_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F884 = GENERATED / "f884_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_interface_target_packet.json"
IN_F877 = GENERATED / "f877_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F880 = GENERATED / "f880_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f889_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_output_schema_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f889_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_output_schema_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P889, IN_F888, IN_F886, IN_F885, IN_F884, IN_F877, IN_F880]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F889",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p889 = load_json(IN_P889)
    f888 = load_json(IN_F888)
    f886 = load_json(IN_F886)
    f885 = load_json(IN_F885)
    f884 = load_json(IN_F884)
    f877 = load_json(IN_F877)
    f880 = load_json(IN_F880)

    p889_theorem = p889.get("theorem_result") or {}
    f888_target = f888.get("target_object") or {}
    f886_target = f886.get("target_object") or {}
    f885_target = f885.get("target_object") or {}
    f884_target = f884.get("target_object") or {}
    f877_target = f877.get("target_object") or {}
    f880_target = f880.get("target_object") or {}

    if (
        p889.get("status")
        == "P889_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        and p889_theorem.get("output_schema_statement_class_candidate_supported_now") is True
        and p889_theorem.get("exact_output_schema_statement_exported_now") is False
        and p889_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
        and p889_theorem.get("next_honest_move_is_freeze_exact_statement_target") is True
        and f888.get("status")
        == "F888_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F889_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F889_REQUIRES_REVIEW"

    artifact = {
        "stage": "F889",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedShiftInterfaceRuleDomainAdmissionExactRequiredFormStatementDomainAdmissionExactRequiredFormStatementDomainAdmissionExactOutputSchemaStatementTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p889_exact_statement_probe": rel(IN_P889),
            "f888_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F888),
            "f886_boundary_target_packet": rel(IN_F886),
            "f885_rule_target_packet": rel(IN_F885),
            "f884_interface_target_packet": rel(IN_F884),
            "f877_neighboring_output_target_packet": rel(IN_F877),
            "f880_neighboring_exact_required_form_statement_target_packet": rel(IN_F880),
        },
        "why_this_packet_exists": [
            "F888 already freezes the deeper lawful refined domain-admission output-schema object and names one exact missing statement field.",
            "P889 shows that neighboring statement slots exist, but the exact statement needed by the new deeper lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1",
            "goal": "Freeze the exact deeper lawful refined domain-admission output-schema statement object still missing for the lawful refined shift-interface-rule T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F888 output-schema target and not silently replace the problem.",
                },
                {
                    "name": "statement_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring statement-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_statement_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring statement slots remain nonidentical to the new-lane statement.",
                },
                {
                    "name": "exact_statement_required_form_ref",
                    "required": True,
                    "hard_limit": "Must state what exact statement form is required for the new lane without claiming that the statement already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring statement slots, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_output_schema_target_ref": f888_target.get("object_id"),
            "statement_class_candidate_support_refs": [
                f886_target.get("object_id"),
                f885_target.get("object_id"),
                f884_target.get("object_id"),
                f877_target.get("object_id"),
                f880_target.get("object_id"),
            ],
            "neighboring_statement_slot_refs": [
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_output_schema_statement",
                "exact_required_form_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves statement-level output slots around the deeper lawful refined domain-admission lane, but only as neighboring target fields.",
            "No current export yet names the exact output-schema statement required by F888.",
            "F889 freezes that exact missing statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_statement_required_form_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact deeper lawful refined output-schema statement already exists.",
            "Does not claim that any neighboring output-schema statement slot silently discharges the new lawful refined lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F889",
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
