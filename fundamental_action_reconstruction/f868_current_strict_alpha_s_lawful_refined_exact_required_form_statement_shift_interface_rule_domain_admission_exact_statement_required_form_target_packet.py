#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P868 = GENERATED / "p868_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_blocked.json"
IN_F867 = GENERATED / "f867_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F866 = GENERATED / "f866_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_packet.json"
IN_F865 = GENERATED / "f865_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_packet.json"
IN_F864 = GENERATED / "f864_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F863 = GENERATED / "f863_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F862 = GENERATED / "f862_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_interface_target_packet.json"
IN_F858 = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f868_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_packet.json"
OUT_SUMMARY = GENERATED / "f868_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P868, IN_F867, IN_F866, IN_F865, IN_F864, IN_F863, IN_F862, IN_F858]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F868",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p868 = load_json(IN_P868)
    f867 = load_json(IN_F867)
    f866 = load_json(IN_F866)
    f865 = load_json(IN_F865)
    f864 = load_json(IN_F864)
    f863 = load_json(IN_F863)
    f862 = load_json(IN_F862)
    f858 = load_json(IN_F858)

    p868_theorem = p868.get("theorem_result") or {}
    f867_target = f867.get("target_object") or {}
    f866_target = f866.get("target_object") or {}
    f865_target = f865.get("target_object") or {}
    f864_target = f864.get("target_object") or {}
    f863_target = f863.get("target_object") or {}
    f862_target = f862.get("target_object") or {}
    f858_target = f858.get("target_object") or {}

    if (
        p868.get("status")
        == "P868_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        and p868_theorem.get("statement_form_class_candidate_supported_now") is True
        and p868_theorem.get("exact_statement_required_form_exported_now") is False
        and p868_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
        and p868_theorem.get("next_honest_move_is_freeze_exact_required_form_target") is True
        and f867.get("status")
        == "F867_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F868_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F868_REQUIRES_REVIEW"

    artifact = {
        "stage": "F868",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedExactRequiredFormStatementShiftInterfaceRuleDomainAdmissionExactStatementRequiredFormTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p868_required_form_probe": rel(IN_P868),
            "f867_lawful_refined_shift_interface_rule_exact_statement_target_packet": rel(IN_F867),
            "f866_lawful_refined_shift_interface_rule_output_schema_target_packet": rel(IN_F866),
            "f865_lawful_refined_shift_interface_rule_admission_target_packet": rel(IN_F865),
            "f864_boundary_target_packet": rel(IN_F864),
            "f863_rule_target_packet": rel(IN_F863),
            "f862_interface_target_packet": rel(IN_F862),
            "f858_neighboring_exact_required_form_statement_target_packet": rel(IN_F858),
        },
        "why_this_packet_exists": [
            "F867 already freezes the lawful refined shift-interface-rule exact statement object and names one exact missing required-form field.",
            "P868 shows that neighboring statement slots and neighboring required-form targets exist, but the exact required form needed by the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_v1",
            "goal": "Freeze the exact statement-required form object still missing for the lawful refined shift-interface-rule T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "exact_output_schema_statement_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F867 statement target and not silently replace the problem.",
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
            "exact_output_schema_statement_target_ref": f867_target.get("object_id"),
            "statement_form_class_candidate_support_refs": [
                f866_target.get("object_id"),
                f865_target.get("object_id"),
                f864_target.get("object_id"),
                f863_target.get("object_id"),
                f862_target.get("object_id"),
                f858_target.get("object_id"),
            ],
            "neighboring_form_slot_refs": [
                "exact_output_schema_statement",
                "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the lawful refined shift-interface-rule lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F867.",
            "F868 freezes that exact missing required-form object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_required_form_statement_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact statement-required form already exists.",
            "Does not claim that any neighboring statement slot or neighboring form support silently discharges the lawful refined new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s lawful refined exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F868",
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
