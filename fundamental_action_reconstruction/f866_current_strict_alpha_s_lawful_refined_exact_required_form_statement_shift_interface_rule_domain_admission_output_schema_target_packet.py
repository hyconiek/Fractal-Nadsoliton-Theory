#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P866 = GENERATED / "p866_current_strict_alpha_s_output_schema_class_candidate_supported_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_blocked.json"
IN_F865 = GENERATED / "f865_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_packet.json"
IN_F864 = GENERATED / "f864_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F863 = GENERATED / "f863_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F862 = GENERATED / "f862_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_lawful_refined_exact_required_form_statement_interface_target_packet.json"
IN_F855 = GENERATED / "f855_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F858 = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f866_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_packet.json"
OUT_SUMMARY = GENERATED / "f866_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P866, IN_F865, IN_F864, IN_F863, IN_F862, IN_F855, IN_F858]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F866",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p866 = load_json(IN_P866)
    f865 = load_json(IN_F865)
    f864 = load_json(IN_F864)
    f863 = load_json(IN_F863)
    f862 = load_json(IN_F862)
    f855 = load_json(IN_F855)
    f858 = load_json(IN_F858)

    p866_theorem = p866.get("theorem_result") or {}
    f865_target = f865.get("target_object") or {}
    f864_target = f864.get("target_object") or {}
    f863_target = f863.get("target_object") or {}
    f862_target = f862.get("target_object") or {}
    f855_target = f855.get("target_object") or {}
    f858_target = f858.get("target_object") or {}

    if (
        p866.get("status")
        == "P866_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OUTPUT_SCHEMA_BLOCKED"
        and p866_theorem.get("output_schema_class_candidate_supported_now") is True
        and p866_theorem.get(
            "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_exported_now"
        )
        is False
        and p866_theorem.get("sharp_blocker_field")
        == "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema"
        and p866_theorem.get("next_honest_move_is_freeze_exact_output_schema_target") is True
        and f865.get("status")
        == "F865_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F866_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F866_REQUIRES_REVIEW"

    artifact = {
        "stage": "F866",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedExactRequiredFormStatementShiftInterfaceRuleDomainAdmissionOutputSchemaTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p866_output_schema_probe": rel(IN_P866),
            "f865_lawful_refined_shift_interface_rule_domain_admission_target_packet": rel(IN_F865),
            "f864_combined_boundary_target_packet": rel(IN_F864),
            "f863_upstream_rule_target_packet": rel(IN_F863),
            "f862_upstream_interface_target_packet": rel(IN_F862),
            "f855_neighboring_output_target_packet": rel(IN_F855),
            "f858_neighboring_exact_required_form_statement_target_packet": rel(IN_F858),
        },
        "why_this_packet_exists": [
            "F865 already freezes the lawful refined shift-interface-rule domain-admission object and names one exact missing output field.",
            "P866 shows that output-schema class support exists on nearby lanes, but the exact output schema for the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1",
            "goal": "Freeze the exact lawful refined exact-required-form-statement shift-interface-rule domain-admission output-schema object still missing for the new lawful T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F865 lawful refined admission target and not silently replace the problem.",
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
                    "hard_limit": "Must keep explicit neighboring output-schema or output-statement classes without silently reusing them as the new lawful refined admission output schema.",
                },
                {
                    "name": "exact_output_schema_statement",
                    "required": True,
                    "hard_limit": "Must state what exact lawful refined admission output would be required for the new lane without claiming that such output already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring output-schema classes, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_target_ref": f865_target.get("object_id"),
            "output_schema_class_candidate_support_refs": [
                f864_target.get("object_id"),
                f863_target.get("object_id"),
                f862_target.get("object_id"),
                f855_target.get("object_id"),
                f858_target.get("object_id"),
            ],
            "upstream_rule_or_interface_output_refs": [
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
            ],
            "neighboring_output_schema_or_statement_refs": [
                "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_v1",
                "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1",
                "lawful_refined_exact_required_form_statement_domain_admission_output_schema",
                "exact_output_schema_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves output-schema class around the new lawful refined shift-interface-rule lane, but only at upstream, boundary, and neighboring target level.",
            "No current export yet names the exact lawful refined shift-interface-rule domain-admission output schema required by F865.",
            "F866 freezes that exact missing output-schema object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_output_schema_statement for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_output_schema_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact lawful refined shift-interface-rule domain-admission output schema already exists.",
            "Does not claim that any upstream or neighboring output-schema class silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s lawful refined exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F866",
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
