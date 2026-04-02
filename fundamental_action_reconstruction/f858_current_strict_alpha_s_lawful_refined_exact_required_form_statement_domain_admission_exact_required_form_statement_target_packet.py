#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P858 = GENERATED / "p858_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_blocked.json"
IN_F857 = GENERATED / "f857_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
IN_F856 = GENERATED / "f856_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F855 = GENERATED / "f855_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F854 = GENERATED / "f854_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_target_packet.json"
IN_F853 = GENERATED / "f853_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F852 = GENERATED / "f852_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F851 = GENERATED / "f851_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_refined_exact_required_form_statement_interface_target_packet.json"
IN_F847 = GENERATED / "f847_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f858_current_strict_alpha_s_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P858, IN_F857, IN_F856, IN_F855, IN_F854, IN_F853, IN_F852, IN_F851, IN_F847]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F858",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p858 = load_json(IN_P858)
    f857 = load_json(IN_F857)
    f856 = load_json(IN_F856)
    f855 = load_json(IN_F855)
    f854 = load_json(IN_F854)
    f853 = load_json(IN_F853)
    f852 = load_json(IN_F852)
    f851 = load_json(IN_F851)
    f847 = load_json(IN_F847)

    p858_theorem = p858.get("theorem_result") or {}
    f857_target = f857.get("target_object") or {}
    f856_target = f856.get("target_object") or {}
    f855_target = f855.get("target_object") or {}
    f854_target = f854.get("target_object") or {}
    f853_target = f853.get("target_object") or {}
    f852_target = f852.get("target_object") or {}
    f851_target = f851.get("target_object") or {}
    f847_target = f847.get("target_object") or {}

    if (
        p858.get("status")
        == "P858_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        and p858_theorem.get("required_form_statement_class_candidate_supported_now") is True
        and p858_theorem.get("exact_required_form_statement_exported_now") is False
        and p858_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
        and p858_theorem.get("next_honest_move_is_freeze_exact_required_form_statement_target") is True
        and f857.get("status")
        == "F857_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F858_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F858_REQUIRES_REVIEW"

    artifact = {
        "stage": "F858",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedExactRequiredFormStatementDomainAdmissionExactRequiredFormStatementTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p858_required_form_statement_probe": rel(IN_P858),
            "f857_required_form_target_packet": rel(IN_F857),
            "f856_refined_exact_statement_target_packet": rel(IN_F856),
            "f855_refined_output_schema_target_packet": rel(IN_F855),
            "f854_refined_lawful_admission_target_packet": rel(IN_F854),
            "f853_boundary_target_packet": rel(IN_F853),
            "f852_rule_target_packet": rel(IN_F852),
            "f851_interface_target_packet": rel(IN_F851),
            "f847_neighboring_exact_required_form_statement_target_packet": rel(IN_F847),
        },
        "why_this_packet_exists": [
            "F857 already freezes the refined exact statement-required-form object and names one exact missing required-form-statement field.",
            "P858 shows that neighboring statement slots and neighboring required-form targets exist, but the exact required-form statement needed by the refined lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1",
            "goal": "Freeze the exact required-form statement object still missing for the refined lawful T213/T216 -> alpha_s exact-required-form-statement domain-admission lane.",
            "required_fields": [
                {
                    "name": "exact_statement_required_form_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F857 required-form target and not silently replace the problem.",
                },
                {
                    "name": "required_form_statement_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring statement/form-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_statement_or_form_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring statement or form slots remain nonidentical to the refined new-lane statement.",
                },
                {
                    "name": "exact_required_form_statement_ref",
                    "required": True,
                    "hard_limit": "Must state what exact required-form statement is needed for the refined new lane without claiming that the statement already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring statement or form slots, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "exact_statement_required_form_target_ref": f857_target.get("object_id"),
            "required_form_statement_class_candidate_support_refs": [
                f856_target.get("object_id"),
                f855_target.get("object_id"),
                f854_target.get("object_id"),
                f853_target.get("object_id"),
                f852_target.get("object_id"),
                f851_target.get("object_id"),
                f847_target.get("object_id"),
            ],
            "neighboring_statement_or_form_slot_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_refined_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the refined new lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F857.",
            "F858 freezes that exact missing required-form statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow same-lane exhaustion-boundary audit testing whether any further passive split remains below alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact required-form statement already exists.",
            "Does not claim that any neighboring statement or form slot silently discharges the refined new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F858",
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
