#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P901 = GENERATED / "p901_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F900 = GENERATED / "f900_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F899 = GENERATED / "f899_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F898 = GENERATED / "f898_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_target_packet.json"
IN_F897 = GENERATED / "f897_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F896 = GENERATED / "f896_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_adapter_or_carrier_rule_target_packet.json"
IN_F895 = GENERATED / "f895_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_interface_target_packet.json"
IN_F891 = GENERATED / "f891_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f901_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
OUT_SUMMARY = GENERATED / "f901_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P901, IN_F900, IN_F899, IN_F898, IN_F897, IN_F896, IN_F895, IN_F891]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F901",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p901 = load_json(IN_P901)
    f900 = load_json(IN_F900)
    f899 = load_json(IN_F899)
    f898 = load_json(IN_F898)
    f897 = load_json(IN_F897)
    f896 = load_json(IN_F896)
    f895 = load_json(IN_F895)
    f891 = load_json(IN_F891)

    p901_theorem = p901.get("theorem_result") or {}
    f900_target = f900.get("target_object") or {}
    f899_target = f899.get("target_object") or {}
    f898_target = f898.get("target_object") or {}
    f897_target = f897.get("target_object") or {}
    f896_target = f896.get("target_object") or {}
    f895_target = f895.get("target_object") or {}
    f891_target = f891.get("target_object") or {}

    if (
        p901.get("status")
        == "P901_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        and p901_theorem.get("statement_form_class_candidate_supported_now") is True
        and p901_theorem.get("exact_statement_required_form_exported_now") is False
        and p901_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
        and p901_theorem.get("next_honest_move_is_freeze_exact_required_form_target") is True
        and f900.get("status")
        == "F900_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F901_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F901_REQUIRES_REVIEW"

    artifact = {
        "stage": "F901",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedSirDaErfsDaErfsDaExactRequiredFormStatementDomainAdmissionExactStatementRequiredFormTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p901_required_form_probe": rel(IN_P901),
            "f900_lawful_refined_deeper_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F900),
            "f899_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F899),
            "f898_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F898),
            "f897_boundary_target_packet": rel(IN_F897),
            "f896_rule_target_packet": rel(IN_F896),
            "f895_interface_target_packet": rel(IN_F895),
            "f891_neighboring_exact_required_form_statement_target_packet": rel(IN_F891),
        },
        "why_this_packet_exists": [
            "F900 already freezes the deeper lawful refined domain-admission exact-output-schema-statement object and names one exact missing form field.",
            "P901 shows that neighboring statement slots and neighboring required-form supports exist, but the exact required form needed by the new deeper lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1",
            "goal": "Freeze the exact deeper lawful refined domain-admission statement-required form object still missing for the lawful refined shift-interface-rule T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "exact_output_schema_statement_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F900 statement target and not silently replace the problem.",
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
            "exact_output_schema_statement_target_ref": f900_target.get("object_id"),
            "statement_form_class_candidate_support_refs": [
                f899_target.get("object_id"),
                f898_target.get("object_id"),
                f897_target.get("object_id"),
                f896_target.get("object_id"),
                f895_target.get("object_id"),
                f891_target.get("object_id"),
            ],
            "neighboring_form_slot_refs": [
                "exact_output_schema_statement",
                "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the deeper lawful refined domain-admission lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F900.",
            "F901 freezes that exact missing required-form object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_required_form_statement_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact deeper lawful refined statement-required form already exists.",
            "Does not claim that any neighboring statement slot or neighboring form support silently discharges the new lawful refined lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F901",
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
