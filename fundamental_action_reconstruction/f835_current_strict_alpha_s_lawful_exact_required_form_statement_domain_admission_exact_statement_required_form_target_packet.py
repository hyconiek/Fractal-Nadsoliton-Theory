#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P835 = GENERATED / "p835_current_strict_alpha_s_statement_form_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_blocked.json"
IN_F834 = GENERATED / "f834_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F833 = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F832 = GENERATED / "f832_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_target_packet.json"
IN_F831 = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F829 = GENERATED / "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"

OUT = GENERATED / "f835_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet.json"
OUT_SUMMARY = GENERATED / "f835_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_P835,
        IN_F834,
        IN_F833,
        IN_F832,
        IN_F831,
        IN_F830,
        IN_F829,
        IN_F825,
        IN_F822,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F835",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p835 = load_json(IN_P835)
    f834 = load_json(IN_F834)
    f833 = load_json(IN_F833)
    f832 = load_json(IN_F832)
    f831 = load_json(IN_F831)
    f830 = load_json(IN_F830)
    f829 = load_json(IN_F829)
    f825 = load_json(IN_F825)
    f822 = load_json(IN_F822)

    p835_theorem = p835.get("theorem_result") or {}
    f834_target = f834.get("target_object") or {}
    f833_target = f833.get("target_object") or {}
    f832_target = f832.get("target_object") or {}
    f831_target = f831.get("target_object") or {}
    f830_target = f830.get("target_object") or {}
    f829_target = f829.get("target_object") or {}
    f825_target = f825.get("target_object") or {}
    f822_target = f822.get("target_object") or {}

    if (
        p835.get("status")
        == "P835_CURRENT_STRICT_ALPHA_S_STATEMENT_FORM_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_BLOCKED"
        and p835_theorem.get("statement_form_class_candidate_supported_now") is True
        and p835_theorem.get("exact_statement_required_form_exported_now") is False
        and p835_theorem.get("sharp_blocker_field") == "exact_statement_required_form_ref"
        and p835_theorem.get("next_honest_move_is_freeze_exact_required_form_target") is True
        and f834.get("status")
        == "F834_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F835_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F835_REQUIRES_REVIEW"

    artifact = {
        "stage": "F835",
        "packet_name": "CurrentStrictAlphaSLawfulExactRequiredFormStatementDomainAdmissionExactStatementRequiredFormTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p835_required_form_probe": rel(IN_P835),
            "f834_exact_statement_target_packet": rel(IN_F834),
            "f833_output_schema_target_packet": rel(IN_F833),
            "f832_lawful_admission_target_packet": rel(IN_F832),
            "f831_boundary_target_packet": rel(IN_F831),
            "f830_rule_target_packet": rel(IN_F830),
            "f829_interface_target_packet": rel(IN_F829),
            "f825_neighboring_exact_required_form_target_packet": rel(IN_F825),
            "f822_neighboring_schema_output_target_packet": rel(IN_F822),
        },
        "why_this_packet_exists": [
            "F834 already freezes the exact statement object and names one exact missing required-form field.",
            "P835 shows that neighboring statement slots and neighboring required-form targets exist, but the exact required form needed by the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1",
            "goal": "Freeze the exact statement-required form object still missing for the new T213/T216 -> alpha_s exact-required-form-statement lane.",
            "required_fields": [
                {
                    "name": "exact_output_schema_statement_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F834 statement target and not silently replace the problem.",
                },
                {
                    "name": "statement_form_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring form-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_form_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring form slots remain nonidentical to the new lane form.",
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
            "exact_output_schema_statement_target_ref": f834_target.get("object_id"),
            "statement_form_class_candidate_support_refs": [
                f833_target.get("object_id"),
                f832_target.get("object_id"),
                f831_target.get("object_id"),
                f830_target.get("object_id"),
                f829_target.get("object_id"),
                f825_target.get("object_id"),
                f822_target.get("object_id"),
            ],
            "neighboring_form_slot_refs": [
                "exact_output_schema_statement",
                "lawful_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "lawful_schema_domain_admission_exact_required_form_statement_target",
            ],
        },
        "current_honest_reading": [
            "The repo preserves form-like structure around the new lane, but only through neighboring target fields and neighboring required-form targets.",
            "No current export yet names the exact statement-required form required by F834.",
            "F835 freezes that exact missing required-form object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_required_form_statement_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_statement_required_form_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact required form already exists.",
            "Does not claim that any neighboring statement slot or neighboring required-form target silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F835",
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
