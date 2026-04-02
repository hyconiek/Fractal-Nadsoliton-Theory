#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P834 = GENERATED / "p834_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_F833 = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F831 = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F829 = GENERATED / "f829_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F822 = GENERATED / "f822_current_strict_alpha_s_lawful_schema_domain_admission_output_schema_target_packet.json"
IN_F825 = GENERATED / "f825_current_strict_alpha_s_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f834_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f834_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P834, IN_F833, IN_F831, IN_F830, IN_F829, IN_F822, IN_F825]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F834",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p834 = load_json(IN_P834)
    f833 = load_json(IN_F833)
    f831 = load_json(IN_F831)
    f830 = load_json(IN_F830)
    f829 = load_json(IN_F829)
    f822 = load_json(IN_F822)
    f825 = load_json(IN_F825)

    p834_theorem = p834.get("theorem_result") or {}
    f833_target = f833.get("target_object") or {}
    f831_target = f831.get("target_object") or {}
    f830_target = f830.get("target_object") or {}
    f829_target = f829.get("target_object") or {}
    f822_target = f822.get("target_object") or {}
    f825_target = f825.get("target_object") or {}

    if (
        p834.get("status")
        == "P834_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        and p834_theorem.get("output_schema_statement_class_candidate_supported_now") is True
        and p834_theorem.get("exact_output_schema_statement_exported_now") is False
        and p834_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
        and p834_theorem.get("next_honest_move_is_freeze_exact_statement_target") is True
        and f833.get("status")
        == "F833_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F834_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F834_REQUIRES_REVIEW"

    artifact = {
        "stage": "F834",
        "packet_name": "CurrentStrictAlphaSLawfulExactRequiredFormStatementDomainAdmissionExactOutputSchemaStatementTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p834_exact_statement_probe": rel(IN_P834),
            "f833_output_schema_target_packet": rel(IN_F833),
            "f831_boundary_target_packet": rel(IN_F831),
            "f830_rule_target_packet": rel(IN_F830),
            "f829_interface_target_packet": rel(IN_F829),
            "f822_neighboring_schema_output_target_packet": rel(IN_F822),
            "f825_downstream_exact_required_form_statement_target_packet": rel(IN_F825),
        },
        "why_this_packet_exists": [
            "F833 already freezes the exact output-schema object and names one exact missing statement field.",
            "P834 shows that neighboring output-schema statement slots exist, but the exact statement needed by the new lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1",
            "goal": "Freeze the exact output-schema statement object still missing for the new T213/T216 -> alpha_s exact-required-form-statement lane.",
            "required_fields": [
                {
                    "name": "lawful_exact_required_form_statement_domain_admission_output_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F833 output-schema target and not silently replace the problem.",
                },
                {
                    "name": "statement_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring statement-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_statement_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring statement slots remain nonidentical to the new lane statement.",
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
            "lawful_exact_required_form_statement_domain_admission_output_schema_target_ref": f833_target.get("object_id"),
            "statement_class_candidate_support_refs": [
                f831_target.get("object_id"),
                f830_target.get("object_id"),
                f829_target.get("object_id"),
                f822_target.get("object_id"),
                f825_target.get("object_id"),
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
            "The repo preserves statement-level output slots around the new lane, but only as neighboring target fields.",
            "No current export yet names the exact output-schema statement required by F833.",
            "F834 freezes that exact missing statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_statement_required_form_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact output-schema statement already exists.",
            "Does not claim that any neighboring output-schema statement slot silently discharges the new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F834",
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
