#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P933 = GENERATED / "p933_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_output_schema_statement_blocked.json"
IN_F932 = GENERATED / "f932_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet.json"
IN_F930 = GENERATED / "f930_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F929 = GENERATED / "f929_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json"
IN_F928 = GENERATED / "f928_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.json"
IN_F921 = GENERATED / "f921_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet.json"
IN_F924 = GENERATED / "f924_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f933_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_output_schema_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f933_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_output_schema_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P933, IN_F932, IN_F930, IN_F929, IN_F928, IN_F921, IN_F924]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F933",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p933 = load_json(IN_P933)
    f932 = load_json(IN_F932)
    f930 = load_json(IN_F930)
    f929 = load_json(IN_F929)
    f928 = load_json(IN_F928)
    f921 = load_json(IN_F921)
    f924 = load_json(IN_F924)

    p933_theorem = p933.get("theorem_result") or {}
    f932_target = f932.get("target_object") or {}
    f930_target = f930.get("target_object") or {}
    f929_target = f929.get("target_object") or {}
    f928_target = f928.get("target_object") or {}
    f921_target = f921.get("target_object") or {}
    f924_target = f924.get("target_object") or {}

    if (
        p933.get("status")
        == "P933_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        and p933_theorem.get("output_schema_statement_class_candidate_supported_now") is True
        and p933_theorem.get("exact_output_schema_statement_exported_now") is False
        and p933_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
        and p933_theorem.get("next_honest_move_is_freeze_exact_statement_target") is True
        and f932.get("status")
        == "F932_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F933_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F933_REQUIRES_REVIEW"

    artifact = {
        "stage": "F933",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedSirDaErfsDaErfsDaDaErfsDomainAdmissionExactOutputSchemaStatementTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p933_exact_statement_probe": rel(IN_P933),
            "f932_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F932),
            "f930_boundary_target_packet": rel(IN_F930),
            "f929_rule_target_packet": rel(IN_F929),
            "f928_interface_target_packet": rel(IN_F928),
            "f921_neighboring_output_target_packet": rel(IN_F921),
            "f924_neighboring_exact_required_form_statement_target_packet": rel(IN_F924),
        },
        "why_this_packet_exists": [
            "F932 already freezes the deeper lawful refined domain-admission output-schema object and names one exact missing statement field.",
            "P933 shows that neighboring statement slots exist, but the exact statement needed by the new deeper lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1",
            "goal": "Freeze the exact deeper lawful refined domain-admission output-schema statement object still missing for the lawful refined shift-interface-rule T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "lawful_refined_deeper_domain_admission_output_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F932 output-schema target and not silently replace the problem.",
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
            "lawful_refined_deeper_domain_admission_output_schema_target_ref": f932_target.get("object_id"),
            "statement_class_candidate_support_refs": [
                f930_target.get("object_id"),
                f929_target.get("object_id"),
                f928_target.get("object_id"),
                f921_target.get("object_id"),
                f924_target.get("object_id"),
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
            "No current export yet names the exact output-schema statement required by F932.",
            "F933 freezes that exact missing statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_statement_required_form_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_v1 without silent domain identification.",
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
        "stage": "F933",
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
