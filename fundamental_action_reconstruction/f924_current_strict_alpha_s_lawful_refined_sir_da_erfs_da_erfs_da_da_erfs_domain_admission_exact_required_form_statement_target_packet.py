#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P924 = GENERATED / "p924_current_strict_alpha_s_required_form_statement_class_candidate_supported_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_blocked.json"
IN_F923 = GENERATED / "f923_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_statement_required_form_target_packet.json"
IN_F922 = GENERATED / "f922_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_output_schema_statement_target_packet.json"
IN_F921 = GENERATED / "f921_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_output_schema_target_packet.json"
IN_F920 = GENERATED / "f920_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_target_packet.json"
IN_F919 = GENERATED / "f919_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F918 = GENERATED / "f918_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_adapter_or_carrier_rule_target_packet.json"
IN_F917 = GENERATED / "f917_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_interface_target_packet.json"
IN_F913 = GENERATED / "f913_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f924_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f924_current_strict_alpha_s_lawful_refined_sir_da_erfs_da_erfs_da_da_erfs_domain_admission_exact_required_form_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P924, IN_F923, IN_F922, IN_F921, IN_F920, IN_F919, IN_F918, IN_F917, IN_F913]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F924",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p924 = load_json(IN_P924)
    f923 = load_json(IN_F923)
    f922 = load_json(IN_F922)
    f921 = load_json(IN_F921)
    f920 = load_json(IN_F920)
    f919 = load_json(IN_F919)
    f918 = load_json(IN_F918)
    f917 = load_json(IN_F917)
    f913 = load_json(IN_F913)

    p924_theorem = p924.get("theorem_result") or {}
    f923_target = f923.get("target_object") or {}
    f922_target = f922.get("target_object") or {}
    f921_target = f921.get("target_object") or {}
    f920_target = f920.get("target_object") or {}
    f919_target = f919.get("target_object") or {}
    f918_target = f918.get("target_object") or {}
    f917_target = f917.get("target_object") or {}
    f913_target = f913.get("target_object") or {}

    if (
        p924.get("status")
        == "P924_CURRENT_STRICT_ALPHA_S_REQUIRED_FORM_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_BLOCKED"
        and p924_theorem.get("required_form_statement_class_candidate_supported_now") is True
        and p924_theorem.get("exact_required_form_statement_exported_now") is False
        and p924_theorem.get("sharp_blocker_field") == "exact_required_form_statement_ref"
        and p924_theorem.get("next_honest_move_is_freeze_exact_required_form_statement_target") is True
        and f923.get("status")
        == "F923_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F924_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SIR_DA_ERFS_DA_ERFS_DA_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F924_REQUIRES_REVIEW"

    artifact = {
        "stage": "F924",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedSirDaErfsDaErfsDaDaErfsDomainAdmissionExactRequiredFormStatementTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p924_required_form_statement_probe": rel(IN_P924),
            "f923_lawful_refined_deeper_domain_admission_exact_statement_required_form_target_packet": rel(IN_F923),
            "f922_lawful_refined_deeper_domain_admission_exact_output_schema_statement_target_packet": rel(IN_F922),
            "f921_lawful_refined_deeper_domain_admission_output_schema_target_packet": rel(IN_F921),
            "f920_lawful_refined_deeper_domain_admission_target_packet": rel(IN_F920),
            "f919_boundary_target_packet": rel(IN_F919),
            "f918_rule_target_packet": rel(IN_F918),
            "f917_interface_target_packet": rel(IN_F917),
            "f913_neighboring_exact_required_form_statement_target_packet": rel(IN_F913),
        },
        "why_this_packet_exists": [
            "F923 already freezes the deeper lawful refined domain-admission exact-statement-required-form object and names one exact missing required-form-statement field.",
            "P924 shows that neighboring statement slots and neighboring required-form supports exist, but the exact required-form statement needed by the new deeper lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1",
            "goal": "Freeze the exact deeper lawful refined domain-admission required-form statement object still missing for the lawful refined shift-interface-rule T213/T216 -> alpha_s lane.",
            "required_fields": [
                {
                    "name": "exact_statement_required_form_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F923 required-form target and not silently replace the problem.",
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
            "exact_statement_required_form_target_ref": f923_target.get("object_id"),
            "required_form_statement_class_candidate_support_refs": [
                f922_target.get("object_id"),
                f921_target.get("object_id"),
                f920_target.get("object_id"),
                f919_target.get("object_id"),
                f918_target.get("object_id"),
                f917_target.get("object_id"),
                f913_target.get("object_id"),
            ],
            "neighboring_statement_or_form_slot_refs": [
                "exact_statement_required_form_ref",
                "exact_output_schema_statement",
                "lawful_refined_deeper_exact_required_form_statement_domain_admission_output_schema",
                "boundary_output_schema",
                "selected_interface_output_schema",
                "exact_interface_output_schema",
                "exact_required_form_statement",
            ],
        },
        "current_honest_reading": [
            "The repo preserves statement and form-like structure around the deeper lawful refined domain-admission lane, but only through neighboring target fields and neighboring old-lane targets.",
            "No current export yet names the exact required-form statement required by F923 for the lawful refined T213/T216 -> alpha_s route.",
            "F924 freezes that exact missing required-form-statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow same-lane exhaustion-boundary audit testing whether any further passive split remains below alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact deeper lawful refined required-form statement already exists.",
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
        "stage": "F924",
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
