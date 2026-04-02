#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P845 = GENERATED / "p845_current_strict_alpha_s_output_schema_statement_class_candidate_supported_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_blocked.json"
IN_F844 = GENERATED / "f844_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F842 = GENERATED / "f842_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
IN_F841 = GENERATED / "f841_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_F840 = GENERATED / "f840_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_interface_target_packet.json"
IN_F833 = GENERATED / "f833_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_output_schema_target_packet.json"
IN_F836 = GENERATED / "f836_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f845_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet.json"
OUT_SUMMARY = GENERATED / "f845_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_output_schema_statement_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P845, IN_F844, IN_F842, IN_F841, IN_F840, IN_F833, IN_F836]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F845",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p845 = load_json(IN_P845)
    f844 = load_json(IN_F844)
    f842 = load_json(IN_F842)
    f841 = load_json(IN_F841)
    f840 = load_json(IN_F840)
    f833 = load_json(IN_F833)
    f836 = load_json(IN_F836)

    p845_theorem = p845.get("theorem_result") or {}
    f844_target = f844.get("target_object") or {}
    f842_target = f842.get("target_object") or {}
    f841_target = f841.get("target_object") or {}
    f840_target = f840.get("target_object") or {}
    f833_target = f833.get("target_object") or {}
    f836_target = f836.get("target_object") or {}

    if (
        p845.get("status")
        == "P845_CURRENT_STRICT_ALPHA_S_OUTPUT_SCHEMA_STATEMENT_CLASS_CANDIDATE_SUPPORTED_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_BLOCKED"
        and p845_theorem.get("output_schema_statement_class_candidate_supported_now") is True
        and p845_theorem.get("exact_output_schema_statement_exported_now") is False
        and p845_theorem.get("sharp_blocker_field") == "exact_output_schema_statement"
        and p845_theorem.get("next_honest_move_is_freeze_exact_statement_target") is True
        and f844.get("status")
        == "F844_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F845_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_OUTPUT_SCHEMA_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F845_REQUIRES_REVIEW"

    artifact = {
        "stage": "F845",
        "packet_name": "CurrentStrictAlphaSLawfulExactRequiredFormStatementDomainAdmissionExactOutputSchemaStatementTargetRefined_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p845_exact_statement_probe": rel(IN_P845),
            "f844_refined_output_schema_target_packet": rel(IN_F844),
            "f842_boundary_target_packet": rel(IN_F842),
            "f841_rule_target_packet": rel(IN_F841),
            "f840_interface_target_packet": rel(IN_F840),
            "f833_old_lawful_output_target_packet": rel(IN_F833),
            "f836_downstream_exact_required_form_statement_target_packet": rel(IN_F836),
        },
        "why_this_packet_exists": [
            "F844 already freezes the refined lawful output-schema object and names one exact missing statement field.",
            "P845 shows that neighboring output-schema statement slots exist, but the exact statement needed by the refined lane is still unexported.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_output_schema_statement_target_v1",
            "goal": "Freeze the exact output-schema statement object still missing for the refined lawful T213/T216 -> alpha_s exact-required-form-statement lane.",
            "required_fields": [
                {
                    "name": "lawful_exact_required_form_statement_domain_admission_output_schema_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F844 refined output-schema target and not silently replace the problem.",
                },
                {
                    "name": "statement_class_candidate_support_refs",
                    "required": True,
                    "hard_limit": "Must preserve only candidate-grade neighboring statement-class support and must not promote it into exact discharge.",
                },
                {
                    "name": "neighboring_statement_slot_refs",
                    "required": True,
                    "hard_limit": "Must keep explicit which neighboring statement slots remain nonidentical to the refined new-lane statement.",
                },
                {
                    "name": "exact_statement_required_form_ref",
                    "required": True,
                    "hard_limit": "Must state what exact statement form is required for the refined new lane without claiming that the statement already exists.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny exact new-lane admission, silent reuse of neighboring statement slots, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "target_refs": {
            "lawful_exact_required_form_statement_domain_admission_output_schema_target_ref": f844_target.get("object_id"),
            "statement_class_candidate_support_refs": [
                f842_target.get("object_id"),
                f841_target.get("object_id"),
                f840_target.get("object_id"),
                f833_target.get("object_id"),
                f836_target.get("object_id"),
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
            "The repo preserves statement-level output slots around the refined lawful lane, but only as neighboring target fields.",
            "No current export yet names the exact output-schema statement required by F844.",
            "F845 freezes that exact missing statement object without pretending that lawful admission already exists.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply exact_statement_required_form_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_refined_exact_output_schema_statement_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that exact output-schema statement already exists.",
            "Does not claim that any neighboring output-schema statement slot silently discharges the refined new lane.",
            "Does not claim that the T213/T216 lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F845",
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
