#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P870 = GENERATED / "p870_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_same_lane_exhaustion_boundary_audit_probe.json"
IN_F868 = GENERATED / "f868_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_statement_required_form_target_packet.json"
IN_F869 = GENERATED / "f869_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_target_packet.json"

OUT = GENERATED / "f870_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_continuation_boundary_packet.json"
OUT_SUMMARY = GENERATED / "f870_current_strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_continuation_boundary_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P870, IN_F868, IN_F869]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F870",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p870 = load_json(IN_P870)
    f868 = load_json(IN_F868)
    f869 = load_json(IN_F869)

    p870_theorem = p870.get("theorem_result") or {}
    f868_target = f868.get("target_object") or {}
    f869_target = f869.get("target_object") or {}

    if (
        p870.get("status")
        == "P870_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_SAME_LANE_EXHAUSTION_BOUNDARY_AUDIT_PROBE"
        and p870_theorem.get("exact_required_form_statement_exported_on_current_repo_state") is False
        and p870_theorem.get("same_lane_passive_groundwork_exhausted") is True
        and p870_theorem.get("next_honest_move_requires_continuation_boundary_export") is True
        and f868.get("status")
        == "F868_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_STATEMENT_REQUIRED_FORM_TARGET_PACKET_NO_FALSE_PASS"
        and f869.get("status")
        == "F869_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F870_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_EXACT_REQUIRED_FORM_STATEMENT_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
    else:
        status = "F870_REQUIRES_REVIEW"

    artifact = {
        "stage": "F870",
        "packet_name": "CurrentStrictAlphaSLawfulRefinedExactRequiredFormStatementShiftInterfaceRuleDomainAdmissionExactRequiredFormStatementContinuationBoundary_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p870_same_lane_exhaustion_boundary_probe": rel(IN_P870),
            "f868_refined_exact_statement_required_form_target_packet": rel(IN_F868),
            "f869_refined_exact_required_form_statement_target_packet": rel(IN_F869),
        },
        "exported_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_exact_required_form_statement_continuation_boundary_v1",
            "goal": "Export the continuation boundary after the current lawful refined shift-interface-rule same-lane passive groundwork under the missing exact required-form statement is exhausted.",
            "same_lane_exhaustion_boundary_ref": p870_theorem.get("boundary_name"),
            "current_exact_required_form_statement_target_ref": f869_target.get("object_id"),
            "upstream_exact_statement_required_form_target_ref": f868_target.get("object_id"),
            "admitted_next_move_classes": [
                "export_one_genuinely_new_same_lane_exact_required_form_statement_source",
                "shift_to_a_different_required_form_statement_provider_class_lane",
            ],
            "forbidden_repetition_clause": "no_further_same_level_passive_split_under_unchanged_exact_required_form_statement_ref_on_current_lawful_refined_shift_interface_rule_lane_as_primary_move",
            "scope": "strict_alpha_s_lawful_refined_exact_required_form_statement_shift_interface_rule_domain_admission_continuation_boundary_only",
        },
        "current_honest_reading": [
            "The repo now exports an explicit continuation boundary for the current lawful refined shift-interface-rule exact-required-form-statement lane.",
            "The current passive chain has been exported as far as the present honest local split goes on that lawful refined lane.",
            "The next honest move is now constrained to a genuinely new exact required-form statement source or a different provider class, not another fabricated passive split.",
        ],
        "recommended_next_move": "Attack one genuinely new same-lane exact required-form statement source candidate for the current lawful refined shift-interface-rule exact-required-form-statement lane before attempting any further promotion on the current lane.",
        "hard_limits": [
            "Does not claim that the exact required-form statement already exists.",
            "Does not claim that the lawful refined shift-interface-rule lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that any neighboring statement or form slot silently discharges the refined new lane.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F870",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "current_exact_required_form_statement_target_ref": artifact["exported_object"]["current_exact_required_form_statement_target_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
