#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-20"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P884 = GENERATED / "p884_current_strict_alpha_s_no_exact_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_interface_target_freeze_required.json"
IN_F883 = GENERATED / "f883_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_different_required_form_statement_provider_class_shift_candidate_reference_packet.json"
IN_F880 = GENERATED / "f880_current_strict_alpha_s_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"

OUT = GENERATED / "f884_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_interface_target_packet.json"
OUT_SUMMARY = GENERATED / "f884_current_strict_alpha_s_pair12_provider_shift_to_lawful_refined_sir_da_exact_required_form_statement_da_exact_required_form_statement_interface_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P884, IN_F883, IN_F880, IN_P764]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F884",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p884 = load_json(IN_P884)
    f883 = load_json(IN_F883)
    f880 = load_json(IN_F880)
    p764 = load_json(IN_P764)

    p884_theorem = p884.get("theorem_result") or {}
    missing_target = p884.get("exact_missing_interface_target_candidate") or {}
    f883_export = f883.get("exported_object") or {}
    f880_target = f880.get("target_object") or {}

    if (
        p884.get("status")
        == "P884_CURRENT_STRICT_ALPHA_S_NO_EXACT_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_EXACT_REQUIRED_FORM_STATEMENT_DA_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_FREEZE_REQUIRED"
        and p884_theorem.get("exact_pair12_provider_shift_to_lawful_refined_domain_admission_exact_required_form_statement_interface_target_exported")
        is False
        and p884_theorem.get("next_honest_move_requires_freezing_exact_shift_to_lawful_refined_domain_admission_exact_required_form_statement_interface_target")
        is True
        and f883.get("status")
        == "F883_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
        and f880.get("status")
        == "F880_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_REFINED_SHIFT_INTERFACE_RULE_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F884_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_PROVIDER_SHIFT_TO_LAWFUL_REFINED_SIR_DA_EXACT_REQUIRED_FORM_STATEMENT_DA_EXACT_REQUIRED_FORM_STATEMENT_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F884_REQUIRES_REVIEW"

    artifact = {
        "stage": "F884",
        "packet_name": "CurrentStrictAlphaSPair12ProviderShiftToLawfulRefinedSIRDAExactRequiredFormStatementDAExactRequiredFormStatementInterfaceTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p884_shift_to_lawful_refined_domain_admission_exact_required_form_statement_interface_audit_probe": rel(IN_P884),
            "f883_domain_admission_different_provider_class_shift_candidate_reference_packet": rel(IN_F883),
            "f880_domain_admission_downstream_exact_required_form_statement_target_packet": rel(IN_F880),
            "p764_own_lane_missing_interface_target_summary": rel(IN_P764),
        },
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_interface_target_v1",
            "goal": "Freeze the exact missing interface target from the admitted T213/T216 provider-class shift candidate lane into the current lawful refined alpha_s shift-interface-rule domain-admission exact-required-form-statement domain-admission exact-required-form-statement problem without claiming interface realization.",
            "required_fields": [
                {
                    "name": "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F883 candidate-reference lane and not silently replace it."
                },
                {
                    "name": "downstream_lawful_refined_domain_admission_exact_required_form_statement_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F880 downstream target and not silently replace the problem."
                },
                {
                    "name": "self_lane_missing_interface_target_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that the candidate lane still has its own unresolved T218 interface target and cannot silently bypass it."
                },
                {
                    "name": "shift_interface_adapter_or_carrier_identification_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the exact adapter or carrier-safe identification rule; silent domain identification is forbidden."
                },
                {
                    "name": "nontransfer_boundary_ref",
                    "required": True,
                    "hard_limit": "Must block silent reuse of earlier-lane alpha_s or foreign-domain interface artifacts."
                },
                {
                    "name": "future_route_grade_ref",
                    "required": True,
                    "hard_limit": "Must keep this interface target at future-route target grade until a real interface is exported."
                },
                {
                    "name": "exact_interface_output_schema",
                    "required": True,
                    "hard_limit": "Must state what exact lawful output would enter the current lawful refined alpha_s domain-admission exact-required-form-statement problem if the interface were exported."
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny interface realization, provider-class shift success, QCD closure, and ToE closure."
                }
            ]
        },
        "target_refs": {
            "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref": f883_export.get("object_id"),
            "downstream_lawful_refined_domain_admission_exact_required_form_statement_target_ref": f880_target.get("object_id"),
            "self_lane_missing_interface_target_ref": p764.get("t218_target_name"),
            "missing_target_candidate_id": missing_target.get("candidate_id"),
        },
        "current_honest_reading": [
            "The repo now exports the exact missing lawful refined interface target between the admitted T213/T216 candidate lane and the current lawful refined alpha_s domain-admission exact-required-form-statement problem.",
            "This sits above the admitted lawful refined different-provider-class candidate lane and below any adapter or rule that would instantiate it.",
            "It does not claim that the interface exists; it only localizes the missing lawful refined object sharply.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current export can lawfully supply shift_interface_adapter_or_carrier_identification_rule_ref for alpha_s_pair12_source_side_branch_selection_provider_lawful_refined_shift_interface_rule_domain_admission_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_interface_target_v1 without silent domain identification.",
        "hard_limits": [
            "Does not claim that the exact required-form statement already exists.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
            "Does not claim that any adapter or carrier-identification rule is already exported.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F884",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref": artifact["target_refs"][
            "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref"
        ],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
