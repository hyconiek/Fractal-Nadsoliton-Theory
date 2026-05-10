#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P831 = GENERATED / "p831_current_strict_alpha_s_action_schema_candidate_supported_exact_required_form_statement_domain_admission_boundary_blocked.json"
IN_F830 = GENERATED / "f830_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_packet.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_F819 = GENERATED / "f819_current_strict_alpha_s_pair12_source_side_branch_selection_provider_shift_to_schema_adapter_or_carrier_identification_rule_target_packet.json"

OUT = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet.json"
OUT_SUMMARY = GENERATED / "f831_current_strict_alpha_s_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P831, IN_F830, IN_P764, IN_P788, IN_F819]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F831",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p831 = load_json(IN_P831)
    f830 = load_json(IN_F830)
    p764 = load_json(IN_P764)
    p788 = load_json(IN_P788)
    f819 = load_json(IN_F819)

    split = p831.get("clause_split_audit") or {}
    f830_target = f830.get("target_object") or {}
    f830_refs = f830.get("target_refs") or {}
    p764_theorem = p764.get("theorem_result") or {}
    f819_target = f819.get("target_object") or {}

    if (
        p831.get("status")
        == "P831_CURRENT_STRICT_ALPHA_S_ACTION_SCHEMA_CANDIDATE_SUPPORTED_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_BOUNDARY_BLOCKED"
        and split.get("sharp_blocker_clause") == "exact_required_form_statement_domain_admission_boundary"
        and f830.get("status")
        == "F830_EXECUTED_CURRENT_STRICT_ALPHA_S_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_SHIFT_TO_EXACT_REQUIRED_FORM_STATEMENT_ADAPTER_OR_CARRIER_IDENTIFICATION_RULE_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F831_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_OR_NONIDENTIFICATION_BOUNDARY_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F831_REQUIRES_REVIEW"

    artifact = {
        "stage": "F831",
        "packet_name": "CurrentStrictAlphaSExactRequiredFormStatementDomainAdmissionOrNonidentificationBoundaryTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p831_clause_split_audit": rel(IN_P831),
            "f830_rule_target_packet": rel(IN_F830),
            "p764_own_lane_interface_target_probe": rel(IN_P764),
            "p788_generic_alpha_s_adapter_probe": rel(IN_P788),
            "f819_old_same_lane_schema_rule_target_packet": rel(IN_F819),
        },
        "why_this_packet_exists": [
            "F830 freezes the larger adapter-or-carrier-identification rule target and includes both action-schema and exact-required-form-statement-domain-admission clauses.",
            "P831 shows that the sharper blocker is lawful exact-required-form-statement domain admission / nonidentification, not generic action-schema shape.",
        ],
        "target_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_domain_admission_or_nonidentification_boundary_target_v1",
            "goal": "Freeze the exact exact-required-form-statement domain admission or nonidentification boundary object still missing before the new lane can move closer to an export-grade rule.",
            "required_fields": [
                {
                    "name": "shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact F830 rule target and not silently replace it.",
                },
                {
                    "name": "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact admitted T213/T216 candidate lane and not silently replace it.",
                },
                {
                    "name": "candidate_action_schema_support_ref",
                    "required": True,
                    "hard_limit": "Must point only to candidate-supported action-schema support, not to a claimed exported action-schema.",
                },
                {
                    "name": "downstream_exact_required_form_statement_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact downstream exact-required-form-statement target and not silently replace it.",
                },
                {
                    "name": "generic_alpha_s_adapter_absence_ref",
                    "required": True,
                    "hard_limit": "Must keep explicit that no generic alpha_s adapter is exported on current repo state.",
                },
                {
                    "name": "old_same_lane_schema_nontransfer_ref",
                    "required": True,
                    "hard_limit": "Must explicitly deny silent reuse of the old F819 same-lane schema-side rule target as if it already discharged the new exact-required-form-statement boundary.",
                },
                {
                    "name": "boundary_output_schema",
                    "required": True,
                    "hard_limit": "Must state what lawful admission or lawful continued nonidentification would output for the F830 lane.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny realized domain entry, provider-class shift success, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The current sharp blocker is no longer generic action-schema shape.",
            "It is the missing exact-required-form-statement domain admission / nonidentification boundary for the new T213/T216 -> alpha_s route.",
            "F831 freezes that exact missing object without promoting the current lane into the alpha_s exact-required-form-statement domain.",
        ],
        "target_refs": {
            "shift_to_exact_required_form_statement_adapter_or_carrier_identification_rule_target_ref": f830_target.get("object_id"),
            "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref": f830_refs.get(
                "different_required_form_statement_provider_class_shift_candidate_reference_lane_ref"
            ),
            "candidate_action_schema_support_ref": p764_theorem.get("t218_target_name"),
            "downstream_exact_required_form_statement_target_ref": f830_refs.get(
                "downstream_exact_required_form_statement_target_ref"
            ),
            "generic_alpha_s_adapter_absence_ref": "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE",
            "old_same_lane_schema_nontransfer_ref": f819_target.get("object_id"),
        },
        "recommended_next_move": "Build one narrow probe testing whether lawful exact-required-form-statement domain admission or continued nonidentification can be justified from current exports without silent domain identification.",
        "hard_limits": [
            "Does not claim that the action-schema is already exported.",
            "Does not claim that the exact-required-form-statement domain admission boundary already exists.",
            "Does not claim that the T213/T216 lane already enters the alpha_s exact-required-form-statement domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F831",
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
