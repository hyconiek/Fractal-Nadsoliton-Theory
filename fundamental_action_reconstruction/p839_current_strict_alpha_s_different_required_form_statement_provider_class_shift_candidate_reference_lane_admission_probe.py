#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F836 = GENERATED / "f836_current_strict_alpha_s_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_packet.json"
IN_F837 = GENERATED / "f837_current_strict_alpha_s_exact_required_form_statement_continuation_boundary_packet.json"
IN_F838 = GENERATED / "f838_current_strict_alpha_s_exact_required_form_statement_provider_class_shift_requirement_packet.json"
IN_P759 = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"
IN_P760 = GENERATED / "p760_current_strict_t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_audit_probe_summary.json"
IN_P761 = GENERATED / "p761_current_strict_t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_audit_probe_summary.json"
IN_P762 = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"

OUT = GENERATED / "p839_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_lane_admission_probe.json"
OUT_SUMMARY = GENERATED / "p839_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_lane_admission_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def json_lacks_alpha_s_tag(obj: dict[str, Any]) -> bool:
    return "alpha_s" not in json.dumps(obj, sort_keys=True)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F836, IN_F837, IN_F838, IN_P759, IN_P760, IN_P761, IN_P762, IN_P763, IN_P764]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P839",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f836 = load_json(IN_F836)
    f837 = load_json(IN_F837)
    f838 = load_json(IN_F838)
    p759 = load_json(IN_P759)
    p760 = load_json(IN_P760)
    p761 = load_json(IN_P761)
    p762 = load_json(IN_P762)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)

    f836_target = f836.get("target_object") or {}
    f837_export = f837.get("exported_object") or {}
    f838_export = f838.get("exported_object") or {}
    p762_attempt = p762.get("first_actual_t213_realization_attempt") or {}

    checks = [
        {
            "id": "f837_allows_shift_to_different_required_form_statement_provider_class_lane",
            "pass": (
                f837.get("status")
                == "F837_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
                and "shift_to_a_different_required_form_statement_provider_class_lane"
                in (f837_export.get("admitted_next_move_classes") or [])
            ),
            "details": "F837 already exports that a different required-form-statement provider-class lane is an admitted continuation class.",
        },
        {
            "id": "f838_freezes_required_form_statement_provider_class_shift_requirement",
            "pass": (
                f838.get("status")
                == "F838_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET_NO_FALSE_PASS"
                and f838_export.get("remaining_admitted_move_class")
                == "shift_to_a_different_required_form_statement_provider_class_lane"
            ),
            "details": "F838 already freezes required-form-statement provider-class shift as the remaining continuation class on the current alpha_s lane.",
        },
        {
            "id": "f836_still_freezes_the_exact_downstream_blocker",
            "pass": (
                f836.get("status")
                == "F836_EXECUTED_CURRENT_STRICT_ALPHA_S_LAWFUL_EXACT_REQUIRED_FORM_STATEMENT_DOMAIN_ADMISSION_EXACT_REQUIRED_FORM_STATEMENT_TARGET_PACKET_NO_FALSE_PASS"
                and f836_target.get("object_id")
                == "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_target_v1"
            ),
            "details": "F836 still freezes the exact downstream blocker for which any lawful provider-class shift would have to supply a route.",
        },
        {
            "id": "p759_exports_genuinely_new_source_side_branch_selection_provider_target",
            "pass": (
                p759.get("status") == "PASS_STRICT_T213_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_TARGET_EXPORTED"
                and p759.get("t213_target_exported_on_current_repo_state") is True
                and p759.get("current_t213_target_is_source_side_observer_free") is True
                and p759.get("current_t213_target_is_pair12_typed_and_branch_sensitive") is True
                and p759.get("current_t213_target_is_chart_sensitive_and_residual_datum_carrier_binding") is True
                and p759.get("current_t213_target_is_external_to_current_exported_p731_continuation_family")
                is True
                and p759.get("current_t213_target_is_future_route_only") is True
            ),
            "details": "P759 exports a genuinely new source-side branch-selection provider target external to the exhausted current family.",
        },
        {
            "id": "p760_keeps_actual_realization_nonexport",
            "pass": (
                p760.get("status")
                == "PASS_STRICT_T214_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
                and p760.get("t214_target_exported_on_current_repo_state") is False
                and p760.get("current_repo_still_does_not_export_actual_realization_of_t213_target") is True
            ),
            "details": "P760 keeps actual realization of the new provider lane explicitly below export grade.",
        },
        {
            "id": "p761_and_p762_place_the_lane_into_an_active_own_attempt_without_success_claim",
            "pass": (
                p761.get("status")
                == "PASS_STRICT_T215_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_DIRECTION_ACTIVATION_BOUNDARY_AUDIT_PROBE"
                or p761.get("status")
                == "PASS_STRICT_T215_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_DIRECTION_ACTIVATION_BOUNDARY_AUDITED"
            ) and (
                p761.get("actual_t213_realization_direction_is_now_active_primary_t173_branch_on_current_repo_state")
                is True
                and p762.get("status")
                == "PASS_STRICT_T216_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
                and p762.get("t216_attempt_exported_on_current_repo_state") is True
                and p762.get("next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt")
                is True
                and p762.get("first_actual_t213_realization_attempt_keeps_success_failure_open") is True
                and bool(p762_attempt.get("immediate_missing_interface"))
            ),
            "details": "P761/P762 make the lane active on its own terms, but still only as an attempt with open success/failure.",
        },
        {
            "id": "p763_and_p764_keep_the_exact_missing_interface_explicit_and_future_only",
            "pass": (
                p763.get("status")
                == "PASS_STRICT_T217_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_INTERFACE_NONEXPORT_AUDITED"
                and p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported") is True
                and p763.get("current_t216_attempt_stalls_exactly_at_the_named_missing_interface") is True
                and p764.get("status")
                == "PASS_STRICT_T218_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_EXPORTED"
                and p764.get("t218_target_exported_on_current_repo_state") is True
                and p764.get("current_t218_target_is_future_route_only") is True
                and p764.get("current_t218_target_freezes_exact_t216_immediate_missing_interface") is True
            ),
            "details": "P763/P764 keep the own-lane missing interface exact and explicit, which blocks any false reading of realized success.",
        },
        {
            "id": "the_t213_t216_lane_remains_self_contained_and_not_yet_alpha_s_interfaced",
            "pass": (
                all(json_lacks_alpha_s_tag(obj) for obj in [p759, p760, p761, p762, p763, p764])
                and p764.get("next_honest_move_is_actual_export_of_frozen_exact_missing_interface_or_attempt_level_failure_boundary")
                is True
            ),
            "details": "The T213/T216 lane remains self-contained on its own provider class and still exports no alpha_s-side shift interface or exact required-form statement.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "different_required_form_statement_provider_class_shift_candidate_reference_lane_admitted": all_pass,
        "different_required_form_statement_provider_class_shift_candidate_reference_lane_grade": "reference_context_candidate_only"
        if all_pass
        else None,
        "candidate_reference_lane_id": "pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_candidate_reference_lane_v1"
        if all_pass
        else None,
        "candidate_reference_lane_has_own_active_attempt": True if all_pass else None,
        "alpha_s_shift_interface_exported": False if all_pass else None,
        "exact_required_form_statement_exported": False if all_pass else None,
        "provider_class_shift_realized": False if all_pass else None,
        "next_honest_move_requires_exact_alpha_s_shift_interface_audit_for_different_required_form_statement_provider_class_lane": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P839_CURRENT_STRICT_ALPHA_S_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_ALPHA_S_SHIFT_INTERFACE_BLOCKED"
        if all_pass
        else "P839_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P839",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f836_exact_required_form_statement_target_packet": rel(IN_F836),
            "f837_continuation_boundary_packet": rel(IN_F837),
            "f838_provider_class_shift_requirement_packet": rel(IN_F838),
            "p759_t213_target_probe_summary": rel(IN_P759),
            "p760_t214_actual_realization_nonexport_summary": rel(IN_P760),
            "p761_t215_activation_boundary_summary": rel(IN_P761),
            "p762_t216_actual_attempt_summary": rel(IN_P762),
            "p763_t217_missing_interface_boundary_summary": rel(IN_P763),
            "p764_t218_missing_interface_target_summary": rel(IN_P764),
        },
        "theorem_result": theorem_result,
        "candidate_reference_lane": {
            "candidate_id": theorem_result["candidate_reference_lane_id"],
            "provider_target_name": p759.get("t213_target_name"),
            "actual_attempt_name": p762.get("t216_attempt_name"),
            "actual_attempt_shape": p762_attempt.get("attempt_shape"),
            "immediate_missing_interface": p763.get("exact_named_missing_interface"),
            "lane_grade": theorem_result["different_required_form_statement_provider_class_shift_candidate_reference_lane_grade"],
        },
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current alpha_s exact required-form statement lane is exhausted enough that a different provider-class shift is admitted.",
            "The T213/T216 pair12 source-side branch-selection provider lane is real, active on its own lane, and external to the exhausted current family.",
            "But that lane still has only its own attempt and own missing interface; it does not yet export an alpha_s-side shift interface or exact required-form statement for the current blocker.",
        ],
        "recommended_next_packet": {
            "id": "F839_CURRENT_STRICT_ALPHA_S_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET",
            "goal": "Export the T213/T216 lane only as a different required-form-statement provider-class shift candidate reference lane for the current alpha_s exact required-form statement problem.",
            "export_object_id": "alpha_s_pair12_source_side_branch_selection_provider_lawful_exact_required_form_statement_domain_admission_exact_required_form_statement_shift_candidate_reference_lane_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P839",
        "status": status,
        "as_of": AS_OF,
        "different_required_form_statement_provider_class_shift_candidate_reference_lane_admitted": theorem_result[
            "different_required_form_statement_provider_class_shift_candidate_reference_lane_admitted"
        ],
        "different_required_form_statement_provider_class_shift_candidate_reference_lane_grade": theorem_result[
            "different_required_form_statement_provider_class_shift_candidate_reference_lane_grade"
        ],
        "alpha_s_shift_interface_exported": theorem_result["alpha_s_shift_interface_exported"],
        "provider_class_shift_realized": theorem_result["provider_class_shift_realized"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
