#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P828 = GENERATED / "p828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_lane_admission_probe.json"
IN_F826 = GENERATED / "f826_current_strict_alpha_s_exact_required_form_statement_continuation_boundary_packet.json"
IN_F827 = GENERATED / "f827_current_strict_alpha_s_exact_required_form_statement_provider_class_shift_requirement_packet.json"
IN_P759 = GENERATED / "p759_current_strict_t213_pair12_source_side_branch_selection_provider_target_probe_summary.json"
IN_P761 = GENERATED / "p761_current_strict_t215_pair12_source_side_branch_selection_provider_actual_realization_direction_activation_boundary_audit_probe_summary.json"
IN_P762 = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"

OUT = GENERATED / "f828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_packet.json"
OUT_SUMMARY = GENERATED / "f828_current_strict_alpha_s_different_required_form_statement_provider_class_shift_candidate_reference_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P828, IN_F826, IN_F827, IN_P759, IN_P761, IN_P762, IN_P763, IN_P764]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F828",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p828 = load_json(IN_P828)
    f826 = load_json(IN_F826)
    f827 = load_json(IN_F827)
    p759 = load_json(IN_P759)
    p761 = load_json(IN_P761)
    p762 = load_json(IN_P762)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)

    p828_theorem = p828.get("theorem_result") or {}
    f826_export = f826.get("exported_object") or {}
    f827_export = f827.get("exported_object") or {}

    if (
        p828.get("status")
        == "P828_CURRENT_STRICT_ALPHA_S_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_ALPHA_S_SHIFT_INTERFACE_BLOCKED"
        and p828_theorem.get("different_required_form_statement_provider_class_shift_candidate_reference_lane_admitted")
        is True
        and p828_theorem.get("different_required_form_statement_provider_class_shift_candidate_reference_lane_grade")
        == "reference_context_candidate_only"
        and p828_theorem.get("alpha_s_shift_interface_exported") is False
        and p828_theorem.get("provider_class_shift_realized") is False
        and f826.get("status")
        == "F826_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
        and f827.get("status")
        == "F827_EXECUTED_CURRENT_STRICT_ALPHA_S_EXACT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET_NO_FALSE_PASS"
    ):
        status = "F828_EXECUTED_CURRENT_STRICT_ALPHA_S_DIFFERENT_REQUIRED_FORM_STATEMENT_PROVIDER_CLASS_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
    else:
        status = "F828_REQUIRES_REVIEW"

    artifact = {
        "stage": "F828",
        "packet_name": "CurrentStrictAlphaSDifferentRequiredFormStatementProviderClassShiftCandidateReference_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p828_shift_candidate_reference_lane_probe": rel(IN_P828),
            "f826_continuation_boundary_packet": rel(IN_F826),
            "f827_provider_class_shift_requirement_packet": rel(IN_F827),
            "p759_t213_target_probe_summary": rel(IN_P759),
            "p761_t215_activation_boundary_summary": rel(IN_P761),
            "p762_t216_actual_attempt_summary": rel(IN_P762),
            "p763_t217_missing_interface_boundary_summary": rel(IN_P763),
            "p764_t218_missing_interface_target_summary": rel(IN_P764),
        },
        "exported_object": {
            "object_id": "alpha_s_pair12_source_side_branch_selection_provider_exact_required_form_statement_shift_candidate_reference_lane_v1",
            "goal": "Export the T213/T216 lane only as a different required-form-statement provider-class shift candidate reference lane for the alpha_s exact required-form statement problem.",
            "continuation_boundary_ref": f826_export.get("object_id"),
            "provider_class_shift_requirement_ref": f827_export.get("object_id"),
            "candidate_reference_lane_ref": p828_theorem.get("candidate_reference_lane_id"),
            "candidate_reference_lane_grade": p828_theorem.get(
                "different_required_form_statement_provider_class_shift_candidate_reference_lane_grade"
            ),
            "provider_target_name": p759.get("t213_target_name"),
            "own_lane_activation_boundary_ref": p761.get("t215_boundary_name"),
            "own_lane_actual_attempt_ref": p762.get("t216_attempt_name"),
            "own_lane_missing_interface_boundary_ref": p763.get("t217_boundary_name"),
            "own_lane_missing_interface_target_ref": p764.get("t218_target_name"),
            "alpha_s_shift_interface_status": "blocked_nonexport",
            "exact_required_form_statement_status": "blocked_nonexport",
            "provider_class_shift_realization_status": "not_realized",
        },
        "current_honest_reading": [
            "The repo now exports a different provider-class shift candidate reference lane for the alpha_s exact required-form statement problem.",
            "This lane is more than a generic hint: it already has its own future-only target, active own-lane attempt, and exact missing interface target.",
            "But it still does not export the exact alpha_s-side shift interface, so the correct reading remains candidate-reference only.",
        ],
        "recommended_next_move": "Audit whether the T213/T216 lane exports any exact alpha_s-side shift interface or adapter target for the exact required-form statement problem without silent domain identification.",
        "hard_limits": [
            "Does not claim that the exact required-form statement already exists.",
            "Does not claim that the T213/T216 lane already enters the alpha_s domain.",
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim that any exact alpha_s source binding already exists.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F828",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "candidate_reference_lane_grade": artifact["exported_object"]["candidate_reference_lane_grade"],
        "alpha_s_shift_interface_status": artifact["exported_object"]["alpha_s_shift_interface_status"],
        "provider_class_shift_realization_status": artifact["exported_object"][
            "provider_class_shift_realization_status"
        ],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
