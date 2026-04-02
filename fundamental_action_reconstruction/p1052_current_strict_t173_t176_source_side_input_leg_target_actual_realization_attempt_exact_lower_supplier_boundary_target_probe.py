#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1050 = GENERATED / "p1050_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_probe_summary.json"
IN_P1051 = GENERATED / "p1051_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_verdict_or_exact_lower_supplier_boundary_nonexport_audit_probe_summary.json"
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_P768 = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"
IN_P769 = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T303 = ROOT / "T303_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1052_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_exact_lower_supplier_boundary_target_probe.json"
OUT_SUMMARY = GENERATED / "p1052_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_exact_lower_supplier_boundary_target_probe_summary.json"

P1050_STATUS = "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
P1051_STATUS = (
    "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_LOWER_SUPPLIER_BOUNDARY_NONEXPORT_AUDITED"
)
P767_STATUS = (
    "PASS_STRICT_T221_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_NONEXPORT_AUDITED"
)
P768_STATUS = (
    "PASS_STRICT_T222_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_EXPORTED"
)
P769_STATUS = (
    "PASS_STRICT_T223_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
)
P770_STATUS = (
    "PASS_STRICT_T224_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
)

T302_ATTEMPT_NAME = (
    "W_strict_t173_t176_inversion_sensitive_pair12_branch_separation_bridge_source_side_input_leg_actual_realization_attempt_v1"
)
T220_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1"
)
T224_ATTEMPT_NAME = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_v1"
)
TARGET_NAME = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_actual_realization_attempt_exact_lower_supplier_boundary_target_v1"
)
EXACT_NAMED_SUBINTERFACE = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1050, IN_P1051, IN_P767, IN_P768, IN_P769, IN_P770, IN_T303]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1052",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1050 = load_json(IN_P1050)
    p1051 = load_json(IN_P1051)
    p767 = load_json(IN_P767)
    p768 = load_json(IN_P768)
    p769 = load_json(IN_P769)
    p770 = load_json(IN_P770)
    t303_text = load_text(IN_T303)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    p1051_nonexport_audit_passed = (
        p1051.get("status") == P1051_STATUS
        and p1051.get("t302_attempt_name") == T302_ATTEMPT_NAME
        and p1051.get("t302_verdict_exported_on_current_repo_state") is False
        and p1051.get("t302_exact_lower_supplier_boundary_exported_on_current_repo_state") is False
        and p1051.get("next_honest_move_is_freeze_exact_lower_supplier_boundary_target_below_t302") is True
    )

    t302_attempt_is_exported = (
        p1050.get("status") == P1050_STATUS
        and p1050.get("attempt_name") == T302_ATTEMPT_NAME
        and p1050.get("attempt_exported_on_current_repo_state") is True
    )

    exact_lower_route_local_support_family_remains_explicit = (
        p767.get("status") == P767_STATUS
        and p767.get("exact_named_missing_subinterface") == EXACT_NAMED_SUBINTERFACE
        and p768.get("status") == P768_STATUS
        and p768.get("t222_target_exported_on_current_repo_state") is True
        and p769.get("status") == P769_STATUS
        and p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target") is True
        and p770.get("status") == P770_STATUS
        and p770.get("t224_attempt_name") == T224_ATTEMPT_NAME
        and p770.get("t224_attempt_exported_on_current_repo_state") is True
    )

    t303_target_shape_is_frozen = all(
        needle in t303_text
        for needle in [
            "T303_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_SPEC_NO_FALSE_PASS",
            TARGET_NAME,
            T302_ATTEMPT_NAME,
            EXACT_NAMED_SUBINTERFACE,
            T220_ATTEMPT_NAME,
            T224_ATTEMPT_NAME,
            "target_must_not_reenter_exhausted_pair12_entry_point_same_lane_as_primary_descent := yes",
        ]
    )

    target_exported_on_current_repo_state = (
        p1051_nonexport_audit_passed
        and t302_attempt_is_exported
        and exact_lower_route_local_support_family_remains_explicit
        and t303_target_shape_is_frozen
    )

    add_check(
        "p1051_nonexport_audit_passed",
        p1051_nonexport_audit_passed,
        True,
        "P1051 already freezes that the exact T302 attempt still has neither lawful verdict nor exact lower supplier-boundary.",
    )
    add_check(
        "t302_attempt_is_exported",
        t302_attempt_is_exported,
        True,
        "P1050 already exports the exact T302 noncyclic source-side input-leg actual-realization attempt.",
    )
    add_check(
        "exact_lower_route_local_support_family_remains_explicit",
        exact_lower_route_local_support_family_remains_explicit,
        True,
        "The older exact T220/T222/T224 lower support family remains explicit around the same named chart-label-retaining pair12 typed seed subinterface.",
    )
    add_check(
        "t303_target_shape_is_frozen",
        t303_target_shape_is_frozen,
        True,
        "T303 freezes one exact lower supplier-boundary target beneath T302 while keeping the noncyclic guardrail explicit.",
    )
    add_check(
        "target_exported_on_current_repo_state",
        target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact lower supplier-boundary target beneath the exact T302 source-side input-leg attempt.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_EXPORTED"
        if not blocking and target_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_LOWER_SUPPLIER_BOUNDARY_TARGET_PROBE"
    )

    artifact = {
        "stage": "P1052",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "t302_attempt_name": T302_ATTEMPT_NAME,
        "target_exported_on_current_repo_state": target_exported_on_current_repo_state,
        "keeps_t302_verdict_open": True,
        "keeps_source_side_input_leg_actual_export_open": True,
        "keeps_full_c_v1_transported_section_lift_open": True,
        "keeps_t176_open": True,
        "uses_noncyclic_guardrail": True,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_name": artifact["target_name"],
        "t302_attempt_name": artifact["t302_attempt_name"],
        "target_exported_on_current_repo_state": artifact["target_exported_on_current_repo_state"],
        "keeps_t302_verdict_open": artifact["keeps_t302_verdict_open"],
        "keeps_source_side_input_leg_actual_export_open": artifact[
            "keeps_source_side_input_leg_actual_export_open"
        ],
        "keeps_full_c_v1_transported_section_lift_open": artifact[
            "keeps_full_c_v1_transported_section_lift_open"
        ],
        "keeps_t176_open": artifact["keeps_t176_open"],
        "uses_noncyclic_guardrail": artifact["uses_noncyclic_guardrail"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
