#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P760 = GENERATED / "p760_current_strict_t214_pair12_source_side_branch_selection_provider_actual_realization_nonexport_audit_probe_summary.json"
IN_P762 = GENERATED / "p762_current_strict_t216_pair12_source_side_branch_selection_provider_actual_realization_attempt_probe_summary.json"

OUT_JSON = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P742, IN_P747, IN_P760, IN_P762]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P763",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p742 = load_json(IN_P742)
    p747 = load_json(IN_P747)
    p760 = load_json(IN_P760)
    p762 = load_json(IN_P762)

    attempt = p762.get("first_actual_t213_realization_attempt") or {}
    immediate_missing_interface = str(attempt.get("immediate_missing_interface") or "")

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

    t216_attempt_already_exported_and_still_open = (
        bool(p762.get("t216_attempt_exported_on_current_repo_state"))
        and bool(p762.get("next_primary_t173_move_reduced_to_one_first_actual_t213_realization_attempt"))
        and bool(p762.get("first_actual_t213_realization_attempt_keeps_success_failure_open"))
        and immediate_missing_interface
        == "chart_sensitive_pair12_typed_descent_from_Sigma_sel_src_target_v1_to_the_surviving_F301_pair12_carrier_without_Q_basis_sel_v1_terminal_collapse_and_without_projector_only_atlas_collapse"
    )

    current_q_basis_terminal_collapse_is_still_the_strongest_exported_codomain_continuation = (
        bool(
            p742.get(
                "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation"
            )
        )
        and bool(
            p742.get(
                "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed"
            )
        )
        and not bool(
            p742.get(
                "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation"
            )
        )
    )

    current_local_atlas_projector_only_collapse_is_still_the_strongest_atlas_side_state = (
        bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported"))
        and bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"))
        and bool(p747.get("current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"))
        and not bool(p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge"))
    )

    current_repo_still_does_not_export_actual_realization_of_t213_target = (
        not bool(p760.get("t214_target_exported_on_current_repo_state"))
        and bool(p760.get("current_repo_still_does_not_export_actual_realization_of_t213_target"))
    )

    current_t216_attempt_immediate_missing_interface_is_still_unexported = (
        t216_attempt_already_exported_and_still_open
        and current_q_basis_terminal_collapse_is_still_the_strongest_exported_codomain_continuation
        and current_local_atlas_projector_only_collapse_is_still_the_strongest_atlas_side_state
        and current_repo_still_does_not_export_actual_realization_of_t213_target
    )

    current_t216_attempt_stalls_exactly_at_the_named_missing_interface = (
        current_t216_attempt_immediate_missing_interface_is_still_unexported
    )

    next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary = (
        current_t216_attempt_immediate_missing_interface_is_still_unexported
    )

    add_check(
        "t216_attempt_already_exported_and_still_open",
        t216_attempt_already_exported_and_still_open,
        True,
        "P762 already exports one exact first actual-realization attempt instance and keeps its success/failure open.",
    )
    add_check(
        "current_q_basis_terminal_collapse_is_still_the_strongest_exported_codomain_continuation",
        current_q_basis_terminal_collapse_is_still_the_strongest_exported_codomain_continuation,
        True,
        "P742 still freezes that the strongest current exported continuation out of Sigma_sel_src_target_v1 ends only at Q_basis_sel_v1 and not at a pair1/pair2 typed residual-datum continuation.",
    )
    add_check(
        "current_local_atlas_projector_only_collapse_is_still_the_strongest_atlas_side_state",
        current_local_atlas_projector_only_collapse_is_still_the_strongest_atlas_side_state,
        True,
        "P747 still freezes that the strongest current atlas-side state remains only the projector-level sign-gauge-safe local pair1/pair2 atlas lane with no positive bridge from Sigma_sel_src_target_v1.",
    )
    add_check(
        "current_repo_still_does_not_export_actual_realization_of_t213_target",
        current_repo_still_does_not_export_actual_realization_of_t213_target,
        True,
        "P760 still freezes that the repo does not yet export one actual realization of the T213 provider target.",
    )
    add_check(
        "current_t216_attempt_immediate_missing_interface_is_still_unexported",
        current_t216_attempt_immediate_missing_interface_is_still_unexported,
        True,
        "Therefore the exact missing interface named in T216 remains unexported on the current repo state.",
    )
    add_check(
        "current_t216_attempt_stalls_exactly_at_the_named_missing_interface",
        current_t216_attempt_stalls_exactly_at_the_named_missing_interface,
        True,
        "So the current first actual-realization attempt stalls exactly at the missing chart-sensitive pair1/pair2 typed descent interface named in T216.",
    )
    add_check(
        "next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary",
        next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary,
        True,
        "Hence the next honest move is now either export that exact interface or freeze an attempt-level failure boundary if the route stalls.",
    )

    status = (
        "PASS_STRICT_T217_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_INTERFACE_NONEXPORT_AUDITED"
        if not blocking and current_t216_attempt_immediate_missing_interface_is_still_unexported
        else "FAIL_STRICT_T217_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_INTERFACE_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P763",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "t217_boundary_name": "StrictT216ActualRealizationAttemptImmediateMissingInterfaceNonexportBoundary_strict_v1",
            "t217_boundary_exported_on_current_repo_state": status
            == "PASS_STRICT_T217_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_INTERFACE_NONEXPORT_AUDITED",
            "current_t216_attempt_immediate_missing_interface_is_still_unexported": current_t216_attempt_immediate_missing_interface_is_still_unexported,
            "current_t216_attempt_stalls_exactly_at_the_named_missing_interface": current_t216_attempt_stalls_exactly_at_the_named_missing_interface,
            "exact_named_missing_interface": immediate_missing_interface,
            "next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary": next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary,
            "no_false_pass": True,
        },
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P763",
        "status": status,
        "as_of": AS_OF,
        "t217_boundary_name": artifact["theorem_result"]["t217_boundary_name"],
        "t217_boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "t217_boundary_exported_on_current_repo_state"
        ],
        "current_t216_attempt_immediate_missing_interface_is_still_unexported": artifact[
            "theorem_result"
        ]["current_t216_attempt_immediate_missing_interface_is_still_unexported"],
        "current_t216_attempt_stalls_exactly_at_the_named_missing_interface": artifact[
            "theorem_result"
        ]["current_t216_attempt_stalls_exactly_at_the_named_missing_interface"],
        "exact_named_missing_interface": artifact["theorem_result"][
            "exact_named_missing_interface"
        ],
        "next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary": artifact[
            "theorem_result"
        ]["next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
