#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1105 = GENERATED / "p1105_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_realization_nonexport_audit_probe_summary.json"
IN_T323 = ROOT / "T323_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe_summary.json"

ATTEMPT_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTargetActualRealizationAttemptExactRouteCoherenceWitnessRefinementTargetActualRealizationAttempt_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1105, IN_T323]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1106",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1105 = load_json(IN_P1105)
    t323_text = load_text(IN_T323)

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

    p1105_already_freezes_target_as_future_only_not_actual = (
        p1105.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(p1105.get("next_honest_move_is_exact_actual_realization_attempt_of_same_t322_target"))
    )

    t323_exact_actual_realization_attempt_frozen = all(
        needle in t323_text
        for needle in [
            ATTEMPT_NAME,
            "attempt_is_over_exact_t322_route_coherence_witness_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_route_coherence_witness_refinement_target := yes",
            "attempt_keeps_minimal_onrd_boundary_as_source_seed := yes",
            "attempt_keeps_active_bridge_as_target_context := yes",
            "attempt_uses_genuinely_new_inversion_sensitive_source_side_provider_class_route := yes",
            "attempt_is_not_within_exported_noncyclic_provider_split_family := yes",
            "attempt_must_not_promote_to_exact_reduction_by_fiat := yes",
            "attempt_must_not_promote_to_lawful_supplier_by_fiat := yes",
            "attempt_must_not_promote_to_solution_or_strict_physical_orientation_datum_by_fiat := yes",
            "attempt_must_remain_below_T183_T176_QW2191_and_ToE_closure := yes",
        ]
    )

    t323_attempt_exported_on_current_repo_state = (
        p1105_already_freezes_target_as_future_only_not_actual
        and t323_exact_actual_realization_attempt_frozen
    )

    t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open = (
        t323_attempt_exported_on_current_repo_state
        and "attempt_must_not_promote_to_exact_reduction_by_fiat := yes" in t323_text
        and "attempt_must_not_promote_to_lawful_supplier_by_fiat := yes" in t323_text
        and "attempt_must_not_promote_to_solution_or_strict_physical_orientation_datum_by_fiat := yes" in t323_text
        and "attempt_must_remain_below_T183_T176_QW2191_and_ToE_closure := yes" in t323_text
    )

    add_check(
        "p1105_already_freezes_target_as_future_only_not_actual",
        p1105_already_freezes_target_as_future_only_not_actual,
        True,
        "P1105 already freezes that the exact T322 target still remains future-only and not actually realized.",
    )
    add_check(
        "t323_exact_actual_realization_attempt_frozen",
        t323_exact_actual_realization_attempt_frozen,
        True,
        "T323 freezes one exact first actual-realization attempt on the same T322 target.",
    )
    add_check(
        "t323_attempt_exported_on_current_repo_state",
        t323_attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt on the same T322 branch.",
    )
    add_check(
        "t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open",
        t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open,
        True,
        "That attempt still keeps reduction, supplier, solution, orientation, and closure questions open.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and t323_attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P1106",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "t323_attempt_exported_on_current_repo_state": t323_attempt_exported_on_current_repo_state,
        "t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open": t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t323_attempt_exported_on_current_repo_state": artifact["t323_attempt_exported_on_current_repo_state"],
        "t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open": artifact["t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
