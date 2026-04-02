#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1101 = GENERATED / "p1101_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_realization_nonexport_audit_probe_summary.json"
IN_T321 = ROOT / "T321_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe_summary.json"

ATTEMPT_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTargetActualRealizationAttemptExactRouteCoherenceRefinementTargetActualRealizationAttempt_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1101, IN_T321]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1102",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1101 = load_json(IN_P1101)
    t321_text = load_text(IN_T321)

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

    p1101_already_freezes_target_as_future_only_not_actual = (
        p1101.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(p1101.get("next_honest_move_is_exact_actual_realization_attempt_of_same_t320_target"))
    )

    t321_exact_actual_realization_attempt_frozen = all(
        needle in t321_text
        for needle in [
            ATTEMPT_NAME,
            "attempt_is_over_exact_t320_route_coherence_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_route_coherence_refinement_target := yes",
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

    t321_attempt_exported_on_current_repo_state = (
        p1101_already_freezes_target_as_future_only_not_actual
        and t321_exact_actual_realization_attempt_frozen
    )

    t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open = (
        t321_attempt_exported_on_current_repo_state
        and "attempt_must_not_promote_to_exact_reduction_by_fiat := yes" in t321_text
        and "attempt_must_not_promote_to_lawful_supplier_by_fiat := yes" in t321_text
        and "attempt_must_not_promote_to_solution_or_strict_physical_orientation_datum_by_fiat := yes" in t321_text
        and "attempt_must_remain_below_T183_T176_QW2191_and_ToE_closure := yes" in t321_text
    )

    add_check(
        "p1101_already_freezes_target_as_future_only_not_actual",
        p1101_already_freezes_target_as_future_only_not_actual,
        True,
        "P1101 already freezes that the exact T320 target still remains future-only and not actually realized.",
    )
    add_check(
        "t321_exact_actual_realization_attempt_frozen",
        t321_exact_actual_realization_attempt_frozen,
        True,
        "T321 freezes one exact first actual-realization attempt on the same T320 target.",
    )
    add_check(
        "t321_attempt_exported_on_current_repo_state",
        t321_attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt on the same T320 branch.",
    )
    add_check(
        "t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open",
        t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open,
        True,
        "That attempt still keeps reduction, supplier, solution, orientation, and closure questions open.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and t321_attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P1102",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "t321_attempt_exported_on_current_repo_state": t321_attempt_exported_on_current_repo_state,
        "t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open": t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t321_attempt_exported_on_current_repo_state": artifact["t321_attempt_exported_on_current_repo_state"],
        "t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open": artifact["t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
