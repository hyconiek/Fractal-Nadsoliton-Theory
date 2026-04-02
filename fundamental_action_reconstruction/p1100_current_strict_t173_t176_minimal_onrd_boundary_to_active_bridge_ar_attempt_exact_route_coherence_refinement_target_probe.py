#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-03-29"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1099 = GENERATED / "p1099_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_verdict_or_exact_route_coherence_refinement_nonexport_audit_probe_summary.json"
IN_T320 = ROOT / "T320_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1100_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_probe.json"
OUT_SUMMARY = GENERATED / "p1100_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_probe_summary.json"

TARGET_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTargetActualRealizationAttemptExactRouteCoherenceRefinementTarget_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1099, IN_T320]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1100",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1099 = load_json(IN_P1099)
    t320_text = load_text(IN_T320)

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

    p1099_already_orders_route_coherence_refinement_target = (
        p1099.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1099.get("next_honest_move_is_freeze_exact_route_coherence_refinement_target_below_same_attempt"))
    )

    t320_exact_route_coherence_refinement_target_frozen = all(
        needle in t320_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t319_actual_reduction_attempt := yes",
            "target_is_exact_refinement_of_same_route_coherence_bundle := yes",
            "target_keeps_minimal_onrd_boundary_as_source_seed := yes",
            "target_keeps_active_bridge_as_target_context := yes",
            "target_uses_genuinely_new_inversion_sensitive_source_side_provider_class_route := yes",
            "target_is_not_within_exported_noncyclic_provider_split_family := yes",
            "target_must_not_promote_to_exact_reduction_by_fiat := yes",
            "target_must_not_promote_to_lawful_supplier_by_fiat := yes",
            "target_must_not_promote_to_solution_or_strict_physical_orientation_datum_by_fiat := yes",
            "target_must_remain_below_T183_T176_QW2191_and_ToE_closure := yes",
        ]
    )

    t320_target_exported_on_current_repo_state = (
        p1099_already_orders_route_coherence_refinement_target
        and t320_exact_route_coherence_refinement_target_frozen
    )

    t320_target_keeps_reduction_supplier_solution_orientation_and_closure_open = (
        t320_target_exported_on_current_repo_state
        and "target_must_not_promote_to_exact_reduction_by_fiat := yes" in t320_text
        and "target_must_not_promote_to_lawful_supplier_by_fiat := yes" in t320_text
        and "target_must_not_promote_to_solution_or_strict_physical_orientation_datum_by_fiat := yes" in t320_text
        and "target_must_remain_below_T183_T176_QW2191_and_ToE_closure := yes" in t320_text
    )

    add_check(
        "p1099_already_orders_route_coherence_refinement_target",
        p1099_already_orders_route_coherence_refinement_target,
        True,
        "P1099 already orders continuation toward the exact route-coherence refinement target below T319.",
    )
    add_check(
        "t320_exact_route_coherence_refinement_target_frozen",
        t320_exact_route_coherence_refinement_target_frozen,
        True,
        "T320 freezes one exact lower route-coherence refinement target below the same T319 attempt.",
    )
    add_check(
        "t320_target_exported_on_current_repo_state",
        t320_target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact lower route-coherence refinement target on the same T319 branch.",
    )
    add_check(
        "t320_target_keeps_reduction_supplier_solution_orientation_and_closure_open",
        t320_target_keeps_reduction_supplier_solution_orientation_and_closure_open,
        True,
        "That target still keeps reduction, supplier, solution, orientation, and closure questions open.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_EXPORTED"
        if not blocking and t320_target_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET"
    )

    artifact = {
        "stage": "P1100",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "t320_target_exported_on_current_repo_state": t320_target_exported_on_current_repo_state,
        "t320_target_keeps_reduction_supplier_solution_orientation_and_closure_open": t320_target_keeps_reduction_supplier_solution_orientation_and_closure_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t320_target_exported_on_current_repo_state": artifact["t320_target_exported_on_current_repo_state"],
        "t320_target_keeps_reduction_supplier_solution_orientation_and_closure_open": artifact["t320_target_keeps_reduction_supplier_solution_orientation_and_closure_open"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
