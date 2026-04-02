#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1103 = GENERATED / "p1103_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_realization_attempt_verdict_or_exact_route_coherence_witness_refinement_nonexport_audit_probe_summary.json"
IN_T322 = ROOT / "T322_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1104_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_probe.json"
OUT_SUMMARY = GENERATED / "p1104_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_probe_summary.json"

TARGET_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTargetActualRealizationAttemptExactRouteCoherenceWitnessRefinementTarget_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1103, IN_T322]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1104",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1103 = load_json(IN_P1103)
    t322_text = load_text(IN_T322)

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

    p1103_already_orders_exact_route_coherence_witness_target = (
        p1103.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1103.get("next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt"))
    )

    t322_exact_route_coherence_witness_target_frozen = all(
        needle in t322_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t321_route_coherence_actual_realization_attempt := yes",
            "target_is_exact_route_coherence_witness_refinement_of_same_provider_class_route := yes",
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

    t322_target_exported_on_current_repo_state = (
        p1103_already_orders_exact_route_coherence_witness_target
        and t322_exact_route_coherence_witness_target_frozen
    )

    t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open = (
        t322_target_exported_on_current_repo_state
        and "target_must_not_promote_to_exact_reduction_by_fiat := yes" in t322_text
        and "target_must_not_promote_to_lawful_supplier_by_fiat := yes" in t322_text
        and "target_must_not_promote_to_solution_or_strict_physical_orientation_datum_by_fiat := yes" in t322_text
        and "target_must_remain_below_T183_T176_QW2191_and_ToE_closure := yes" in t322_text
    )

    add_check(
        "p1103_already_orders_exact_route_coherence_witness_target",
        p1103_already_orders_exact_route_coherence_witness_target,
        True,
        "P1103 already freezes that the next honest move is one exact route-coherence-witness refinement target below T321.",
    )
    add_check(
        "t322_exact_route_coherence_witness_target_frozen",
        t322_exact_route_coherence_witness_target_frozen,
        True,
        "T322 freezes one exact route-coherence-witness refinement target on the same provider-class route.",
    )
    add_check(
        "t322_target_exported_on_current_repo_state",
        t322_target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact route-coherence-witness refinement target below T321.",
    )
    add_check(
        "t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open",
        t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open,
        True,
        "That target still keeps reduction, supplier, solution, orientation, and closure questions open.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        if not blocking and t322_target_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET"
    )

    artifact = {
        "stage": "P1104",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "t322_target_exported_on_current_repo_state": t322_target_exported_on_current_repo_state,
        "t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open": t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t322_target_exported_on_current_repo_state": artifact["t322_target_exported_on_current_repo_state"],
        "t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open": artifact["t322_target_keeps_reduction_supplier_solution_orientation_and_closure_open"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
