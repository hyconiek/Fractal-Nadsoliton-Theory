#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1108 = GENERATED / "p1108_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_target_probe_summary.json"
IN_P1107 = GENERATED / "p1107_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_audit_probe_summary.json"
IN_T324 = ROOT / "T324_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_SPEC.md"

# Shortened generated basenames avoid filesystem filename-length limits.
OUT_JSON = GENERATED / "p1109_current_strict_t173_t176_minimal_onrd_ssl_rcw_actual_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1109_current_strict_t173_t176_minimal_onrd_ssl_rcw_actual_audit_probe_summary.json"

TARGET_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTargetActualRealizationAttemptExactRouteCoherenceWitnessRefinementTargetActualRealizationAttemptSharperSameLaneRouteCoherenceWitnessRefinementTarget_v1"
CURRENT_THEOREM_FILE = "N944_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_actual_realization_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1109_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "T324_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_SPEC.md",
        "N943_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_THEOREM.md",
        "p1108_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_sharper_same_lane_route_coherence_witness_refinement_target_probe.py",
        "p1108_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_target_probe.json",
        "p1108_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_target_probe_summary.json",
        "T325_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N945_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1110_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_sharper_same_lane_route_coherence_witness_refinement_target_actual_realization_attempt_probe.py",
        "p1110_current_strict_t173_t176_minimal_onrd_ssl_rcw_actual_attempt_probe.json",
        "p1110_current_strict_t173_t176_minimal_onrd_ssl_rcw_actual_attempt_probe_summary.json",
        OUT_JSON.name,
        OUT_SUMMARY.name,
    }
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded_names:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if TARGET_NAME in text and "actual_realization_attempt" in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1108, IN_P1107, IN_T324]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1109",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1108 = load_json(IN_P1108)
    p1107 = load_json(IN_P1107)
    t324_text = load_text(IN_T324)
    positive_candidates = scan_positive_actual_realization_candidates()

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

    t324_target_already_exported_only_at_future_only_strength = (
        p1108.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and bool(p1108.get("t324_target_exported_on_current_repo_state"))
        and bool(p1108.get("t324_target_keeps_reduction_supplier_solution_orientation_and_closure_open"))
    )

    p1107_branch_ordering_already_prefers_exact_sharper_same_lane_route_coherence_witness_first = (
        p1107.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1107.get("next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt"))
    )

    t324_same_exact_sharper_same_lane_route_coherence_witness_route_still_frozen = all(
        needle in t324_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t323_route_coherence_witness_actual_realization_attempt := yes",
            "target_is_exact_sharper_same_lane_route_coherence_witness_refinement_of_same_provider_class_route := yes",
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

    current_repo_has_exported_actual_realization_of_t324_target = len(positive_candidates) > 0

    t324_target_still_remains_future_only_not_actual_export = (
        t324_target_already_exported_only_at_future_only_strength
        and p1107_branch_ordering_already_prefers_exact_sharper_same_lane_route_coherence_witness_first
        and t324_same_exact_sharper_same_lane_route_coherence_witness_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t324_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t324_target = (
        t324_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "t324_target_already_exported_only_at_future_only_strength",
        t324_target_already_exported_only_at_future_only_strength,
        True,
        "P1108 already exports the exact T324 sharper same-lane route-coherence-witness target only at future-only strength.",
    )
    add_check(
        "p1107_branch_ordering_already_prefers_exact_sharper_same_lane_route_coherence_witness_first",
        p1107_branch_ordering_already_prefers_exact_sharper_same_lane_route_coherence_witness_first,
        True,
        "P1107 already orders continuation toward the exact sharper same-lane route-coherence-witness branch first.",
    )
    add_check(
        "t324_same_exact_sharper_same_lane_route_coherence_witness_route_still_frozen",
        t324_same_exact_sharper_same_lane_route_coherence_witness_route_still_frozen,
        True,
        "T324 still freezes the same exact sharper same-lane route-coherence-witness refinement route below the fixed T323 attempt.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t324_target",
        current_repo_has_exported_actual_realization_of_t324_target,
        False,
        "No stronger actual-realization artifact for this exact T324 target is exported on the current repo state.",
    )
    add_check(
        "t324_target_still_remains_future_only_not_actual_export",
        t324_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T324 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t324_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t324_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T324 target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t324_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1109",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "current_repo_has_exported_actual_realization_of_t324_target": current_repo_has_exported_actual_realization_of_t324_target,
        "t324_target_still_remains_future_only_not_actual_export": t324_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t324_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t324_target,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_name": artifact["target_name"],
        "current_repo_has_exported_actual_realization_of_t324_target": artifact["current_repo_has_exported_actual_realization_of_t324_target"],
        "t324_target_still_remains_future_only_not_actual_export": artifact["t324_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t324_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t324_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
