#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1106_JSON = GENERATED / "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe.json"
IN_T323 = ROOT / "T323_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

# Shortened generated basenames avoid filesystem filename-length limits.
OUT_JSON = GENERATED / "p1107_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1107_current_strict_t173_t176_minimal_onrd_sharper_same_lane_rcw_audit_probe_summary.json"

ATTEMPT_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTargetActualRealizationAttemptExactRouteCoherenceWitnessRefinementTargetActualRealizationAttempt_v1"
SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STEM = "sharper_same_lane_route_coherence_witness_refinement"
CURRENT_THEOREM_FILE = "N942_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_THEOREM.md"
FUTURE_TARGET_FILE = "T324_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_SPEC.md"
FUTURE_TARGET_THEOREM = "N943_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_THEOREM.md"
FUTURE_TARGET_PROBE = "p1108_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_sharper_same_lane_route_coherence_witness_refinement_target_probe.py"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_verdict_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1107_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
        "T323_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N941_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_realization_attempt_probe.py",
        "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe.json",
        "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe_summary.json",
        "p1107_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_same_lane_route_coherence_witness_refinement_nonexport_audit_probe.json",
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
            if ATTEMPT_NAME not in text:
                continue
            if any(
                needle in text
                for needle in (
                    "actual_realization_attempt_verdict",
                    "lawful_realization_verdict",
                    "lawful_verdict",
                    "explicit_success_verdict",
                    "success_verdict",
                )
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def scan_positive_sharper_same_lane_route_coherence_witness_refinement_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1107_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
        "T323_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N941_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_realization_attempt_probe.py",
        "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe.json",
        "p1106_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_attempt_probe_summary.json",
        "p1107_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_same_lane_route_coherence_witness_refinement_nonexport_audit_probe.json",
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
            if ATTEMPT_NAME not in text:
                continue
            if SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STEM in text and f"future_{SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STEM}" not in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1106_JSON, IN_T323]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1107",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1106 = load_json(IN_P1106_JSON)
    t323_text = load_text(IN_T323)

    verdict_candidates = scan_positive_verdict_candidates()
    sharper_candidates = scan_positive_sharper_same_lane_route_coherence_witness_refinement_candidates()

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

    t323_exact_attempt_already_exported_and_still_open = (
        p1106.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p1106.get("t323_attempt_exported_on_current_repo_state"))
        and bool(p1106.get("t323_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open"))
        and str(p1106.get("attempt_name") or "") == ATTEMPT_NAME
    )

    t323_attempt_explicitly_keeps_open_bundle_below_realization = all(
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

    current_repo_has_lawful_verdict_for_exact_t323_attempt = len(verdict_candidates) > 0
    current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt = len(sharper_candidates) > 0

    next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt = (
        t323_exact_attempt_already_exported_and_still_open
        and t323_attempt_explicitly_keeps_open_bundle_below_realization
        and not current_repo_has_lawful_verdict_for_exact_t323_attempt
        and not current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt
    )

    add_check(
        "t323_exact_attempt_already_exported_and_still_open",
        t323_exact_attempt_already_exported_and_still_open,
        True,
        "P1106 already exports one exact first actual-realization attempt on the T322 branch and keeps it open.",
    )
    add_check(
        "t323_attempt_explicitly_keeps_open_bundle_below_realization",
        t323_attempt_explicitly_keeps_open_bundle_below_realization,
        True,
        "T323 explicitly keeps reduction, supplier, solution, orientation, and closure questions below the attempt.",
    )
    add_check(
        "current_repo_has_lawful_verdict_for_exact_t323_attempt",
        current_repo_has_lawful_verdict_for_exact_t323_attempt,
        False,
        "The current repo should not already export a lawful verdict stronger than the frozen T323 attempt state.",
    )
    add_check(
        "current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt",
        current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt,
        False,
        "The current repo should not already export the sharper same-lane route-coherence-witness refinement below the same T323 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt",
        next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt,
        True,
        "If both stronger artifacts are absent, the next honest move is to freeze one exact sharper same-lane route-coherence-witness refinement target below the same attempt.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking and next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1107",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "current_repo_has_lawful_verdict_for_exact_t323_attempt": current_repo_has_lawful_verdict_for_exact_t323_attempt,
        "current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt": current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt,
        "next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt": next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "current_repo_has_lawful_verdict_for_exact_t323_attempt": artifact["current_repo_has_lawful_verdict_for_exact_t323_attempt"],
        "current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt": artifact["current_repo_has_exact_sharper_same_lane_route_coherence_witness_refinement_below_t323_attempt"],
        "next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt": artifact["next_honest_move_is_freeze_exact_sharper_same_lane_route_coherence_witness_refinement_target_below_same_attempt"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
