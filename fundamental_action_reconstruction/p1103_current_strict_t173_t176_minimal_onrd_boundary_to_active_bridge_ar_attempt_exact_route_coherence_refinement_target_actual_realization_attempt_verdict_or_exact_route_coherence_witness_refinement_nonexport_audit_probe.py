#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1102_JSON = GENERATED / "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe.json"
IN_T321 = ROOT / "T321_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1103_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_realization_attempt_verdict_or_exact_route_coherence_witness_refinement_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1103_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_realization_attempt_verdict_or_exact_route_coherence_witness_refinement_nonexport_audit_probe_summary.json"

ATTEMPT_NAME = "MinimalONRDBoundaryToActiveBridgeExactReductionTargetActualRealizationAttemptExactRouteCoherenceRefinementTargetActualRealizationAttempt_v1"
ROUTE_COHERENCE_WITNESS_REFINEMENT_STEM = "route_coherence_witness_refinement"
CURRENT_THEOREM_FILE = "N938_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_THEOREM.md"
FUTURE_TARGET_FILE = "T322_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_SPEC.md"
FUTURE_TARGET_THEOREM = "N939_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_TARGET_THEOREM.md"
FUTURE_TARGET_PROBE = "p1104_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_witness_refinement_target_probe.py"


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
        "P1103_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
        "T321_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N937_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_realization_attempt_probe.py",
        "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe.json",
        "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe_summary.json",
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


def scan_positive_route_coherence_witness_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1103_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
        "T321_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N937_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_realization_attempt_probe.py",
        "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe.json",
        "p1102_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_ar_attempt_exact_route_coherence_refinement_target_actual_attempt_probe_summary.json",
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
            if ROUTE_COHERENCE_WITNESS_REFINEMENT_STEM in text and f"future_{ROUTE_COHERENCE_WITNESS_REFINEMENT_STEM}" not in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1102_JSON, IN_T321]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1103",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1102 = load_json(IN_P1102_JSON)
    t321_text = load_text(IN_T321)

    verdict_candidates = scan_positive_verdict_candidates()
    witness_candidates = scan_positive_route_coherence_witness_candidates()

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

    t321_exact_attempt_already_exported_and_still_open = (
        p1102.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p1102.get("t321_attempt_exported_on_current_repo_state"))
        and bool(p1102.get("t321_attempt_keeps_reduction_supplier_solution_orientation_and_closure_open"))
        and str(p1102.get("attempt_name") or "") == ATTEMPT_NAME
    )

    t321_attempt_explicitly_keeps_open_bundle_below_realization = all(
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

    current_repo_has_lawful_verdict_for_exact_t321_attempt = len(verdict_candidates) > 0
    current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt = len(witness_candidates) > 0

    next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt = (
        t321_exact_attempt_already_exported_and_still_open
        and t321_attempt_explicitly_keeps_open_bundle_below_realization
        and not current_repo_has_lawful_verdict_for_exact_t321_attempt
        and not current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt
    )

    add_check(
        "t321_exact_attempt_already_exported_and_still_open",
        t321_exact_attempt_already_exported_and_still_open,
        True,
        "P1102 already exports one exact first actual-realization attempt on the T320 branch and keeps it open.",
    )
    add_check(
        "t321_attempt_explicitly_keeps_open_bundle_below_realization",
        t321_attempt_explicitly_keeps_open_bundle_below_realization,
        True,
        "T321 explicitly keeps reduction, supplier, solution, orientation, and closure questions below the attempt.",
    )
    add_check(
        "current_repo_has_lawful_verdict_for_exact_t321_attempt",
        current_repo_has_lawful_verdict_for_exact_t321_attempt,
        False,
        "The current repo should not already export a lawful verdict stronger than the frozen T321 attempt state.",
    )
    add_check(
        "current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt",
        current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt,
        False,
        "The current repo should not already export the sharper route-coherence-witness refinement below the same T321 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt",
        next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt,
        True,
        "If both stronger artifacts are absent, the next honest move is to freeze one exact route-coherence-witness refinement target below the same attempt.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking and next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt
        else "FAIL_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_AR_ATTEMPT_EXACT_ROUTE_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_ROUTE_COHERENCE_WITNESS_REFINEMENT_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1103",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "current_repo_has_lawful_verdict_for_exact_t321_attempt": current_repo_has_lawful_verdict_for_exact_t321_attempt,
        "current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt": current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt,
        "next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt": next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "current_repo_has_lawful_verdict_for_exact_t321_attempt": artifact["current_repo_has_lawful_verdict_for_exact_t321_attempt"],
        "current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt": artifact["current_repo_has_exact_route_coherence_witness_refinement_below_t321_attempt"],
        "next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt": artifact["next_honest_move_is_freeze_exact_route_coherence_witness_refinement_target_below_same_attempt"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
