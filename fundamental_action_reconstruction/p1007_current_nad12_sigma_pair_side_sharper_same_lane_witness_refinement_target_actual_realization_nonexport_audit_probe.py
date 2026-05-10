#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1006 = GENERATED / "p1006_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe.json"
IN_P1005 = GENERATED / "p1005_current_nad12_sigma_pair_side_sharper_same_lane_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe.json"
IN_T284 = ROOT / "T284_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1007_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p1007_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe_summary.json"

TARGET_NAME = "Lambda_nad12_sigma_residual_shannon_pair_side_candidate_factor_coherence_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_v1"
EXPECTED_ATTEMPT_NAME = "Lambda_nad12_sigma_residual_shannon_pair_side_candidate_factor_coherence_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_actual_realization_attempt_v1"
CURRENT_THEOREM_FILE = "N840_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"


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
        "P1007_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "T284_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_SPEC.md",
        "N839_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_THEOREM.md",
        "p1006_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_target_actual_realization_attempt_exact_sharper_same_lane_witness_refinement_target_probe.py",
        "p1006_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe.json",
        "p1006_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe_summary.json",
        "T285_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N841_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1008_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_target_actual_realization_attempt_probe.py",
        "p1008_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe.json",
        "p1008_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe_summary.json",
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
            if EXPECTED_ATTEMPT_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1006, IN_P1005, IN_T284]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1007",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1006 = load_json(IN_P1006)
    p1005 = load_json(IN_P1005)
    t284_text = load_text(IN_T284)
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

    p1006_already_exports_exact_t284_target_only_at_future_only_strength = (
        p1006.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and bool(p1006.get("t284_target_exported_on_current_repo_state"))
        and bool(p1006.get("t284_target_keeps_open_bundle_below_actual_support"))
    )

    p1005_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first = (
        p1005.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1005.get("next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt"))
    )

    t284_same_exact_sharper_same_lane_route_still_frozen = all(
        needle in t284_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t283_sharper_same_lane_witness_actual_realization_attempt := yes",
            "target_is_exact_sharper_same_lane_witness_refinement_of_same_lambda_branch := yes",
            "target_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "target_retains_pair_indexed_noncyclic_provider_split_scope := yes",
            "target_preserves_same_lambda_witness_branch_scope := yes",
            "target_uses_existing_candidate_factor_and_witness_strata_without_promoting_them_by_fiat := yes",
            "target_must_not_promote_to_actual_pair_realization_side_provider_support_by_fiat := yes",
            "target_must_not_promote_to_actual_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "target_must_remain_below_actual_pair_realization_side_provider_support := yes",
            "target_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
            "target_must_remain_below_actual_QW2191_discharge := yes",
        ]
    )

    current_repo_has_exported_actual_realization_of_t284_target = len(positive_candidates) > 0

    t284_target_still_remains_future_only_not_actual_export = (
        p1006_already_exports_exact_t284_target_only_at_future_only_strength
        and p1005_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first
        and t284_same_exact_sharper_same_lane_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t284_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t284_target = (
        t284_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "p1006_already_exports_exact_t284_target_only_at_future_only_strength",
        p1006_already_exports_exact_t284_target_only_at_future_only_strength,
        True,
        "P1006 already exports the exact T284 sharper same-lane witness target only at future-only strength.",
    )
    add_check(
        "p1005_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first",
        p1005_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first,
        True,
        "P1005 already orders continuation toward the exact sharper same-lane witness branch first.",
    )
    add_check(
        "t284_same_exact_sharper_same_lane_route_still_frozen",
        t284_same_exact_sharper_same_lane_route_still_frozen,
        True,
        "T284 still freezes the same exact sharper same-lane witness route below the fixed T283 attempt.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t284_target",
        current_repo_has_exported_actual_realization_of_t284_target,
        False,
        "No stronger actual-realization artifact for this exact T284 target is exported on the current repo state.",
    )
    add_check(
        "t284_target_still_remains_future_only_not_actual_export",
        t284_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T284 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t284_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t284_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T284 target.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t284_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1007",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "expected_attempt_name": EXPECTED_ATTEMPT_NAME,
        "current_repo_has_exported_actual_realization_of_t284_target": current_repo_has_exported_actual_realization_of_t284_target,
        "t284_target_still_remains_future_only_not_actual_export": t284_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t284_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t284_target,
        "actual_realization_candidates": positive_candidates,
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
        "current_repo_has_exported_actual_realization_of_t284_target": artifact["current_repo_has_exported_actual_realization_of_t284_target"],
        "t284_target_still_remains_future_only_not_actual_export": artifact["t284_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t284_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t284_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
