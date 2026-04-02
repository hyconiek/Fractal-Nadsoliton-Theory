#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1004 = GENERATED / "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe.json"
IN_T283 = ROOT / "T283_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1005_current_nad12_sigma_pair_side_sharper_same_lane_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p1005_current_nad12_sigma_pair_side_sharper_same_lane_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe_summary.json"

ATTEMPT_NAME = "Lambda_nad12_sigma_residual_shannon_pair_side_candidate_factor_coherence_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_actual_realization_attempt_v1"
SHARPER_SAME_LANE_WITNESS_REFINEMENT_STEM = "sharper_same_lane_witness_refinement"
CURRENT_THEOREM_FILE = "N838_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_THEOREM.md"
FUTURE_TARGET_FILE = "T284_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_SPEC.md"
FUTURE_TARGET_THEOREM = "N839_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_THEOREM.md"
FUTURE_TARGET_PROBE = "p1006_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_target_actual_realization_attempt_exact_sharper_same_lane_witness_refinement_target_probe.py"


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
        "P1005_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "P1003_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N836_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "p1003_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_target_actual_realization_nonexport_audit_probe.py",
        "p1003_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe.json",
        "p1003_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe_summary.json",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
        "T283_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N837_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_target_actual_realization_attempt_probe.py",
        "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe.json",
        "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe_summary.json",
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


def scan_positive_sharper_same_lane_witness_refinement_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1005_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        "P1003_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "N836_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
        "p1003_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_target_actual_realization_nonexport_audit_probe.py",
        "p1003_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe.json",
        "p1003_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_nonexport_probe_summary.json",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
        "T283_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N837_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_target_actual_realization_attempt_probe.py",
        "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe.json",
        "p1004_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe_summary.json",
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
            if SHARPER_SAME_LANE_WITNESS_REFINEMENT_STEM in text and f"future_{SHARPER_SAME_LANE_WITNESS_REFINEMENT_STEM}" not in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1004, IN_T283]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1005",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1004 = load_json(IN_P1004)
    t283_text = load_text(IN_T283)

    verdict_candidates = scan_positive_verdict_candidates()
    sharper_candidates = scan_positive_sharper_same_lane_witness_refinement_candidates()

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

    t283_exact_attempt_already_exported_and_still_open = (
        p1004.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p1004.get("t283_attempt_exported_on_current_repo_state"))
        and bool(p1004.get("t283_attempt_keeps_open_bundle_below_actual_support"))
        and str(p1004.get("attempt_name") or "") == ATTEMPT_NAME
    )

    t283_attempt_explicitly_keeps_open_bundle_below_actual_support = all(
        needle in t283_text
        for needle in [
            ATTEMPT_NAME,
            "attempt_is_over_exact_t282_sharper_same_lane_witness_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_sharper_same_lane_witness_refinement_target := yes",
            "attempt_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "attempt_retains_pair_indexed_noncyclic_provider_split_scope := yes",
            "attempt_preserves_same_lambda_witness_branch_scope := yes",
            "attempt_uses_existing_candidate_factor_and_witness_strata_without_promoting_them_by_fiat := yes",
            "attempt_must_not_promote_to_actual_pair_realization_side_provider_support_by_fiat := yes",
            "attempt_must_not_promote_to_actual_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "attempt_must_remain_below_actual_pair_realization_side_provider_support := yes",
            "attempt_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
        ]
    )

    # Preserve exact repo logic even if the prior attempt stays open only at a bundle level.
    if not t283_attempt_explicitly_keeps_open_bundle_below_actual_support:
        t283_attempt_explicitly_keeps_open_bundle_below_actual_support = (
            ATTEMPT_NAME in t283_text
            and "attempt_must_remain_below_actual_pair_realization_side_provider_support := yes" in t283_text
            and "attempt_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t283_text
        )

    current_repo_has_lawful_verdict_for_exact_t283_attempt = len(verdict_candidates) > 0
    current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt = len(sharper_candidates) > 0

    current_repo_still_lacks_lawful_verdict_for_exact_t283_attempt = (
        t283_exact_attempt_already_exported_and_still_open
        and t283_attempt_explicitly_keeps_open_bundle_below_actual_support
        and not current_repo_has_lawful_verdict_for_exact_t283_attempt
    )

    current_repo_still_lacks_exact_sharper_same_lane_witness_refinement_below_t283_attempt = (
        current_repo_still_lacks_lawful_verdict_for_exact_t283_attempt
        and not current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt
    )

    next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt = (
        current_repo_still_lacks_lawful_verdict_for_exact_t283_attempt
        and current_repo_still_lacks_exact_sharper_same_lane_witness_refinement_below_t283_attempt
    )

    add_check(
        "t283_exact_attempt_already_exported_and_still_open",
        t283_exact_attempt_already_exported_and_still_open,
        True,
        "P1004 already exports one exact first actual-realization attempt on the sharpened T282 branch and keeps it open.",
    )
    add_check(
        "t283_attempt_explicitly_keeps_open_bundle_below_actual_support",
        t283_attempt_explicitly_keeps_open_bundle_below_actual_support,
        True,
        "T283 explicitly keeps pair-support, theta, population, residual-bridge-support, and loop-break questions below the attempt.",
    )
    add_check(
        "current_repo_has_lawful_verdict_for_exact_t283_attempt",
        current_repo_has_lawful_verdict_for_exact_t283_attempt,
        False,
        "The current repo should not already export a lawful verdict stronger than the frozen T283 attempt state.",
    )
    add_check(
        "current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt",
        current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt,
        False,
        "The current repo should not already export the sharper same-lane witness refinement below the same T283 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt",
        next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt,
        True,
        "If both stronger artifacts are absent, the next honest move is to freeze one exact sharper same-lane witness refinement target below the same attempt.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking and next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_AUDIT"
    )

    artifact = {
        "stage": "P1005",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "current_repo_has_lawful_verdict_for_exact_t283_attempt": current_repo_has_lawful_verdict_for_exact_t283_attempt,
        "current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt": current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt,
        "verdict_candidates": verdict_candidates,
        "sharper_same_lane_witness_candidates": sharper_candidates,
        "next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt": next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "current_repo_has_lawful_verdict_for_exact_t283_attempt": artifact["current_repo_has_lawful_verdict_for_exact_t283_attempt"],
        "current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt": artifact["current_repo_has_exact_sharper_same_lane_witness_refinement_below_t283_attempt"],
        "next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt": artifact["next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
