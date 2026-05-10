#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P995 = GENERATED / "p995_current_nad12_sigma_pair_side_candidate_factor_coherence_target_actual_nonexport_probe.json"
IN_T279 = ROOT / "T279_CURRENT_NAD12_SIGMA_PAIR_SIDE_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p996_current_nad12_sigma_pair_side_candidate_factor_coherence_target_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p996_current_nad12_sigma_pair_side_candidate_factor_coherence_target_actual_attempt_probe_summary.json"

ATTEMPT_NAME = "Lambda_nad12_sigma_residual_shannon_pair_side_candidate_factor_coherence_refinement_target_actual_realization_attempt_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    try:
        p995 = load_json(IN_P995)
        t279_text = load_text(IN_T279)
    except FileNotFoundError as exc:
        missing_path = Path(exc.filename)
        artifact = {
            "stage": "P996",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(missing_path.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

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

    p995_already_freezes_target_as_future_only_not_actual = (
        p995.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(p995.get("next_honest_move_is_exact_actual_realization_attempt_of_same_t278_target"))
    )

    t279_exact_actual_realization_attempt_frozen = all(
        needle in t279_text
        for needle in [
            ATTEMPT_NAME,
            "attempt_is_over_exact_t278_candidate_factor_coherence_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_candidate_factor_coherence_refinement_target := yes",
            "attempt_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "attempt_retains_pair_indexed_noncyclic_provider_split_scope := yes",
            "attempt_preserves_same_lambda_witness_branch_scope := yes",
            "attempt_uses_existing_candidate_factor_strata_without_promoting_them_by_fiat := yes",
            "attempt_must_not_promote_to_actual_pair_realization_side_provider_support_by_fiat := yes",
            "attempt_must_not_promote_to_actual_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "attempt_must_remain_below_actual_pair_realization_side_provider_support := yes",
            "attempt_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
            "attempt_must_remain_below_actual_QW2191_discharge := yes",
        ]
    )

    t279_attempt_exported_on_current_repo_state = (
        p995_already_freezes_target_as_future_only_not_actual
        and t279_exact_actual_realization_attempt_frozen
    )

    t279_attempt_keeps_open_bundle_below_actual_support = (
        t279_attempt_exported_on_current_repo_state
        and "attempt_must_remain_below_actual_pair_realization_side_provider_support := yes" in t279_text
        and "attempt_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t279_text
    )

    add_check(
        "p995_already_freezes_target_as_future_only_not_actual",
        p995_already_freezes_target_as_future_only_not_actual,
        True,
        "P995 already freezes that the exact T278 target still remains future-only and not actually realized.",
    )
    add_check(
        "t279_exact_actual_realization_attempt_frozen",
        t279_exact_actual_realization_attempt_frozen,
        True,
        "T279 freezes one exact first actual-realization attempt on the same T278 target.",
    )
    add_check(
        "t279_attempt_exported_on_current_repo_state",
        t279_attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt on the same T278 branch.",
    )
    add_check(
        "t279_attempt_keeps_open_bundle_below_actual_support",
        t279_attempt_keeps_open_bundle_below_actual_support,
        True,
        "That attempt still keeps pair-support, theta, population, residual-bridge-support, and loop-break questions open.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and t279_attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P996",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "t279_attempt_exported_on_current_repo_state": t279_attempt_exported_on_current_repo_state,
        "t279_attempt_keeps_open_bundle_below_actual_support": t279_attempt_keeps_open_bundle_below_actual_support,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t279_attempt_exported_on_current_repo_state": artifact["t279_attempt_exported_on_current_repo_state"],
        "t279_attempt_keeps_open_bundle_below_actual_support": artifact["t279_attempt_keeps_open_bundle_below_actual_support"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
