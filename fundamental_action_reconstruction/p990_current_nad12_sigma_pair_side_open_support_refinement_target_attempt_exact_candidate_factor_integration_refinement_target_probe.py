#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P989 = GENERATED / "p989_current_nad12_sigma_pair_side_open_support_refinement_target_attempt_verdict_or_candidate_factor_integration_nonexport_probe.json"
IN_T276 = ROOT / "T276_CURRENT_NAD12_SIGMA_PAIR_SIDE_OPEN_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p990_current_nad12_sigma_pair_side_open_support_refinement_target_attempt_candidate_factor_integration_target_probe.json"
OUT_SUMMARY = GENERATED / "p990_current_nad12_sigma_pair_side_open_support_refinement_target_attempt_candidate_factor_integration_target_probe_summary.json"

TARGET_NAME = "Lambda_nad12_sigma_residual_shannon_pair_side_open_support_refinement_attempt_candidate_factor_integration_refinement_target_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P989, IN_T276]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P990",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p989 = load_json(IN_P989)
    t276_text = load_text(IN_T276)

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

    p989_already_freezes_absence_of_verdict_and_factor_integration_export = (
        p989.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_OPEN_SUPPORT_REFINEMENT_TARGET_ATTEMPT_VERDICT_OR_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p989.get("next_honest_move_is_freeze_exact_candidate_factor_integration_refinement_target_below_same_attempt"))
    )

    t276_exact_candidate_factor_integration_target_frozen = all(
        needle in t276_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t275_open_support_refinement_actual_realization_attempt := yes",
            "target_is_exact_candidate_factor_integration_refinement_of_open_bundle := yes",
            "target_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "target_retains_pair_indexed_noncyclic_provider_split_scope := yes",
            "target_preserves_same_lambda_witness_branch_scope := yes",
            "target_uses_existing_omega_feeder_theta_population_loop_candidate_strata_without_promoting_them_by_fiat := yes",
            "target_must_not_promote_to_actual_pair_realization_side_provider_support_by_fiat := yes",
            "target_must_not_promote_to_actual_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "target_must_remain_below_actual_pair_realization_side_provider_support := yes",
            "target_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
            "target_must_remain_below_actual_QW2191_discharge := yes",
        ]
    )

    t276_target_exported_on_current_repo_state = (
        p989_already_freezes_absence_of_verdict_and_factor_integration_export
        and t276_exact_candidate_factor_integration_target_frozen
    )

    t276_target_keeps_open_bundle_below_actual_support = (
        t276_target_exported_on_current_repo_state
        and "target_must_remain_below_actual_pair_realization_side_provider_support := yes" in t276_text
        and "target_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t276_text
    )

    add_check(
        "p989_already_freezes_absence_of_verdict_and_factor_integration_export",
        p989_already_freezes_absence_of_verdict_and_factor_integration_export,
        True,
        "P989 already freezes that the exact T275 attempt still lacks both a lawful verdict and a sharper candidate-factor integration refinement export.",
    )
    add_check(
        "t276_exact_candidate_factor_integration_target_frozen",
        t276_exact_candidate_factor_integration_target_frozen,
        True,
        "T276 freezes one exact candidate-factor integration refinement target below the same exact T275 attempt.",
    )
    add_check(
        "t276_target_exported_on_current_repo_state",
        t276_target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact candidate-factor integration refinement target on the same Lambda witness branch.",
    )
    add_check(
        "t276_target_keeps_open_bundle_below_actual_support",
        t276_target_keeps_open_bundle_below_actual_support,
        True,
        "That target still keeps pair-support, theta, population, residual-bridge-support, and loop-break questions below actual support.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_OPEN_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_EXPORTED"
        if not blocking and t276_target_exported_on_current_repo_state
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_OPEN_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET"
    )

    artifact = {
        "stage": "P990",
        "status": status,
        "as_of": AS_OF,
        "t276_target_name": TARGET_NAME,
        "t276_target_exported_on_current_repo_state": t276_target_exported_on_current_repo_state,
        "t276_target_keeps_open_bundle_below_actual_support": t276_target_keeps_open_bundle_below_actual_support,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t276_target_exported_on_current_repo_state": artifact["t276_target_exported_on_current_repo_state"],
        "t276_target_keeps_open_bundle_below_actual_support": artifact["t276_target_keeps_open_bundle_below_actual_support"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
