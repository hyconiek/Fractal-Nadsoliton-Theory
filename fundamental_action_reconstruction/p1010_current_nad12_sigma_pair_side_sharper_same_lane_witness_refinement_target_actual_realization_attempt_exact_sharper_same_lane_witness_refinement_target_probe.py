#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1009 = GENERATED / "p1009_current_nad12_sigma_pair_side_sharper_same_lane_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe.json"
IN_T286 = ROOT / "T286_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1010_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe.json"
OUT_SUMMARY = GENERATED / "p1010_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe_summary.json"

TARGET_NAME = "Lambda_nad12_sigma_residual_shannon_pair_side_candidate_factor_coherence_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    try:
        p1009 = load_json(IN_P1009)
        t286_text = load_text(IN_T286)
    except FileNotFoundError as exc:
        missing_path = Path(exc.filename)
        artifact = {
            "stage": "P1010",
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

    p1009_already_freezes_absence_of_verdict_and_sharper_same_lane_export = (
        p1009.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1009.get("next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt"))
    )

    t286_exact_sharper_same_lane_witness_target_frozen = all(
        needle in t286_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t285_sharper_same_lane_witness_actual_realization_attempt := yes",
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

    t286_target_exported_on_current_repo_state = (
        p1009_already_freezes_absence_of_verdict_and_sharper_same_lane_export
        and t286_exact_sharper_same_lane_witness_target_frozen
    )

    t286_target_keeps_open_bundle_below_actual_support = (
        t286_target_exported_on_current_repo_state
        and "target_must_remain_below_actual_pair_realization_side_provider_support := yes" in t286_text
        and "target_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t286_text
    )

    add_check(
        "p1009_already_freezes_absence_of_verdict_and_sharper_same_lane_export",
        p1009_already_freezes_absence_of_verdict_and_sharper_same_lane_export,
        True,
        "P1009 already freezes that the exact T285 attempt still lacks both a lawful verdict and a sharper same-lane witness refinement export.",
    )
    add_check(
        "t286_exact_sharper_same_lane_witness_target_frozen",
        t286_exact_sharper_same_lane_witness_target_frozen,
        True,
        "T286 freezes one exact sharper same-lane witness refinement target below the same exact T285 attempt.",
    )
    add_check(
        "t286_target_exported_on_current_repo_state",
        t286_target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact sharper same-lane witness refinement target on the same Lambda witness branch.",
    )
    add_check(
        "t286_target_keeps_open_bundle_below_actual_support",
        t286_target_keeps_open_bundle_below_actual_support,
        True,
        "That target still keeps pair-support, theta, population, residual-bridge-support, and loop-break questions below actual support.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        if not blocking and t286_target_exported_on_current_repo_state
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET"
    )

    artifact = {
        "stage": "P1010",
        "status": status,
        "as_of": AS_OF,
        "t286_target_name": TARGET_NAME,
        "t286_target_exported_on_current_repo_state": t286_target_exported_on_current_repo_state,
        "t286_target_keeps_open_bundle_below_actual_support": t286_target_keeps_open_bundle_below_actual_support,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t286_target_exported_on_current_repo_state": artifact["t286_target_exported_on_current_repo_state"],
        "t286_target_keeps_open_bundle_below_actual_support": artifact["t286_target_keeps_open_bundle_below_actual_support"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
