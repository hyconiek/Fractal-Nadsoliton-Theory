#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P985 = GENERATED / "p985_current_nad12_sigma_pair_side_provider_support_witness_attempt_verdict_or_open_support_refinement_nonexport_probe_summary.json"
IN_T274 = ROOT / "T274_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_EXACT_OPEN_SUPPORT_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p986_current_nad12_sigma_pair_side_provider_support_witness_attempt_exact_open_support_refinement_target_probe.json"
OUT_SUMMARY = GENERATED / "p986_current_nad12_sigma_pair_side_provider_support_witness_attempt_exact_open_support_refinement_target_probe_summary.json"

OPEN_REFINEMENT_TARGET_NAME = (
    "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_"
    "actual_realization_attempt_open_pair_support_theta_population_residual_bridge_"
    "loop_break_refinement_target_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P985, IN_T274]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P986",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p985 = load_json(IN_P985)
    t274_text = load_text(IN_T274)

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

    p985_already_freezes_absence_of_verdict_and_open_refinement_export = (
        p985.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_VERDICT_OR_EXACT_OPEN_SUPPORT_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p985.get("next_honest_move_is_freeze_exact_open_refinement_target_below_same_attempt"))
    )

    t274_exact_open_refinement_target_frozen = all(
        needle in t274_text
        for needle in [
            OPEN_REFINEMENT_TARGET_NAME,
            "target_is_below_exact_lambda_pair_realization_side_provider_support_witness_actual_realization_attempt := yes",
            "target_is_exact_refinement_of_open_pair_support_theta_population_residual_bridge_loop_break_bundle := yes",
            "target_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "target_retains_pair_indexed_noncyclic_provider_split_scope := yes",
            "target_preserves_same_lambda_witness_branch_scope := yes",
            "target_must_not_promote_to_actual_pair_realization_side_provider_support_by_fiat := yes",
            "target_must_not_promote_to_actual_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "target_must_remain_below_actual_pair_realization_side_provider_support := yes",
            "target_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
            "target_must_remain_below_actual_QW2191_discharge := yes",
        ]
    )

    t274_target_exported_on_current_repo_state = (
        p985_already_freezes_absence_of_verdict_and_open_refinement_export
        and t274_exact_open_refinement_target_frozen
    )

    t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open = (
        t274_target_exported_on_current_repo_state
        and "target_must_remain_below_actual_pair_realization_side_provider_support := yes" in t274_text
        and "target_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t274_text
    )

    add_check(
        "p985_already_freezes_absence_of_verdict_and_open_refinement_export",
        p985_already_freezes_absence_of_verdict_and_open_refinement_export,
        True,
        "P985 already freezes that the exact T273 attempt still lacks both a lawful verdict and a sharper exact open refinement export.",
    )
    add_check(
        "t274_exact_open_refinement_target_frozen",
        t274_exact_open_refinement_target_frozen,
        True,
        "T274 freezes one exact open refinement target below the same exact T273 attempt.",
    )
    add_check(
        "t274_target_exported_on_current_repo_state",
        t274_target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact open refinement target on the same Lambda witness branch.",
    )
    add_check(
        "t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open",
        t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open,
        True,
        "That target still keeps pair-support, theta, population, residual-bridge-support, and loop-break questions open.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_EXACT_OPEN_SUPPORT_REFINEMENT_TARGET_EXPORTED"
        if not blocking and t274_target_exported_on_current_repo_state
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_EXACT_OPEN_SUPPORT_REFINEMENT_TARGET"
    )

    artifact = {
        "stage": "P986",
        "status": status,
        "as_of": AS_OF,
        "t274_target_name": OPEN_REFINEMENT_TARGET_NAME,
        "t274_target_exported_on_current_repo_state": t274_target_exported_on_current_repo_state,
        "t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open": t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t274_target_exported_on_current_repo_state": artifact["t274_target_exported_on_current_repo_state"],
        "t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open": artifact["t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
