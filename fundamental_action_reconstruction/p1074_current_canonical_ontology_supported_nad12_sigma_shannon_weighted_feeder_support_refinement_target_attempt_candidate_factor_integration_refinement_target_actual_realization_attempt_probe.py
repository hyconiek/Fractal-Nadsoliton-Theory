#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1073 = GENERATED / "p1073_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_refinement_target_actual_realization_nonexport_audit_probe_summary.json"
IN_T311 = ROOT / "T311_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1074_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1074_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_actual_attempt_probe_summary.json"

ATTEMPT_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_side_open_feeder_support_refinement_attempt_candidate_factor_integration_refinement_target_actual_realization_attempt_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    try:
        p1073 = load_json(IN_P1073)
        t311_text = load_text(IN_T311)
    except FileNotFoundError as exc:
        missing_path = Path(exc.filename)
        artifact = {
            "stage": "P1074",
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

    p1073_already_freezes_target_as_future_only_not_actual = (
        p1073.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and bool(p1073.get("next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target"))
    )

    t311_exact_actual_realization_attempt_frozen = all(
        needle in t311_text
        for needle in [
            ATTEMPT_NAME,
            "attempt_is_over_exact_t310_candidate_factor_integration_refinement_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_candidate_factor_integration_refinement_target := yes",
            "attempt_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "attempt_retains_feeder_side_noncyclic_provider_split_scope := yes",
            "attempt_preserves_same_sigma_witness_branch_scope := yes",
            "attempt_uses_existing_omega_feeder_theta_population_loop_candidate_strata_without_promoting_them_by_fiat := yes",
            "attempt_must_not_promote_to_actual_feeder_support_side_provider_support_by_fiat := yes",
            "attempt_must_not_promote_to_actual_feeder_support_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "attempt_must_remain_below_actual_feeder_support_side_provider_support := yes",
            "attempt_must_remain_below_actual_feeder_support_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
            "attempt_must_remain_below_actual_QW2191_discharge := yes",
        ]
    )

    t311_attempt_exported_on_current_repo_state = (
        p1073_already_freezes_target_as_future_only_not_actual
        and t311_exact_actual_realization_attempt_frozen
    )

    t311_attempt_keeps_open_bundle_below_actual_support = (
        t311_attempt_exported_on_current_repo_state
        and "attempt_must_remain_below_actual_feeder_support_side_provider_support := yes" in t311_text
        and "attempt_must_remain_below_actual_feeder_support_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t311_text
    )

    add_check(
        "p1073_already_freezes_target_as_future_only_not_actual",
        p1073_already_freezes_target_as_future_only_not_actual,
        True,
        "P1073 already freezes that the exact T310 target still remains future-only and not actually realized.",
    )
    add_check(
        "t311_exact_actual_realization_attempt_frozen",
        t311_exact_actual_realization_attempt_frozen,
        True,
        "T311 freezes one exact first actual-realization attempt on the same T310 target.",
    )
    add_check(
        "t311_attempt_exported_on_current_repo_state",
        t311_attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt on the same T310 branch.",
    )
    add_check(
        "t311_attempt_keeps_open_bundle_below_actual_support",
        t311_attempt_keeps_open_bundle_below_actual_support,
        True,
        "That attempt still keeps feeder-support, theta, population, residual-bridge-support, and loop-break questions open.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and t311_attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P1074",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "t311_attempt_exported_on_current_repo_state": t311_attempt_exported_on_current_repo_state,
        "t311_attempt_keeps_open_bundle_below_actual_support": t311_attempt_keeps_open_bundle_below_actual_support,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t311_attempt_exported_on_current_repo_state": artifact["t311_attempt_exported_on_current_repo_state"],
        "t311_attempt_keeps_open_bundle_below_actual_support": artifact["t311_attempt_keeps_open_bundle_below_actual_support"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
