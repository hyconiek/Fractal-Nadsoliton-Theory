#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1071 = GENERATED / "p1071_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_verdict_or_exact_candidate_factor_integration_refinement_nonexport_audit_probe_summary.json"
IN_T310 = ROOT / "T310_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1072_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_probe.json"
OUT_SUMMARY = GENERATED / "p1072_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_probe_summary.json"

TARGET_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_side_open_feeder_support_refinement_attempt_candidate_factor_integration_refinement_target_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1071, IN_T310]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1072",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1071 = load_json(IN_P1071)
    t310_text = load_text(IN_T310)

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

    p1071_already_freezes_absence_of_verdict_and_factor_integration_export = (
        p1071.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_VERDICT_OR_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1071.get("next_honest_move_is_freeze_exact_candidate_factor_integration_refinement_target_below_same_attempt"))
    )

    t310_exact_candidate_factor_integration_target_frozen = all(
        needle in t310_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t309_open_feeder_support_refinement_actual_realization_attempt := yes",
            "target_is_exact_candidate_factor_integration_refinement_of_open_bundle := yes",
            "target_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "target_retains_feeder_side_noncyclic_provider_split_scope := yes",
            "target_preserves_same_sigma_witness_branch_scope := yes",
            "target_uses_existing_omega_feeder_theta_population_loop_candidate_strata_without_promoting_them_by_fiat := yes",
            "target_must_not_promote_to_actual_feeder_support_side_provider_support_by_fiat := yes",
            "target_must_not_promote_to_actual_feeder_support_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "target_must_remain_below_actual_feeder_support_side_provider_support := yes",
            "target_must_remain_below_actual_feeder_support_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
            "target_must_remain_below_actual_QW2191_discharge := yes",
        ]
    )

    t310_target_exported_on_current_repo_state = (
        p1071_already_freezes_absence_of_verdict_and_factor_integration_export
        and t310_exact_candidate_factor_integration_target_frozen
    )

    t310_target_keeps_open_bundle_below_actual_support = (
        t310_target_exported_on_current_repo_state
        and "target_must_remain_below_actual_feeder_support_side_provider_support := yes" in t310_text
        and "target_must_remain_below_actual_feeder_support_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t310_text
    )

    add_check(
        "p1071_already_freezes_absence_of_verdict_and_factor_integration_export",
        p1071_already_freezes_absence_of_verdict_and_factor_integration_export,
        True,
        "P1071 already freezes that the exact T309 attempt still lacks both a lawful verdict and a sharper candidate-factor integration refinement export.",
    )
    add_check(
        "t310_exact_candidate_factor_integration_target_frozen",
        t310_exact_candidate_factor_integration_target_frozen,
        True,
        "T310 freezes one exact candidate-factor integration refinement target below the same exact T309 attempt.",
    )
    add_check(
        "t310_target_exported_on_current_repo_state",
        t310_target_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact candidate-factor integration refinement target on the same Sigma feeder-side branch.",
    )
    add_check(
        "t310_target_keeps_open_bundle_below_actual_support",
        t310_target_keeps_open_bundle_below_actual_support,
        True,
        "That target still keeps feeder-support, theta, population, residual-bridge-support, and loop-break questions below actual support.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_EXPORTED"
        if not blocking and t310_target_exported_on_current_repo_state
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET"
    )

    artifact = {
        "stage": "P1072",
        "status": status,
        "as_of": AS_OF,
        "t310_target_name": TARGET_NAME,
        "t310_target_exported_on_current_repo_state": t310_target_exported_on_current_repo_state,
        "t310_target_keeps_open_bundle_below_actual_support": t310_target_keeps_open_bundle_below_actual_support,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t310_target_exported_on_current_repo_state": artifact["t310_target_exported_on_current_repo_state"],
        "t310_target_keeps_open_bundle_below_actual_support": artifact["t310_target_keeps_open_bundle_below_actual_support"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
