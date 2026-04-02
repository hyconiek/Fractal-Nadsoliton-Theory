#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1072 = GENERATED / "p1072_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_probe_summary.json"
IN_P1071 = GENERATED / "p1071_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_verdict_or_exact_candidate_factor_integration_refinement_nonexport_audit_probe_summary.json"
IN_T310 = ROOT / "T310_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1073_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_refinement_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1073_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_refinement_target_actual_realization_nonexport_audit_probe_summary.json"

TARGET_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_side_open_feeder_support_refinement_attempt_candidate_factor_integration_refinement_target_v1"
CURRENT_THEOREM_FILE = "N908_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"


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
        "P1073_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "T310_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_SPEC.md",
        "N907_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_THEOREM.md",
        "p1072_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_exact_candidate_factor_integration_refinement_target_probe.py",
        "p1072_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_probe.json",
        "p1072_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_probe_summary.json",
        "T311_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N909_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1074_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_refinement_target_actual_realization_attempt_probe.py",
        "p1074_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_actual_attempt_probe.json",
        "p1074_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_actual_attempt_probe_summary.json",
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
            if TARGET_NAME in text and "actual_realization_attempt" in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1072, IN_P1071, IN_T310]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1073",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1072 = load_json(IN_P1072)
    p1071 = load_json(IN_P1071)
    t310_text = load_text(IN_T310)
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

    t310_target_already_exported_only_at_future_only_strength = (
        p1072.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_EXPORTED"
        and bool(p1072.get("t310_target_exported_on_current_repo_state"))
        and bool(p1072.get("t310_target_keeps_open_bundle_below_actual_support"))
    )

    p1071_branch_ordering_already_prefers_exact_candidate_factor_integration_first = (
        p1071.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_VERDICT_OR_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1071.get("next_honest_move_is_freeze_exact_candidate_factor_integration_refinement_target_below_same_attempt"))
    )

    t310_same_exact_candidate_factor_route_still_frozen = all(
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

    current_repo_has_exported_actual_realization_of_t310_target = len(positive_candidates) > 0

    t310_target_still_remains_future_only_not_actual_export = (
        t310_target_already_exported_only_at_future_only_strength
        and p1071_branch_ordering_already_prefers_exact_candidate_factor_integration_first
        and t310_same_exact_candidate_factor_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t310_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target = (
        t310_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "t310_target_already_exported_only_at_future_only_strength",
        t310_target_already_exported_only_at_future_only_strength,
        True,
        "P1072 already exports the exact T310 candidate-factor integration target only at future-only strength.",
    )
    add_check(
        "p1071_branch_ordering_already_prefers_exact_candidate_factor_integration_first",
        p1071_branch_ordering_already_prefers_exact_candidate_factor_integration_first,
        True,
        "P1071 already orders continuation toward the exact candidate-factor integration branch first.",
    )
    add_check(
        "t310_same_exact_candidate_factor_route_still_frozen",
        t310_same_exact_candidate_factor_route_still_frozen,
        True,
        "T310 still freezes the same exact candidate-factor integration route below the fixed T309 attempt.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t310_target",
        current_repo_has_exported_actual_realization_of_t310_target,
        False,
        "No stronger actual-realization artifact for this exact T310 target is exported on the current repo state.",
    )
    add_check(
        "t310_target_still_remains_future_only_not_actual_export",
        t310_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T310 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T310 target.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t310_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1073",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "current_repo_has_exported_actual_realization_of_t310_target": current_repo_has_exported_actual_realization_of_t310_target,
        "t310_target_still_remains_future_only_not_actual_export": t310_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target,
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
        "current_repo_has_exported_actual_realization_of_t310_target": artifact["current_repo_has_exported_actual_realization_of_t310_target"],
        "t310_target_still_remains_future_only_not_actual_export": artifact["t310_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t310_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
