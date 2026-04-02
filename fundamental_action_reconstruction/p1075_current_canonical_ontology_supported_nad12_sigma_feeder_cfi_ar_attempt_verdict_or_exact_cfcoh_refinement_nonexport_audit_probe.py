#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1074_JSON = GENERATED / "p1074_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_refinement_target_attempt_candidate_factor_integration_target_actual_attempt_probe.json"
IN_T311 = ROOT / "T311_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1075_current_canonical_ontology_supported_nad12_sigma_feeder_cfi_ar_attempt_verdict_or_exact_cfcoh_refinement_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1075_current_canonical_ontology_supported_nad12_sigma_feeder_cfi_ar_attempt_verdict_or_exact_cfcoh_refinement_nonexport_audit_probe_summary.json"

ATTEMPT_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_side_open_feeder_support_refinement_attempt_candidate_factor_integration_refinement_target_actual_realization_attempt_v1"
COHERENCE_REFINEMENT_STEM = "candidate_factor_coherence_refinement"
CURRENT_THEOREM_FILE = "N910_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_CFI_AR_ATTEMPT_VERDICT_OR_EXACT_CFCOH_REFINEMENT_NONEXPORT_AUDIT_THEOREM.md"
FUTURE_TARGET_FILE = "T312_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_TARGET_SPEC.md"
FUTURE_TARGET_THEOREM = "N911_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_TARGET_THEOREM.md"
FUTURE_TARGET_PROBE = "p1076_current_canonical_ontology_supported_nad12_sigma_feeder_cfi_ar_attempt_exact_cfcoh_refinement_target_probe.py"


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
        "P1075_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_CFI_AR_ATTEMPT_VERDICT_OR_EXACT_CFCOH_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
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


def scan_positive_candidate_factor_coherence_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py", "generated/*.json")
    excluded_names = {
        Path(__file__).name,
        CURRENT_THEOREM_FILE,
        "P1075_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_CFI_AR_ATTEMPT_VERDICT_OR_EXACT_CFCOH_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
        FUTURE_TARGET_FILE,
        FUTURE_TARGET_THEOREM,
        FUTURE_TARGET_PROBE,
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
            if ATTEMPT_NAME not in text:
                continue
            if COHERENCE_REFINEMENT_STEM in text and f"future_{COHERENCE_REFINEMENT_STEM}" not in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1074_JSON, IN_T311]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1075",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1074 = load_json(IN_P1074_JSON)
    t311_text = load_text(IN_T311)

    verdict_candidates = scan_positive_verdict_candidates()
    coherence_candidates = scan_positive_candidate_factor_coherence_candidates()

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

    t311_exact_attempt_already_exported_and_still_open = (
        p1074.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and bool(p1074.get("t311_attempt_exported_on_current_repo_state"))
        and bool(p1074.get("t311_attempt_keeps_open_bundle_below_actual_support"))
        and str(p1074.get("attempt_name") or "") == ATTEMPT_NAME
    )

    t311_attempt_explicitly_keeps_open_bundle_below_actual_support = all(
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
        ]
    )

    current_repo_has_lawful_verdict_for_exact_t311_attempt = len(verdict_candidates) > 0
    current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt = len(coherence_candidates) > 0

    current_repo_still_lacks_lawful_verdict_for_exact_t311_attempt = (
        t311_exact_attempt_already_exported_and_still_open
        and t311_attempt_explicitly_keeps_open_bundle_below_actual_support
        and not current_repo_has_lawful_verdict_for_exact_t311_attempt
    )

    current_repo_still_lacks_exact_candidate_factor_coherence_refinement_below_t311_attempt = (
        current_repo_still_lacks_lawful_verdict_for_exact_t311_attempt
        and not current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt
    )

    next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt = (
        current_repo_still_lacks_lawful_verdict_for_exact_t311_attempt
        and current_repo_still_lacks_exact_candidate_factor_coherence_refinement_below_t311_attempt
    )

    add_check(
        "t311_exact_attempt_already_exported_and_still_open",
        t311_exact_attempt_already_exported_and_still_open,
        True,
        "P1074 already exports one exact first actual-realization attempt on the sharpened T310 branch and keeps it open.",
    )
    add_check(
        "t311_attempt_explicitly_keeps_open_bundle_below_actual_support",
        t311_attempt_explicitly_keeps_open_bundle_below_actual_support,
        True,
        "T311 explicitly keeps feeder-support, theta, population, residual-bridge-support, and loop-break questions below the attempt.",
    )
    add_check(
        "current_repo_has_lawful_verdict_for_exact_t311_attempt",
        current_repo_has_lawful_verdict_for_exact_t311_attempt,
        False,
        "The current repo should not already export a lawful verdict stronger than the frozen T311 attempt state.",
    )
    add_check(
        "current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt",
        current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt,
        False,
        "The current repo should not already export the sharper candidate-factor coherence refinement below the same T311 attempt.",
    )
    add_check(
        "next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt",
        next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt,
        True,
        "If both stronger artifacts are absent, the next honest move is to freeze one exact candidate-factor coherence refinement target below the same attempt.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_NONEXPORT_AUDITED"
        if not blocking and next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_REFINEMENT_TARGET_ATTEMPT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_CANDIDATE_FACTOR_COHERENCE_REFINEMENT_AUDIT"
    )

    artifact = {
        "stage": "P1075",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "current_repo_has_lawful_verdict_for_exact_t311_attempt": current_repo_has_lawful_verdict_for_exact_t311_attempt,
        "current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt": current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt,
        "verdict_candidates": verdict_candidates,
        "candidate_factor_coherence_candidates": coherence_candidates,
        "next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt": next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "current_repo_has_lawful_verdict_for_exact_t311_attempt": artifact["current_repo_has_lawful_verdict_for_exact_t311_attempt"],
        "current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt": artifact["current_repo_has_exact_candidate_factor_coherence_refinement_below_t311_attempt"],
        "next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt": artifact["next_honest_move_is_freeze_exact_candidate_factor_coherence_refinement_target_below_same_attempt"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
