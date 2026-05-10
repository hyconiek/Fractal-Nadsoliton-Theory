#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-28"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1084 = GENERATED / "p1084_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_probe.json"
IN_P1083 = GENERATED / "p1083_current_canonical_ontology_supported_nad12_sigma_feeder_cfwit_attempt_verdict_or_exact_sharper_same_lane_witness_refinement_nonexport_audit_probe.json"
IN_T316 = ROOT / "T316_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_CFWIT_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p1085_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p1085_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_actual_nonexport_probe_summary.json"

TARGET_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_side_candidate_factor_coherence_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_v1"
EXPECTED_ATTEMPT_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_side_candidate_factor_coherence_witness_refinement_target_actual_realization_attempt_sharper_same_lane_witness_refinement_target_actual_realization_attempt_v1"
CURRENT_THEOREM_FILE = "N920_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"


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
        "P1085_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "T316_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_CFWIT_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_SPEC.md",
        "N919_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_CFWIT_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_THEOREM.md",
        "p1084_current_canonical_ontology_supported_nad12_sigma_feeder_cfwit_attempt_exact_sharper_same_lane_witness_refinement_target_probe.py",
        "p1084_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_probe.json",
        "p1084_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_probe_summary.json",
        "T317_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N921_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p1086_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_refinement_target_actual_realization_attempt_probe.py",
        "p1086_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_actual_attempt_probe.json",
        "p1086_current_canonical_ontology_supported_nad12_sigma_feeder_sharper_same_lane_witness_target_actual_attempt_probe_summary.json",
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

    prerequisites = [IN_P1084, IN_P1083, IN_T316]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1085",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1084 = load_json(IN_P1084)
    p1083 = load_json(IN_P1083)
    t316_text = load_text(IN_T316)
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

    p1084_already_exports_exact_t316_target_only_at_future_only_strength = (
        p1084.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_EXPORTED"
        and bool(p1084.get("t316_target_exported_on_current_repo_state"))
        and bool(p1084.get("t316_target_keeps_open_bundle_below_actual_support"))
    )

    p1083_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first = (
        p1083.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p1083.get("next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt"))
    )

    t316_same_exact_sharper_same_lane_route_still_frozen = all(
        needle in t316_text
        for needle in [
            TARGET_NAME,
            "target_is_below_exact_t315_candidate_factor_coherence_witness_actual_realization_attempt := yes",
            "target_is_exact_sharper_same_lane_witness_refinement_of_same_sigma_branch := yes",
            "target_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "target_retains_feeder_side_noncyclic_provider_split_scope := yes",
            "target_preserves_same_sigma_witness_branch_scope := yes",
            "target_uses_existing_candidate_factor_and_witness_strata_without_promoting_them_by_fiat := yes",
            "target_must_not_promote_to_actual_feeder_support_side_provider_support_by_fiat := yes",
            "target_must_not_promote_to_actual_feeder_support_theta_export_pair_population_residual_bridge_support_or_loop_break_by_fiat := yes",
            "target_must_remain_below_actual_feeder_support_side_provider_support := yes",
            "target_must_remain_below_actual_feeder_support_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
        ]
    )

    current_repo_has_exported_actual_realization_of_t316_target = len(positive_candidates) > 0

    t316_target_still_remains_future_only_not_actual_export = (
        p1084_already_exports_exact_t316_target_only_at_future_only_strength
        and p1083_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first
        and t316_same_exact_sharper_same_lane_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t316_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t316_target = (
        t316_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "p1084_already_exports_exact_t316_target_only_at_future_only_strength",
        p1084_already_exports_exact_t316_target_only_at_future_only_strength,
        True,
        "P1084 already exports the exact T316 sharper same-lane witness target only at future-only strength.",
    )
    add_check(
        "p1083_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first",
        p1083_branch_ordering_already_prefers_exact_sharper_same_lane_witness_first,
        True,
        "P1083 already orders continuation toward the exact sharper same-lane witness branch first.",
    )
    add_check(
        "t316_same_exact_sharper_same_lane_route_still_frozen",
        t316_same_exact_sharper_same_lane_route_still_frozen,
        True,
        "T316 still freezes the same exact sharper same-lane witness route below the fixed T315 attempt.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t316_target",
        current_repo_has_exported_actual_realization_of_t316_target,
        False,
        "No stronger actual-realization artifact for this exact T316 target is exported on the current repo state.",
    )
    add_check(
        "t316_target_still_remains_future_only_not_actual_export",
        t316_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T316 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t316_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t316_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T316 target.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t316_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1085",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "expected_attempt_name": EXPECTED_ATTEMPT_NAME,
        "current_repo_has_exported_actual_realization_of_t316_target": current_repo_has_exported_actual_realization_of_t316_target,
        "t316_target_still_remains_future_only_not_actual_export": t316_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t316_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t316_target,
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
        "current_repo_has_exported_actual_realization_of_t316_target": artifact["current_repo_has_exported_actual_realization_of_t316_target"],
        "t316_target_still_remains_future_only_not_actual_export": artifact["t316_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t316_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t316_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
