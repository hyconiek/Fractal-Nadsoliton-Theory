#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P986 = GENERATED / "p986_current_nad12_sigma_pair_side_provider_support_witness_attempt_exact_open_support_refinement_target_probe_summary.json"
IN_P985 = GENERATED / "p985_current_nad12_sigma_pair_side_provider_support_witness_attempt_verdict_or_open_support_refinement_nonexport_probe_summary.json"
IN_T274 = ROOT / "T274_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_EXACT_OPEN_SUPPORT_REFINEMENT_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p987_current_nad12_sigma_pair_side_provider_support_witness_attempt_open_support_refinement_target_actual_nonexport_probe.json"
OUT_SUMMARY = GENERATED / "p987_current_nad12_sigma_pair_side_provider_support_witness_attempt_open_support_refinement_target_actual_nonexport_probe_summary.json"

TARGET_NAME = (
    "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_"
    "actual_realization_attempt_open_pair_support_theta_population_residual_bridge_"
    "loop_break_refinement_target_v1"
)
CURRENT_THEOREM_FILE = "N820_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_OPEN_SUPPORT_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md"


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
        "P987_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_OPEN_SUPPORT_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "T274_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_EXACT_OPEN_SUPPORT_REFINEMENT_TARGET_SPEC.md",
        "N819_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_EXACT_OPEN_SUPPORT_REFINEMENT_TARGET_THEOREM.md",
        "p986_current_nad12_sigma_pair_side_provider_support_witness_attempt_exact_open_support_refinement_target_probe.py",
        "p986_current_nad12_sigma_pair_side_provider_support_witness_attempt_exact_open_support_refinement_target_probe.json",
        "p986_current_nad12_sigma_pair_side_provider_support_witness_attempt_exact_open_support_refinement_target_probe_summary.json",
        "T275_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_OPEN_SUPPORT_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "N821_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_OPEN_SUPPORT_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
        "p988_current_nad12_sigma_pair_side_provider_support_witness_attempt_open_support_refinement_target_actual_realization_attempt_probe.py",
        "p988_current_nad12_sigma_pair_side_provider_support_witness_attempt_open_support_refinement_target_actual_attempt_probe.json",
        "p988_current_nad12_sigma_pair_side_provider_support_witness_attempt_open_support_refinement_target_actual_attempt_probe_summary.json",
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

    prerequisites = [IN_P986, IN_P985, IN_T274]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P987",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p986 = load_json(IN_P986)
    p985 = load_json(IN_P985)
    t274_text = load_text(IN_T274)
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

    t274_target_already_exported_only_at_future_only_strength = (
        p986.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_EXACT_OPEN_SUPPORT_REFINEMENT_TARGET_EXPORTED"
        and bool(p986.get("t274_target_exported_on_current_repo_state"))
        and bool(p986.get("t274_target_keeps_pair_support_theta_population_residual_bridge_loop_break_open"))
    )

    p985_branch_ordering_already_prefers_exact_open_refinement_first = (
        p985.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_VERDICT_OR_EXACT_OPEN_SUPPORT_REFINEMENT_NONEXPORT_AUDITED"
        and bool(p985.get("next_honest_move_is_freeze_exact_open_refinement_target_below_same_attempt"))
    )

    t274_same_exact_open_refinement_route_still_frozen = all(
        needle in t274_text
        for needle in [
            TARGET_NAME,
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

    current_repo_has_exported_actual_realization_of_t274_target = len(positive_candidates) > 0

    t274_target_still_remains_future_only_not_actual_export = (
        t274_target_already_exported_only_at_future_only_strength
        and p985_branch_ordering_already_prefers_exact_open_refinement_first
        and t274_same_exact_open_refinement_route_still_frozen
        and not current_repo_has_exported_actual_realization_of_t274_target
    )

    next_honest_move_is_exact_actual_realization_attempt_of_same_t274_target = (
        t274_target_still_remains_future_only_not_actual_export
    )

    add_check(
        "t274_target_already_exported_only_at_future_only_strength",
        t274_target_already_exported_only_at_future_only_strength,
        True,
        "P986 already exports the exact T274 open-support refinement target only at future-only strength.",
    )
    add_check(
        "p985_branch_ordering_already_prefers_exact_open_refinement_first",
        p985_branch_ordering_already_prefers_exact_open_refinement_first,
        True,
        "P985 already orders continuation toward the exact open-support refinement branch first.",
    )
    add_check(
        "t274_same_exact_open_refinement_route_still_frozen",
        t274_same_exact_open_refinement_route_still_frozen,
        True,
        "T274 still freezes the same exact open-support refinement route below the fixed T273 attempt.",
    )
    add_check(
        "current_repo_has_exported_actual_realization_of_t274_target",
        current_repo_has_exported_actual_realization_of_t274_target,
        False,
        "No stronger actual-realization artifact for this exact T274 target is exported on the current repo state.",
    )
    add_check(
        "t274_target_still_remains_future_only_not_actual_export",
        t274_target_still_remains_future_only_not_actual_export,
        True,
        "Therefore the exact T274 target still remains future-only and not actually realized.",
    )
    add_check(
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t274_target",
        next_honest_move_is_exact_actual_realization_attempt_of_same_t274_target,
        True,
        "The next honest move is now one exact actual-realization attempt on the same T274 target.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_OPEN_SUPPORT_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and t274_target_still_remains_future_only_not_actual_export
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_SIDE_PROVIDER_SUPPORT_WITNESS_ATTEMPT_OPEN_SUPPORT_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P987",
        "status": status,
        "as_of": AS_OF,
        "target_name": TARGET_NAME,
        "current_repo_has_exported_actual_realization_of_t274_target": current_repo_has_exported_actual_realization_of_t274_target,
        "t274_target_still_remains_future_only_not_actual_export": t274_target_still_remains_future_only_not_actual_export,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t274_target": next_honest_move_is_exact_actual_realization_attempt_of_same_t274_target,
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
        "current_repo_has_exported_actual_realization_of_t274_target": artifact["current_repo_has_exported_actual_realization_of_t274_target"],
        "t274_target_still_remains_future_only_not_actual_export": artifact["t274_target_still_remains_future_only_not_actual_export"],
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t274_target": artifact["next_honest_move_is_exact_actual_realization_attempt_of_same_t274_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
