#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"
IN_P983 = GENERATED / "p983_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_nonexport_probe_summary.json"
IN_T273 = ROOT / "T273_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p984_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p984_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_attempt_probe_summary.json"

ATTEMPT_NAME = "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_actual_realization_attempt_v1"
WITNESS_NAME = "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F949, IN_P983, IN_T273]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P984",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f949 = load_json(IN_F949)
    p983 = load_json(IN_P983)
    t273_text = load_text(IN_T273)

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

    preferred_first_noncyclic_pivot_branch_already_selected = (
        f949.get("preferred_first_pivot_branch") == WITNESS_NAME
        and bool(f949.get("packet_exported_on_current_repo_state"))
    )

    witness_nonexport_boundary_already_exported = (
        p983.get("status")
        == "PASS_CURRENT_NAD12_SIGMA_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p983.get("witness_name") == WITNESS_NAME
        and bool(p983.get("witness_still_remains_future_only_not_actual_export"))
        and bool(p983.get("next_honest_move_is_exact_actual_realization_attempt_of_same_witness"))
    )

    t273_attempt_shape_frozen = all(
        needle in t273_text
        for needle in [
            ATTEMPT_NAME,
            "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1",
            "attempt_is_over_exact_lambda_pair_realization_side_provider_support_witness := yes",
            "attempt_is_actual_realization_attempt_of_exact_pair_realization_side_provider_support_witness := yes",
            "attempt_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "attempt_retains_pair_indexed_noncyclic_provider_split_scope := yes",
            "attempt_must_not_promote_to_actual_pair_realization_side_provider_support_by_fiat := yes",
            "attempt_must_not_promote_to_actual_theta_export_pair_population_or_loop_break_by_fiat := yes",
            "attempt_must_remain_below_actual_pair_realization_side_provider_support := yes",
            "attempt_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
        ]
    )

    attempt_exported_on_current_repo_state = (
        preferred_first_noncyclic_pivot_branch_already_selected
        and witness_nonexport_boundary_already_exported
        and t273_attempt_shape_frozen
    )

    attempt_keeps_support_theta_population_loopbreak_open = (
        attempt_exported_on_current_repo_state
        and "attempt_must_remain_below_actual_pair_realization_side_provider_support := yes" in t273_text
        and "attempt_must_remain_below_actual_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t273_text
    )

    add_check(
        "preferred_first_noncyclic_pivot_branch_already_selected",
        preferred_first_noncyclic_pivot_branch_already_selected,
        True,
        "F949 already selects the Lambda witness branch as the preferred first pivot branch.",
    )
    add_check(
        "witness_nonexport_boundary_already_exported",
        witness_nonexport_boundary_already_exported,
        True,
        "P983 already freezes that the Lambda witness remains future-only and not actually realized.",
    )
    add_check(
        "t273_attempt_shape_frozen",
        t273_attempt_shape_frozen,
        True,
        "T273 exports one exact first actual-realization attempt shape on the same witness branch.",
    )
    add_check(
        "attempt_exported_on_current_repo_state",
        attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt on the chosen pivot branch.",
    )
    add_check(
        "attempt_keeps_support_theta_population_loopbreak_open",
        attempt_keeps_support_theta_population_loopbreak_open,
        True,
        "That attempt still keeps support/theta/population/loop-break questions open.",
    )

    status = (
        "PASS_CURRENT_NAD12_SIGMA_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_NAD12_SIGMA_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P984",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "attempt_exported_on_current_repo_state": attempt_exported_on_current_repo_state,
        "attempt_keeps_support_theta_population_loopbreak_open": attempt_keeps_support_theta_population_loopbreak_open,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "attempt_exported_on_current_repo_state": artifact["attempt_exported_on_current_repo_state"],
        "attempt_keeps_support_theta_population_loopbreak_open": artifact["attempt_keeps_support_theta_population_loopbreak_open"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
