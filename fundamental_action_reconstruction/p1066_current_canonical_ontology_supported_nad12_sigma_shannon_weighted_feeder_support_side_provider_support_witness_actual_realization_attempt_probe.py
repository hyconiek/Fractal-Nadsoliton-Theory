#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F963 = GENERATED / "f963_current_nad12_sigma_pair_side_sharper_same_lane_witness_refinement_exhaustion_and_feeder_pivot_packet_summary.json"
IN_P1065 = GENERATED / "p1065_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_side_provider_support_witness_actual_realization_nonexport_audit_probe_summary.json"
IN_T307 = ROOT / "T307_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1066_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_side_provider_support_witness_actual_realization_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1066_current_canonical_ontology_supported_nad12_sigma_shannon_weighted_feeder_support_side_provider_support_witness_actual_realization_attempt_probe_summary.json"

ATTEMPT_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_support_side_provider_support_witness_actual_realization_attempt_v1"
WITNESS_NAME = "Sigma_nad12_sigma_residual_shannon_feeder_support_side_provider_support_witness_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F963, IN_P1065, IN_T307]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1066",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f963 = load_json(IN_F963)
    p1065 = load_json(IN_P1065)
    t307_text = load_text(IN_T307)

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
        f963.get("status")
        == "F963_EXECUTED_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_EXHAUSTION_AND_FEEDER_PIVOT_PACKET_NO_FALSE_PASS"
        and f963.get("preferred_first_pivot_branch") == WITNESS_NAME
    )

    witness_nonexport_boundary_already_exported = (
        p1065.get("status")
        == "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1065.get("witness_name") == WITNESS_NAME
        and bool(p1065.get("witness_still_remains_future_only_not_actual_export"))
        and bool(p1065.get("next_honest_move_is_exact_actual_realization_attempt_of_same_witness"))
    )

    t307_attempt_shape_frozen = all(
        needle in t307_text
        for needle in [
            ATTEMPT_NAME,
            WITNESS_NAME,
            "attempt_is_over_exact_sigma_feeder_support_side_provider_support_witness := yes",
            "attempt_is_actual_realization_attempt_of_exact_feeder_support_side_provider_support_witness := yes",
            "attempt_retains_nad12_sigma_residual_shannon_weighted_provenance := yes",
            "attempt_retains_feeder_side_noncyclic_provider_split_scope := yes",
            "attempt_must_not_promote_to_actual_feeder_support_side_provider_support_by_fiat := yes",
            "attempt_must_not_promote_to_actual_feeder_support_theta_export_pair_population_or_loop_break_by_fiat := yes",
            "attempt_must_remain_below_actual_feeder_support_side_provider_support := yes",
            "attempt_must_remain_below_actual_feeder_support_theta_export_pair_population_residual_bridge_support_and_loop_break := yes",
        ]
    )

    attempt_exported_on_current_repo_state = (
        preferred_first_noncyclic_pivot_branch_already_selected
        and witness_nonexport_boundary_already_exported
        and t307_attempt_shape_frozen
    )

    attempt_keeps_support_theta_population_loopbreak_open = (
        attempt_exported_on_current_repo_state
        and "attempt_must_remain_below_actual_feeder_support_side_provider_support := yes" in t307_text
        and "attempt_must_remain_below_actual_feeder_support_theta_export_pair_population_residual_bridge_support_and_loop_break := yes" in t307_text
    )

    add_check(
        "preferred_first_noncyclic_pivot_branch_already_selected",
        preferred_first_noncyclic_pivot_branch_already_selected,
        True,
        "F963 already selects the Sigma feeder-side witness branch as the preferred first pivot branch.",
    )
    add_check(
        "witness_nonexport_boundary_already_exported",
        witness_nonexport_boundary_already_exported,
        True,
        "P1065 already freezes that the Sigma feeder-side witness remains future-only and not actually realized.",
    )
    add_check(
        "t307_attempt_shape_frozen",
        t307_attempt_shape_frozen,
        True,
        "T307 exports one exact first actual-realization attempt shape on the same feeder-side witness branch.",
    )
    add_check(
        "attempt_exported_on_current_repo_state",
        attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt on the chosen feeder-side pivot branch.",
    )
    add_check(
        "attempt_keeps_support_theta_population_loopbreak_open",
        attempt_keeps_support_theta_population_loopbreak_open,
        True,
        "That attempt still keeps support/theta/population/loop-break questions open.",
    )

    status = (
        "PASS_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_CANONICAL_ONTOLOGY_SUPPORTED_NAD12_SIGMA_SHANNON_WEIGHTED_FEEDER_SUPPORT_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P1066",
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
