#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-05-10"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1137 = GENERATED / "p1137_current_strict_t173_t176_existing_t339_actual_realization_nonexport_audit_probe_summary.json"
IN_T340 = ROOT / "T340_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1138_current_strict_t173_t176_existing_t339_actual_realization_attempt_export_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1138_current_strict_t173_t176_existing_t339_actual_realization_attempt_export_audit_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T179_TARGET = "PositiveCorridorOddEvenLaneSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"
TARGET_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_v1"
)
PARENT_ATTEMPT_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_v1"
)
ATTEMPT_NAME = (
    "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_exact_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_exact_even_further_lower_boundary_target_actual_realization_attempt_exact_yet_further_lower_boundary_target_actual_realization_attempt_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1137, IN_T340]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1138",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1137 = load_json(IN_P1137)
    t340_text = load_text(IN_T340)

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

    p1137_nonexport_audit_passed = (
        p1137.get("status")
        == "PASS_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1137.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1137.get("target_candidate") == T178_TARGET
        and p1137.get("current_best_lower_support_family_target") == T180_TARGET
        and p1137.get("target_name") == TARGET_NAME
        and p1137.get("attempt_object_id") == PARENT_ATTEMPT_NAME
        and p1137.get("current_repo_has_exported_actual_realization_of_t339_target") is False
        and p1137.get("next_honest_move_is_exact_actual_realization_attempt_of_same_t339_target") is True
    )

    t340_attempt_shape_frozen = all(
        needle in t340_text
        for needle in [
            "T340_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_SPEC_NO_FALSE_PASS",
            ATTEMPT_NAME,
            TARGET_NAME,
            PARENT_ATTEMPT_NAME,
            T180_TARGET,
            T179_TARGET,
            "PARTIAL_POSITIVE_SOURCE_POLARITY_REDUCES_ATLAS_ENTRY_CORRIDOR_ONLY",
            ACTIVE_BRIDGE,
            T178_TARGET,
            "attempt_is_over_exact_t339_yet_further_lower_boundary_target := yes",
            "attempt_is_actual_realization_attempt_of_exact_yet_further_lower_boundary_target := yes",
            "attempt_must_not_promote_to_lawful_verdict_for_t329_by_fiat := yes",
            "attempt_must_not_promote_to_actual_t180_export_by_fiat := yes",
            "attempt_must_not_promote_to_actual_t179_export_by_fiat := yes",
            "attempt_must_not_promote_to_actual_t178_export_by_fiat := yes",
            "attempt_must_not_promote_to_lawful_supplier_by_fiat := yes",
            "attempt_must_not_promote_to_solution_by_fiat := yes",
            "attempt_must_not_promote_to_strict_physical_orientation_datum_by_fiat := yes",
            "attempt_must_remain_below_actual_T183_discharge := yes",
            "attempt_must_remain_below_actual_T176_discharge := yes",
            "attempt_must_remain_below_actual_QW2191_discharge := yes",
        ]
    )

    attempt_exported_on_current_repo_state = (
        p1137_nonexport_audit_passed and t340_attempt_shape_frozen
    )

    add_check(
        "p1137_nonexport_audit_passed",
        p1137_nonexport_audit_passed,
        True,
        "P1137 already freezes that the exact T339 target still lacks actual realization on the current repo state.",
    )
    add_check(
        "t340_attempt_shape_frozen",
        t340_attempt_shape_frozen,
        True,
        "T340 freezes one exact first actual-realization attempt instance over the exact T339 target.",
    )
    add_check(
        "attempt_exported_on_current_repo_state",
        attempt_exported_on_current_repo_state,
        True,
        "Therefore the current repo exports one exact first actual-realization attempt over the exact T339 yet-further lower-boundary target.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_T339_ACTUAL_REALIZATION_ATTEMPT_PROBE"
    )

    artifact = {
        "stage": "P1138",
        "status": status,
        "as_of": AS_OF,
        "attempt_name": ATTEMPT_NAME,
        "target_name": TARGET_NAME,
        "parent_attempt_name": PARENT_ATTEMPT_NAME,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "target_candidate": T178_TARGET,
        "current_best_lower_support_family_target": T180_TARGET,
        "intermediate_lane_support_target": T179_TARGET,
        "attempt_exported_on_current_repo_state": attempt_exported_on_current_repo_state,
        "counts_as_lawful_verdict_for_t329": False,
        "counts_as_actual_t180_export": False,
        "counts_as_actual_t179_export": False,
        "counts_as_actual_t178_export": False,
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "counts_as_t183_discharge": False,
        "counts_as_t176_discharge": False,
        "counts_as_qw2191_discharge": False,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "attempt_name": artifact["attempt_name"],
        "target_name": artifact["target_name"],
        "parent_attempt_name": artifact["parent_attempt_name"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "target_candidate": artifact["target_candidate"],
        "current_best_lower_support_family_target": artifact["current_best_lower_support_family_target"],
        "intermediate_lane_support_target": artifact["intermediate_lane_support_target"],
        "attempt_exported_on_current_repo_state": artifact["attempt_exported_on_current_repo_state"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
