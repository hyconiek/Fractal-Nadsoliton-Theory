#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-25"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1120 = GENERATED / "p1120_current_strict_t173_t176_f975_f960_t178_actual_attempt_probe_summary.json"
IN_P726 = GENERATED / "p726_current_strict_t180_positive_corridor_outer_interior_chart_selection_bridge_nonexport_audit_probe_summary.json"
IN_T327 = ROOT / "T327_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_SOURCE_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1121_current_strict_t173_t176_f975_f960_t178_actual_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1121_current_strict_t173_t176_f975_f960_t178_actual_attempt_verdict_or_exact_further_lower_boundary_nonexport_audit_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"
ATTEMPT_OBJECT = "SourceTopologyToAtlasChartSeedSelectionBridge_current_best_actual_realization_attempt_v1"


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

    prerequisites = [IN_P1120, IN_P726, IN_T327]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1121",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1120 = load_json(IN_P1120)
    p726 = load_json(IN_P726)
    t327_text = load_text(IN_T327)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append({
            "id": check_id,
            "actual": actual,
            "expected": expected,
            "pass": passed,
            "meaning": meaning,
        })
        if not passed:
            blocking.append(check_id)

    p1120_attempt_exported = (
        p1120.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        and p1120.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1120.get("target_candidate") == T178_TARGET
        and p1120.get("current_best_lower_refinement_target") == T180_TARGET
        and p1120.get("attempt_exported_on_current_repo_state") is True
    )

    p726_lower_support_known = (
        p726.get("status") == "PASS_POSITIVE_CORRIDOR_OUTER_INTERIOR_CHART_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p726.get("t180_target_name") == T180_TARGET
        and p726.get("t180_target_exported_on_current_repo_state") is False
    )

    t327_attempt_shape_frozen = all(
        needle in t327_text
        for needle in [
            ATTEMPT_OBJECT,
            f"target_candidate := {T178_TARGET}",
            f"active_missing_bridge := {ACTIVE_BRIDGE}",
            f"current_best_lower_refinement_target := {T180_TARGET}",
            "counts_as_actual_t178_export := no",
            "counts_as_lawful_supplier := no",
            "counts_as_solution := no",
            "counts_as_strict_physical_orientation_datum := no",
        ]
    )

    self_names = {
        Path(__file__).name,
        "P1121_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT_PROBE.md",
    }
    verdict_artifacts = sorted(
        path.name
        for path in ROOT.iterdir()
        if "CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_VERDICT" in path.name
        and path.name not in self_names
    )
    further_lower_boundary_target_artifacts = sorted(
        path.name
        for path in ROOT.iterdir()
        if "CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXACT_FURTHER_LOWER_BOUNDARY_TARGET" in path.name
    )

    lawful_verdict_present_for_same_attempt = len(verdict_artifacts) > 0
    exact_further_lower_boundary_target_below_same_attempt_present = (
        len(further_lower_boundary_target_artifacts) > 0
    )

    add_check(
        "p1120_attempt_exported",
        p1120_attempt_exported,
        True,
        "P1120 already freezes that the first actual-realization attempt beneath the admitted T178 candidate is exported.",
    )
    add_check(
        "p726_lower_support_known",
        p726_lower_support_known,
        True,
        "P726 already keeps the sharper lower support family explicit at T180 while still nonexported.",
    )
    add_check(
        "t327_attempt_shape_frozen",
        t327_attempt_shape_frozen,
        True,
        "T327 already freezes the scoped current-best actual-realization attempt with T180 as lower refinement guide only.",
    )
    add_check(
        "lawful_verdict_present_for_same_attempt",
        lawful_verdict_present_for_same_attempt,
        False,
        "No current export yet upgrades the T327 attempt into one lawful success/failure verdict.",
    )
    add_check(
        "exact_further_lower_boundary_target_below_same_attempt_present",
        exact_further_lower_boundary_target_below_same_attempt_present,
        False,
        "No current export yet freezes one exact further lower-boundary target explicitly beneath the same T327 attempt.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDITED"
        if not blocking
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_FURTHER_LOWER_BOUNDARY_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1121",
        "status": status,
        "as_of": AS_OF,
        "attempt_object_id": ATTEMPT_OBJECT,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "target_candidate": T178_TARGET,
        "current_best_lower_support_family_target": T180_TARGET,
        "lawful_verdict_present_for_same_attempt": lawful_verdict_present_for_same_attempt,
        "exact_further_lower_boundary_target_below_same_attempt_present": exact_further_lower_boundary_target_below_same_attempt_present,
        "known_verdict_artifacts": verdict_artifacts,
        "known_further_lower_boundary_target_artifacts": further_lower_boundary_target_artifacts,
        "checks": checks,
        "blocking_checks": blocking,
        "next_honest_move": "freeze_one_exact_further_lower_boundary_target_beneath_the_current_best_actual_realization_attempt",
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "attempt_object_id": artifact["attempt_object_id"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "target_candidate": artifact["target_candidate"],
        "current_best_lower_support_family_target": artifact["current_best_lower_support_family_target"],
        "lawful_verdict_present_for_same_attempt": artifact["lawful_verdict_present_for_same_attempt"],
        "exact_further_lower_boundary_target_below_same_attempt_present": artifact["exact_further_lower_boundary_target_below_same_attempt_present"],
        "next_honest_move": artifact["next_honest_move"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
