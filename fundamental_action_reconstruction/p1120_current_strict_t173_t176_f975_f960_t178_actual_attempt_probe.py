#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-02"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1119 = GENERATED / "p1119_current_strict_t173_t176_f975_f960_t178_candidate_actual_nonexport_audit_probe_summary.json"
IN_T327 = ROOT / "T327_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_SOURCE_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"

OUT_JSON = GENERATED / "p1120_current_strict_t173_t176_f975_f960_t178_actual_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1120_current_strict_t173_t176_f975_f960_t178_actual_attempt_probe_summary.json"

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

    prerequisites = [IN_P1119, IN_T327]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1120",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1119 = load_json(IN_P1119)
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

    p1119_nonexport_audit_passed = (
        p1119.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and p1119.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1119.get("current_best_new_narrow_provider_class_seed_candidate") == T178_TARGET
        and p1119.get("candidate_actual_realization_exported_on_current_repo_state") is False
        and p1119.get("current_best_exact_lower_refinement_target") == T180_TARGET
    )

    t327_attempt_shape_frozen = all(
        needle in t327_text for needle in [
            ATTEMPT_OBJECT,
            f"target_candidate := {T178_TARGET}",
            f"active_missing_bridge := {ACTIVE_BRIDGE}",
            "current_source_positive_polarity_corridor := [pair1, pair2, pair3, pair5]",
            "unique_chart_seed_selected_on_current_repo_state := no",
            f"current_best_lower_refinement_target := {T180_TARGET}",
            "PositiveCorridorOddEvenLaneSelectionBridge_global_C_v1_strict_v1",
            "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1",
            "counts_as_actual_t178_export := no",
            "counts_as_lawful_supplier := no",
            "counts_as_solution := no",
            "counts_as_strict_physical_orientation_datum := no",
        ]
    )

    attempt_exported_on_current_repo_state = p1119_nonexport_audit_passed and t327_attempt_shape_frozen

    add_check("p1119_nonexport_audit_passed", p1119_nonexport_audit_passed, True, "P1119 already freezes that the admitted T178 candidate is still not actually realized and points to T180 as the sharpest lower refinement target.")
    add_check("t327_attempt_shape_frozen", t327_attempt_shape_frozen, True, "T327 freezes one first actual-realization attempt object beneath the admitted T178 candidate.")
    add_check("attempt_exported_on_current_repo_state", attempt_exported_on_current_repo_state, True, "Therefore the repo now exports one first actual-realization attempt beneath the admitted T178 seed candidate, while keeping all discharge claims open.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if not blocking and attempt_exported_on_current_repo_state else
        "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_ACTUAL_REALIZATION_ATTEMPT"
    )

    artifact = {
        "stage": "P1120",
        "status": status,
        "as_of": AS_OF,
        "attempt_object_id": ATTEMPT_OBJECT,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "target_candidate": T178_TARGET,
        "current_best_lower_refinement_target": T180_TARGET,
        "attempt_exported_on_current_repo_state": attempt_exported_on_current_repo_state,
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "checks": checks,
        "blocking_checks": blocking,
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
        "current_best_lower_refinement_target": artifact["current_best_lower_refinement_target"],
        "attempt_exported_on_current_repo_state": artifact["attempt_exported_on_current_repo_state"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
