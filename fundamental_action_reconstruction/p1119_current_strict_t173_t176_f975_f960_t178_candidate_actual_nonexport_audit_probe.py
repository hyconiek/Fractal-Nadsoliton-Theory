#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-02"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1118 = GENERATED / "p1118_current_strict_t173_t176_f974_f960_t178_chart_seed_candidate_admission_probe_summary.json"
IN_F975 = GENERATED / "f975_current_strict_t173_t176_f974_f960_t178_chart_seed_candidate_packet_summary.json"
IN_P723 = GENERATED / "p723_current_strict_t178_source_topology_to_atlas_chart_seed_selection_bridge_nonexport_audit_probe_summary.json"
IN_P724 = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe_summary.json"
IN_P725 = GENERATED / "p725_current_strict_t179_positive_corridor_odd_even_lane_selection_bridge_nonexport_audit_probe_summary.json"
IN_P726 = GENERATED / "p726_current_strict_t180_positive_corridor_outer_interior_chart_selection_bridge_nonexport_audit_probe_summary.json"

OUT_JSON = GENERATED / "p1119_current_strict_t173_t176_f975_f960_t178_candidate_actual_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1119_current_strict_t173_t176_f975_f960_t178_candidate_actual_nonexport_audit_probe_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
T180_TARGET = "PositiveCorridorOuterInteriorChartSelectionBridge_global_C_v1_strict_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1118, IN_F975, IN_P723, IN_P724, IN_P725, IN_P726]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1119",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1118 = load_json(IN_P1118)
    f975 = load_json(IN_F975)
    p723 = load_json(IN_P723)
    p724 = load_json(IN_P724)
    p725 = load_json(IN_P725)
    p726 = load_json(IN_P726)

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

    t178_candidate_admitted = (
        p1118.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_ADMITTED"
        and p1118.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1118.get("current_best_new_narrow_provider_class_seed_candidate") == T178_TARGET
        and p1118.get("candidate_admissible_against_existing_f960_bridge_target") is True
        and f975.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_PACKET_EXPORTED"
        and f975.get("packet_exported_on_current_repo_state") is True
    )

    t178_itself_still_unexported = (
        p723.get("status") == "PASS_SOURCE_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p723.get("t178_target_name") == T178_TARGET
        and p723.get("t178_target_exported_on_current_repo_state") is False
    )

    positive_corridor_only_partial = (
        p724.get("status") == "PARTIAL_POSITIVE_SOURCE_POLARITY_REDUCES_ATLAS_ENTRY_CORRIDOR_ONLY"
        and p724.get("atlas_entry_roots_compatible_with_current_positive_source_polarity") == ["pair1", "pair2", "pair3", "pair5"]
        and p724.get("unique_chart_seed_selected") is False
    )

    sharper_lower_family_still_unexported = (
        p725.get("status") == "PASS_POSITIVE_CORRIDOR_ODD_EVEN_LANE_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p725.get("t179_target_exported_on_current_repo_state") is False
        and p726.get("status") == "PASS_POSITIVE_CORRIDOR_OUTER_INTERIOR_CHART_SELECTION_BRIDGE_NONEXPORT_AUDITED"
        and p726.get("t180_target_name") == T180_TARGET
        and p726.get("t180_target_exported_on_current_repo_state") is False
    )

    candidate_actual_realization_exported_on_current_repo_state = False

    add_check("t178_candidate_admitted", t178_candidate_admitted, True, "P1118/F975 already admit and freeze T178 as the current best narrow seed candidate beneath F960.")
    add_check("t178_itself_still_unexported", t178_itself_still_unexported, True, "P723 still keeps the T178 bridge target unexported.")
    add_check("positive_corridor_only_partial", positive_corridor_only_partial, True, "P724 still reduces source-side entry only to a positive corridor and does not select a unique chart seed.")
    add_check("sharper_lower_family_still_unexported", sharper_lower_family_still_unexported, True, "P725/P726 still keep the sharper lower lane-selection family unexported.")
    add_check("candidate_actual_realization_exported_on_current_repo_state", candidate_actual_realization_exported_on_current_repo_state, False, "Therefore the current-best T178 seed candidate is still not actually realized on the present repo state.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking else
        "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F975_F960_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1119",
        "status": status,
        "as_of": AS_OF,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "current_best_new_narrow_provider_class_seed_candidate": T178_TARGET,
        "candidate_actual_realization_exported_on_current_repo_state": candidate_actual_realization_exported_on_current_repo_state,
        "current_best_exact_lower_refinement_target": T180_TARGET,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_best_new_narrow_provider_class_seed_candidate": artifact["current_best_new_narrow_provider_class_seed_candidate"],
        "candidate_actual_realization_exported_on_current_repo_state": artifact["candidate_actual_realization_exported_on_current_repo_state"],
        "current_best_exact_lower_refinement_target": artifact["current_best_exact_lower_refinement_target"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
