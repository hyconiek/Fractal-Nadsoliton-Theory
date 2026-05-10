#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1118 = GENERATED / "p1118_current_strict_t173_t176_f974_f960_t178_chart_seed_candidate_admission_probe_summary.json"
IN_F975 = ROOT / "F975_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_SOURCE_TO_ATLAS_CHART_SEED_SELECTION_BRIDGE_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_PACKET.md"

OUT_JSON = GENERATED / "f975_current_strict_t173_t176_f974_f960_t178_chart_seed_candidate_packet.json"
OUT_SUMMARY = GENERATED / "f975_current_strict_t173_t176_f974_f960_t178_chart_seed_candidate_packet_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
T178_TARGET = "SourceTopologyToAtlasChartSeedSelectionBridge_global_C_v1_strict_v1"
ACTIVE_CONTRACT = "attack_existing_t178_chart_seed_candidate_beneath_existing_f960_bridge_target"


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

    prerequisites = [IN_P1118, IN_F975]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F975",
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
    f975_text = load_text(IN_F975)

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

    p1118_candidate_admission_passed = (
        p1118.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_ADMITTED"
        and p1118.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1118.get("current_best_new_narrow_provider_class_seed_candidate") == T178_TARGET
        and p1118.get("candidate_admissible_against_existing_f960_bridge_target") is True
        and p1118.get("candidate_grade") == "current_best_new_narrow_provider_class_seed_only"
        and p1118.get("active_contract") == ACTIVE_CONTRACT
    )

    f975_packet_shape_frozen = all(
        needle in f975_text for needle in [
            "Xi_current_strict_t173_t176_existing_f974_f960_t178_chart_seed_candidate_packet_v1",
            f"active_missing_bridge := {ACTIVE_BRIDGE}",
            f"existing_t178_target := {T178_TARGET}",
            "source_topology_lane_is_physics_facing_but_chart_blind := yes",
            "existing_t178_target_exported_on_current_repo_state := no",
            "unique_chart_seed_selected_on_current_repo_state := no",
            "source_side_input_leg_same_lane_reentry_disallowed := yes",
            f"current_best_new_narrow_provider_class_seed_candidate := {T178_TARGET}",
            "candidate_grade := current_best_new_narrow_provider_class_seed_only",
            f"active_contract := {ACTIVE_CONTRACT}",
            "counts_as_lawful_supplier := no",
            "counts_as_solution := no",
            "counts_as_strict_physical_orientation_datum := no",
        ]
    )

    packet_exported_on_current_repo_state = p1118_candidate_admission_passed and f975_packet_shape_frozen

    add_check("p1118_candidate_admission_passed", p1118_candidate_admission_passed, True, "P1118 already freezes the admission of the existing T178 target as the current best narrow seed candidate beneath F960.")
    add_check("f975_packet_shape_frozen", f975_packet_shape_frozen, True, "F975 freezes the packet shape for that seed-candidate admission explicitly.")
    add_check("packet_exported_on_current_repo_state", packet_exported_on_current_repo_state, True, "Therefore the repo now exports one honest packet freezing the existing T178 target as the current best new narrow provider-class seed candidate beneath F960.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state else
        "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F974_F960_EXISTING_T178_CURRENT_BEST_NEW_NARROW_PROVIDER_CLASS_SEED_CANDIDATE_PACKET"
    )

    artifact = {
        "stage": "F975",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "current_best_new_narrow_provider_class_seed_candidate": T178_TARGET,
        "candidate_grade": "current_best_new_narrow_provider_class_seed_only",
        "active_contract": ACTIVE_CONTRACT,
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
        "packet_exported_on_current_repo_state": artifact["packet_exported_on_current_repo_state"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_best_new_narrow_provider_class_seed_candidate": artifact["current_best_new_narrow_provider_class_seed_candidate"],
        "candidate_grade": artifact["candidate_grade"],
        "active_contract": artifact["active_contract"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
