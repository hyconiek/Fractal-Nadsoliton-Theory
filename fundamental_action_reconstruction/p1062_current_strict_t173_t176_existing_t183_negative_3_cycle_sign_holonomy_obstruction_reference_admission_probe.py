#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P730 = GENERATED / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P687 = GENERATED / "p687_current_strict_t173_global_edge_sign_coherence_solvability_audit_probe_summary.json"
IN_F691 = GENERATED / "f691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_export_packet_summary.json"
IN_F960 = GENERATED / "f960_current_strict_t173_t176_existing_t183_residual_datum_pair12_orbit_direction_selection_bridge_target_packet_summary.json"
IN_N456 = ROOT / "N456_CURRENT_FIRST_STRICT_PHASE_12_AUT_Z12_CANONICAL_12_CYCLE_SUCCESSOR_MAP_NONEXISTENCE_THEOREM.md"
IN_N460 = ROOT / "N460_CURRENT_FIRST_STRICT_AX20_Z12_Z2_DENSITY_OPERATOR_BERRY_HOLONOMY_THETA_SUPPLY_NONDERIVATION_CLOSURE_THEOREM.md"

OUT_JSON = GENERATED / "p1062_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1062_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_admission_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P729, IN_P730, IN_P687, IN_F691, IN_F960, IN_N456, IN_N460]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P1062",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p729 = load_json(IN_P729)
    p730 = load_json(IN_P730)
    p687 = load_json(IN_P687)
    f691 = load_json(IN_F691)
    f960 = load_json(IN_F960)
    n456_text = IN_N456.read_text(encoding="utf-8")
    n460_text = IN_N460.read_text(encoding="utf-8")

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

    add_check(
        "t183_split_localized_as_delta_k_vs_delta_minus_k",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving strict split as delta_k versus delta_-k.",
    )
    add_check(
        "current_direction_free_shannon_lane_still_negative",
        bool(p730.get("current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch")),
        False,
        "P730/N726 already keep the current direction-free Shannon lane nonselective on the surviving branches.",
    )
    add_check(
        "negative_3_cycle_triangle_witness_present",
        bool(p687.get("triangle_witness_present")),
        True,
        "P687 already exports an explicit negative 3-cycle witness.",
    )
    add_check(
        "global_chart_sign_system_unsolved",
        bool(p687.get("sign_system_solvable")),
        False,
        "P687 already proves the global per-chart sign system is not solvable under current fixed representatives.",
    )
    add_check(
        "current_t174_f691_layer_is_not_strict_physical_orientation",
        bool(f691.get("counts_as_strict_physical_orientation_datum")),
        False,
        "F691 keeps the current oriented edge sign lift inside convention scope only.",
    )
    add_check(
        "n456_blocks_canonical_oriented_z12_successor_map",
        "There exists **no** Aut-invariant 12-cycle successor map" in n456_text,
        True,
        "N456 already blocks a canonical oriented Z12 successor map under Aut(Z12) gauge.",
    )
    add_check(
        "n460_blocks_berry_holonomy_theta_overpromotion",
        "Berry/holonomy language cannot be used as strict theta supply" in n460_text,
        True,
        "N460 already blocks naive Berry/holonomy overpromotion into strict-core selector supply.",
    )
    add_check(
        "f960_recommends_new_inversion_sensitive_provider_class",
        "genuinely new inversion-sensitive source-side provider class"
        in str(f960.get("recommended_next_move", "")),
        True,
        "F960 already freezes the active next need as a genuinely new inversion-sensitive source-side provider class.",
    )

    discharged = len(blocking) == 0
    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_T183_NEGATIVE_3_CYCLE_SIGN_HOLONOMY_OBSTRUCTION_REFERENCE_ADMISSION_AUDITED"
        if discharged
        else "P1062_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FRONTIER_STATE"
    )

    artifact = {
        "stage": "P1062",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "P729": str(IN_P729.relative_to(REPO)),
            "P730": str(IN_P730.relative_to(REPO)),
            "P687": str(IN_P687.relative_to(REPO)),
            "F691": str(IN_F691.relative_to(REPO)),
            "F960": str(IN_F960.relative_to(REPO)),
            "N456": str(IN_N456.relative_to(REPO)),
            "N460": str(IN_N460.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "reference_object_id": "Negative3CycleSignHolonomyObstructionReference_global_C_v1_strict_v1",
        "supports_search_for": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
        "triangle_witness_present": True,
        "counts_as_strict_physical_orientation_datum": False,
        "counts_as_strict_berry_holonomy_primitive": False,
        "admissible_as_reference_only": discharged,
        "admissible_as_lawful_supplier": False,
        "reference_grade": "obstruction_reference_only",
        "audit_conclusion": {
            "current_negative_3_cycle_witness_is_real": True,
            "current_negative_3_cycle_witness_is_frontier_relevant": True,
            "current_negative_3_cycle_witness_lawfully_supplies_t183": False,
            "current_negative_3_cycle_witness_is_admissible_as_reference_only": discharged,
            "next_honest_move": (
                "search_one_genuinely_new_inversion_sensitive_source_side_provider_class_"
                "using_the_negative_3_cycle_witness_only_as_an_obstruction_reference_clue"
            ),
        },
        "hard_limits": [
            "No T183 discharge.",
            "No T176 discharge.",
            "No strict physical orientation datum.",
            "No strict Berry/holonomy primitive.",
            "No QW-2191 discharge.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "reference_object_id": artifact["reference_object_id"],
        "supports_search_for": artifact["supports_search_for"],
        "triangle_witness_present": True,
        "counts_as_strict_physical_orientation_datum": False,
        "counts_as_strict_berry_holonomy_primitive": False,
        "admissible_as_reference_only": discharged,
        "admissible_as_lawful_supplier": False,
        "reference_grade": "obstruction_reference_only",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
