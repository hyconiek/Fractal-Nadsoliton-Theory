#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"
TOL = 1.0e-12

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P737 = GENERATED / "p737_current_strict_t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_F469 = GENERATED / "selector_transition_global_c_v1_strict_v1.json"
IN_F470 = GENERATED / "selector_state_global_c_v1_projective_strict_v1.json"
IN_P474 = GENERATED / "p474_current_strict_global_projective_selector_state_gluing_consistency_audit_probe_summary.json"

OUT_JSON = GENERATED / "p738_current_strict_t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p738_current_strict_t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def small_enough(value: Any, tol: float = TOL) -> bool:
    return isinstance(value, (int, float)) and abs(float(value)) <= tol


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P737, IN_F469, IN_F470, IN_P474]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P738",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p731 = load_json(IN_P731)
    p737 = load_json(IN_P737)
    f469 = load_json(IN_F469)
    f470 = load_json(IN_F470)
    p474 = load_json(IN_P474)

    transition_cocycle = f469.get("cocycle_discipline") or {}
    transition_hard_limits = f469.get("hard_limits") or []
    transition_non_promotion = transition_cocycle.get("explicit_non_promotion") or []

    state_type = f470.get("state_type") or {}
    state_gluing = f470.get("gluing") or {}
    state_hard_limits = f470.get("hard_limits") or []

    current_global_selector_transition_object_exported = (
        f469.get("stage") == "F469"
        and f469.get("object") == "SelectorTransition_global_C_v1_strict_v1"
        and f469.get("charts") == ["pair1", "pair2", "pair3", "pair4", "pair5"]
        and isinstance(f469.get("gluing_laws_ref"), str)
    )
    current_global_selector_transition_object_is_projector_section_only = (
        current_global_selector_transition_object_exported
        and transition_cocycle.get("level") == "projector_section_level"
        and contains_token(transition_hard_limits, "Projector-section cocycle only")
        and any("sign-sensitive physical orientation datum" in str(entry) for entry in transition_non_promotion)
    )
    current_global_projective_selector_state_object_exported = (
        f470.get("stage") == "F470"
        and f470.get("object") == "SelectorState_global_C_v1_projective_strict_v1"
        and state_type.get("level") == "projective_ray_state"
        and f470.get("global_transition_object") == "SelectorTransition_global_C_v1_strict_v1"
        and sorted((f470.get("charts") or {}).keys()) == ["pair1", "pair2", "pair3", "pair4", "pair5"]
    )
    current_global_projective_selector_state_is_projective_ray_level_sign_gauge_safe = (
        current_global_projective_selector_state_object_exported
        and "u and -u identified" in str(state_type.get("sign_gauge") or "")
        and bool(state_gluing.get("projector_section_level"))
        and contains_token(state_hard_limits, "Projective/ray-level object only")
    )
    current_global_projective_selector_state_gluing_consistent = (
        bool(p474.get("overall_pass"))
        and small_enough(p474.get("max_projector_transport_max_abs_residual"))
        and small_enough(p474.get("max_projector_cocycle_max_abs_residual"))
    )
    current_global_projective_selector_state_lane_exported = (
        current_global_selector_transition_object_exported
        and current_global_projective_selector_state_object_exported
        and current_global_projective_selector_state_gluing_consistent
    )
    current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe = (
        current_global_projective_selector_state_lane_exported
        and current_global_selector_transition_object_is_projector_section_only
        and current_global_projective_selector_state_is_projective_ray_level_sign_gauge_safe
    )
    current_global_projective_selector_state_lane_is_pair12_orbit_direction_typed = (
        current_global_projective_selector_state_lane_exported
        and not current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe
    )
    p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_global_projective_selector_state_lane_is_pair12_orbit_direction_typed
    )

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
        "p731_pair12_witness_split_already_opposite_and_unpromoted",
        {
            "split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "pair1_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_sign": p731.get("pair2_w_break_branch_score_sign"),
            "t185_exported": bool(p731.get("t185_target_exported_on_current_repo_state")),
        },
        {
            "split_separated": True,
            "pair1_sign": "negative",
            "pair2_sign": "positive",
            "t185_exported": False,
        },
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the typed promotion bridge remains unexported.",
    )
    add_check(
        "p737_local_pair12_projector_lane_already_exported_but_sign_gauge_safe",
        {
            "lane_exported": bool(p737.get("current_local_pair12_projector_atlas_glue_lane_exported")),
            "lane_sign_gauge_safe": bool(p737.get("current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe")),
            "split_descends": bool(p737.get("p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane")),
            "t191_exported": bool(p737.get("t191_target_exported_on_current_repo_state")),
        },
        {
            "lane_exported": True,
            "lane_sign_gauge_safe": True,
            "split_descends": False,
            "t191_exported": False,
        },
        "P737 already proves that the strongest current local pair1/pair2 atlas/glue lane exists, but remains projector-level and sign-gauge-safe.",
    )
    add_check(
        "current_global_selector_transition_object_exported_and_projector_section_only",
        {
            "exported": current_global_selector_transition_object_exported,
            "projector_section_only": current_global_selector_transition_object_is_projector_section_only,
        },
        {
            "exported": True,
            "projector_section_only": True,
        },
        "The current repo already exports one explicit global selector transition object, but its cocycle/gluing discipline remains projector-section-level only.",
    )
    add_check(
        "current_global_projective_selector_state_object_exported_and_projective_sign_gauge_safe",
        {
            "exported": current_global_projective_selector_state_object_exported,
            "projective_ray_level": current_global_projective_selector_state_is_projective_ray_level_sign_gauge_safe,
        },
        {
            "exported": True,
            "projective_ray_level": True,
        },
        "The current repo already exports one explicit global projective selector state object on C_v1, but only at projective/ray level with sign-gauge-safe semantics.",
    )
    add_check(
        "current_global_projective_selector_state_gluing_consistent",
        current_global_projective_selector_state_gluing_consistent,
        True,
        "The exported global projective selector state remains globally glued consistently on current artifacts.",
    )
    add_check(
        "current_global_projective_selector_state_lane_exported",
        current_global_projective_selector_state_lane_exported,
        True,
        "So the current repo already exports one real global projective selector transition/state lane on C_v1.",
    )
    add_check(
        "current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe",
        current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe,
        True,
        "That current global projective selector transition/state lane remains projective/ray-level and sign-gauge-safe on current exports.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane",
        p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current global projective selector transition/state lane as one typed branch distinction.",
    )
    add_check(
        "t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_exported",
        False,
        False,
        "The repo still does not export the global projective selector state pair1/pair2 orbit-direction descent bridge on current repo state.",
    )

    status = (
        "PASS_GLOBAL_PROJECTIVE_SELECTOR_STATE_PAIR12_ORBIT_DIRECTION_DESCENT_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P738_REQUIRES_REVIEW_CHANGED_GLOBAL_PROJECTIVE_SELECTOR_STATE_LANE_STATE"
    )

    artifact = {
        "stage": "P738",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_nonexport_boundary_only",
        "inputs": {
            "p731_summary_ref": str(IN_P731.relative_to(REPO)),
            "p737_summary_ref": str(IN_P737.relative_to(REPO)),
            "global_selector_transition_ref": str(IN_F469.relative_to(REPO)),
            "global_projective_selector_state_ref": str(IN_F470.relative_to(REPO)),
            "p474_summary_ref": str(IN_P474.relative_to(REPO)),
        },
        "checks": checks,
        "current_global_selector_transition_object_exported": current_global_selector_transition_object_exported,
        "current_global_selector_transition_object_is_projector_section_only": current_global_selector_transition_object_is_projector_section_only,
        "current_global_projective_selector_state_object_exported": current_global_projective_selector_state_object_exported,
        "current_global_projective_selector_state_is_projective_ray_level_sign_gauge_safe": current_global_projective_selector_state_is_projective_ray_level_sign_gauge_safe,
        "current_global_projective_selector_state_gluing_consistent": current_global_projective_selector_state_gluing_consistent,
        "current_global_projective_selector_state_lane_exported": current_global_projective_selector_state_lane_exported,
        "current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe": current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe,
        "p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane": p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane,
        "t192_target_exported_on_current_repo_state": False,
        "blocking_mismatches": blocking,
        "no_false_pass": True,
    }
    summary = {
        "stage": "P738",
        "status": status,
        "as_of": AS_OF,
        "current_global_selector_transition_object_exported": current_global_selector_transition_object_exported,
        "current_global_selector_transition_object_is_projector_section_only": current_global_selector_transition_object_is_projector_section_only,
        "current_global_projective_selector_state_object_exported": current_global_projective_selector_state_object_exported,
        "current_global_projective_selector_state_is_projective_ray_level_sign_gauge_safe": current_global_projective_selector_state_is_projective_ray_level_sign_gauge_safe,
        "current_global_projective_selector_state_gluing_consistent": current_global_projective_selector_state_gluing_consistent,
        "current_global_projective_selector_state_lane_exported": current_global_projective_selector_state_lane_exported,
        "current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe": current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe,
        "p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane": p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane,
        "t192_target_exported_on_current_repo_state": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
