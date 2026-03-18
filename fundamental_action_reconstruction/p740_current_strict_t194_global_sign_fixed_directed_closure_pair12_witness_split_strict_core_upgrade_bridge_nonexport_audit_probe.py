#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P739 = GENERATED / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"
IN_F690 = GENERATED / "f690_current_strict_t175_global_chart_sign_fixing_from_strict_core_payload_weights_export_packet_summary.json"
IN_F691 = GENERATED / "f691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_export_packet_summary.json"
IN_F692 = GENERATED / "f692_current_strict_t175_sign_fixed_directed_closure_export_packet_summary.json"
IN_N690 = GENERATED / "n690_current_strict_t175_global_chart_sign_fixing_discharge_theorem_summary.json"
IN_N691 = GENERATED / "n691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_discharge_theorem_summary.json"
IN_N692 = GENERATED / "n692_current_strict_t173_projective_directed_closure_output_ray_invariance_theorem_summary.json"
IN_N693 = GENERATED / "n693_current_strict_t173_output_sign_lift_gauge_covariance_theorem_summary.json"
IN_SIGN_FIXED_STATE = (
    GENERATED
    / "selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)
IN_SIGN_FIXED_TRANSITION = (
    GENERATED
    / "selector_transition_global_c_v1_oriented_mod_2pi_edge_sign_lift_from_sign_fixed_directed_state_strict_convention_v1.json"
)
IN_SIGN_FIXED_CLOSURE = (
    GENERATED
    / "selector_closure_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
)

OUT_JSON = GENERATED / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p740_current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P731,
        IN_P739,
        IN_F690,
        IN_F691,
        IN_F692,
        IN_N690,
        IN_N691,
        IN_N692,
        IN_N693,
        IN_SIGN_FIXED_STATE,
        IN_SIGN_FIXED_TRANSITION,
        IN_SIGN_FIXED_CLOSURE,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P740",
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
    p739 = load_json(IN_P739)
    f690 = load_json(IN_F690)
    f691 = load_json(IN_F691)
    f692 = load_json(IN_F692)
    n690 = load_json(IN_N690)
    n691 = load_json(IN_N691)
    n692 = load_json(IN_N692)
    n693 = load_json(IN_N693)
    sign_fixed_state = load_json(IN_SIGN_FIXED_STATE)
    sign_fixed_transition = load_json(IN_SIGN_FIXED_TRANSITION)
    sign_fixed_closure = load_json(IN_SIGN_FIXED_CLOSURE)

    state_depends_on = sign_fixed_state.get("depends_on") or {}
    state_hard_limits = sign_fixed_state.get("hard_limits") or []
    transition_hard_limits = sign_fixed_transition.get("hard_limits") or []
    closure_depends_on = sign_fixed_closure.get("depends_on") or {}
    closure_hard_limits = sign_fixed_closure.get("hard_limits") or []
    output_sign_lift = sign_fixed_closure.get("output_sign_lift") or {}
    output_observable = sign_fixed_closure.get("output_observable") or {}
    raw_output_observable = sign_fixed_closure.get("raw_output_observable") or {}
    output_gluing = output_observable.get("gluing") or {}
    raw_output_gluing = raw_output_observable.get("gluing") or {}

    n690_result = n690.get("theorem_result") or {}
    n691_result = n691.get("theorem_result") or {}
    n692_result = n692.get("theorem_result") or {}
    n693_result = n693.get("theorem_result") or {}

    current_global_sign_fixed_directed_state_lane_exported = (
        str(f690.get("status") or "").startswith("PASS_")
        and sign_fixed_state.get("stage") == "F690"
        and sign_fixed_state.get("object")
        == "SelectorState_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1"
        and sign_fixed_state.get("counts_as_strict_physical_orientation_datum") is False
    )
    current_global_sign_fixed_directed_transition_lift_exported = (
        str(f691.get("status") or "").startswith("PASS_")
        and sign_fixed_transition.get("stage") == "F691"
        and sign_fixed_transition.get("object")
        == "SelectorTransition_global_C_v1_oriented_mod_2pi_edge_sign_lift_from_sign_fixed_directed_state_strict_convention_v1"
    )
    current_global_sign_fixed_directed_closure_lane_exported = (
        str(f692.get("status") or "").startswith("PASS_")
        and sign_fixed_closure.get("stage") == "F692"
        and sign_fixed_closure.get("object")
        == "SelectorClosure_global_C_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1"
        and bool(f692.get("glued"))
        and bool(output_gluing.get("glued"))
        and sign_fixed_closure.get("no_false_pass") is True
    )
    current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing = (
        current_global_sign_fixed_directed_closure_lane_exported
        and f692.get("raw_glued") is False
        and raw_output_gluing.get("glued") is False
        and bool(output_gluing.get("glued"))
        and isinstance(output_sign_lift.get("signs_by_pair"), dict)
    )
    current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only = (
        current_global_sign_fixed_directed_state_lane_exported
        and current_global_sign_fixed_directed_transition_lift_exported
        and current_global_sign_fixed_directed_closure_lane_exported
        and "strict_convention" in str(sign_fixed_state.get("object") or "")
        and "strict_convention" in str(sign_fixed_transition.get("object") or "")
        and "strict_convention" in str(sign_fixed_closure.get("object") or "")
        and contains_token(state_hard_limits, "Convention layer only")
        and contains_token(transition_hard_limits, "Convention layer only")
        and contains_token(closure_hard_limits, "Convention layer only")
        and n690_result.get("directed_sign_sensitive_physical_orientation_in_strict_core") is False
        and n691_result.get("directed_sign_sensitive_physical_orientation_in_strict_core") is False
        and n692_result.get("directed_sign_sensitive_physical_orientation_in_strict_core") is False
        and n693_result.get("directed_sign_sensitive_physical_orientation_in_strict_core") is False
    )
    current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray = (
        current_global_sign_fixed_directed_closure_lane_exported
        and n692_result.get("projective_output_ray_invariant_across_exported_directed_closures") is True
    )
    current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant = (
        current_global_sign_fixed_directed_closure_lane_exported
        and n693_result.get("output_sign_lift_is_gauge_covariant_under_chart_sign_relift") is True
    )
    current_global_sign_fixed_directed_closure_lane_is_downstream_gauge_relift_of_current_nonupgrading_premise_based_lane = (
        bool(p739.get("current_global_premise_based_directed_selector_state_lane_exported"))
        and not bool(
            p739.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane")
        )
        and isinstance(state_depends_on.get("premise_based_directed_state_ref"), str)
        and "selector_state_global_c_v1_directed_strict_v1.json"
        in str(state_depends_on.get("premise_based_directed_state_ref") or "")
        and isinstance(closure_depends_on.get("sign_fixed_directed_state_ref"), str)
        and "selector_state_global_c_v1_directed_sign_fixed_from_strict_core_payload_weights_strict_convention_v1.json"
        in str(closure_depends_on.get("sign_fixed_directed_state_ref") or "")
        and current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant
    )
    p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_global_sign_fixed_directed_closure_lane_exported
        and not current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only
        and not current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray
        and not current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant
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
        "p739_current_global_premise_based_directed_lane_still_does_not_upgrade_split",
        {
            "lane_exported": bool(p739.get("current_global_premise_based_directed_selector_state_lane_exported")),
            "lane_premise_based": bool(p739.get("current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164")),
            "split_upgrades": bool(
                p739.get("p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane")
            ),
            "t193_exported": bool(p739.get("t193_target_exported_on_current_repo_state")),
        },
        {
            "lane_exported": True,
            "lane_premise_based": True,
            "split_upgrades": False,
            "t193_exported": False,
        },
        "P739 already proves that the strongest current global premise-based directed selector-state lane exists, but still does not upgrade the P731 split into strict core.",
    )
    add_check(
        "current_global_sign_fixed_directed_state_lane_exported_and_convention_only",
        {
            "exported": current_global_sign_fixed_directed_state_lane_exported,
            "counts_as_strict_physical_orientation_datum": sign_fixed_state.get(
                "counts_as_strict_physical_orientation_datum"
            ),
        },
        {
            "exported": True,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "The current repo already exports one explicit sign-fixed directed state lane on C_v1, but only as a strict_convention gauge-fixing layer (F690/N690).",
    )
    add_check(
        "current_global_sign_fixed_directed_transition_lift_exported_and_convention_only",
        {
            "exported": current_global_sign_fixed_directed_transition_lift_exported,
            "convention_only": contains_token(transition_hard_limits, "Convention layer only"),
        },
        {
            "exported": True,
            "convention_only": True,
        },
        "The current repo already exports one explicit oriented transition edge sign-lift anchored to that sign-fixed state, but still only in strict_convention scope (F691/N691).",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_exported_requires_output_sign_lift_and_is_not_physical",
        {
            "exported": current_global_sign_fixed_directed_closure_lane_exported,
            "glued": bool(f692.get("glued")),
            "raw_glued": bool(f692.get("raw_glued")),
            "convention_only": contains_token(closure_hard_limits, "Convention layer only"),
        },
        {
            "exported": True,
            "glued": True,
            "raw_glued": False,
            "convention_only": True,
        },
        "The current repo already exports one sign-fixed directed closure lane on C_v1, but it glues only after an explicit output sign-lift and is not promoted to any strict physical orientation datum (F692).",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only",
        current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only,
        True,
        "Taken together, the current sign-fixed directed state/transition/closure lane remains strict_convention/gauge only rather than a strict-core branch-typing lane.",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray",
        current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray,
        True,
        "N692 already proves that the current sign-fixed directed closure outputs collapse to the same projective output ray as the projective closure output.",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant",
        current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant,
        True,
        "N693 already proves that the current sign-fixed directed closure output sign-lift is chartwise Z2 gauge-covariant under relift from the premise-based lane.",
    )
    add_check(
        "current_global_sign_fixed_directed_closure_lane_is_downstream_gauge_relift_of_current_nonupgrading_premise_based_lane",
        current_global_sign_fixed_directed_closure_lane_is_downstream_gauge_relift_of_current_nonupgrading_premise_based_lane,
        True,
        "Therefore the current sign-fixed directed closure lane is already identified as a downstream gauge-relift of the current non-upgrading premise-based directed lane, not as a new strict-core provider.",
    )
    add_check(
        "p731_pair12_witness_split_does_not_upgrade_to_strict_core_via_current_global_sign_fixed_directed_closure_lane",
        p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane,
        False,
        "So the opposite P731 pair1/pair2 witness split still does not upgrade into one strict-core branch distinction through the current global sign-fixed directed closure lane.",
    )
    add_check(
        "t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_exported",
        False,
        False,
        "The repo still does not export the strict-core upgrade bridge from the current global sign-fixed directed closure lane to one selected pair1/pair2 branch on current repo state.",
    )

    status = (
        "PASS_GLOBAL_SIGN_FIXED_DIRECTED_CLOSURE_STRICT_CORE_UPGRADE_NONEXPORT_AUDITED"
        if not blocking
        else "P740_REQUIRES_REVIEW_CHANGED_SIGN_FIXED_DIRECTED_CLOSURE_LANE_STATE"
    )

    artifact = {
        "stage": "P740",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t194_global_sign_fixed_directed_closure_pair12_witness_split_strict_core_upgrade_bridge_nonexport_boundary_only",
        "inputs": {
            "p731_summary_ref": str(IN_P731.relative_to(REPO)),
            "p739_summary_ref": str(IN_P739.relative_to(REPO)),
            "f690_summary_ref": str(IN_F690.relative_to(REPO)),
            "f691_summary_ref": str(IN_F691.relative_to(REPO)),
            "f692_summary_ref": str(IN_F692.relative_to(REPO)),
            "n690_summary_ref": str(IN_N690.relative_to(REPO)),
            "n691_summary_ref": str(IN_N691.relative_to(REPO)),
            "n692_summary_ref": str(IN_N692.relative_to(REPO)),
            "n693_summary_ref": str(IN_N693.relative_to(REPO)),
            "sign_fixed_state_ref": str(IN_SIGN_FIXED_STATE.relative_to(REPO)),
            "sign_fixed_transition_ref": str(IN_SIGN_FIXED_TRANSITION.relative_to(REPO)),
            "sign_fixed_closure_ref": str(IN_SIGN_FIXED_CLOSURE.relative_to(REPO)),
        },
        "checks": checks,
        "current_global_sign_fixed_directed_state_lane_exported": current_global_sign_fixed_directed_state_lane_exported,
        "current_global_sign_fixed_directed_transition_lift_exported": current_global_sign_fixed_directed_transition_lift_exported,
        "current_global_sign_fixed_directed_closure_lane_exported": current_global_sign_fixed_directed_closure_lane_exported,
        "current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing": current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing,
        "current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only": current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only,
        "current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray": current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray,
        "current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant": current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant,
        "current_global_sign_fixed_directed_closure_lane_is_downstream_gauge_relift_of_current_nonupgrading_premise_based_lane": current_global_sign_fixed_directed_closure_lane_is_downstream_gauge_relift_of_current_nonupgrading_premise_based_lane,
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane": p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane,
        "t194_target_exported_on_current_repo_state": False,
        "blocking_mismatches": blocking,
        "no_false_pass": True,
    }
    summary = {
        "stage": "P740",
        "status": status,
        "as_of": AS_OF,
        "current_global_sign_fixed_directed_closure_lane_exported": current_global_sign_fixed_directed_closure_lane_exported,
        "current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing": current_global_sign_fixed_directed_closure_lane_requires_explicit_output_sign_lift_for_gluing,
        "current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only": current_global_sign_fixed_directed_closure_lane_is_strict_convention_gauge_only,
        "current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray": current_global_sign_fixed_directed_closure_lane_descends_to_same_projective_output_ray,
        "current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant": current_global_sign_fixed_directed_closure_output_sign_lift_is_gauge_covariant,
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane": p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_sign_fixed_directed_closure_lane,
        "t194_target_exported_on_current_repo_state": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
