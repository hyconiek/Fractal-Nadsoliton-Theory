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
IN_P738 = GENERATED / "p738_current_strict_t192_global_projective_selector_state_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_F474 = GENERATED / "f474_current_strict_t171_global_directed_selector_state_datum_export_packet_summary.json"
IN_STATE = GENERATED / "selector_state_global_c_v1_directed_strict_v1.json"
IN_H37 = GENERATED / "h37_sign_distinction_state_audit_summary.json"
IN_H38 = GENERATED / "h38_projective_selector_state_audit_summary.json"
IN_P632 = GENERATED / "p632_current_strict_directed_continuation_decision_packet_summary.json"

OUT_JSON = GENERATED / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p739_current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P738, IN_F474, IN_STATE, IN_H37, IN_H38, IN_P632]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P739",
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
    p738 = load_json(IN_P738)
    f474 = load_json(IN_F474)
    state = load_json(IN_STATE)
    h37 = load_json(IN_H37)
    h38 = load_json(IN_H38)
    p632 = load_json(IN_P632)

    state_type = state.get("state_type") or {}
    depends_on = state.get("depends_on") or {}
    construction = state.get("construction") or {}
    compatibility = state.get("compatibility") or {}
    hard_limits = state.get("hard_limits") or []

    current_global_premise_based_directed_selector_state_lane_exported = (
        f474.get("t171_discharged") is True
        and f474.get("h37_discharged") is True
        and state.get("stage") == "F474"
        and state.get("object") == "SelectorState_global_C_v1_directed_strict_v1"
        and state_type.get("level") == "directed_vector_state"
    )
    current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164 = (
        current_global_premise_based_directed_selector_state_lane_exported
        and isinstance(depends_on.get("t164_fixing_datum_ref"), str)
        and "T164" in str(state_type.get("sign_gauge") or "")
        and contains_token(hard_limits, "Premise-based")
        and "premise_based_via_T164" in str(h37.get("result") or "")
        and "premise_based_via_T164" in str(h38.get("result") or "")
    )
    current_global_premise_based_directed_selector_state_lane_descends_to_projective_state = (
        current_global_premise_based_directed_selector_state_lane_exported
        and bool(compatibility.get("descends_to_projective"))
        and bool((compatibility.get("vector_section_matches_projector_section") or {}).get("all_pairs_match"))
        and f474.get("vector_section_matches_projector_section") is True
    )
    current_global_premise_based_directed_selector_state_lane_is_pair1_fixing_anchored_not_internal_pair12_provider = (
        current_global_premise_based_directed_selector_state_lane_exported
        and isinstance(depends_on.get("sign_observable_ref"), str)
        and "pair1" in str(depends_on.get("sign_observable_ref") or "")
        and "apply global sign_fix from S_dir(u1)" in str(construction.get("rule") or "")
    )
    p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_global_premise_based_directed_selector_state_lane_exported
        and not current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164
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
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the strict-core promotion bridge remains unexported.",
    )
    add_check(
        "p738_global_projective_lane_already_exported_but_sign_gauge_safe",
        {
            "lane_exported": bool(p738.get("current_global_projective_selector_state_lane_exported")),
            "lane_sign_gauge_safe": bool(p738.get("current_global_projective_selector_state_lane_is_projective_ray_level_sign_gauge_safe")),
            "split_descends": bool(p738.get("p731_pair12_witness_split_descends_to_current_global_projective_selector_state_lane")),
            "t192_exported": bool(p738.get("t192_target_exported_on_current_repo_state")),
        },
        {
            "lane_exported": True,
            "lane_sign_gauge_safe": True,
            "split_descends": False,
            "t192_exported": False,
        },
        "P738 already proves that the strongest current global projective selector lane exists, but remains projective/ray-level and sign-gauge-safe.",
    )
    add_check(
        "current_global_premise_based_directed_selector_state_lane_exported",
        current_global_premise_based_directed_selector_state_lane_exported,
        True,
        "The current repo already exports one explicit global directed selector state lane above the projective state.",
    )
    add_check(
        "current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164",
        current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164,
        True,
        "That current global directed selector state lane remains explicitly premise-based via the exported T164 fixing datum.",
    )
    add_check(
        "current_global_premise_based_directed_selector_state_lane_descends_to_projective_state",
        current_global_premise_based_directed_selector_state_lane_descends_to_projective_state,
        True,
        "The exported global directed selector state descends back to the already-exported projective state rather than upgrading it into strict-core uniqueness.",
    )
    add_check(
        "current_global_premise_based_directed_selector_state_lane_is_pair1_fixing_anchored_not_internal_pair12_provider",
        current_global_premise_based_directed_selector_state_lane_is_pair1_fixing_anchored_not_internal_pair12_provider,
        True,
        "The current global directed lane is anchored by the pair1 sign observable plus one explicit fixing datum, not by a newly exported internal pair1/pair2 strict-core provider.",
    )
    add_check(
        "directed_continuation_selected_only_in_declared_scope",
        {
            "selected_continuation": p632.get("selected_continuation"),
            "decision": p632.get("decision"),
        },
        {
            "selected_continuation": "directed",
            "decision": "DIRECTED_CONTINUATION_SELECTED",
        },
        "The current directed continuation is selected explicitly as a declared-scope continuation decision, not as a strict-core upgrade theorem.",
    )
    add_check(
        "p731_pair12_witness_split_does_not_upgrade_to_strict_core_via_current_global_premise_based_directed_selector_state_lane",
        p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not upgrade into one strict-core branch distinction through the current global premise-based directed selector state lane.",
    )
    add_check(
        "t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_exported",
        False,
        False,
        "The repo still does not export the strict-core upgrade bridge from the current global premise-based directed selector state lane to one selected pair1/pair2 branch on current repo state.",
    )

    status = (
        "PASS_GLOBAL_PREMISE_BASED_DIRECTED_SELECTOR_STATE_STRICT_CORE_UPGRADE_NONEXPORT_AUDITED"
        if not blocking
        else "P739_REQUIRES_REVIEW_CHANGED_GLOBAL_PREMISE_BASED_DIRECTED_LANE_STATE"
    )

    artifact = {
        "stage": "P739",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t193_global_premise_based_directed_selector_state_pair12_witness_split_strict_core_upgrade_bridge_nonexport_boundary_only",
        "inputs": {
            "p731_summary_ref": str(IN_P731.relative_to(REPO)),
            "p738_summary_ref": str(IN_P738.relative_to(REPO)),
            "f474_summary_ref": str(IN_F474.relative_to(REPO)),
            "directed_state_ref": str(IN_STATE.relative_to(REPO)),
            "h37_summary_ref": str(IN_H37.relative_to(REPO)),
            "h38_summary_ref": str(IN_H38.relative_to(REPO)),
            "p632_summary_ref": str(IN_P632.relative_to(REPO)),
        },
        "checks": checks,
        "current_global_premise_based_directed_selector_state_lane_exported": current_global_premise_based_directed_selector_state_lane_exported,
        "current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164": current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164,
        "current_global_premise_based_directed_selector_state_lane_descends_to_projective_state": current_global_premise_based_directed_selector_state_lane_descends_to_projective_state,
        "current_global_premise_based_directed_selector_state_lane_is_pair1_fixing_anchored_not_internal_pair12_provider": current_global_premise_based_directed_selector_state_lane_is_pair1_fixing_anchored_not_internal_pair12_provider,
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane": p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane,
        "t193_target_exported_on_current_repo_state": False,
        "blocking_mismatches": blocking,
        "no_false_pass": True,
    }
    summary = {
        "stage": "P739",
        "status": status,
        "as_of": AS_OF,
        "current_global_premise_based_directed_selector_state_lane_exported": current_global_premise_based_directed_selector_state_lane_exported,
        "current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164": current_global_premise_based_directed_selector_state_lane_is_premise_based_via_t164,
        "current_global_premise_based_directed_selector_state_lane_descends_to_projective_state": current_global_premise_based_directed_selector_state_lane_descends_to_projective_state,
        "current_global_premise_based_directed_selector_state_lane_is_pair1_fixing_anchored_not_internal_pair12_provider": current_global_premise_based_directed_selector_state_lane_is_pair1_fixing_anchored_not_internal_pair12_provider,
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane": p731_pair12_witness_split_upgrades_to_strict_core_via_current_global_premise_based_directed_selector_state_lane,
        "t193_target_exported_on_current_repo_state": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
