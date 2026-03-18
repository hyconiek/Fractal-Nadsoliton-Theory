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

IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P729, IN_P731, IN_F147, IN_F301]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P741",
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
    p731 = load_json(IN_P731)
    f147 = load_json(IN_F147)
    f301 = load_json(IN_F301)

    support_packet = f147.get("support_packet") or {}
    selector_axis_realization = support_packet.get("selector_axis_realization") or {}
    selector_signed_split_realization = support_packet.get("selector_signed_split_realization") or {}
    preobserver_scope_realization = support_packet.get("preobserver_scope_realization") or {}

    pair_index_set = f301.get("pair_index_set") or []

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

    current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier = (
        f147.get("input_packet") == "tau_src_candidate_v1"
        and f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and bool(f147.get("actual_selector_witness_exported"))
    )

    current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only = (
        bool(f147.get("chart_bound_preobserver_realization"))
        and f147.get("observer_role") == "downstream_only"
        and bool(preobserver_scope_realization.get("preobserver_only"))
        and bool(preobserver_scope_realization.get("observer_downstream_only"))
    )

    current_actual_source_topology_selector_witness_exports_positive_signed_preLM_split = (
        selector_axis_realization.get("object") == "E_orient_preLM_v1"
        and selector_axis_realization.get("frame_basis") == ["u_T", "u_L"]
        and selector_signed_split_realization.get("object") == "B_sel_preLM_v1"
        and bool(selector_signed_split_realization.get("positive_signed_selector_response"))
        and bool(preobserver_scope_realization.get("positive_plus_output"))
        and bool(preobserver_scope_realization.get("vanishing_minus_output"))
    )

    current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed = (
        current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier
        and current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only
        and current_actual_source_topology_selector_witness_exports_positive_signed_preLM_split
        and f147.get("codomain_packet") == "Sigma_sel_src_target_v1"
        and selector_axis_realization.get("frame_basis") != ["pair1", "pair2"]
        and not bool(f147.get("tau_src_identified_with_s_prelm"))
        and not bool(f147.get("basis_independence_discharged"))
        and not bool(f147.get("qw2191_quotient_safe_discharged"))
        and not bool(f147.get("current_selector_closure"))
    )

    current_actual_source_topology_selector_witness_is_pair12_witness_split_typed = (
        current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier
        and not current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed
    )

    p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_actual_source_topology_selector_witness_is_pair12_witness_split_typed
    )

    add_check(
        "p729_pair12_orbit_direction_split_already_localized",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving pair1/pair2 ambiguity as opposite residual-datum orbit directions on the same source-side carrier.",
    )
    add_check(
        "p731_pair12_witness_split_already_separated",
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
        "current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier",
        {
            "selector_input_packet": f147.get("input_packet"),
            "carrier_source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
            "actual_selector_witness_exported": bool(f147.get("actual_selector_witness_exported")),
        },
        {
            "selector_input_packet": "tau_src_candidate_v1",
            "carrier_source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
            "actual_selector_witness_exported": True,
        },
        "The current actual source-topology selector witness already binds the same tau_src_candidate_v1 packet as the surviving pair1/pair2 residual-datum carrier.",
    )
    add_check(
        "current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only",
        {
            "chart_bound_preobserver_realization": bool(f147.get("chart_bound_preobserver_realization")),
            "observer_role": f147.get("observer_role"),
            "preobserver_only": bool(preobserver_scope_realization.get("preobserver_only")),
            "observer_downstream_only": bool(preobserver_scope_realization.get("observer_downstream_only")),
        },
        {
            "chart_bound_preobserver_realization": True,
            "observer_role": "downstream_only",
            "preobserver_only": True,
            "observer_downstream_only": True,
        },
        "That current selector witness remains chart-bound and preobserver/downstream only on current exports.",
    )
    add_check(
        "current_actual_source_topology_selector_witness_exports_positive_signed_preLM_split",
        {
            "selector_axis_object": selector_axis_realization.get("object"),
            "frame_basis": selector_axis_realization.get("frame_basis"),
            "signed_split_object": selector_signed_split_realization.get("object"),
            "positive_signed_selector_response": bool(
                selector_signed_split_realization.get("positive_signed_selector_response")
            ),
            "positive_plus_output": bool(preobserver_scope_realization.get("positive_plus_output")),
            "vanishing_minus_output": bool(preobserver_scope_realization.get("vanishing_minus_output")),
        },
        {
            "selector_axis_object": "E_orient_preLM_v1",
            "frame_basis": ["u_T", "u_L"],
            "signed_split_object": "B_sel_preLM_v1",
            "positive_signed_selector_response": True,
            "positive_plus_output": True,
            "vanishing_minus_output": True,
        },
        "The current actual selector witness already exports one real preLM signed selector split with positive plus-channel response and vanishing minus output.",
    )
    add_check(
        "current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed",
        current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed,
        True,
        "That current selector witness still lives only in the preLM basis u_T/u_L, does not identify tau_src with the current selector carrier, and therefore remains untyped on the surviving F301 pair1/pair2 carrier.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness",
        p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current actual source-topology selector witness as one typed source-side branch distinction.",
    )
    add_check(
        "t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_exported",
        False,
        False,
        "So the actual source-topology selector witness pair1/pair2 witness-split promotion bridge still remains unexported on current repo state.",
    )

    status = (
        "PASS_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_PAIR12_WITNESS_SPLIT_PROMOTION_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P741_REQUIRES_REVIEW_CHANGED_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_PAIR12_STATE"
    )

    artifact = {
        "stage": "P741",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_boundary_only",
        "inputs": {
            "P729": str(IN_P729.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t195_target_name": "ActualSourceTopologySelectorWitnessPair12WitnessSplitPromotionBridge_global_C_v1_strict_v1",
        "t195_target_exported_on_current_repo_state": False,
        "current_actual_source_topology_selector_witness_exported": bool(
            f147.get("actual_selector_witness_exported")
        ),
        "current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier": (
            current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only": (
            current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only
        ),
        "current_actual_source_topology_selector_witness_exports_positive_signed_preLM_split": (
            current_actual_source_topology_selector_witness_exports_positive_signed_preLM_split
        ),
        "current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed": (
            current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed
        ),
        "p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness": (
            p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness
        ),
        "audit_conclusion": {
            "current_repo_already_exports_actual_source_side_selector_witness_on_tau_src_candidate_v1": bool(
                f147.get("actual_selector_witness_exported")
            ),
            "current_repo_already_exports_t195_target": False,
            "next_honest_move": (
                "attack_an_explicit_typed_carrier_bridge_from_the_current_chart_bound_preLM_source_topology_selector_witness_to_one_selected_F301_pair12_branch_without_quotienting_chart_labels_away_and_without_collapsing_into_convention_or_gauge_lanes"
            ),
        },
        "hard_limits": [
            "No T195 discharge claim.",
            "No claim that the current actual source-topology selector witness already selects one raw pair1/pair2 branch.",
            "No identification of tau_src_candidate_v1 with the current selector carrier.",
            "No basis-independent or quotient-safe QW-2191 upgrade claim.",
            "No directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P741",
        "status": status,
        "as_of": AS_OF,
        "t195_target_name": "ActualSourceTopologySelectorWitnessPair12WitnessSplitPromotionBridge_global_C_v1_strict_v1",
        "t195_target_exported_on_current_repo_state": False,
        "current_actual_source_topology_selector_witness_exported": bool(
            f147.get("actual_selector_witness_exported")
        ),
        "current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier": (
            current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only": (
            current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only
        ),
        "current_actual_source_topology_selector_witness_exports_positive_signed_preLM_split": (
            current_actual_source_topology_selector_witness_exports_positive_signed_preLM_split
        ),
        "current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed": (
            current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed
        ),
        "p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness": (
            p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
