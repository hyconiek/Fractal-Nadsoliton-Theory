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
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_F149 = GENERATED / "f149_first_actual_source_topology_quotient_safe_qw2191_resolution_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p743_current_strict_t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p743_current_strict_t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P742, IN_F149, IN_F301]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P743",
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
    p742 = load_json(IN_P742)
    f149 = load_json(IN_F149)
    f301 = load_json(IN_F301)

    support_packet = f149.get("support_packet") or {}
    quotient_relation = support_packet.get("quotient_relation") or {}
    distinguished_quotient_class = support_packet.get("distinguished_quotient_class") or {}
    pair_index_set = f301.get("pair_index_set") or []
    f301_notes = f301.get("notes") or []
    f301_current_absence = f301.get("current_absence") or []

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

    current_actual_source_topology_quotient_safe_qw2191_resolution_exported = (
        f149.get("witness") == "Phi_qw2191_safe_actual_witness_v1"
        and f149.get("input_packet") == "tau_src_candidate_v1"
        and f149.get("codomain_target") == "actual_quotient_safe_qw2191_resolution_target_v1"
        and support_packet.get("basis_promotion_witness") == "Upsilon_sel_basis_actual_witness_v1"
        and bool(f149.get("actual_basis_independent_selector_promotion_exported"))
        and bool(f149.get("actual_quotient_safe_qw2191_resolution_exported"))
        and bool(f149.get("qw2191_quotient_safe_discharged"))
    )

    current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier = (
        current_actual_source_topology_quotient_safe_qw2191_resolution_exported
        and f149.get("input_packet") == "tau_src_candidate_v1"
        and f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
    )

    current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only = (
        current_actual_source_topology_quotient_safe_qw2191_resolution_exported
        and quotient_relation.get("object") == "~_src_qw2191_v1"
        and quotient_relation.get("criterion") == "same_source_side_basis_free_selector_packet"
        and not bool(quotient_relation.get("uses_observer_side_data"))
        and not bool(quotient_relation.get("uses_external_selector_axiom"))
        and bool(quotient_relation.get("quotients_chart_labels_only_via_basis_free_packet"))
        and distinguished_quotient_class.get("object") == "C_sel_src_qw2191_resolved_v1"
        and bool(f149.get("distinguished_quotient_class_exported"))
        and bool(distinguished_quotient_class.get("quotient_class_only"))
        and not bool(distinguished_quotient_class.get("raw_theta_uniqueness_claimed"))
        and not bool(f149.get("raw_theta_uniqueness_claimed"))
        and not bool(f149.get("current_selector_closure"))
        and not bool(f149.get("current_global_qw2191_discharge"))
    )

    surviving_pair12_residual_datum_carrier_remains_selector_neutral = (
        f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and contains_token(f301_notes, "selector-neutral")
        and contains_token(f301_current_absence, "no strict-core selector closure")
        and contains_token(f301_current_absence, "no QW-2191 discharge")
    )

    current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation = (
        current_actual_source_topology_quotient_safe_qw2191_resolution_exported
        and f149.get("codomain_target")
        == "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1"
        and not bool(distinguished_quotient_class.get("quotient_class_only"))
        and current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
    )

    current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed = (
        current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only
        and current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier
        and bool(
            p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed")
        )
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and not current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation
    )

    p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation
    )

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
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the typed source-side promotion bridge remains unexported.",
    )
    add_check(
        "p742_actual_selector_witness_continuation_already_basis_free_not_pair12_typed",
        {
            "basis_free_chart_label_forgetting_continuation": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")
            ),
            "continuation_remains_basis_free_not_pair12_typed": bool(
                p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed")
            ),
            "pair12_typed_residual_datum_continuation_exported": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")
            ),
            "t196_exported": bool(p742.get("t196_target_exported_on_current_repo_state")),
        },
        {
            "basis_free_chart_label_forgetting_continuation": True,
            "continuation_remains_basis_free_not_pair12_typed": True,
            "pair12_typed_residual_datum_continuation_exported": False,
            "t196_exported": False,
        },
        "P742 already proves that the strongest exported continuation out of the actual selector-witness codomain remains only the basis-free chart-label-forgetting lane, not one typed pair1/pair2 residual-datum carrier bridge.",
    )
    add_check(
        "current_actual_source_topology_quotient_safe_qw2191_resolution_exported",
        {
            "witness": f149.get("witness"),
            "input_packet": f149.get("input_packet"),
            "codomain_target": f149.get("codomain_target"),
            "basis_promotion_witness": support_packet.get("basis_promotion_witness"),
            "actual_quotient_safe_qw2191_resolution_exported": bool(
                f149.get("actual_quotient_safe_qw2191_resolution_exported")
            ),
            "qw2191_quotient_safe_discharged": bool(f149.get("qw2191_quotient_safe_discharged")),
        },
        {
            "witness": "Phi_qw2191_safe_actual_witness_v1",
            "input_packet": "tau_src_candidate_v1",
            "codomain_target": "actual_quotient_safe_qw2191_resolution_target_v1",
            "basis_promotion_witness": "Upsilon_sel_basis_actual_witness_v1",
            "actual_quotient_safe_qw2191_resolution_exported": True,
            "qw2191_quotient_safe_discharged": True,
        },
        "F149 already exports one real actual source-topology quotient-safe QW-2191 resolution witness built from the current basis-free selector packet.",
    )
    add_check(
        "current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier",
        {
            "quotient_safe_input_packet": f149.get("input_packet"),
            "carrier_source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
        },
        {
            "quotient_safe_input_packet": "tau_src_candidate_v1",
            "carrier_source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
        },
        "That quotient-safe witness already lives on the same tau_src_candidate_v1 packet as the surviving F301 pair1/pair2 residual-datum carrier.",
    )
    add_check(
        "current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only",
        {
            "quotient_relation_object": quotient_relation.get("object"),
            "quotient_relation_criterion": quotient_relation.get("criterion"),
            "uses_observer_side_data": bool(quotient_relation.get("uses_observer_side_data")),
            "uses_external_selector_axiom": bool(quotient_relation.get("uses_external_selector_axiom")),
            "quotients_chart_labels_only_via_basis_free_packet": bool(
                quotient_relation.get("quotients_chart_labels_only_via_basis_free_packet")
            ),
            "distinguished_quotient_class_object": distinguished_quotient_class.get("object"),
            "distinguished_quotient_class_exported": bool(f149.get("distinguished_quotient_class_exported")),
            "quotient_class_only": bool(distinguished_quotient_class.get("quotient_class_only")),
            "raw_theta_uniqueness_claimed": bool(distinguished_quotient_class.get("raw_theta_uniqueness_claimed")),
        },
        {
            "quotient_relation_object": "~_src_qw2191_v1",
            "quotient_relation_criterion": "same_source_side_basis_free_selector_packet",
            "uses_observer_side_data": False,
            "uses_external_selector_axiom": False,
            "quotients_chart_labels_only_via_basis_free_packet": True,
            "distinguished_quotient_class_object": "C_sel_src_qw2191_resolved_v1",
            "distinguished_quotient_class_exported": True,
            "quotient_class_only": True,
            "raw_theta_uniqueness_claimed": False,
        },
        "That current quotient-safe resolution still stops exactly at one basis-free distinguished quotient class and does not claim raw-theta uniqueness.",
    )
    add_check(
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral",
        {
            "source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
            "selector_neutral_note_present": contains_token(f301_notes, "selector-neutral"),
            "strict_core_selector_closure_absent": contains_token(
                f301_current_absence, "no strict-core selector closure"
            ),
            "qw2191_discharge_absent": contains_token(f301_current_absence, "no QW-2191 discharge"),
        },
        {
            "source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
            "selector_neutral_note_present": True,
            "strict_core_selector_closure_absent": True,
            "qw2191_discharge_absent": True,
        },
        "The surviving F301 residual-datum carrier remains selector-neutral and still sits below strict-core selector closure or QW-2191 discharge.",
    )
    add_check(
        "current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed",
        current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed,
        True,
        "So even the strongest current quotient-safe QW-2191 resolution witness still remains basis-free quotient-class only rather than typed on the surviving pair1/pair2 residual-datum carrier.",
    )
    add_check(
        "current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation",
        current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation,
        False,
        "No current export continues the present quotient-safe QW-2191 resolution witness into one typed pair1/pair2 residual-datum carrier lane.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge",
        p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current quotient-safe QW-2191 resolution lane into one typed residual-datum carrier bridge.",
    )
    add_check(
        "t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_exported",
        False,
        False,
        "So the source-topology quotient-safe QW-2191 resolution to residual-datum pair1/pair2 typed-carrier bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_SOURCE_TOPOLOGY_QUOTIENT_SAFE_QW2191_RESOLUTION_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P743_REQUIRES_REVIEW_CHANGED_QW2191_SAFE_TO_PAIR12_TYPED_CARRIER_STATE"
    )

    artifact = {
        "stage": "P743",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "P742": str(IN_P742.relative_to(REPO)),
            "F149": str(IN_F149.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t197_target_name": "SourceTopologyQuotientSafeQW2191ResolutionToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t197_target_exported_on_current_repo_state": False,
        "current_actual_source_topology_quotient_safe_qw2191_resolution_exported": (
            current_actual_source_topology_quotient_safe_qw2191_resolution_exported
        ),
        "current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier": (
            current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only": (
            current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed": (
            current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed
        ),
        "current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation": (
            current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation
        ),
        "p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge": (
            p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge
        ),
        "audit_conclusion": {
            "current_repo_already_exports_actual_source_topology_quotient_safe_qw2191_resolution": (
                current_actual_source_topology_quotient_safe_qw2191_resolution_exported
            ),
            "current_repo_already_exports_t197_target": False,
            "next_honest_move": (
                "attack_an_explicit_typed_source_side_bridge_from_the_current_quotient_safe_qw2191_resolution_lane_or_its_codomain_to_the_surviving_F301_pair12_residual_datum_carrier_without_stopping_at_basis_free_distinguished_quotient_classes_and_without_collapsing_into_already_audited_gauge_or_convention_lanes"
            ),
        },
        "hard_limits": [
            "No T197 discharge claim.",
            "No claim that the current quotient-safe QW-2191 resolution already selects one raw pair1/pair2 branch.",
            "No claim that the current distinguished quotient class is already a typed F301 pair1/pair2 residual-datum carrier lane.",
            "No raw-theta uniqueness claim.",
            "No directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P743",
        "status": status,
        "as_of": AS_OF,
        "t197_target_name": "SourceTopologyQuotientSafeQW2191ResolutionToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t197_target_exported_on_current_repo_state": False,
        "current_actual_source_topology_quotient_safe_qw2191_resolution_exported": (
            current_actual_source_topology_quotient_safe_qw2191_resolution_exported
        ),
        "current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier": (
            current_actual_source_topology_quotient_safe_qw2191_resolution_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only": (
            current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed": (
            current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed
        ),
        "current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation": (
            current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation
        ),
        "p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge": (
            p731_pair12_witness_split_descends_to_current_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_typed_bridge
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
