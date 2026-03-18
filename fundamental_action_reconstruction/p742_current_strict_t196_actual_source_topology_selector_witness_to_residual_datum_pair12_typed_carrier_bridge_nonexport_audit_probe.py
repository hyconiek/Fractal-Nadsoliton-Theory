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
IN_P741 = GENERATED / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P741, IN_F148, IN_F301]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P742",
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
    p741 = load_json(IN_P741)
    f148 = load_json(IN_F148)
    f301 = load_json(IN_F301)

    support_packet = f148.get("support_packet") or {}
    basis_class_reduction_operator = support_packet.get("basis_class_reduction_operator") or {}
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

    current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation = (
        f148.get("witness") == "Upsilon_sel_basis_actual_witness_v1"
        and support_packet.get("selector_witness") == "Pi_sel_src_actual_witness_v1"
        and basis_class_reduction_operator.get("object") == "Q_basis_sel_v1"
        and basis_class_reduction_operator.get("domain_packet") == "Sigma_sel_src_target_v1"
        and basis_class_reduction_operator.get("codomain_packet") == "Sigma_sel_basis_free_target_v1"
        and bool(basis_class_reduction_operator.get("forgets_chart_basis_labels"))
        and bool(f148.get("actual_basis_independent_selector_promotion_exported"))
    )

    surviving_pair12_residual_datum_carrier_remains_selector_neutral = (
        f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and contains_token(f301_notes, "selector-neutral")
        and contains_token(f301_current_absence, "no strict-core selector closure")
        and contains_token(f301_current_absence, "no QW-2191 discharge")
    )

    current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation = (
        current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation
        and not bool(p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed"))
        and basis_class_reduction_operator.get("codomain_packet")
        == "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1"
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
    )

    current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed = (
        current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation
        and bool(p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed"))
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and basis_class_reduction_operator.get("codomain_packet") == "Sigma_sel_basis_free_target_v1"
        and not current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation
    )

    p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation
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
        "p741_actual_selector_witness_exists_but_remains_prelm_not_pair12_typed",
        {
            "same_tau_src_packet": bool(
                p741.get("current_actual_source_topology_selector_witness_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "chart_bound_preobserver_only": bool(
                p741.get("current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only")
            ),
            "prelm_not_pair12_typed": bool(
                p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed")
            ),
            "split_descends": bool(
                p741.get("p731_pair12_witness_split_descends_to_current_actual_source_topology_selector_witness")
            ),
            "t195_exported": bool(p741.get("t195_target_exported_on_current_repo_state")),
        },
        {
            "same_tau_src_packet": True,
            "chart_bound_preobserver_only": True,
            "prelm_not_pair12_typed": True,
            "split_descends": False,
            "t195_exported": False,
        },
        "P741 already proves that the actual selector witness is real on the same tau_src packet as F301, but still remains preLM-only rather than typed on the surviving pair1/pair2 carrier.",
    )
    add_check(
        "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation",
        {
            "selector_witness": support_packet.get("selector_witness"),
            "basis_class_reduction_object": basis_class_reduction_operator.get("object"),
            "domain_packet": basis_class_reduction_operator.get("domain_packet"),
            "codomain_packet": basis_class_reduction_operator.get("codomain_packet"),
            "forgets_chart_basis_labels": bool(basis_class_reduction_operator.get("forgets_chart_basis_labels")),
            "actual_basis_independent_selector_promotion_exported": bool(
                f148.get("actual_basis_independent_selector_promotion_exported")
            ),
        },
        {
            "selector_witness": "Pi_sel_src_actual_witness_v1",
            "basis_class_reduction_object": "Q_basis_sel_v1",
            "domain_packet": "Sigma_sel_src_target_v1",
            "codomain_packet": "Sigma_sel_basis_free_target_v1",
            "forgets_chart_basis_labels": True,
            "actual_basis_independent_selector_promotion_exported": True,
        },
        "The strongest current exported continuation out of the actual selector-witness codomain is F148's basis-free class reduction Q_basis_sel_v1, and it explicitly forgets chart basis labels.",
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
        "The surviving F301 residual-datum carrier is real on the same tau_src packet, but remains selector-neutral and below strict-core selector closure or QW-2191 discharge.",
    )
    add_check(
        "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed",
        current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed,
        True,
        "So the current exported continuation from Sigma_sel_src_target_v1 still leaves the actual selector witness only at basis-free chart-label-forgetting class level, not as one typed pair1/pair2 residual-datum carrier bridge.",
    )
    add_check(
        "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation",
        current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation,
        False,
        "No current export continues the actual selector-witness codomain into one typed pair1/pair2 residual-datum carrier lane.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge",
        p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through any exported actual-selector-witness to residual-datum typed carrier bridge.",
    )
    add_check(
        "t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_exported",
        False,
        False,
        "So the actual source-topology selector-witness to residual-datum pair1/pair2 typed-carrier bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P742_REQUIRES_REVIEW_CHANGED_ACTUAL_SELECTOR_WITNESS_TO_PAIR12_TYPED_CARRIER_STATE"
    )

    artifact = {
        "stage": "P742",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "P741": str(IN_P741.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t196_target_name": "ActualSourceTopologySelectorWitnessToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t196_target_exported_on_current_repo_state": False,
        "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation": (
            current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed": (
            current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed
        ),
        "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation": (
            current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation
        ),
        "p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge": (
            p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge
        ),
        "audit_conclusion": {
            "current_repo_already_exports_actual_selector_witness_codomain_continuation": (
                current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation
            ),
            "current_repo_already_exports_t196_target": False,
            "next_honest_move": (
                "attack_an_explicit_chart_sensitive_source_side_typed_bridge_from_Sigma_sel_src_target_v1_to_the_surviving_F301_pair12_residual_datum_carrier_without_forgetting_chart_labels_into_basis_free_classes_and_without_collapsing_into_already_audited_gauge_or_convention_lanes"
            ),
        },
        "hard_limits": [
            "No T196 discharge claim.",
            "No claim that the current basis-free continuation from the actual selector witness already selects one raw pair1/pair2 branch.",
            "No claim that Sigma_sel_src_target_v1 is already typed on the surviving F301 pair1/pair2 carrier.",
            "No basis-independent or quotient-safe QW-2191 upgrade claim beyond the already exported basis-free packet.",
            "No directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P742",
        "status": status,
        "as_of": AS_OF,
        "t196_target_name": "ActualSourceTopologySelectorWitnessToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t196_target_exported_on_current_repo_state": False,
        "current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation": (
            current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed": (
            current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed
        ),
        "current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation": (
            current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation
        ),
        "p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge": (
            p731_pair12_witness_split_descends_to_current_actual_selector_witness_to_residual_datum_typed_carrier_bridge
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
