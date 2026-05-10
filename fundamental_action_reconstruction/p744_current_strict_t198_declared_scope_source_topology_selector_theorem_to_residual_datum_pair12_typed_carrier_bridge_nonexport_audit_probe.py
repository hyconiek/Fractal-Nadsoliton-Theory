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
IN_P743 = GENERATED / "p743_current_strict_t197_source_topology_quotient_safe_qw2191_resolution_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p744_current_strict_t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p744_current_strict_t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P743, IN_F150, IN_F301]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P744",
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
    p743 = load_json(IN_P743)
    f150 = load_json(IN_F150)
    f301 = load_json(IN_F301)

    support_packet = f150.get("support_packet") or {}
    l5_boundaries = support_packet.get("l5_boundaries") or {}
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

    current_declared_scope_source_topology_selector_theorem_exported = (
        f150.get("witness") == "T14_src_selector_declared_scope_actual_witness_v1"
        and f150.get("input_packet") == "tau_src_candidate_v1"
        and f150.get("codomain_target") == "declared_scope_source_topology_selector_theorem_target_v1"
        and support_packet.get("t14_theorem_spec") == "T14_SourceTopologySelector_Theorem"
        and bool(f150.get("declared_scope_source_topology_selector_theorem_exported"))
    )

    current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier = (
        current_declared_scope_source_topology_selector_theorem_exported
        and f150.get("input_packet") == "tau_src_candidate_v1"
        and f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
    )

    current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only = (
        current_declared_scope_source_topology_selector_theorem_exported
        and bool(f150.get("qw2191_quotient_safe_discharged"))
        and bool(f150.get("declared_scope_only"))
        and not bool(f150.get("raw_theta_uniqueness_claimed"))
        and not bool(f150.get("tau_src_identified_with_s_prelm"))
        and not bool(f150.get("admissible_strict_core_internal_selector_source_object_claimed"))
        and bool(l5_boundaries.get("observer_downstream_only"))
        and bool(p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only"))
    )

    surviving_pair12_residual_datum_carrier_remains_selector_neutral = (
        f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and contains_token(f301_notes, "selector-neutral")
        and contains_token(f301_current_absence, "no strict-core selector closure")
        and contains_token(f301_current_absence, "no QW-2191 discharge")
    )

    current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation = (
        current_declared_scope_source_topology_selector_theorem_exported
        and f150.get("codomain_target") == "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1"
        and not bool(f150.get("declared_scope_only"))
        and bool(f150.get("raw_theta_uniqueness_claimed"))
        and current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
    )

    current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed = (
        current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only
        and current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier
        and bool(
            p743.get("current_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only_not_pair12_typed")
        )
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and not current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation
    )

    p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation
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
        "p743_quotient_safe_qw2191_resolution_already_real_but_not_pair12_typed",
        {
            "resolution_exported": bool(
                p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_exported")
            ),
            "remains_quotient_class_only": bool(
                p743.get("current_actual_source_topology_quotient_safe_qw2191_resolution_remains_quotient_class_only")
            ),
            "pair12_typed_residual_datum_continuation_exported": bool(
                p743.get("current_source_topology_quotient_safe_qw2191_resolution_has_exported_pair12_typed_residual_datum_continuation")
            ),
            "t197_exported": bool(p743.get("t197_target_exported_on_current_repo_state")),
        },
        {
            "resolution_exported": True,
            "remains_quotient_class_only": True,
            "pair12_typed_residual_datum_continuation_exported": False,
            "t197_exported": False,
        },
        "P743 already proves that the strongest current quotient-safe QW-2191 resolution lane is real, but still stops at quotient-class level rather than one typed pair1/pair2 residual-datum carrier bridge.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_exported",
        {
            "witness": f150.get("witness"),
            "input_packet": f150.get("input_packet"),
            "codomain_target": f150.get("codomain_target"),
            "t14_theorem_spec": support_packet.get("t14_theorem_spec"),
            "declared_scope_source_topology_selector_theorem_exported": bool(
                f150.get("declared_scope_source_topology_selector_theorem_exported")
            ),
        },
        {
            "witness": "T14_src_selector_declared_scope_actual_witness_v1",
            "input_packet": "tau_src_candidate_v1",
            "codomain_target": "declared_scope_source_topology_selector_theorem_target_v1",
            "t14_theorem_spec": "T14_SourceTopologySelector_Theorem",
            "declared_scope_source_topology_selector_theorem_exported": True,
        },
        "F150 already exports one real declared-scope source-topology selector theorem witness on tau_src_candidate_v1.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier",
        {
            "theorem_input_packet": f150.get("input_packet"),
            "carrier_source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
        },
        {
            "theorem_input_packet": "tau_src_candidate_v1",
            "carrier_source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
        },
        "That declared-scope theorem already lives on the same tau_src_candidate_v1 packet as the surviving F301 pair1/pair2 residual-datum carrier.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only",
        {
            "qw2191_quotient_safe_discharged": bool(f150.get("qw2191_quotient_safe_discharged")),
            "declared_scope_only": bool(f150.get("declared_scope_only")),
            "raw_theta_uniqueness_claimed": bool(f150.get("raw_theta_uniqueness_claimed")),
            "tau_src_identified_with_s_prelm": bool(f150.get("tau_src_identified_with_s_prelm")),
            "admissible_strict_core_internal_selector_source_object_claimed": bool(
                f150.get("admissible_strict_core_internal_selector_source_object_claimed")
            ),
            "observer_downstream_only": bool(l5_boundaries.get("observer_downstream_only")),
        },
        {
            "qw2191_quotient_safe_discharged": True,
            "declared_scope_only": True,
            "raw_theta_uniqueness_claimed": False,
            "tau_src_identified_with_s_prelm": False,
            "admissible_strict_core_internal_selector_source_object_claimed": False,
            "observer_downstream_only": True,
        },
        "That theorem-level packaging still remains declared-scope only, quotient-safe only, and does not claim raw-theta uniqueness or a strict-core selector source object.",
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
        "current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed",
        current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed,
        True,
        "So even the strongest current theorem-level source-topology selector lane still remains declared-scope quotient-class only rather than typed on the surviving pair1/pair2 residual-datum carrier.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation",
        current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation,
        False,
        "No current export continues the present declared-scope source-topology selector theorem into one typed pair1/pair2 residual-datum carrier lane.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge",
        p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through the current declared-scope theorem lane into one typed residual-datum carrier bridge.",
    )
    add_check(
        "t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_exported",
        False,
        False,
        "So the declared-scope source-topology selector theorem to residual-datum pair1/pair2 typed-carrier bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P744_REQUIRES_REVIEW_CHANGED_DECLARED_SCOPE_THEOREM_TO_PAIR12_TYPED_CARRIER_STATE"
    )

    artifact = {
        "stage": "P744",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "P743": str(IN_P743.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t198_target_name": "DeclaredScopeSourceTopologySelectorTheoremToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t198_target_exported_on_current_repo_state": False,
        "current_declared_scope_source_topology_selector_theorem_exported": (
            current_declared_scope_source_topology_selector_theorem_exported
        ),
        "current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier": (
            current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only": (
            current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed": (
            current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed
        ),
        "current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation": (
            current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation
        ),
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge": (
            p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge
        ),
        "audit_conclusion": {
            "current_repo_already_exports_declared_scope_source_topology_selector_theorem": (
                current_declared_scope_source_topology_selector_theorem_exported
            ),
            "current_repo_already_exports_t198_target": False,
            "next_honest_move": (
                "attack_an_explicit_typed_source_side_bridge_from_the_current_declared_scope_source_topology_selector_theorem_lane_or_its_codomain_to_the_surviving_F301_pair12_residual_datum_carrier_without_stopping_at_declared_scope_basis_free_quotient_classes_and_without_collapsing_into_already_audited_gauge_or_convention_lanes"
            ),
        },
        "hard_limits": [
            "No T198 discharge claim.",
            "No claim that the current declared-scope source-topology selector theorem already selects one raw pair1/pair2 branch.",
            "No claim that the current declared-scope theorem codomain is already a typed F301 pair1/pair2 residual-datum carrier lane.",
            "No raw-theta uniqueness claim.",
            "No directed/sign-sensitive physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P744",
        "status": status,
        "as_of": AS_OF,
        "t198_target_name": "DeclaredScopeSourceTopologySelectorTheoremToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t198_target_exported_on_current_repo_state": False,
        "current_declared_scope_source_topology_selector_theorem_exported": (
            current_declared_scope_source_topology_selector_theorem_exported
        ),
        "current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier": (
            current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only": (
            current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed": (
            current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed
        ),
        "current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation": (
            current_declared_scope_source_topology_selector_theorem_has_exported_pair12_typed_residual_datum_continuation
        ),
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge": (
            p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_to_residual_datum_typed_bridge
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
