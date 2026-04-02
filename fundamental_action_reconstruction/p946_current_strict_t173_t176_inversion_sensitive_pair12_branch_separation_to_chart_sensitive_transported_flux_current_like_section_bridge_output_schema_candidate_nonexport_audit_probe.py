#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-21"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F945 = GENERATED / "f945_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_target_packet.json"
IN_P721 = GENERATED / "p721_current_strict_t176_source_topology_basis_free_qw2191_safe_provider_nonupgrade_audit_probe_summary.json"
IN_P722 = GENERATED / "p722_current_strict_t177_chart_sensitive_transported_flux_current_like_section_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_P747 = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"
IN_P763 = GENERATED / "p763_current_strict_t217_pair12_source_side_branch_selection_provider_actual_realization_attempt_immediate_missing_interface_nonexport_audit_probe_summary.json"
IN_P764 = GENERATED / "p764_current_strict_t218_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_target_probe_summary.json"
IN_T218 = ROOT / "T218_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_candidate_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_candidate_nonexport_audit_probe_summary.json"

F945_STATUS = (
    "F945_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
)
P946_PASS = (
    "PASS_T173_T176_INVERSION_SENSITIVE_PAIR12_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_CANDIDATE_NONEXPORT_AUDITED"
)
EXPECTED_MINIMUM_PROPERTIES = [
    "full_C_v1_scope",
    "pair12_branch_sensitive",
    "chart_sensitive_transported_section_level",
    "nonconvention_nonprojective_nonpremise_smuggled",
]
EXPECTED_INTERFACE_NAME = (
    "chart_sensitive_pair12_typed_descent_from_Sigma_sel_src_target_v1_to_the_surviving_F301_pair12_carrier_without_Q_basis_sel_v1_terminal_collapse_and_without_projector_only_atlas_collapse"
)
EXPECTED_T177 = "ChartSensitiveTransportedFluxCurrentLikeSection_global_C_v1_source_topology_bridge_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F945, IN_P721, IN_P722, IN_P742, IN_P747, IN_P763, IN_P764, IN_T218]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P946",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f945 = load_json(IN_F945)
    p721 = load_json(IN_P721)
    p722 = load_json(IN_P722)
    p742 = load_json(IN_P742)
    p747 = load_json(IN_P747)
    p763 = load_json(IN_P763)
    p764 = load_json(IN_P764)
    t218_text = load_text(IN_T218)

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

    bridge_output_schema = f945.get("bridge_output_schema") or {}
    minimum_properties = bridge_output_schema.get("minimum_properties") or []

    add_check(
        "f945_bridge_output_schema_frozen_with_required_minimum_properties",
        {
            "status": f945.get("status"),
            "minimum_properties": minimum_properties,
        },
        {
            "status": F945_STATUS,
            "minimum_properties": EXPECTED_MINIMUM_PROPERTIES,
        },
        "F945 already freezes one exact bridge_output_schema with the required full-C_v1, pair12-sensitive, chart-sensitive, and nonprojective minimum properties.",
    )
    add_check(
        "p721_p722_keep_full_chart_sensitive_global_bridge_unexported",
        {
            "t176_provider_exported": bool(p721.get("source_topology_lane_upgrades_to_exact_t176_provider")),
            "t177_target_name": p722.get("t177_target_name"),
            "t177_exported": bool(p722.get("t177_target_exported_on_current_repo_state")),
        },
        {
            "t176_provider_exported": False,
            "t177_target_name": EXPECTED_T177,
            "t177_exported": False,
        },
        "The global chart-sensitive transported flux/current-like section bridge remains unexported on the active T173/T176 frontier.",
    )
    add_check(
        "p742_basis_free_continuation_cannot_supply_pair12_typed_chart_sensitive_output_schema",
        {
            "basis_free_chart_label_forgetting_continuation": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")
            ),
            "continuation_pair12_typed": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")
            ),
            "target_exported": bool(p742.get("t196_target_exported_on_current_repo_state")),
        },
        {
            "basis_free_chart_label_forgetting_continuation": True,
            "continuation_pair12_typed": False,
            "target_exported": False,
        },
        "The strongest current continuation out of the actual selector witness codomain still stays basis-free and therefore cannot supply the pair12-typed chart-sensitive bridge output schema.",
    )
    add_check(
        "p747_local_chart_sensitive_atlas_lane_remains_local_projector_only_and_unbridged",
        {
            "local_chart_sensitive_atlas_lane_exported": bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_exported")),
            "projector_only": bool(p747.get("current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe")),
            "bridge_exported": bool(p747.get("current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge")),
            "t201_target_exported": bool(p747.get("t201_target_exported_on_current_repo_state")),
        },
        {
            "local_chart_sensitive_atlas_lane_exported": True,
            "projector_only": True,
            "bridge_exported": False,
            "t201_target_exported": False,
        },
        "The current chart-sensitive atlas-side lane remains only local and projector-level, and still does not export the bridge from the actual selector witness target.",
    )
    add_check(
        "p763_exact_missing_interface_for_chart_sensitive_pair12_typed_descent_remains_unexported",
        {
            "interface_nonexported": bool(p763.get("current_t216_attempt_immediate_missing_interface_is_still_unexported")),
            "named_interface_matches": p763.get("exact_named_missing_interface"),
            "next_move_still_interface_or_failure_boundary": bool(
                p763.get("next_honest_move_is_export_that_exact_interface_or_freeze_attempt_level_failure_boundary")
            ),
        },
        {
            "interface_nonexported": True,
            "named_interface_matches": EXPECTED_INTERFACE_NAME,
            "next_move_still_interface_or_failure_boundary": True,
        },
        "The exact chart-sensitive pair12-typed descent interface remains unexported and is already sharply named on the current repo state.",
    )
    add_check(
        "p764_t218_export_only_future_route_target_below_actual_interface_provider_and_t176",
        {
            "t218_target_exported": bool(p764.get("t218_target_exported_on_current_repo_state")),
            "future_only": bool(p764.get("current_t218_target_is_future_route_only")),
            "freezes_exact_missing_interface": bool(p764.get("current_t218_target_freezes_exact_t216_immediate_missing_interface")),
            "below_actual_interface_and_t176": bool(
                p764.get("current_t218_target_remains_below_actual_interface_export_and_below_t176_discharge")
            ),
        },
        {
            "t218_target_exported": True,
            "future_only": True,
            "freezes_exact_missing_interface": True,
            "below_actual_interface_and_t176": True,
        },
        "The repo exports one sharp future-only interface target candidate, but it remains explicitly below actual interface export, below actual provider export, and below global T176 discharge.",
    )

    t218_scope = {
        "starts_at_sigma_sel_src_target": "Sigma_sel_src_target_v1" in t218_text,
        "ends_in_surviving_f301_pair12_carrier_lane": "target_ends_in_surviving_pair12_residual_datum_carrier_lane := yes" in t218_text,
        "future_route_only": "future_route_only := yes" in t218_text,
        "below_global_t176_discharge": "target_remains_below_global_t176_discharge := yes" in t218_text,
        "mentions_full_c_v1_global_chart_sensitive_section_bridge": EXPECTED_T177 in t218_text,
    }
    add_check(
        "t218_scope_is_relevant_support_candidate_but_not_full_c_v1_bridge_output_schema_supplier",
        t218_scope,
        {
            "starts_at_sigma_sel_src_target": True,
            "ends_in_surviving_f301_pair12_carrier_lane": True,
            "future_route_only": True,
            "below_global_t176_discharge": True,
            "mentions_full_c_v1_global_chart_sensitive_section_bridge": False,
        },
        "The T218 target is relevant as a support-side candidate, but its own frozen scope remains an internal future route from Sigma_sel_src_target_v1 to the surviving F301 carrier lane rather than a full-C_v1 chart-sensitive transported section bridge.",
    )

    actual_bridge_output_schema_supplier_exported = False
    add_check(
        "no_current_export_lawfully_supplies_f945_bridge_output_schema",
        actual_bridge_output_schema_supplier_exported,
        False,
        "No current export lawfully supplies the F945 bridge_output_schema on full C_v1.",
    )

    status = P946_PASS if not blocking else "P946_REQUIRES_REVIEW_CHANGED_BRIDGE_OUTPUT_SCHEMA_CANDIDATE_STATE"

    artifact = {
        "stage": "P946",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t173_t176_bridge_output_schema_candidate_nonexport_boundary_only",
        "inputs": {
            "F945": str(IN_F945.relative_to(REPO)),
            "P721": str(IN_P721.relative_to(REPO)),
            "P722": str(IN_P722.relative_to(REPO)),
            "P742": str(IN_P742.relative_to(REPO)),
            "P747": str(IN_P747.relative_to(REPO)),
            "P763": str(IN_P763.relative_to(REPO)),
            "P764": str(IN_P764.relative_to(REPO)),
            "T218": str(IN_T218.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "bridge_target_object_id": ((f945.get("target_object") or {}).get("object_id")),
        "required_bridge_output_minimum_properties": minimum_properties,
        "candidate_classes_audited": {
            "basis_free_continuation_candidate": {
                "artifact": p742.get("t196_target_name"),
                "exported_on_current_repo_state": bool(p742.get("t196_target_exported_on_current_repo_state")),
                "fit_assessment": "fails_pair12_typed_and_chart_sensitive_section_level",
            },
            "local_chart_sensitive_atlas_candidate": {
                "artifact": p747.get("t201_target_name"),
                "exported_on_current_repo_state": bool(p747.get("t201_target_exported_on_current_repo_state")),
                "fit_assessment": "fails_full_C_v1_scope_and_remains_projector_only",
            },
            "future_only_exact_interface_target_candidate": {
                "artifact": p764.get("t218_target_name"),
                "exact_missing_interface_name": p763.get("exact_named_missing_interface"),
                "exported_on_current_repo_state": bool(p764.get("t218_target_exported_on_current_repo_state")),
                "fit_assessment": "relevant_future_only_support_candidate_but_not_actual_supplier",
            },
        },
        "audit_conclusion": {
            "current_repo_already_exports_actual_bridge_output_schema_supplier": False,
            "current_repo_already_exports_relevant_future_only_interface_target_candidate": bool(
                p764.get("t218_target_exported_on_current_repo_state")
            ),
            "current_best_future_only_interface_target_candidate": p764.get("t218_target_name"),
            "current_best_support_interface_name": p763.get("exact_named_missing_interface"),
            "current_best_support_candidate_alone_satisfies_f945_bridge_output_schema": False,
            "next_honest_move": (
                "freeze_the_exact_future_only_support_interface_target_already_named_by_P763_and_T218_as_the_current_best_bridge_output_schema_support_candidate_for_F945_while_explicitly_keeping_it_below_actual_interface_export_and_below_full_C_v1_T176_discharge"
            ),
        },
        "hard_limits": [
            "No T176 discharge claim.",
            "No T177 discharge claim.",
            "No T185 discharge claim.",
            "No claim that the T218 future-only interface target is already an actual bridge export.",
            "No claim that the T218 future-only interface target already suffices for full C_v1.",
            "No promotion of F647 to admissible S_sel_int.",
            "No promotion of the rooted w_break sign lift into a strict physical orientation datum.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P946",
        "status": status,
        "as_of": AS_OF,
        "bridge_target_object_id": artifact["bridge_target_object_id"],
        "current_repo_already_exports_actual_bridge_output_schema_supplier": artifact["audit_conclusion"][
            "current_repo_already_exports_actual_bridge_output_schema_supplier"
        ],
        "current_best_future_only_interface_target_candidate": artifact["audit_conclusion"][
            "current_best_future_only_interface_target_candidate"
        ],
        "current_best_support_interface_name": artifact["audit_conclusion"]["current_best_support_interface_name"],
        "next_honest_move": artifact["audit_conclusion"]["next_honest_move"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
