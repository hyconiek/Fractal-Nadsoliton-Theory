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
IN_F946 = GENERATED / "f946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_packet_summary.json"
IN_P946 = GENERATED / "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_existing_export_or_near_miss_candidate_audit_probe_summary.json"
IN_P765 = GENERATED / "p765_current_strict_t219_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_nonexport_audit_probe_summary.json"
IN_P766 = GENERATED / "p766_current_strict_t220_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_probe_summary.json"
IN_P767 = GENERATED / "p767_current_strict_t221_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_nonexport_audit_probe_summary.json"
IN_P768 = GENERATED / "p768_current_strict_t222_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_target_probe_summary.json"
IN_P769 = GENERATED / "p769_current_strict_t223_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_nonexport_audit_probe_summary.json"
IN_P770 = GENERATED / "p770_current_strict_t224_pair12_source_side_branch_selection_provider_actual_realization_attempt_chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_immediate_missing_subinterface_actual_realization_attempt_probe_summary.json"
IN_T220 = ROOT / "T220_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md"
IN_T222 = ROOT / "T222_CURRENT_STRICT_PAIR12_SOURCE_SIDE_BRANCH_SELECTION_PROVIDER_ACTUAL_REALIZATION_ATTEMPT_CHART_SENSITIVE_PAIR12_TYPED_DESCENT_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_IMMEDIATE_MISSING_SUBINTERFACE_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p947_current_strict_t173_t176_t220_t222_chart_label_retaining_pair12_typed_seed_subinterface_source_side_input_leg_sufficiency_or_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p947_current_strict_t173_t176_t220_t222_chart_label_retaining_pair12_typed_seed_subinterface_source_side_input_leg_sufficiency_or_nonexport_audit_probe_summary.json"

F945_STATUS = (
    "F945_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_TARGET_PACKET_NO_FALSE_PASS"
)
F946_STATUS = (
    "F946_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_OUTPUT_SCHEMA_TARGET_PACKET_NO_FALSE_PASS"
)
P947_PASS = (
    "PASS_T173_T176_T220_T222_CHART_LABEL_RETAINING_PAIR12_TYPED_SEED_SUBINTERFACE_SOURCE_SIDE_INPUT_LEG_SUFFICIENCY_OR_NONEXPORT_AUDITED"
)
EXPECTED_MINIMUM_PROPERTIES = [
    "full_C_v1_scope",
    "pair12_branch_sensitive",
    "chart_sensitive_transported_section_level",
    "nonconvention_nonprojective_nonpremise_smuggled",
]
EXPECTED_RECOMMENDED_NEXT_MOVE = (
    "Build one narrow probe testing whether the frozen T218 future-only support interface, even if exported, would already suffice for inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_v1 or would still require one additional full_C_v1 transported-section lift."
)
EXPECTED_SUBINTERFACE_NAME = (
    "chart_label_retaining_pair12_typed_seed_from_Sigma_sel_src_target_v1_"
    "toward_the_surviving_F301_pair12_carrier_prior_to_Q_basis_sel_v1_terminal_"
    "collapse_and_prior_to_projector_only_local_pair12_atlas_collapse"
)
EXPECTED_T220_ATTEMPT = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_"
    "chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_v1"
)
EXPECTED_T222_TARGET = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_"
    "chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_"
    "immediate_missing_subinterface_target_v1"
)
EXPECTED_T224_ATTEMPT = (
    "W_strict_t173_pair12_source_side_branch_selection_provider_actual_realization_attempt_"
    "chart_sensitive_pair12_typed_descent_interface_actual_realization_attempt_"
    "immediate_missing_subinterface_actual_realization_attempt_v1"
)
BRIDGE_TARGET_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1"
)
BRIDGE_OUTPUT_SCHEMA_TARGET_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_named_source_side_input_leg_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P947_CURRENT_STRICT_T173_T176_T220_T222_CHART_LABEL_RETAINING_PAIR12_TYPED_SEED_SUBINTERFACE_SOURCE_SIDE_INPUT_LEG_SUFFICIENCY_OR_NONEXPORT_AUDIT_PROBE.md",
        "F947_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_SOURCE_SIDE_INPUT_LEG_TARGET_PACKET.md",
        "f947_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_packet.py",
        "P946_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_EXISTING_EXPORT_OR_NEAR_MISS_CANDIDATE_AUDIT_PROBE.md",
        "p946_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_existing_export_or_near_miss_candidate_audit_probe.py",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if "source_side_input_leg" in text and (
                BRIDGE_TARGET_ID in text
                or BRIDGE_OUTPUT_SCHEMA_TARGET_ID in text
                or EXPECTED_SUBINTERFACE_NAME in text
            ):
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_F945,
        IN_F946,
        IN_P946,
        IN_P765,
        IN_P766,
        IN_P767,
        IN_P768,
        IN_P769,
        IN_P770,
        IN_T220,
        IN_T222,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P947",
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
    f946 = load_json(IN_F946)
    p946 = load_json(IN_P946)
    p765 = load_json(IN_P765)
    p766 = load_json(IN_P766)
    p767 = load_json(IN_P767)
    p768 = load_json(IN_P768)
    p769 = load_json(IN_P769)
    p770 = load_json(IN_P770)
    t220_text = load_text(IN_T220)
    t222_text = load_text(IN_T222)
    named_source_side_input_leg_candidates = scan_named_source_side_input_leg_candidates()

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

    f945_f946_full_c_v1_target_and_probe_recommendation_are_frozen = (
        f945.get("status") == F945_STATUS
        and ((f945.get("target_object") or {}).get("object_id") == BRIDGE_TARGET_ID)
        and ((f945.get("bridge_output_schema") or {}).get("minimum_properties") == EXPECTED_MINIMUM_PROPERTIES)
        and f946.get("status") == F946_STATUS
        and f946.get("target_object_id") == BRIDGE_OUTPUT_SCHEMA_TARGET_ID
        and f946.get("recommended_next_move") == EXPECTED_RECOMMENDED_NEXT_MOVE
    )

    p946_already_localizes_the_question_at_source_side_input_leg_level = (
        p946.get("status")
        == "PASS_T173_T176_BRIDGE_EXISTING_EXPORT_OR_NEAR_MISS_CANDIDATES_AUDITED_NO_LAWFUL_SUPPLIER_FOUND"
        and p946.get("bridge_target_object_id") == BRIDGE_TARGET_ID
        and p946.get("no_current_lawful_bridge_supplier_found") is True
        and "source_side_input_leg" in str(p946.get("narrowest_honest_next_probe_question") or "")
        and "T220/T222" in str(p946.get("narrowest_honest_next_probe_question") or "")
    )

    current_t218_interface_is_still_not_actually_exported = (
        p765.get("current_repo_still_does_not_export_actual_realization_of_t218_target") is True
        and p765.get("current_actual_selector_witness_codomain_still_lacks_actual_chart_sensitive_pair12_typed_descent_interface")
        is True
    )

    current_t220_attempt_is_exported_but_still_open = (
        p766.get("t220_attempt_exported_on_current_repo_state") is True
        and p766.get("t220_attempt_name") == EXPECTED_T220_ATTEMPT
        and p766.get("first_actual_t218_interface_realization_attempt_keeps_success_failure_open") is True
        and ((p766.get("first_actual_t218_interface_realization_attempt") or {}).get("immediate_missing_subinterface"))
        == EXPECTED_SUBINTERFACE_NAME
    )

    current_t220_missing_subinterface_remains_unexported = (
        p767.get("current_repo_still_does_not_export_actual_realization_of_t220_immediate_missing_subinterface")
        is True
        and p767.get("current_t220_attempt_stalls_exactly_at_the_named_missing_subinterface") is True
        and p767.get("exact_named_missing_subinterface") == EXPECTED_SUBINTERFACE_NAME
    )

    current_t222_subinterface_target_is_future_only = (
        p768.get("t222_target_exported_on_current_repo_state") is True
        and p768.get("t222_target_name") == EXPECTED_T222_TARGET
        and p768.get("current_t222_target_is_future_route_only") is True
        and p768.get("current_t222_target_freezes_exact_t220_immediate_missing_subinterface") is True
        and p768.get("current_t222_target_remains_below_actual_subinterface_export_interface_export_and_t176_discharge")
        is True
    )

    current_actual_t222_subinterface_realization_is_still_unexported = (
        p769.get("current_repo_still_does_not_export_actual_realization_of_t222_target") is True
        and p769.get("current_actual_selector_witness_codomain_still_lacks_actual_chart_label_retaining_pair12_typed_seed_subinterface")
        is True
        and p769.get("current_exact_t220_missing_subinterface_still_only_future_target_not_actual_export")
        is True
    )

    current_t224_attempt_is_exported_but_still_open = (
        p770.get("t224_attempt_exported_on_current_repo_state") is True
        and p770.get("t224_attempt_name") == EXPECTED_T224_ATTEMPT
        and p770.get("first_actual_t222_subinterface_realization_attempt_keeps_success_failure_open") is True
        and ((p770.get("first_actual_t222_subinterface_realization_attempt") or {}).get("targeted_subinterface"))
        == EXPECTED_SUBINTERFACE_NAME
    )

    t220_t222_route_scope_is_source_side_local_and_precollapse = all(
        needle in t220_text
        for needle in [
            EXPECTED_T220_ATTEMPT,
            EXPECTED_SUBINTERFACE_NAME,
            "Sigma_sel_src_target_v1",
            "attempt_must_preserve_chart_labels_prior_to_Q_basis_sel_v1_terminal_collapse := yes",
            "attempt_must_not_collapse_to_projector_only_local_pair12_atlas := yes",
        ]
    ) and all(
        needle in t222_text
        for needle in [
            EXPECTED_T222_TARGET,
            EXPECTED_SUBINTERFACE_NAME,
            "target_starts_at_current_actual_selector_witness_codomain := yes",
            "target_points_toward_surviving_pair12_residual_datum_carrier_lane := yes",
            "target_precedes_Q_basis_sel_v1_terminal_collapse := yes",
            "target_precedes_projector_only_local_pair12_atlas_collapse := yes",
            "target_remains_below_global_t176_discharge := yes",
            "future_route_only := yes",
        ]
    )

    current_repo_already_exports_actual_source_side_input_leg_supplier = (
        current_t218_interface_is_still_not_actually_exported is False
        or current_t220_missing_subinterface_remains_unexported is False
        or current_actual_t222_subinterface_realization_is_still_unexported is False
        or len(named_source_side_input_leg_candidates) > 0
    )

    route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift = (
        ((f945.get("bridge_output_schema") or {}).get("minimum_properties") == EXPECTED_MINIMUM_PROPERTIES)
        and current_t222_subinterface_target_is_future_only
        and t220_t222_route_scope_is_source_side_local_and_precollapse
    )

    add_check(
        "f945_f946_full_c_v1_target_and_probe_recommendation_are_frozen",
        f945_f946_full_c_v1_target_and_probe_recommendation_are_frozen,
        True,
        "F945 and F946 already freeze the full-C_v1 bridge target plus the narrower output-schema target and explicitly call for this exact sufficiency-or-nonexport probe.",
    )
    add_check(
        "p946_already_localizes_the_question_at_source_side_input_leg_level",
        p946_already_localizes_the_question_at_source_side_input_leg_level,
        True,
        "The existing P946 near-miss audit already localizes the next honest question exactly at whether the T220/T222 route can lawfully supply the source_side_input_leg.",
    )
    add_check(
        "current_t218_interface_is_still_not_actually_exported",
        current_t218_interface_is_still_not_actually_exported,
        True,
        "P765 already freezes that the stronger chart-sensitive pair12 typed descent interface itself is still not actually exported.",
    )
    add_check(
        "current_t220_attempt_is_exported_but_still_open",
        current_t220_attempt_is_exported_but_still_open,
        True,
        "P766 exports only the first actual T220 attempt and keeps success/failure open.",
    )
    add_check(
        "current_t220_missing_subinterface_remains_unexported",
        current_t220_missing_subinterface_remains_unexported,
        True,
        "P767 freezes that the exact T220 seed/subinterface remains unexported.",
    )
    add_check(
        "current_t222_subinterface_target_is_future_only",
        current_t222_subinterface_target_is_future_only,
        True,
        "P768 exports only a future-only T222 target below actual subinterface export, below actual interface export, and below T176.",
    )
    add_check(
        "current_actual_t222_subinterface_realization_is_still_unexported",
        current_actual_t222_subinterface_realization_is_still_unexported,
        True,
        "P769 freezes that the actual realization of the T222 target is still unexported, so the route still lacks an actual chart-label-retaining pair12 typed seed/subinterface.",
    )
    add_check(
        "current_t224_attempt_is_exported_but_still_open",
        current_t224_attempt_is_exported_but_still_open,
        True,
        "P770 exports only the first actual T222 subinterface realization attempt and still keeps success/failure open.",
    )
    add_check(
        "t220_t222_route_scope_is_source_side_local_and_precollapse",
        t220_t222_route_scope_is_source_side_local_and_precollapse,
        True,
        "T220 and T222 still describe only one route-local, source-side, precollapse seed/subinterface lane from Sigma_sel_src_target_v1 toward the surviving F301 pair12 carrier.",
    )
    add_check(
        "current_repo_already_exports_actual_source_side_input_leg_supplier",
        current_repo_already_exports_actual_source_side_input_leg_supplier,
        False,
        "No current export lawfully supplies the actual source_side_input_leg for the frozen bridge family.",
    )
    add_check(
        "route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift",
        route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift,
        True,
        "Even a future actual export of the route-local seed/subinterface would still remain only a source-side leg and would not by itself equal the full-C_v1 chart-sensitive transported-section bridge output schema.",
    )

    status = P947_PASS if not blocking else "P947_REQUIRES_REVIEW_CHANGED_SOURCE_SIDE_INPUT_LEG_STATE"

    artifact = {
        "stage": "P947",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t173_t176_t220_t222_source_side_input_leg_sufficiency_or_nonexport_boundary_only",
        "inputs": {
            "F945": rel(IN_F945),
            "F946": rel(IN_F946),
            "P946": rel(IN_P946),
            "P765": rel(IN_P765),
            "P766": rel(IN_P766),
            "P767": rel(IN_P767),
            "P768": rel(IN_P768),
            "P769": rel(IN_P769),
            "P770": rel(IN_P770),
            "T220": rel(IN_T220),
            "T222": rel(IN_T222),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "named_source_side_input_leg_candidates_found_elsewhere": named_source_side_input_leg_candidates,
        "audit_conclusion": {
            "bridge_target_object_id": BRIDGE_TARGET_ID,
            "bridge_output_schema_target_id": BRIDGE_OUTPUT_SCHEMA_TARGET_ID,
            "current_repo_already_exports_actual_source_side_input_leg_supplier": current_repo_already_exports_actual_source_side_input_leg_supplier,
            "current_best_route_local_attempt_name": EXPECTED_T220_ATTEMPT,
            "current_best_future_only_subinterface_target_name": EXPECTED_T222_TARGET,
            "current_best_open_actual_subinterface_attempt_name": EXPECTED_T224_ATTEMPT,
            "exact_route_local_seed_subinterface_name": EXPECTED_SUBINTERFACE_NAME,
            "route_local_seed_subinterface_scope": {
                "starts_at": "Sigma_sel_src_target_v1",
                "points_toward": "surviving_F301_pair12_residual_datum_carrier_lane",
                "precedes_q_basis_terminal_collapse": True,
                "precedes_projector_only_local_atlas_collapse": True,
                "remains_below_actual_interface_export": True,
                "remains_below_t176_discharge": True,
            },
            "route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift": (
                route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift
            ),
            "next_honest_move": "Freeze the exact missing source_side_input_leg target for the F946 bridge-output-schema family without claiming that the leg already exists or that the remaining full_C_v1 transported-section lift is solved.",
        },
        "current_honest_reading": [
            "The current T220/T222 family exports only route-local attempts and future-only targets around the chart-label-retaining pair12 typed seed/subinterface.",
            "No current artifact yet exports one actual source_side_input_leg that could lawfully feed the frozen F946 bridge-output-schema target.",
            "Even a future actual export of that route-local source-side leg would still remain below the full-C_v1 chart-sensitive transported-section completion required by F945.",
        ],
        "hard_limits": [
            "Does not claim that an actual source_side_input_leg is already exported.",
            "Does not claim success of the frozen T220 attempt.",
            "Does not claim success of the frozen T222 actual-realization attempt.",
            "Does not claim that the route-local seed/subinterface already suffices for full C_v1.",
            "Does not claim that the additional transported-section lift already exists.",
            "Does not claim that T176 is discharged.",
            "Does not claim that T177 is discharged.",
            "Does not claim that T185 is discharged.",
            "Does not claim that F647 is admissible S_sel_int.",
            "Does not claim that the rooted w_break sign lift is a strict physical orientation datum.",
            "Does not claim kernel-alone/global QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P947",
        "status": status,
        "as_of": AS_OF,
        "bridge_target_object_id": BRIDGE_TARGET_ID,
        "bridge_output_schema_target_id": BRIDGE_OUTPUT_SCHEMA_TARGET_ID,
        "no_current_source_side_input_leg_supplier_found": not current_repo_already_exports_actual_source_side_input_leg_supplier,
        "additional_full_c_v1_transported_section_lift_still_required": (
            route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift
        ),
        "recommended_next_move": artifact["audit_conclusion"]["next_honest_move"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
