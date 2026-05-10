#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1091 = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe_summary.json"
IN_F966 = GENERATED / "f966_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_packet_summary.json"
IN_P1114 = GENERATED / "p1114_current_strict_t173_t176_post_f970_onrd_to_f966_rejoin_audit_probe_summary.json"
IN_F971 = GENERATED / "f971_current_strict_t173_t176_post_f970_onrd_to_f966_route_packet_summary.json"
IN_P1011 = GENERATED / "p1011_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admission_probe_summary.json"
IN_P1034 = GENERATED / "p1034_current_strict_qw2191_nadsoliton_neural_character_information_primary_selector_support_reference_admission_probe_summary.json"
IN_P1047 = GENERATED / "p1047_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_audit_probe_summary.json"
IN_F956 = GENERATED / "f956_current_strict_qw2191_post_stop_neural_bridge_lane_to_strict_int_selector_source_frontier_route_decision_packet_summary.json"
IN_P1048 = GENERATED / "p1048_current_strict_qw2191_post_stop_route_rejoin_to_existing_t173_frontier_audit_probe_summary.json"

OUT_JSON = GENERATED / "p1115_current_strict_t173_t176_f966_selector_source_subroute_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1115_current_strict_t173_t176_f966_selector_source_subroute_audit_probe_summary.json"

ACTIVE_MISSING_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
PRIMARY_SUBROUTE = "explicit_strict_internal_selector_source_derivation_frontier"
PRIMARY_TARGET = "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P1091,
        IN_F966,
        IN_P1114,
        IN_F971,
        IN_P1011,
        IN_P1034,
        IN_P1047,
        IN_F956,
        IN_P1048,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1115",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1091 = load_json(IN_P1091)
    f966 = load_json(IN_F966)
    p1114 = load_json(IN_P1114)
    f971 = load_json(IN_F971)
    p1011 = load_json(IN_P1011)
    p1034 = load_json(IN_P1034)
    p1047 = load_json(IN_P1047)
    f956 = load_json(IN_F956)
    p1048 = load_json(IN_P1048)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append({
            "id": check_id,
            "actual": actual,
            "expected": expected,
            "pass": passed,
            "meaning": meaning,
        })
        if not passed:
            blocking.append(check_id)

    existing_f966_live_contract_confirmed = (
        p1091.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_DECISION_AUDITED"
        and p1091.get("preferred_search_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
        and p1091.get("active_missing_bridge") == ACTIVE_MISSING_BRIDGE
        and f966.get("status") == "F966_EXECUTED_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_PACKET_NO_FALSE_PASS"
        and f966.get("route_object_id") == "PostF965FailureMapConstrainedToExportedNoncyclicProviderSplitNonSameLaneUpgradeRoute_v1"
        and f966.get("preferred_search_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
        and f966.get("counts_as_lawful_supplier") is False
    )

    post_f970_rejoin_to_f966_confirmed = (
        p1114.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_REJOIN_AUDITED"
        and p1114.get("rejoin_existing_f966_non_same_lane_route") is True
        and f971.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_PACKET_EXPORTED"
        and f971.get("packet_exported_on_current_repo_state") is True
        and f971.get("primary_continuation_route") == "existing_f966_non_same_lane_upgrade_route_contract"
    )

    scpc_like_lane_reference_only_confirmed = (
        p1011.get("status") == "P1011_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_SELECTOR_INTERFACE_BLOCKED"
        and p1011.get("info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted") is True
        and p1011.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
        and p1011.get("strict_selector_interface_exported") is False
        and p1011.get("provider_class_shift_realized") is False
    )

    neural_support_lane_reference_only_confirmed = (
        p1034.get("status") == "P1034_CURRENT_STRICT_QW2191_NADSOLITON_NEURAL_CHARACTER_INFORMATION_PRIMARY_SELECTOR_SUPPORT_REFERENCE_ADMITTED_INTERFACE_STILL_BLOCKED"
        and p1034.get("nadsoliton_neural_character_support_reference_admitted") is True
        and p1034.get("support_reference_grade") == "cross_repo_support_reference_only"
        and p1034.get("strict_selector_interface_exported") is False
        and p1034.get("strict_selector_source_exported") is False
    )

    exported_selector_source_frontier_confirmed = (
        p1047.get("status") == "PASS_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1047.get("strict_int_selector_source_frontier_open") is True
        and p1047.get("primary_continuation_route") == PRIMARY_SUBROUTE
        and f956.get("status") == "PASS_CURRENT_STRICT_QW2191_POST_STOP_NEURAL_BRIDGE_LANE_TO_STRICT_INT_SELECTOR_SOURCE_FRONTIER_ROUTE_DECISION_PACKET_EXPORTED"
        and f956.get("packet_exported_on_current_repo_state") is True
        and f956.get("primary_continuation_route") == PRIMARY_SUBROUTE
    )

    rejoin_to_existing_t173_frontier_confirmed = (
        p1048.get("status") == "PASS_CURRENT_STRICT_QW2191_POST_STOP_ROUTE_REJOIN_TO_EXISTING_T173_FRONTIER_AUDITED"
        and p1048.get("primary_continuation_target") == PRIMARY_TARGET
    )

    strongest_existing_exported_subroute_is_selector_source_frontier = (
        existing_f966_live_contract_confirmed
        and post_f970_rejoin_to_f966_confirmed
        and scpc_like_lane_reference_only_confirmed
        and neural_support_lane_reference_only_confirmed
        and exported_selector_source_frontier_confirmed
        and rejoin_to_existing_t173_frontier_confirmed
    )

    add_check("existing_f966_live_contract_confirmed", existing_f966_live_contract_confirmed, True, "P1091/F966 already freeze one live non-same-lane contract.")
    add_check("post_f970_rejoin_to_f966_confirmed", post_f970_rejoin_to_f966_confirmed, True, "P1114/F971 already rejoin current work to the live F966 contract after the ONRD stop.")
    add_check("scpc_like_lane_reference_only_confirmed", scpc_like_lane_reference_only_confirmed, True, "The SCPC-like lane remains reference-only and does not export a strict selector interface.")
    add_check("neural_support_lane_reference_only_confirmed", neural_support_lane_reference_only_confirmed, True, "The neural-character lane remains support-reference-only and does not export a strict selector source.")
    add_check("exported_selector_source_frontier_confirmed", exported_selector_source_frontier_confirmed, True, "P1047/F956 already export the strict selector-source frontier as the strongest live strict continuation after the neural bridge stop.")
    add_check("rejoin_to_existing_t173_frontier_confirmed", rejoin_to_existing_t173_frontier_confirmed, True, "P1048 already rejoins that selector-source frontier to existing T173.")
    add_check("strongest_existing_exported_subroute_is_selector_source_frontier", strongest_existing_exported_subroute_is_selector_source_frontier, True, "Therefore the strongest already exported subroute inside the live F966 contract is the strict selector-source frontier rejoining T173.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_DECISION_AUDITED"
        if not blocking and strongest_existing_exported_subroute_is_selector_source_frontier
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_DECISION_AUDIT"
    )

    artifact = {
        "stage": "P1115",
        "status": status,
        "as_of": AS_OF,
        "strongest_existing_exported_subroute": PRIMARY_SUBROUTE,
        "primary_existing_subroute_target": PRIMARY_TARGET,
        "active_missing_bridge": ACTIVE_MISSING_BRIDGE,
        "non_selected_existing_lanes": {
            "lambda_branch_info_primary_scpc_like": "reference_context_candidate_only",
            "nadsoliton_neural_character_support_reference": "cross_repo_support_reference_only",
            "minimal_onrd_same_lane": "stopped_as_primary_move",
        },
        "current_primary_work_contract": "anchor_existing_f966_non_same_lane_search_on_exported_strict_selector_source_frontier_subroute",
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "strongest_existing_exported_subroute": artifact["strongest_existing_exported_subroute"],
        "primary_existing_subroute_target": artifact["primary_existing_subroute_target"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_primary_work_contract": artifact["current_primary_work_contract"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
