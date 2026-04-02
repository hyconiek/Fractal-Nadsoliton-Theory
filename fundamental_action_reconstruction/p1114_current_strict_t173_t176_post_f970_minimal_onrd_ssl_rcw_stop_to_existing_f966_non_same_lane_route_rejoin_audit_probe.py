#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1113 = GENERATED / "p1113_current_strict_t173_t176_minimal_onrd_ssl_rcw_stagnation_audit_probe_summary.json"
IN_F970 = GENERATED / "f970_current_strict_t173_t176_minimal_onrd_ssl_rcw_stop_packet_summary.json"
IN_P1091 = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe_summary.json"
IN_F966 = GENERATED / "f966_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_packet_summary.json"

OUT_JSON = GENERATED / "p1114_current_strict_t173_t176_post_f970_onrd_to_f966_rejoin_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1114_current_strict_t173_t176_post_f970_onrd_to_f966_rejoin_audit_probe_summary.json"

ACTIVE_MISSING_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1113, IN_F970, IN_P1091, IN_F966]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1114",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1113 = load_json(IN_P1113)
    f970 = load_json(IN_F970)
    p1091 = load_json(IN_P1091)
    f966 = load_json(IN_F966)

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

    f970_stop_packet_already_exported = (
        p1113.get("status") == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_AND_STOP_AUDITED"
        and p1113.get("same_lane_stagnation_boundary_reached") is True
        and f970.get("status") == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_SHARPER_SAME_LANE_ROUTE_COHERENCE_WITNESS_REFINEMENT_STAGNATION_STOP_PACKET_EXPORTED"
        and f970.get("packet_exported_on_current_repo_state") is True
        and f970.get("same_lane_deeper_route_coherence_witness_descent_disallowed_as_primary_move") is True
    )

    existing_f966_non_same_lane_route_packet_already_exported = (
        p1091.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_DECISION_AUDITED"
        and f966.get("status") == "F966_EXECUTED_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_PACKET_NO_FALSE_PASS"
        and f966.get("route_object_id") == "PostF965FailureMapConstrainedToExportedNoncyclicProviderSplitNonSameLaneUpgradeRoute_v1"
        and f966.get("preferred_search_family") == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
    )

    active_missing_bridge_unchanged = (
        p1091.get("active_missing_bridge") == ACTIVE_MISSING_BRIDGE
    )

    strongest_honest_post_f970_continuation_rejoins_existing_f966_non_same_lane_route = (
        f970_stop_packet_already_exported
        and existing_f966_non_same_lane_route_packet_already_exported
        and active_missing_bridge_unchanged
        and f970.get("restart_requires_new_blocker_cut_or_non_same_lane_upgrade_route") is True
        and p1091.get("allowed_next_move_contract")
        == "search_one_genuinely_new_non_same_lane_upgrade_route_within_exported_noncyclic_provider_split_family_or_one_genuinely_new_inversion_sensitive_source_side_provider_class"
    )

    add_check("f970_stop_packet_already_exported", f970_stop_packet_already_exported, True, "P1113/F970 already freeze that the ONRD same-lane route is stopped as a primary strategy.")
    add_check("existing_f966_non_same_lane_route_packet_already_exported", existing_f966_non_same_lane_route_packet_already_exported, True, "P1091/F966 already freeze the admissible non-same-lane continuation contract.")
    add_check("active_missing_bridge_unchanged", active_missing_bridge_unchanged, True, "The active missing bridge remains unchanged across the stop and rejoin context.")
    add_check("strongest_honest_post_f970_continuation_rejoins_existing_f966_non_same_lane_route", strongest_honest_post_f970_continuation_rejoins_existing_f966_non_same_lane_route, True, "Therefore the strongest honest post-F970 continuation is to rejoin the existing F966 non-same-lane route contract.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_REJOIN_AUDITED"
        if not blocking and strongest_honest_post_f970_continuation_rejoins_existing_f966_non_same_lane_route
        else "FAIL_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_REJOIN_AUDIT"
    )

    artifact = {
        "stage": "P1114",
        "status": status,
        "as_of": AS_OF,
        "rejoin_existing_f966_non_same_lane_route": strongest_honest_post_f970_continuation_rejoins_existing_f966_non_same_lane_route,
        "active_missing_bridge": ACTIVE_MISSING_BRIDGE,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "rejoin_existing_f966_non_same_lane_route": artifact["rejoin_existing_f966_non_same_lane_route"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
