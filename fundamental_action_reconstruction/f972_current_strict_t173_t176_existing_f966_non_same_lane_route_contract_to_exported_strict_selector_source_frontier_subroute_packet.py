#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1115 = GENERATED / "p1115_current_strict_t173_t176_f966_selector_source_subroute_audit_probe_summary.json"
IN_F972 = ROOT / "F972_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_PACKET.md"

OUT_JSON = GENERATED / "f972_current_strict_t173_t176_f966_selector_source_subroute_packet.json"
OUT_SUMMARY = GENERATED / "f972_current_strict_t173_t176_f966_selector_source_subroute_packet_summary.json"

PRIMARY_SUBROUTE = "explicit_strict_internal_selector_source_derivation_frontier"
PRIMARY_TARGET = "T173_CURRENT_STRICT_CORE_SELECTOR_CLOSURE_AND_KERNEL_ALONE_QW2191_DISCHARGE_TARGET_SPEC"
ACTIVE_MISSING_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P1115, IN_F972]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F972",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1115 = load_json(IN_P1115)
    f972_text = load_text(IN_F972)

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

    p1115_subroute_decision_passed = (
        p1115.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_DECISION_AUDITED"
        and p1115.get("strongest_existing_exported_subroute") == PRIMARY_SUBROUTE
        and p1115.get("primary_existing_subroute_target") == PRIMARY_TARGET
        and p1115.get("active_missing_bridge") == ACTIVE_MISSING_BRIDGE
    )

    f972_packet_shape_frozen = all(
        needle in f972_text for needle in [
            "Xi_current_strict_t173_t176_existing_f966_non_same_lane_route_contract_to_exported_strict_selector_source_frontier_subroute_packet_v1",
            "existing_f966_live_contract_confirmed := yes",
            "post_f970_rejoin_to_f966_confirmed := yes",
            "scpc_like_lane_reference_only_confirmed := yes",
            "neural_support_lane_reference_only_confirmed := yes",
            "exported_selector_source_frontier_confirmed := yes",
            "rejoin_to_existing_t173_frontier_confirmed := yes",
            f"strongest_existing_exported_subroute := {PRIMARY_SUBROUTE}",
            f"primary_existing_subroute_target := {PRIMARY_TARGET}",
            f"active_missing_bridge := {ACTIVE_MISSING_BRIDGE}",
            "current_primary_work_contract := anchor_existing_f966_non_same_lane_search_on_exported_strict_selector_source_frontier_subroute",
        ]
    )

    packet_exported_on_current_repo_state = p1115_subroute_decision_passed and f972_packet_shape_frozen

    add_check("p1115_subroute_decision_passed", p1115_subroute_decision_passed, True, "P1115 already freezes the strongest existing exported subroute decision positively.")
    add_check("f972_packet_shape_frozen", f972_packet_shape_frozen, True, "F972 freezes the selector-source-subroute packet shape explicitly.")
    add_check("packet_exported_on_current_repo_state", packet_exported_on_current_repo_state, True, "Therefore the current repo exports one honest packet anchoring the live F966 contract on the already exported strict selector-source frontier subroute.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F966_NON_SAME_LANE_ROUTE_CONTRACT_TO_EXPORTED_STRICT_SELECTOR_SOURCE_FRONTIER_SUBROUTE_PACKET"
    )

    artifact = {
        "stage": "F972",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "strongest_existing_exported_subroute": PRIMARY_SUBROUTE,
        "primary_existing_subroute_target": PRIMARY_TARGET,
        "active_missing_bridge": ACTIVE_MISSING_BRIDGE,
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
        "packet_exported_on_current_repo_state": artifact["packet_exported_on_current_repo_state"],
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
