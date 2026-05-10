#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1116 = GENERATED / "p1116_current_strict_t173_t176_f972_to_f960_active_bridge_route_audit_probe_summary.json"
IN_F973 = ROOT / "F973_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_PACKET.md"

OUT_JSON = GENERATED / "f973_current_strict_t173_t176_f972_to_f960_active_bridge_route_packet.json"
OUT_SUMMARY = GENERATED / "f973_current_strict_t173_t176_f972_to_f960_active_bridge_route_packet_summary.json"

ACTIVE_FRONTIER = "existing_t183_residual_datum_pair12_orbit_direction_selection_frontier"
ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
CURRENT_WORK_CONTRACT = "build_one_narrow_probe_for_a_genuinely_new_inversion_sensitive_source_side_provider_class_against_the_existing_f960_bridge_target"


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

    prerequisites = [IN_P1116, IN_F973]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F973",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1116 = load_json(IN_P1116)
    f973_text = load_text(IN_F973)

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

    p1116_route_decision_passed = (
        p1116.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_DECISION_AUDITED"
        and p1116.get("strongest_exact_active_frontier") == ACTIVE_FRONTIER
        and p1116.get("strongest_exact_active_frontier_target") == ACTIVE_BRIDGE
        and p1116.get("current_primary_work_contract") == CURRENT_WORK_CONTRACT
    )

    f973_packet_shape_frozen = all(
        needle in f973_text for needle in [
            "Xi_current_strict_t173_t176_existing_f972_selector_source_subroute_to_existing_f960_active_bridge_frontier_route_packet_v1",
            "existing_f972_selector_source_subroute_confirmed := yes",
            "t173_residual_pair12_split_localized := yes",
            "existing_t183_frontier_confirmed := yes",
            "no_current_lawful_supplier_confirmed := yes",
            "exact_active_bridge_target_exported := yes",
            f"strongest_exact_active_frontier := {ACTIVE_FRONTIER}",
            f"strongest_exact_active_frontier_target := {ACTIVE_BRIDGE}",
            f"current_primary_work_contract := {CURRENT_WORK_CONTRACT}",
            f"active_missing_bridge := {ACTIVE_BRIDGE}",
        ]
    )

    packet_exported_on_current_repo_state = p1116_route_decision_passed and f973_packet_shape_frozen

    add_check("p1116_route_decision_passed", p1116_route_decision_passed, True, "P1116 already freezes the strongest exact active frontier decision positively.")
    add_check("f973_packet_shape_frozen", f973_packet_shape_frozen, True, "F973 freezes the active-bridge-frontier packet shape explicitly.")
    add_check("packet_exported_on_current_repo_state", packet_exported_on_current_repo_state, True, "Therefore the repo exports one honest packet anchoring the live F972 selector-source subroute on the existing F960/T183 exact active bridge frontier.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F972_SELECTOR_SOURCE_SUBROUTE_TO_EXISTING_F960_ACTIVE_BRIDGE_FRONTIER_ROUTE_PACKET"
    )

    artifact = {
        "stage": "F973",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "strongest_exact_active_frontier": ACTIVE_FRONTIER,
        "strongest_exact_active_frontier_target": ACTIVE_BRIDGE,
        "current_primary_work_contract": CURRENT_WORK_CONTRACT,
        "no_false_pass": True,
        "checks": checks,
        "blocking_checks": blocking,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "packet_exported_on_current_repo_state": artifact["packet_exported_on_current_repo_state"],
        "strongest_exact_active_frontier": artifact["strongest_exact_active_frontier"],
        "strongest_exact_active_frontier_target": artifact["strongest_exact_active_frontier_target"],
        "current_primary_work_contract": artifact["current_primary_work_contract"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
