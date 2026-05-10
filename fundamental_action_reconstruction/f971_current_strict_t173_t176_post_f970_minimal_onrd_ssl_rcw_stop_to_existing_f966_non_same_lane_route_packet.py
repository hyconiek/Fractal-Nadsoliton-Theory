#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1114 = GENERATED / "p1114_current_strict_t173_t176_post_f970_onrd_to_f966_rejoin_audit_probe_summary.json"
IN_F971 = ROOT / "F971_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_PACKET.md"

OUT_JSON = GENERATED / "f971_current_strict_t173_t176_post_f970_onrd_to_f966_route_packet.json"
OUT_SUMMARY = GENERATED / "f971_current_strict_t173_t176_post_f970_onrd_to_f966_route_packet_summary.json"


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

    prerequisites = [IN_P1114, IN_F971]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F971",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1114 = load_json(IN_P1114)
    f971_text = load_text(IN_F971)

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

    p1114_rejoin_audit_passed = (
        p1114.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_REJOIN_AUDITED"
        and p1114.get("rejoin_existing_f966_non_same_lane_route") is True
    )

    f971_packet_shape_frozen = all(
        needle in f971_text for needle in [
            "Xi_current_strict_t173_t176_post_f970_minimal_onrd_ssl_rcw_stop_to_existing_f966_non_same_lane_route_packet_v1",
            "stopped_onrd_same_lane_confirmed := yes",
            "existing_f966_non_same_lane_route_packet_confirmed := yes",
            "active_missing_bridge_unchanged := yes",
            "primary_continuation_route := existing_f966_non_same_lane_upgrade_route_contract",
            "current_primary_work_contract := rejoin_existing_f966_non_same_lane_route_do_not_spawn_competing_post_f970_same_lane",
        ]
    )

    packet_exported_on_current_repo_state = p1114_rejoin_audit_passed and f971_packet_shape_frozen

    add_check("p1114_rejoin_audit_passed", p1114_rejoin_audit_passed, True, "P1114 already freezes the post-F970 rejoin audit positively.")
    add_check("f971_packet_shape_frozen", f971_packet_shape_frozen, True, "F971 freezes the route-rejoin packet shape explicitly.")
    add_check("packet_exported_on_current_repo_state", packet_exported_on_current_repo_state, True, "Therefore the current repo exports one honest post-F970 rejoin packet to the existing F966 non-same-lane route.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_POST_F970_MINIMAL_ONRD_SSL_RCW_STOP_TO_EXISTING_F966_NON_SAME_LANE_ROUTE_PACKET"
    )

    artifact = {
        "stage": "F971",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "primary_continuation_route": "existing_f966_non_same_lane_upgrade_route_contract",
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
        "primary_continuation_route": artifact["primary_continuation_route"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
