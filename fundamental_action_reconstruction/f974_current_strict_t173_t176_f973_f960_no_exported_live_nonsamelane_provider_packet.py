#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

AS_OF = "2026-04-01"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1117 = GENERATED / "p1117_current_strict_t173_t176_f973_f960_no_exported_live_nonsamelane_provider_audit_probe_summary.json"
IN_F974 = ROOT / "F974_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_PACKET.md"

OUT_JSON = GENERATED / "f974_current_strict_t173_t176_f973_f960_no_exported_live_nonsamelane_provider_packet.json"
OUT_SUMMARY = GENERATED / "f974_current_strict_t173_t176_f973_f960_no_exported_live_nonsamelane_provider_packet_summary.json"

ACTIVE_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"
NEW_WORK_CONTRACT = "build_one_genuinely_new_narrow_probe_against_existing_f960_bridge_target"


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

    prerequisites = [IN_P1117, IN_F974]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "F974",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1117 = load_json(IN_P1117)
    f974_text = load_text(IN_F974)

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

    p1117_absence_audit_passed = (
        p1117.get("status") == "PASS_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_AUDITED"
        and p1117.get("already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present") is False
        and p1117.get("active_missing_bridge") == ACTIVE_BRIDGE
        and p1117.get("current_primary_work_contract") == NEW_WORK_CONTRACT
    )

    f974_packet_shape_frozen = all(
        needle in f974_text for needle in [
            "Xi_current_strict_t173_t176_existing_f973_f960_no_already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_packet_v1",
            "active_f960_bridge_frontier_confirmed := yes",
            "old_pair12_provider_frontier_exists_but_is_nonprimary := yes",
            "pair12_same_lane_reentry_disallowed := yes",
            "live_contract_requires_genuinely_new_non_same_lane_or_new_provider_class := yes",
            "already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present := no",
            f"current_primary_work_contract := {NEW_WORK_CONTRACT}",
            f"active_missing_bridge := {ACTIVE_BRIDGE}",
        ]
    )

    packet_exported_on_current_repo_state = p1117_absence_audit_passed and f974_packet_shape_frozen

    add_check("p1117_absence_audit_passed", p1117_absence_audit_passed, True, "P1117 already freezes the absence of an already-live non-same-lane provider candidate positively.")
    add_check("f974_packet_shape_frozen", f974_packet_shape_frozen, True, "F974 freezes the negative decision packet shape explicitly.")
    add_check("packet_exported_on_current_repo_state", packet_exported_on_current_repo_state, True, "Therefore the repo exports one honest packet freezing that no already-live non-same-lane inversion-sensitive source-side provider candidate is currently available beneath F960.")

    status = (
        "PASS_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_PACKET_EXPORTED"
        if not blocking and packet_exported_on_current_repo_state
        else "FAIL_CURRENT_STRICT_T173_T176_EXISTING_F973_F960_NO_ALREADY_EXPORTED_LIVE_NON_SAME_LANE_INVERSION_SENSITIVE_SOURCE_SIDE_PROVIDER_CANDIDATE_PACKET"
    )

    artifact = {
        "stage": "F974",
        "status": status,
        "as_of": AS_OF,
        "packet_exported_on_current_repo_state": packet_exported_on_current_repo_state,
        "already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present": False,
        "active_missing_bridge": ACTIVE_BRIDGE,
        "current_primary_work_contract": NEW_WORK_CONTRACT,
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
        "already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present": artifact["already_exported_live_non_same_lane_inversion_sensitive_source_side_provider_candidate_present"],
        "active_missing_bridge": artifact["active_missing_bridge"],
        "current_primary_work_contract": artifact["current_primary_work_contract"],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
