#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2139 = GEN / "p2139_s1089_strict_cmp2_nonsynthetic_readiness_freeze_report.json"
OUT = GEN / "p2140_s1090_strict_cmp2_blocked_path_resolution_router.json"
MD = GEN / "p2140_s1090_strict_cmp2_blocked_path_resolution_router.md"

SCHEMA_VERSION = "p2140_s1090_v1"
TIMESTAMP_UTC = "2026-05-26T00:00:00+00:00"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path)}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2139 = load(IN_2139)

    report = p2139.get("nonsynthetic_readiness_freeze_report", {}) or {}
    blocked = report.get("missing_or_blocked_stages", []) or []
    fully_unblocked = bool(report.get("fully_unblocked", False))

    if fully_unblocked:
        route = "ADVANCE_TO_NEXT_FORMAL_OBSTRUCTION_FAMILY"
        actions = [
            "freeze_non_synthetic_cmp2_ci_inflation_report",
            "open_next_strict_obstruction_workstream",
        ]
    else:
        route = "DATA_ACQUISITION_AND_RERUN_REQUIRED"
        actions = [
            "deliver_real_cmp2_backend_rows_extension_v1_json",
            "rerun_p2138_orchestrator",
            "confirm_p2132_to_p2135_pass",
        ]

    payload = {
        "schema_version": SCHEMA_VERSION,
        "packet_id": "P2140",
        "stage_id": "S1090",
        "produced_by": Path(__file__).name,
        "timestamp_utc": TIMESTAMP_UTC,
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CMP2_BLOCKED_PATH_RESOLUTION_ROUTER_WITH_TRACE",
        "blocked_path_resolution_router": {
            "source_report": str(IN_2139.relative_to(ROOT)),
            "fully_unblocked": fully_unblocked,
            "blocked_stages": blocked,
            "selected_route": route,
            "action_queue": actions,
            "scope_limit": "routing/decision support only; not theorem-grade closure",
        },
        "recommended_next_honest_step": {
            "id": "P2141_candidate",
            "goal": "execute action_queue and re-evaluate router after P2138/P2139 refresh",
        },
        "gatekeeper_checks": {
            "router_exported": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": False,
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": False,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join([
            "# P2140 S1090: strict CMP2 blocked-path resolution router",
            "",
            f"- Status: `{payload['status']}`",
            f"- Result kind: `{payload['result_kind']}`",
            f"- Selected route: `{route}`",
            "",
            "This stage converts readiness-freeze state into explicit action routing.",
            "No theorem-grade global closure claim is made.",
        ]) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
