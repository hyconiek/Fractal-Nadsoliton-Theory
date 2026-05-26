#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2183 = GEN / "p2183_s1133_strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge.json"
OUT = GEN / "p2184_s1134_strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard.json"
MD = GEN / "p2184_s1134_strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2183 = load(IN_2183)
    bridge = p2183.get("strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge", {}) or {}
    rows = bridge.get("ledger_refresh_bridge_rows", []) or []

    regressions = []
    for row in rows:
        if not row.get("scope_delta_checked", False):
            regressions.append({"packet": row.get("packet"), "issue": "scope_delta_not_checked"})
        if row.get("scope_delta_blocking_present", False):
            regressions.append({"packet": row.get("packet"), "issue": "scope_delta_blocking_present"})

    dashboard = {
        "dashboard_id": "QW2191_EXTERNAL_AUDIT_REPLAY_DASHBOARD_V1",
        "source_governance_bridge_packet": str(IN_2183.relative_to(ROOT)),
        "monitored_packets": [r.get("packet") for r in rows],
        "row_count": len(rows),
        "scope_delta_regression_count": len(regressions),
        "scope_delta_regressions": regressions,
        "status": "PASS_NO_SCOPE_DELTA_REGRESSIONS" if not regressions else "OPEN_SCOPE_DELTA_REGRESSION",
    }

    result_kind = (
        "PASS_STRICT_QW2191_EXTERNAL_AUDIT_REPLAY_DASHBOARD_SCOPE_DELTA_REGRESSION_GUARD"
        if not regressions
        else "OPEN_STRICT_QW2191_SCOPE_DELTA_REGRESSION_GUARD"
    )

    payload = {
        "schema_version": "p2184_s1134_v1",
        "packet_id": "P2184",
        "stage_id": "S1134",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard": dashboard,
        "recommended_next_honest_step": {
            "id": "P2185_candidate",
            "goal": "promote this dashboard into periodic CI sentinel checks with explicit non-bridge scope tagging",
        },
        "gatekeeper_checks": {
            "dashboard_exported": True,
            "scope_delta_regressions_absent": not regressions,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2183.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2183.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2184 S1134: strict QW-2191 external audit replay dashboard + scope-delta regression guard",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Monitored packets: `{len(rows)}`",
                f"- Scope-delta regressions: `{len(regressions)}`",
                "",
                "No selector closure or global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
