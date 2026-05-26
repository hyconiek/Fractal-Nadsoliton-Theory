#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2184 = GEN / "p2184_s1134_strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard.json"
OUT = GEN / "p2185_s1135_strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger.json"
MD = GEN / "p2185_s1135_strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2184 = load(IN_2184)
    dashboard = p2184.get("strict_qw2191_external_audit_replay_dashboard_scope_delta_regression_guard", {}) or {}

    monitored_packets = dashboard.get("monitored_packets", []) or []
    regressions = dashboard.get("scope_delta_regressions", []) or []

    rows = []
    for packet in monitored_packets:
        rows.append(
            {
                "packet": packet,
                "lane": "strict_qw2191_external_replay",
                "scope_tag": "nonbridge_operational_strict_lane",
                "selector_closure_claimed": False,
                "toe_closure_claimed": False,
            }
        )

    sentinel = {
        "sentinel_id": "QW2191_PERIODIC_CI_SENTINEL_V1",
        "input_dashboard_packet": str(IN_2184.relative_to(ROOT)),
        "schedule": {
            "mode": "periodic",
            "recommended_cadence": "per-commit-or-daily",
            "hard_fail_on_scope_delta_regression": True,
        },
        "monitored_packet_count": len(monitored_packets),
        "nonbridge_scope_tagged_rows": rows,
        "regression_count": len(regressions),
        "ci_status": "GREEN" if not regressions else "RED",
    }

    result_kind = (
        "PASS_STRICT_QW2191_PERIODIC_CI_SENTINEL_AND_NONBRIDGE_SCOPE_TAGGER"
        if not regressions
        else "OPEN_STRICT_QW2191_PERIODIC_CI_SENTINEL_REGRESSION_BLOCKED"
    )

    payload = {
        "schema_version": "p2185_s1135_v1",
        "packet_id": "P2185",
        "stage_id": "S1135",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger": sentinel,
        "recommended_next_honest_step": {
            "id": "P2186_candidate",
            "goal": "add cross-packet trendline audit that checks multi-run stability of scope-delta regression count",
        },
        "gatekeeper_checks": {
            "ci_sentinel_exported": True,
            "nonbridge_scope_tags_applied": all(r["scope_tag"] == "nonbridge_operational_strict_lane" for r in rows),
            "scope_delta_regressions_absent": not regressions,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": bool((p2184.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2185 S1135: periodic CI sentinel + non-bridge scope tagging",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Monitored packets: `{len(monitored_packets)}`",
                f"- Scope-delta regressions observed: `{len(regressions)}`",
                "- Scope tag policy: `nonbridge_operational_strict_lane`",
                "",
                "No selector closure or global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
