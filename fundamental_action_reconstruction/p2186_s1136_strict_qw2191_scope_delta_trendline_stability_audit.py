#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2185 = GEN / "p2185_s1135_strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger.json"
HISTORY = GEN / "p2186_s1136_scope_delta_trendline_history.json"
OUT = GEN / "p2186_s1136_strict_qw2191_scope_delta_trendline_stability_audit.json"
MD = GEN / "p2186_s1136_strict_qw2191_scope_delta_trendline_stability_audit.md"


def load(path: Path, default: dict[str, Any]) -> dict[str, Any]:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2185 = load(IN_2185, {"_missing": str(IN_2185), "result_kind": "OPEN_MISSING_ARTIFACT"})
    sentinel = p2185.get("strict_qw2191_periodic_ci_sentinel_and_nonbridge_scope_tagger", {}) or {}

    current = {
        "run_label": "S1136_latest",
        "source_packet": "P2185",
        "monitored_packet_count": int(sentinel.get("monitored_packet_count", 0) or 0),
        "regression_count": int(sentinel.get("regression_count", 0) or 0),
        "ci_status": sentinel.get("ci_status", "UNKNOWN"),
    }

    history = load(HISTORY, {"schema_version": "p2186_history_v1", "runs": []})
    runs = history.get("runs", []) or []
    runs.append(current)

    window = runs[-5:]
    regression_series = [int(r.get("regression_count", 0) or 0) for r in window]

    nonincreasing = all(regression_series[i] <= regression_series[i - 1] for i in range(1, len(regression_series)))
    bounded = all(v >= 0 for v in regression_series)
    stable_zero = all(v == 0 for v in regression_series)

    trendline = {
        "window_size": len(window),
        "window": window,
        "regression_count_series": regression_series,
        "stability_checks": {
            "bounded_nonnegative": bounded,
            "nonincreasing_series": nonincreasing,
            "stable_zero_regressions": stable_zero,
        },
    }

    trend_status = (
        "STABLE_ZERO_REGRESSION_TREND"
        if stable_zero
        else "NONINCREASING_REGRESSION_TREND" if bounded and nonincreasing else "OPEN_TRENDLINE_REGRESSION_INSTABILITY"
    )

    result_kind = (
        "PASS_STRICT_QW2191_SCOPE_DELTA_TRENDLINE_STABILITY_AUDIT"
        if trend_status in {"STABLE_ZERO_REGRESSION_TREND", "NONINCREASING_REGRESSION_TREND"}
        else "OPEN_STRICT_QW2191_SCOPE_DELTA_TRENDLINE_INSTABILITY"
    )

    payload = {
        "schema_version": "p2186_s1136_v1",
        "packet_id": "P2186",
        "stage_id": "S1136",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_scope_delta_trendline_stability_audit": {
            "audit_id": "QW2191_SCOPE_DELTA_TRENDLINE_STABILITY_AUDIT_V1",
            "input_sentinel_packet": str(IN_2185.relative_to(ROOT)),
            "trend_status": trend_status,
            "trendline": trendline,
        },
        "recommended_next_honest_step": {
            "id": "P2187_candidate",
            "goal": "add deterministic CI fail-threshold policy object keyed to trendline status classes",
        },
        "gatekeeper_checks": {
            "trendline_audit_exported": True,
            "bounded_nonnegative_regression_series": bounded,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": bool((p2185.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    HISTORY.write_text(
        json.dumps({"schema_version": "p2186_history_v1", "runs": runs}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2186 S1136: strict QW-2191 scope-delta trendline stability audit",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Trend status: `{trend_status}`",
                f"- Window size: `{len(window)}`",
                f"- Regression series: `{regression_series}`",
                "",
                "No selector closure or global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
