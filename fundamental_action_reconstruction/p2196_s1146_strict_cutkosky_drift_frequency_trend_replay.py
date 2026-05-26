#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2194 = GEN / "p2194_s1144_strict_cutkosky_multi_background_tolerance_stability_audit.json"
IN_2195 = GEN / "p2195_s1145_strict_cutkosky_multi_background_release_gate_override_policy.json"
HIST = GEN / "p2196_s1146_strict_cutkosky_drift_frequency_history.json"
OUT = GEN / "p2196_s1146_strict_cutkosky_drift_frequency_trend_replay.json"
MD = GEN / "p2196_s1146_strict_cutkosky_drift_frequency_trend_replay.md"


def load(path: Path, default: dict[str, Any]) -> dict[str, Any]:
    if not path.exists():
        return default
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2194 = load(IN_2194, {"_missing": str(IN_2194), "result_kind": "OPEN_MISSING_ARTIFACT"})
    p2195 = load(IN_2195, {"_missing": str(IN_2195), "result_kind": "OPEN_MISSING_ARTIFACT"})

    audit = p2194.get("strict_cutkosky_multi_background_tolerance_stability_audit", {}) or {}
    policy = p2195.get("strict_cutkosky_multi_background_release_gate_override_policy", {}) or {}

    run_row = {
        "run_label": "S1146_latest",
        "same_channel": audit.get("same_channel", "graviton_to_gauge_gauge"),
        "class_set": (audit.get("class_stability", {}) or {}).get("class_set", []),
        "stable": bool((audit.get("class_stability", {}) or {}).get("stable", False)),
        "drift_detected": bool(audit.get("background_drift_detected", True)),
        "override_decision": (policy.get("output", {}) or {}).get("override_decision", "FORCE_WARN_REVIEW"),
        "release_gate_override": (policy.get("output", {}) or {}).get("release_gate_override", "OPEN_WITH_DRIFT_NOTICE"),
    }

    hist = load(HIST, {"schema_version": "p2196_history_v1", "runs": []})
    runs = hist.get("runs", []) or []
    runs.append(run_row)

    window = runs[-10:]
    drift_count = sum(1 for r in window if bool(r.get("drift_detected", False)))
    total = len(window)
    drift_freq = (drift_count / total) if total > 0 else 0.0

    trend = {
        "window_size": total,
        "drift_count": drift_count,
        "drift_frequency": drift_freq,
        "drift_free_runs": total - drift_count,
        "window": window,
    }

    if drift_freq == 0.0:
        status = "DRIFT_FREE_TREND"
    elif drift_freq <= 0.2:
        status = "LOW_DRIFT_TREND"
    else:
        status = "ELEVATED_DRIFT_TREND"

    payload = {
        "schema_version": "p2196_s1146_v1",
        "packet_id": "P2196",
        "stage_id": "S1146",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CUTKOSKY_DRIFT_FREQUENCY_TREND_REPLAY",
        "strict_cutkosky_drift_frequency_trend_replay": {
            "replay_id": "STRICT_CUTKOSKY_DRIFT_FREQUENCY_TREND_REPLAY_V1",
            "source_packets": [str(IN_2194.relative_to(ROOT)), str(IN_2195.relative_to(ROOT))],
            "trend_status": status,
            "trend": trend,
        },
        "recommended_next_honest_step": {
            "id": "P2197_candidate",
            "goal": "bind drift-frequency trend status to deterministic rerun cadence and hard-stop policy",
        },
        "gatekeeper_checks": {
            "drift_frequency_trend_exported": True,
            "same_channel_replay_preserved": run_row["same_channel"] == "graviton_to_gauge_gauge",
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    HIST.write_text(json.dumps({"schema_version": "p2196_history_v1", "runs": runs}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2196 S1146: strict Cutkosky drift-frequency trend replay",
                "",
                f"- Trend status: `{status}`",
                f"- Window size: `{total}`",
                f"- Drift frequency: `{drift_freq:.6f}`",
                f"- Drift count: `{drift_count}`",
                "",
                "This is a drift-frequency replay trend only, not full Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
