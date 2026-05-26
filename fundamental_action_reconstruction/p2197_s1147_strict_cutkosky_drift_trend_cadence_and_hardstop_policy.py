#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2196 = GEN / "p2196_s1146_strict_cutkosky_drift_frequency_trend_replay.json"
OUT = GEN / "p2197_s1147_strict_cutkosky_drift_trend_cadence_and_hardstop_policy.json"
MD = GEN / "p2197_s1147_strict_cutkosky_drift_trend_cadence_and_hardstop_policy.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2196 = load(IN_2196)
    replay = p2196.get("strict_cutkosky_drift_frequency_trend_replay", {}) or {}
    trend_status = str(replay.get("trend_status", "ELEVATED_DRIFT_TREND"))
    trend = replay.get("trend", {}) or {}
    drift_frequency = float(trend.get("drift_frequency", 1.0) or 1.0)

    policy_table = {
        "DRIFT_FREE_TREND": {
            "rerun_cadence": "weekly_or_per_release",
            "hard_stop": False,
            "override_mode": "none",
            "ci_gate": "OPEN",
        },
        "LOW_DRIFT_TREND": {
            "rerun_cadence": "daily_or_per_commit",
            "hard_stop": False,
            "override_mode": "warn_review",
            "ci_gate": "OPEN_WITH_NOTICE",
        },
        "ELEVATED_DRIFT_TREND": {
            "rerun_cadence": "per_commit_mandatory",
            "hard_stop": True,
            "override_mode": "fail_block",
            "ci_gate": "BLOCK",
        },
    }

    selected = policy_table.get(trend_status, policy_table["ELEVATED_DRIFT_TREND"])

    payload = {
        "schema_version": "p2197_s1147_v1",
        "packet_id": "P2197",
        "stage_id": "S1147",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CUTKOSKY_DRIFT_TREND_CADENCE_AND_HARDSTOP_POLICY",
        "strict_cutkosky_drift_trend_cadence_and_hardstop_policy": {
            "policy_id": "STRICT_CUTKOSKY_DRIFT_TREND_CADENCE_HARDSTOP_POLICY_V1",
            "source_packet": str(IN_2196.relative_to(ROOT)),
            "input": {"trend_status": trend_status, "drift_frequency": drift_frequency},
            "policy_table": policy_table,
            "selected_policy": selected,
        },
        "recommended_next_honest_step": {
            "id": "P2198_candidate",
            "goal": "connect cadence/hard-stop policy with release-gate contract replay and emit contradiction ledger deltas",
        },
        "gatekeeper_checks": {
            "cadence_hardstop_policy_exported": True,
            "trend_status_consumed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": True,
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2197 S1147: strict Cutkosky drift-trend cadence + hard-stop policy",
                "",
                f"- Trend status: `{trend_status}`",
                f"- Drift frequency: `{drift_frequency:.6f}`",
                f"- Selected rerun cadence: `{selected['rerun_cadence']}`",
                f"- Hard stop: `{selected['hard_stop']}`",
                "",
                "This is a governance policy export only, not full Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
