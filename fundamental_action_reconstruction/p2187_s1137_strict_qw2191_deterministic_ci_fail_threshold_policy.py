#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2186 = GEN / "p2186_s1136_strict_qw2191_scope_delta_trendline_stability_audit.json"
OUT = GEN / "p2187_s1137_strict_qw2191_deterministic_ci_fail_threshold_policy.json"
MD = GEN / "p2187_s1137_strict_qw2191_deterministic_ci_fail_threshold_policy.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2186 = load(IN_2186)
    audit = p2186.get("strict_qw2191_scope_delta_trendline_stability_audit", {}) or {}
    trend_status = audit.get("trend_status", "OPEN_TRENDLINE_REGRESSION_INSTABILITY")
    trendline = audit.get("trendline", {}) or {}
    series = [int(x) for x in (trendline.get("regression_count_series", []) or [])]
    latest = series[-1] if series else 0

    policy_table = {
        "STABLE_ZERO_REGRESSION_TREND": {
            "ci_outcome": "PASS",
            "severity": "none",
            "deterministic_action": "allow_merge",
            "max_allowed_latest_regression_count": 0,
        },
        "NONINCREASING_REGRESSION_TREND": {
            "ci_outcome": "WARN",
            "severity": "medium",
            "deterministic_action": "allow_merge_with_notice",
            "max_allowed_latest_regression_count": 1,
        },
        "OPEN_TRENDLINE_REGRESSION_INSTABILITY": {
            "ci_outcome": "FAIL",
            "severity": "high",
            "deterministic_action": "block_merge",
            "max_allowed_latest_regression_count": 0,
        },
    }

    selected = policy_table.get(trend_status, policy_table["OPEN_TRENDLINE_REGRESSION_INSTABILITY"])

    threshold_respected = latest <= int(selected["max_allowed_latest_regression_count"])
    final_ci_decision = selected["ci_outcome"]
    if selected["ci_outcome"] in {"PASS", "WARN"} and not threshold_respected:
        final_ci_decision = "FAIL"

    result_kind = (
        "PASS_STRICT_QW2191_DETERMINISTIC_CI_FAIL_THRESHOLD_POLICY"
        if final_ci_decision in {"PASS", "WARN"}
        else "OPEN_STRICT_QW2191_CI_FAIL_THRESHOLD_POLICY_BLOCKED"
    )

    payload = {
        "schema_version": "p2187_s1137_v1",
        "packet_id": "P2187",
        "stage_id": "S1137",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_deterministic_ci_fail_threshold_policy": {
            "policy_id": "QW2191_DETERMINISTIC_CI_FAIL_THRESHOLD_POLICY_V1",
            "input_trendline_audit_packet": str(IN_2186.relative_to(ROOT)),
            "trend_status": trend_status,
            "latest_regression_count": latest,
            "policy_table": policy_table,
            "selected_policy": selected,
            "threshold_respected": threshold_respected,
            "final_ci_decision": final_ci_decision,
        },
        "recommended_next_honest_step": {
            "id": "P2188_candidate",
            "goal": "wire deterministic policy decision into a compact machine-checkable release gate contract",
        },
        "gatekeeper_checks": {
            "deterministic_policy_exported": True,
            "threshold_check_executed": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": bool((p2186.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2187 S1137: strict QW-2191 deterministic CI fail-threshold policy",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Trend status: `{trend_status}`",
                f"- Latest regression count: `{latest}`",
                f"- Final CI decision: `{final_ci_decision}`",
                "",
                "No selector closure or global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
