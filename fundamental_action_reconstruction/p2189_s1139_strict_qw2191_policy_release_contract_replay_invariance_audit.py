#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2187 = GEN / "p2187_s1137_strict_qw2191_deterministic_ci_fail_threshold_policy.json"
IN_2188 = GEN / "p2188_s1138_strict_qw2191_release_gate_contract_from_ci_policy.json"
OUT = GEN / "p2189_s1139_strict_qw2191_policy_release_contract_replay_invariance_audit.json"
MD = GEN / "p2189_s1139_strict_qw2191_policy_release_contract_replay_invariance_audit.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def expected_gate(ci_decision: str) -> str:
    table = {
        "PASS": "RELEASE_GATE_OPEN",
        "WARN": "RELEASE_GATE_OPEN_WITH_NOTICE",
        "FAIL": "RELEASE_GATE_BLOCKED",
    }
    return table.get(ci_decision, "RELEASE_GATE_BLOCKED")


def main() -> None:
    GEN.mkdir(exist_ok=True)

    p2187 = load(IN_2187)
    p2188 = load(IN_2188)

    policy = p2187.get("strict_qw2191_deterministic_ci_fail_threshold_policy", {}) or {}
    contract = p2188.get("strict_qw2191_release_gate_contract_from_ci_policy", {}) or {}

    ci_decision = str(policy.get("final_ci_decision", "FAIL"))
    trend_status = str(policy.get("trend_status", "OPEN_TRENDLINE_REGRESSION_INSTABILITY"))
    latest_regression_count = int(policy.get("latest_regression_count", 0) or 0)

    expected = {
        "release_gate_decision": expected_gate(ci_decision),
        "input_final_ci_decision": ci_decision,
        "input_trend_status": trend_status,
        "input_latest_regression_count": latest_regression_count,
    }

    actual_inputs = contract.get("inputs", {}) or {}
    actual_output = contract.get("output", {}) or {}

    checks = [
        {
            "id": "contract_input_ci_decision_matches_policy",
            "expected": ci_decision,
            "actual": actual_inputs.get("final_ci_decision"),
            "pass": actual_inputs.get("final_ci_decision") == ci_decision,
        },
        {
            "id": "contract_input_trend_status_matches_policy",
            "expected": trend_status,
            "actual": actual_inputs.get("trend_status"),
            "pass": actual_inputs.get("trend_status") == trend_status,
        },
        {
            "id": "contract_input_latest_regression_matches_policy",
            "expected": latest_regression_count,
            "actual": actual_inputs.get("latest_regression_count"),
            "pass": actual_inputs.get("latest_regression_count") == latest_regression_count,
        },
        {
            "id": "contract_output_release_gate_decision_is_deterministic",
            "expected": expected["release_gate_decision"],
            "actual": actual_output.get("release_gate_decision"),
            "pass": actual_output.get("release_gate_decision") == expected["release_gate_decision"],
        },
        {
            "id": "contract_scope_tag_preserved",
            "expected": "nonbridge_operational_strict_lane",
            "actual": (contract.get("scope_policy", {}) or {}).get("scope_tag"),
            "pass": (contract.get("scope_policy", {}) or {}).get("scope_tag") == "nonbridge_operational_strict_lane",
        },
    ]

    all_pass = all(c["pass"] for c in checks)
    result_kind = (
        "PASS_STRICT_QW2191_POLICY_RELEASE_CONTRACT_REPLAY_INVARIANCE_AUDIT"
        if all_pass
        else "OPEN_STRICT_QW2191_POLICY_RELEASE_CONTRACT_REPLAY_INVARIANCE_MISMATCH"
    )

    payload = {
        "schema_version": "p2189_s1139_v1",
        "packet_id": "P2189",
        "stage_id": "S1139",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_policy_release_contract_replay_invariance_audit": {
            "audit_id": "QW2191_POLICY_RELEASE_CONTRACT_REPLAY_INVARIANCE_AUDIT_V1",
            "input_packets": [str(IN_2187.relative_to(ROOT)), str(IN_2188.relative_to(ROOT))],
            "expected_contract_slice": expected,
            "checks": checks,
            "all_checks_pass": all_pass,
        },
        "recommended_next_honest_step": {
            "id": "P2190_candidate",
            "goal": "export first strict-lane real-valued discontinuity integral packet on one fixed graviton->gauge_gauge channel/background",
        },
        "gatekeeper_checks": {
            "replay_invariance_audit_exported": True,
            "all_invariance_checks_pass": all_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": bool((p2188.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2189 S1139: policy->release contract replay invariance audit",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Invariance checks: `{len(checks)}`",
                f"- All checks pass: `{all_pass}`",
                "",
                "No selector closure or global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
