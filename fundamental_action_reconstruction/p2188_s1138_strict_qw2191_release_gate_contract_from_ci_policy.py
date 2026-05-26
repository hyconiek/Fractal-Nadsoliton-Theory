#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2187 = GEN / "p2187_s1137_strict_qw2191_deterministic_ci_fail_threshold_policy.json"
OUT = GEN / "p2188_s1138_strict_qw2191_release_gate_contract_from_ci_policy.json"
MD = GEN / "p2188_s1138_strict_qw2191_release_gate_contract_from_ci_policy.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_release_decision(ci_decision: str) -> str:
    table = {
        "PASS": "RELEASE_GATE_OPEN",
        "WARN": "RELEASE_GATE_OPEN_WITH_NOTICE",
        "FAIL": "RELEASE_GATE_BLOCKED",
    }
    return table.get(ci_decision, "RELEASE_GATE_BLOCKED")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2187 = load(IN_2187)
    policy = p2187.get("strict_qw2191_deterministic_ci_fail_threshold_policy", {}) or {}

    ci_decision = str(policy.get("final_ci_decision", "FAIL"))
    trend_status = str(policy.get("trend_status", "OPEN_TRENDLINE_REGRESSION_INSTABILITY"))
    latest_regression_count = int(policy.get("latest_regression_count", 0) or 0)
    gate_decision = normalize_release_decision(ci_decision)

    release_gate_contract = {
        "contract_id": "QW2191_RELEASE_GATE_CONTRACT_V1",
        "contract_kind": "machine_checkable_release_gate",
        "input_policy_packet": str(IN_2187.relative_to(ROOT)),
        "inputs": {
            "final_ci_decision": ci_decision,
            "trend_status": trend_status,
            "latest_regression_count": latest_regression_count,
        },
        "deterministic_mapping": {
            "PASS": "RELEASE_GATE_OPEN",
            "WARN": "RELEASE_GATE_OPEN_WITH_NOTICE",
            "FAIL": "RELEASE_GATE_BLOCKED",
        },
        "output": {
            "release_gate_decision": gate_decision,
            "release_gate_open": gate_decision in {"RELEASE_GATE_OPEN", "RELEASE_GATE_OPEN_WITH_NOTICE"},
            "requires_notice": gate_decision == "RELEASE_GATE_OPEN_WITH_NOTICE",
            "release_blocked": gate_decision == "RELEASE_GATE_BLOCKED",
        },
        "scope_policy": {
            "lane": "strict_qw2191_external_replay",
            "scope_tag": "nonbridge_operational_strict_lane",
            "selector_closure_claimed": False,
            "toe_closure_claimed": False,
        },
    }

    result_kind = (
        "PASS_STRICT_QW2191_RELEASE_GATE_CONTRACT_FROM_CI_POLICY"
        if release_gate_contract["output"]["release_gate_open"]
        else "OPEN_STRICT_QW2191_RELEASE_GATE_CONTRACT_BLOCKED"
    )

    payload = {
        "schema_version": "p2188_s1138_v1",
        "packet_id": "P2188",
        "stage_id": "S1138",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_release_gate_contract_from_ci_policy": release_gate_contract,
        "recommended_next_honest_step": {
            "id": "P2189_candidate",
            "goal": "add policy-contract replay verifier against historical packets to check contract invariance across reruns",
        },
        "gatekeeper_checks": {
            "release_gate_contract_exported": True,
            "deterministic_mapping_explicit": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_cutkosky_closure_proven": False,
            "full_d3_covariance_transport_proven": bool((p2187.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2188 S1138: strict QW-2191 release-gate contract from deterministic CI policy",
                "",
                f"- Result kind: `{result_kind}`",
                f"- CI decision input: `{ci_decision}`",
                f"- Release gate decision: `{gate_decision}`",
                f"- Trend status input: `{trend_status}`",
                "",
                "No selector closure or global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
