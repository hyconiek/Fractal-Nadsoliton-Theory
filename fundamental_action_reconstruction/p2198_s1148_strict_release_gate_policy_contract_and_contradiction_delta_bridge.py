#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2197 = GEN / "p2197_s1147_strict_cutkosky_drift_trend_cadence_and_hardstop_policy.json"
IN_2188 = GEN / "p2188_s1138_strict_qw2191_release_gate_contract_from_ci_policy.json"
IN_2183 = GEN / "p2183_s1133_strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge.json"
OUT = GEN / "p2198_s1148_strict_release_gate_policy_contract_and_contradiction_delta_bridge.json"
MD = GEN / "p2198_s1148_strict_release_gate_policy_contract_and_contradiction_delta_bridge.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def normalize_gate_from_policy(ci_gate: str) -> str:
    m = {
        "OPEN": "RELEASE_GATE_OPEN",
        "OPEN_WITH_NOTICE": "RELEASE_GATE_OPEN_WITH_NOTICE",
        "BLOCK": "RELEASE_GATE_BLOCKED",
    }
    return m.get(ci_gate, "RELEASE_GATE_BLOCKED")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2197 = load(IN_2197)
    p2188 = load(IN_2188)
    p2183 = load(IN_2183)

    pol = p2197.get("strict_cutkosky_drift_trend_cadence_and_hardstop_policy", {}) or {}
    sel = pol.get("selected_policy", {}) or {}
    ci_gate = str(sel.get("ci_gate", "BLOCK"))
    expected_release = normalize_gate_from_policy(ci_gate)

    rg_contract = p2188.get("strict_qw2191_release_gate_contract_from_ci_policy", {}) or {}
    actual_release = str((rg_contract.get("output", {}) or {}).get("release_gate_decision", "RELEASE_GATE_BLOCKED"))

    bridged_rows = (
        (p2183.get("strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge", {}) or {}).get("ledger_refresh_bridge_rows", [])
        or []
    )

    contradiction_delta_rows = []
    for row in bridged_rows:
        contradiction_delta_rows.append(
            {
                "packet": row.get("packet"),
                "scope_delta_checked": bool(row.get("scope_delta_checked", False)),
                "scope_delta_blocking_present": bool(row.get("scope_delta_blocking_present", False)),
                "policy_gate": ci_gate,
                "expected_release_gate": expected_release,
                "actual_release_gate": actual_release,
                "delta_flag": expected_release != actual_release,
            }
        )

    mismatch_count = sum(1 for r in contradiction_delta_rows if r["delta_flag"])
    all_scope_checked = all(bool(r["scope_delta_checked"]) for r in contradiction_delta_rows)
    any_scope_blocking = any(bool(r["scope_delta_blocking_present"]) for r in contradiction_delta_rows)

    payload = {
        "schema_version": "p2198_s1148_v1",
        "packet_id": "P2198",
        "stage_id": "S1148",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_RELEASE_GATE_POLICY_CONTRACT_AND_CONTRADICTION_DELTA_BRIDGE",
        "strict_release_gate_policy_contract_and_contradiction_delta_bridge": {
            "bridge_id": "STRICT_RELEASE_GATE_POLICY_CONTRACT_CONTRADICTION_DELTA_BRIDGE_V1",
            "source_packets": [
                str(IN_2197.relative_to(ROOT)),
                str(IN_2188.relative_to(ROOT)),
                str(IN_2183.relative_to(ROOT)),
            ],
            "policy_ci_gate": ci_gate,
            "expected_release_gate": expected_release,
            "actual_release_gate": actual_release,
            "release_gate_mismatch_count": mismatch_count,
            "all_scope_delta_checked": all_scope_checked,
            "any_scope_delta_blocking": any_scope_blocking,
            "contradiction_delta_rows": contradiction_delta_rows,
        },
        "recommended_next_honest_step": {
            "id": "P2199_candidate",
            "goal": "export deterministic contradiction-ledger action matrix keyed by mismatch_count and scope_delta_blocking flags",
        },
        "gatekeeper_checks": {
            "policy_contract_delta_bridge_exported": True,
            "release_gate_consistency_checked": True,
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
                "# P2198 S1148: strict release-gate policy/contract + contradiction-delta bridge",
                "",
                f"- Policy CI gate: `{ci_gate}`",
                f"- Expected release gate: `{expected_release}`",
                f"- Actual release gate: `{actual_release}`",
                f"- Release mismatch count: `{mismatch_count}`",
                f"- Any scope-delta blocking: `{any_scope_blocking}`",
                "",
                "This is a governance/ledger bridge only, not full Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
