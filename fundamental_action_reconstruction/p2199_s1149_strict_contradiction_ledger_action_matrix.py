#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2198 = GEN / "p2198_s1148_strict_release_gate_policy_contract_and_contradiction_delta_bridge.json"
OUT = GEN / "p2199_s1149_strict_contradiction_ledger_action_matrix.json"
MD = GEN / "p2199_s1149_strict_contradiction_ledger_action_matrix.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2198 = load(IN_2198)
    bridge = p2198.get("strict_release_gate_policy_contract_and_contradiction_delta_bridge", {}) or {}

    mismatch_count = int(bridge.get("release_gate_mismatch_count", 0) or 0)
    all_scope_checked = bool(bridge.get("all_scope_delta_checked", False))
    any_scope_blocking = bool(bridge.get("any_scope_delta_blocking", False))

    matrix = {
        "A0_clean_open": {
            "condition": "mismatch_count==0 and all_scope_delta_checked and not any_scope_blocking",
            "action": "keep_release_gate_state_and_continue_periodic_replay",
            "severity": "none",
        },
        "A1_scope_coverage_gap": {
            "condition": "mismatch_count==0 and not all_scope_delta_checked",
            "action": "force_scope_delta_recheck_and_hold_promotion",
            "severity": "medium",
        },
        "A2_scope_blocking_present": {
            "condition": "mismatch_count==0 and any_scope_blocking",
            "action": "open_contradiction_ticket_and_pin_warn",
            "severity": "high",
        },
        "A3_policy_contract_mismatch": {
            "condition": "mismatch_count>0",
            "action": "block_release_gate_and_replay_policy_contract_consistency",
            "severity": "critical",
        },
    }

    if mismatch_count > 0:
        action_id = "A3_policy_contract_mismatch"
    elif not all_scope_checked:
        action_id = "A1_scope_coverage_gap"
    elif any_scope_blocking:
        action_id = "A2_scope_blocking_present"
    else:
        action_id = "A0_clean_open"

    selected = matrix[action_id]

    payload = {
        "schema_version": "p2199_s1149_v1",
        "packet_id": "P2199",
        "stage_id": "S1149",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": "PASS_STRICT_CONTRADICTION_LEDGER_ACTION_MATRIX",
        "strict_contradiction_ledger_action_matrix": {
            "matrix_id": "STRICT_CONTRADICTION_LEDGER_ACTION_MATRIX_V1",
            "source_packet": str(IN_2198.relative_to(ROOT)),
            "inputs": {
                "release_gate_mismatch_count": mismatch_count,
                "all_scope_delta_checked": all_scope_checked,
                "any_scope_delta_blocking": any_scope_blocking,
            },
            "action_matrix": matrix,
            "selected_action_id": action_id,
            "selected_action": selected,
        },
        "recommended_next_honest_step": {
            "id": "P2200_candidate",
            "goal": "bind selected ledger action to deterministic issue-ticket schema and replay scheduling contract",
        },
        "gatekeeper_checks": {
            "action_matrix_exported": True,
            "selected_action_resolved": True,
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
                "# P2199 S1149: strict contradiction-ledger action matrix",
                "",
                f"- Selected action id: `{action_id}`",
                f"- Action: `{selected['action']}`",
                f"- Severity: `{selected['severity']}`",
                "",
                "This is an action-matrix governance export only, not full Cutkosky closure.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
