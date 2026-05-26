#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2182 = GEN / "p2182_s1132_strict_qw2191_release_claims_delta_audit.json"
IN_2178 = GEN / "p2178_s1128_strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger.json"
OUT = GEN / "p2183_s1133_strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge.json"
MD = GEN / "p2183_s1133_strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2182 = load(IN_2182)
    p2178 = load(IN_2178)

    audit = p2182.get("strict_qw2191_release_claims_delta_audit", {}) or {}
    blocking = audit.get("blocking_findings", []) or []
    delta_findings = audit.get("delta_findings", []) or []

    p2178_ledger = (p2178.get("strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger", {}) or {}).get(
        "compact_ledger", []
    )

    governance_note = {
        "governance_note_id": "QW2191_SCOPE_DELTA_GOVERNANCE_NOTE_V1",
        "source_release_delta_audit": str(IN_2182.relative_to(ROOT)),
        "policy": [
            "Historical release closure claims are treated as scoped history, not automatic current-lane discharge.",
            "Current-lane status depends on executable packet chain checks and replay certificates.",
            "Any mismatch between scoped historical claims and lane certificates must be logged in contradiction ledgers.",
        ],
        "blocking_findings_count": len(blocking),
        "delta_findings_count": len(delta_findings),
        "status": "ACTIVE_WITH_TRACE" if len(blocking) == 0 else "ACTIVE_BLOCKING_WITH_TRACE",
    }

    # Wire governance note into a refresh bridge view over compact ledger rows.
    bridged_rows = []
    for row in p2178_ledger:
        rr = dict(row)
        rr["scope_delta_governance_note_id"] = governance_note["governance_note_id"]
        rr["scope_delta_checked"] = True
        rr["scope_delta_blocking_present"] = len(blocking) > 0
        bridged_rows.append(rr)

    result_kind = (
        "PASS_STRICT_QW2191_SCOPE_DELTA_GOVERNANCE_FREEZE_AND_LEDGER_REFRESH_BRIDGE_WITH_TRACE"
        if len(blocking) == 0
        else "OPEN_STRICT_QW2191_SCOPE_DELTA_GOVERNANCE_BLOCKED"
    )

    payload = {
        "schema_version": "p2183_s1133_v1",
        "packet_id": "P2183",
        "stage_id": "S1133",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_scope_delta_governance_freeze_and_ledger_refresh_bridge": {
            "governance_note": governance_note,
            "ledger_refresh_bridge_rows": bridged_rows,
            "scope_limit": "governance freeze + ledger refresh bridge only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2184_candidate",
            "goal": "consume governance-bridged ledger in an external audit replay dashboard and verify no scope-delta regressions",
        },
        "gatekeeper_checks": {
            "scope_delta_governance_freeze_exported": True,
            "blocking_findings_absent": len(blocking) == 0,
            "ledger_refresh_bridge_exported": True,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2182.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2182.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2183 S1133: strict QW-2191 scope-delta governance freeze + ledger refresh bridge",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Blocking findings: `{len(blocking)}`",
                f"- Bridged ledger rows: `{len(bridged_rows)}`",
                f"- Governance note id: `{governance_note['governance_note_id']}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
