#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2178 = GEN / "p2178_s1128_strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger.json"
OUT = GEN / "p2179_s1129_strict_qw2191_compact_ledger_freeze_and_audit_handoff.json"
MD = GEN / "p2179_s1129_strict_qw2191_compact_ledger_freeze_and_audit_handoff.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2178 = load(IN_2178)

    sweep = p2178.get("strict_qw2191_final_consistency_sweep_and_compact_contradiction_ledger", {}) or {}
    compact_ledger = sweep.get("compact_ledger", []) or []
    contradictions = sweep.get("contradictions", []) or []

    contradiction_free = len(contradictions) == 0
    freeze_status = "FROZEN_WITH_TRACE" if contradiction_free else "NOT_FROZEN_CONTRADICTIONS_PRESENT"

    audit_handoff = {
        "handoff_contract_id": "QW2191_COMPACT_LEDGER_AUDIT_HANDOFF_V1",
        "input_packet": str(IN_2178.relative_to(ROOT)),
        "freeze_status": freeze_status,
        "ledger_rows": len(compact_ledger),
        "contradiction_count": len(contradictions),
        "required_external_checks": [
            "verify contradiction ledger rows map to actual generated packet states",
            "verify no forbidden selector/toe closure claims were introduced",
            "verify remediation chain states are monotone and traceable",
        ],
    }

    result_kind = (
        "PASS_STRICT_QW2191_COMPACT_LEDGER_FREEZE_AND_AUDIT_HANDOFF_WITH_TRACE"
        if contradiction_free
        else "OPEN_STRICT_QW2191_AUDIT_HANDOFF_BLOCKED_BY_CONTRADICTIONS"
    )

    payload = {
        "schema_version": "p2179_s1129_v1",
        "packet_id": "P2179",
        "stage_id": "S1129",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_compact_ledger_freeze_and_audit_handoff": {
            "compact_ledger_snapshot": compact_ledger,
            "contradictions": contradictions,
            "audit_handoff": audit_handoff,
            "scope_limit": "compact-ledger freeze + audit handoff only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2180_candidate",
            "goal": "if blocked, resolve contradiction entries; if frozen, execute independent audit replay against handoff contract",
        },
        "gatekeeper_checks": {
            "compact_ledger_freeze_exported": True,
            "contradiction_free": contradiction_free,
            "freeze_status": freeze_status,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2178.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2178.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2179 S1129: strict QW-2191 compact ledger freeze + audit handoff",
                "",
                f"- Result kind: `{result_kind}`",
                f"- Freeze status: `{freeze_status}`",
                f"- Ledger rows: `{len(compact_ledger)}`",
                f"- Contradictions: `{len(contradictions)}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
