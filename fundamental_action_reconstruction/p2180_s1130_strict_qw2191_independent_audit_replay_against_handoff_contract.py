#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_2179 = GEN / "p2179_s1129_strict_qw2191_compact_ledger_freeze_and_audit_handoff.json"
OUT = GEN / "p2180_s1130_strict_qw2191_independent_audit_replay_against_handoff_contract.json"
MD = GEN / "p2180_s1130_strict_qw2191_independent_audit_replay_against_handoff_contract.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "result_kind": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GEN.mkdir(exist_ok=True)
    p2179 = load(IN_2179)

    pack = p2179.get("strict_qw2191_compact_ledger_freeze_and_audit_handoff", {}) or {}
    handoff = pack.get("audit_handoff", {}) or {}
    ledger = pack.get("compact_ledger_snapshot", []) or []

    expected_contract = "QW2191_COMPACT_LEDGER_AUDIT_HANDOFF_V1"
    contract_ok = handoff.get("handoff_contract_id") == expected_contract

    # Replay checks (independent from previous packet logic, only against handoff contract data)
    check_rows = len(ledger) == int(handoff.get("ledger_rows", -1))
    check_contradictions = len(pack.get("contradictions", []) or []) == int(handoff.get("contradiction_count", -1))
    check_guardrails = all(bool(row.get("no_selector_closure_claimed", False)) and bool(row.get("no_toe_closure_claimed", False)) for row in ledger)

    checks = {
        "contract_id_matches": contract_ok,
        "ledger_row_count_matches": check_rows,
        "contradiction_count_matches": check_contradictions,
        "guardrail_rows_hold": check_guardrails,
    }
    all_checks_pass = all(checks.values())

    result_kind = (
        "PASS_STRICT_QW2191_INDEPENDENT_AUDIT_REPLAY_WITH_TRACE"
        if all_checks_pass
        else "OPEN_STRICT_QW2191_INDEPENDENT_AUDIT_REPLAY_MISMATCH"
    )

    payload = {
        "schema_version": "p2180_s1130_v1",
        "packet_id": "P2180",
        "stage_id": "S1130",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-26T00:00:00+00:00",
        "status": "OPEN_PARTIAL_PROGRESS_WITH_TRACE",
        "result_kind": result_kind,
        "strict_qw2191_independent_audit_replay_against_handoff_contract": {
            "source_handoff_packet": str(IN_2179.relative_to(ROOT)),
            "expected_contract_id": expected_contract,
            "checks": checks,
            "all_checks_pass": all_checks_pass,
            "scope_limit": "independent replay against handoff contract only; no selector closure claim",
        },
        "recommended_next_honest_step": {
            "id": "P2181_candidate",
            "goal": "if mismatch appears, patch handoff contract producer and rerun replay; if pass, freeze replay certificate",
        },
        "gatekeeper_checks": {
            "independent_audit_replay_exported": True,
            "all_checks_pass": all_checks_pass,
            "no_selector_closure_claimed": True,
            "no_toe_closure_claimed": True,
            "full_d3_covariance_transport_proven": bool((p2179.get("gatekeeper_checks", {}) or {}).get("full_d3_covariance_transport_proven", False)),
            "full_cutkosky_closure_proven": False,
            "c3_theorem_proven": bool((p2179.get("gatekeeper_checks", {}) or {}).get("c3_theorem_proven", False)),
        },
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    MD.write_text(
        "\n".join(
            [
                "# P2180 S1130: strict QW-2191 independent audit replay against handoff contract",
                "",
                f"- Result kind: `{result_kind}`",
                f"- all checks pass: `{all_checks_pass}`",
                f"- checks: `{checks}`",
                "",
                "No global ToE closure claim is made.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
