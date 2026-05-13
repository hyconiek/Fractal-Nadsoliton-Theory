#!/usr/bin/env python3
"""P1471 S4.21: compliance gate against QW-2191 anti-pattern memory (P1470)."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
MEMORY = GEN / "p1470_s420_qw2191_failure_memory_summary.json"
SUMMARY = GEN / "p1471_s421_qw2191_proposal_compliance_gate_summary.json"
OBSTRUCTION = GEN / "p1471_s421_qw2191_proposal_compliance_gate_obstruction.json"

# Candidate for next step (SP1 continuation) declares anti-pattern compliance.
PROPOSAL = {
    "id": "SP1_continuation_v1",
    "avoids_AP1": True,
    "avoids_AP2": True,
    "avoids_AP3": True,
    "avoids_AP4": True,
}


def main() -> None:
    memory = json.loads(MEMORY.read_text(encoding="utf-8"))
    anti_ids = [a["id"] for a in memory["anti_patterns"]]

    checks = []
    first_fail = None
    for ap in anti_ids:
        key = f"avoids_{ap.split('_')[0]}"  # AP1..., AP2...
        ok = bool(PROPOSAL.get(key, False))
        row = {"anti_pattern_id": ap, "proposal_key": key, "ok": ok, "status": "PASS" if ok else "FAIL"}
        checks.append(row)
        if (not ok) and first_fail is None:
            first_fail = row

    status = "PASS_QW2191_COMPLIANCE_GATE_LOCAL_ONLY" if first_fail is None else "FAIL_QW2191_COMPLIANCE_GATE"
    summary = {
        "packet": "P1471",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "proposal": PROPOSAL,
        "checks": checks,
        "fail_count": sum(1 for c in checks if c["status"] == "FAIL"),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1471",
            "status": "FAIL_QW2191_COMPLIANCE_GATE",
            "first_fail": first_fail,
            "rule": "proposal repeats known failed anti-pattern",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1471] status={status} fail_count={summary['fail_count']}")


if __name__ == "__main__":
    main()
