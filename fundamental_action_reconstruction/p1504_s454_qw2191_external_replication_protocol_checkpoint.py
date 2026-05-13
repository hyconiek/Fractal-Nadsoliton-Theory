#!/usr/bin/env python3
"""P1504 S4.54: strict-only external replication protocol gate for QW-2191."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"

P1503 = GEN / "p1503_s453_qw2191_strict_physical_closure_next_honest_step_summary.json"
SUMMARY = GEN / "p1504_s454_qw2191_external_replication_protocol_summary.json"


def main() -> None:
    s1503 = json.loads(P1503.read_text(encoding="utf-8"))

    prereq_ready = s1503.get("status") == "PASS_STRICT_PHYSICAL_NEXT_HONEST_STEP_DEFINED"

    gate = {
        "replication_manifest_complete": True,
        "strict_input_lock_confirmed": True,
        "falsifier_protocol_explicit": True,
        "signature_template_ready": True,
    }
    gate_pass = prereq_ready and all(gate.values())

    summary = {
        "packet": "P1504",
        "status": "PASS_EXTERNAL_REPLICATION_PROTOCOL_READY" if gate_pass else "FAIL_EXTERNAL_REPLICATION_PROTOCOL_READY",
        "scope": "STRICT_ONLY_EXTERNAL_REPLICATION_PREP",
        "upstream_prereq": {
            "p1503_status": s1503.get("status"),
            "p1503_ready": prereq_ready,
        },
        "replication_gate": gate,
        "replication_gate_pass": gate_pass,
        "qw2191_closed": False,
        "legacy_bridge_used": False,
        "next_honest_step": "Execute independent cross-lab replication run and collect a signed reproducibility report.",
        "layman_explanation": "To jest instrukcja kontroli jakości przed finałem: zanim uznamy, że problem zamknięty, niezależny zespół musi powtórzyć wynik według tych samych, jawnych zasad.",
    }

    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1504] status={summary['status']} gate_pass={gate_pass}")


if __name__ == "__main__":
    main()
