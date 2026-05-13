#!/usr/bin/env python3
"""P1438: strict-route precheck that requires P1435 PASS before new p14xx checkpoints."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"

    gate_path = gen / "p1435_pass_scope_enforcement_gate_summary.json"
    if not gate_path.exists():
        status = "FAIL_PRECHECK_MISSING_GATE_SUMMARY"
        gate_status = "MISSING"
    else:
        gate = json.loads(gate_path.read_text(encoding="utf-8"))
        gate_status = gate.get("status", "UNKNOWN")
        status = "PASS_PRECHECK_GATE_OPEN" if gate_status == "PASS_ENFORCEMENT_GATE" else "FAIL_PRECHECK_GATE_NOT_PASS"

    summary = {
        "packet": "P1438",
        "status": status,
        "route": "F-Nadsoliton=>L_SM+L_GR",
        "required_precheck": "p1435_pass_scope_enforcement_gate_summary.status == PASS_ENFORCEMENT_GATE",
        "observed_gate_status": gate_status,
        "legacy_import_used": False,
        "policy": "block_new_p14xx_checkpoint_publication_when_precheck_fails",
        "next_honest_step": "if PASS continue strict selector-source work; if FAIL remediate gate first",
        "scope_of_pass": "local_contract",
        "strict_core_qw2191_closed": False,
    }

    out = gen / "p1438_strict_route_precheck_gate_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[P1438] status={status} gate={gate_status}")


if __name__ == "__main__":
    main()
