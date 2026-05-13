#!/usr/bin/env python3
"""P1465 S4.15: run P1463/P1464 after shared gate-core refactor and export consistency summary."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
SUMMARY = GEN / "p1465_s415_policy_gate_core_refactor_summary.json"

RUN = [
    ["python3", str(ROOT / "p1463_s413_policy_aware_rerun_gate_checkpoint.py")],
    ["python3", str(ROOT / "p1464_s414_gate_integration_dryrun_checkpoint.py")],
]


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    logs = []
    for cmd in RUN:
        proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logs.append({"cmd": " ".join(cmd), "stdout": proc.stdout.strip()})

    s1463 = readj(GEN / "p1463_s413_policy_aware_rerun_gate_summary.json")
    s1464 = readj(GEN / "p1464_s414_gate_integration_dryrun_summary.json")

    status_ok = (
        s1463.get("status") == "PASS_POLICY_GATE_ACTIVE"
        and s1464.get("status") == "PASS_GATE_INTEGRATION_DRYRUN"
        and int(s1463.get("blocked_count", 0)) >= 1
        and int(s1464.get("blocked_count", 0)) >= 1
    )

    summary = {
        "packet": "P1465",
        "status": "PASS_POLICY_CORE_REFACTOR_LOCAL_ONLY" if status_ok else "FAIL_POLICY_CORE_REFACTOR",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "checks": {
            "p1463_status": s1463.get("status"),
            "p1463_blocked_count": s1463.get("blocked_count"),
            "p1464_status": s1464.get("status"),
            "p1464_blocked_count": s1464.get("blocked_count"),
        },
        "rerun_logs": logs,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1465] status={summary['status']}")


if __name__ == "__main__":
    main()
