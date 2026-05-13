#!/usr/bin/env python3
"""P1406 checkpoint: export PC4 drift obstruction and trigger PC5."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1405_strict_only_pc4_edge_robustness_sweep_summary.json"
OUT_PATH = ROOT / "generated" / "p1406_strict_only_pc4_drift_local_obstruction_export_and_pc5_trigger_summary.json"


def main() -> None:
    p1405 = json.loads(IN_PATH.read_text(encoding="utf-8"))
    robust_fail = p1405.get("pc4_edge_robustness_verdict") == "ROBUST_FAIL"

    summary = {
        "packet_id": "P1406",
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_packet": "P1405",
        "input_verdict": p1405.get("pc4_edge_robustness_verdict", "UNKNOWN"),
        "pc4_local_obstruction_id": "PC4-DRIFT-v1" if robust_fail else "NOT_APPLICABLE",
        "pc4_local_obstruction_status": "EXPORTED" if robust_fail else "NOT_EXPORTED",
        "pc4_loop_status": "CLOSED_NONCYCLIC" if robust_fail else "OPEN",
        "pc5_trigger": "ACTIVATED" if robust_fail else "NOT_ACTIVATED",
        "l_b1_03_export_status": "NOT_EXPORTED",
        "b1_status": "OPEN",
        "next_packet": "P1407_STRICT_ONLY_PC5_PROVIDER_BASELINE_DESIGN",
        "no_false_pass": True,
    }

    OUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1406] wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
