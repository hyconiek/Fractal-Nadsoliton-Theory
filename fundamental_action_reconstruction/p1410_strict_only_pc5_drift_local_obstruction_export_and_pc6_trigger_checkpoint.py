#!/usr/bin/env python3
"""P1410 checkpoint: export PC5 drift obstruction and trigger PC6."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1409_strict_only_pc5_edge_robustness_sweep_summary.json"
OUT_PATH = ROOT / "generated" / "p1410_strict_only_pc5_drift_local_obstruction_export_and_pc6_trigger_summary.json"


def main() -> None:
    p1409 = json.loads(IN_PATH.read_text(encoding="utf-8"))
    robust_fail = p1409.get("pc5_edge_robustness_verdict") == "ROBUST_FAIL"

    summary = {
        "packet_id": "P1410",
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_packet": "P1409",
        "input_verdict": p1409.get("pc5_edge_robustness_verdict", "UNKNOWN"),
        "pc5_local_obstruction_id": "PC5-DRIFT-v1" if robust_fail else "NOT_APPLICABLE",
        "pc5_local_obstruction_status": "EXPORTED" if robust_fail else "NOT_EXPORTED",
        "pc5_loop_status": "CLOSED_NONCYCLIC" if robust_fail else "OPEN",
        "pc6_trigger": "ACTIVATED" if robust_fail else "NOT_ACTIVATED",
        "l_b1_03_export_status": "NOT_EXPORTED",
        "b1_status": "OPEN",
        "next_packet": "P1411_STRICT_ONLY_PC6_PROVIDER_BASELINE_DESIGN",
        "no_false_pass": True,
    }

    OUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1410] wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
