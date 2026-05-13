#!/usr/bin/env python3
"""P1396 checkpoint: export PC2 local obstruction and trigger PC3 design."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1395_strict_only_pc2_drift_epsilon_edge_robustness_run_summary.json"
OUT_PATH = ROOT / "generated" / "p1396_strict_only_pc2_drift_local_obstruction_export_and_pc3_design_trigger_summary.json"


def main() -> None:
    p1395 = json.loads(IN_PATH.read_text(encoding="utf-8"))

    robust_verdict = p1395.get("pc2_drift_edge_verdict", "UNKNOWN")
    obstruction_exported = robust_verdict == "ROBUST_FAIL"

    summary = {
        "packet_id": "P1396",
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_packet": "P1395",
        "input_robust_verdict": robust_verdict,
        "pc2_local_obstruction_id": "PC2-DRIFT-v1" if obstruction_exported else "NOT_APPLICABLE",
        "pc2_local_obstruction_status": "EXPORTED" if obstruction_exported else "NOT_EXPORTED",
        "pc2_loop_status": "CLOSED_NONCYCLIC" if obstruction_exported else "OPEN",
        "l_b1_03_export_status": "NOT_EXPORTED",
        "b1_status": "OPEN",
        "pc3_trigger": "ACTIVATED" if obstruction_exported else "NOT_ACTIVATED",
        "pc3_policy": "NEW_PROVIDER_CLASS_REQUIRED" if obstruction_exported else "UNSET",
        "next_packet": "P1397_STRICT_ONLY_PC3_PROVIDER_BASELINE_DESIGN",
        "no_false_pass": True,
    }

    OUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1396] wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
