#!/usr/bin/env python3
"""P1397 checkpoint: freeze strict-only PC3 baseline after P1396 trigger."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1396_strict_only_pc2_drift_local_obstruction_export_and_pc3_design_trigger_summary.json"
OUT_PATH = ROOT / "generated" / "p1397_strict_only_pc3_provider_baseline_design_summary.json"


def main() -> None:
    p1396 = json.loads(IN_PATH.read_text(encoding="utf-8"))
    triggered = p1396.get("pc3_trigger") == "ACTIVATED"

    summary = {
        "packet_id": "P1397",
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_packet": "P1396",
        "pc3_trigger_seen": p1396.get("pc3_trigger", "UNKNOWN"),
        "provider_class_id": "PC3_strict_selector_anchor_v1" if triggered else "UNSET",
        "anchor_family": "A3_noncyclic_selector_stabilizers" if triggered else "UNSET",
        "inherits_from_pc2_loop": False,
        "epsilon_sign_v1": 0.05 if triggered else None,
        "epsilon_drift_v1": 0.04 if triggered else None,
        "pc3_baseline_status": "FROZEN_READY_FOR_FIRST_RUN" if triggered else "BLOCKED_NO_TRIGGER",
        "l_b1_03_export_status": "NOT_EXPORTED",
        "b1_status": "OPEN",
        "next_packet": "P1398_STRICT_ONLY_PC3_FIRST_DUAL_METRIC_RUN",
        "no_false_pass": True,
    }

    OUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1397] wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
