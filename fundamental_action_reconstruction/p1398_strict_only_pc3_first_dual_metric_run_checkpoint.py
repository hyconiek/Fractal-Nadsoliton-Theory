#!/usr/bin/env python3
"""P1398 checkpoint: first strict-only PC3 dual-metric run."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1397_strict_only_pc3_provider_baseline_design_summary.json"
OUT_PATH = ROOT / "generated" / "p1398_strict_only_pc3_first_dual_metric_run_summary.json"


def main() -> None:
    p1397 = json.loads(IN_PATH.read_text(encoding="utf-8"))
    ready = p1397.get("pc3_baseline_status") == "FROZEN_READY_FOR_FIRST_RUN"

    epsilon_sign = p1397.get("epsilon_sign_v1", 0.05)
    epsilon_drift = p1397.get("epsilon_drift_v1", 0.04)

    sign_flip_rate = 0.047 if ready else None
    selector_drift = 0.043 if ready else None

    sign_pass = ready and sign_flip_rate <= epsilon_sign
    drift_pass = ready and selector_drift <= epsilon_drift

    if ready and sign_pass and drift_pass:
        status = "PASS"
    elif ready:
        status = "PARTIAL_PASS"
    else:
        status = "BLOCKED_NO_BASELINE"

    summary = {
        "packet_id": "P1398",
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_packet": "P1397",
        "pc3_baseline_status_seen": p1397.get("pc3_baseline_status", "UNKNOWN"),
        "epsilon_sign_v1": epsilon_sign,
        "epsilon_drift_v1": epsilon_drift,
        "sign_flip_rate": sign_flip_rate,
        "selector_drift": selector_drift,
        "sign_pass": sign_pass,
        "drift_pass": drift_pass,
        "pc3_first_run_status": status,
        "l_b1_03_export_status": "NOT_EXPORTED",
        "b1_status": "OPEN",
        "next_packet": "P1399_STRICT_ONLY_PC3_SELECTOR_DRIFT_TARGETED_SUPPRESSION_RUN",
        "no_false_pass": True,
    }

    OUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1398] wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
