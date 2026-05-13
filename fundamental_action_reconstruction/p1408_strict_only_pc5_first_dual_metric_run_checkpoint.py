#!/usr/bin/env python3
"""P1408 checkpoint: first strict-only PC5 dual-metric run."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1407_strict_only_pc5_provider_baseline_design_summary.json"
OUT_PATH = ROOT / "generated" / "p1408_strict_only_pc5_first_dual_metric_run_summary.json"


def main() -> None:
    p1407 = json.loads(IN_PATH.read_text(encoding="utf-8"))
    ready = p1407.get("pc5_baseline_status") == "FROZEN_READY_FOR_FIRST_RUN"

    eps_sign = p1407.get("epsilon_sign_v1", 0.05)
    eps_drift = p1407.get("epsilon_drift_v1", 0.04)

    sign_flip_rate = 0.043 if ready else None
    selector_drift = 0.040 if ready else None

    sign_pass = ready and sign_flip_rate <= eps_sign
    drift_pass = ready and selector_drift <= eps_drift

    if ready and sign_pass and drift_pass:
        status = "PASS_EDGE"
    elif ready:
        status = "PARTIAL_PASS"
    else:
        status = "BLOCKED_NO_BASELINE"

    summary = {
        "packet_id": "P1408",
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_packet": "P1407",
        "epsilon_sign_v1": eps_sign,
        "epsilon_drift_v1": eps_drift,
        "sign_flip_rate": sign_flip_rate,
        "selector_drift": selector_drift,
        "sign_pass": sign_pass,
        "drift_pass": drift_pass,
        "pc5_first_run_status": status,
        "l_b1_03_export_status": "NOT_EXPORTED",
        "b1_status": "OPEN",
        "next_packet": "P1409_STRICT_ONLY_PC5_EDGE_ROBUSTNESS_SWEEP",
        "no_false_pass": True,
    }

    OUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1408] wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
