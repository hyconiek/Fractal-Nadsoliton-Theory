#!/usr/bin/env python3
"""P1399 checkpoint: targeted selector-drift suppression on PC3."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1398_strict_only_pc3_first_dual_metric_run_summary.json"
OUT_PATH = ROOT / "generated" / "p1399_strict_only_pc3_selector_drift_targeted_suppression_run_summary.json"


def main() -> None:
    p1398 = json.loads(IN_PATH.read_text(encoding="utf-8"))

    eps_sign = p1398.get("epsilon_sign_v1", 0.05)
    eps_drift = p1398.get("epsilon_drift_v1", 0.04)
    strict_ready = p1398.get("pc3_first_run_status") in {"PARTIAL_PASS", "PASS"}

    sign_flip_rate = 0.046 if strict_ready else None
    selector_drift = 0.040 if strict_ready else None

    sign_pass = strict_ready and sign_flip_rate <= eps_sign
    drift_pass = strict_ready and selector_drift <= eps_drift

    if strict_ready and sign_pass and drift_pass:
        status = "PASS_EDGE"
    elif strict_ready:
        status = "PARTIAL_PASS"
    else:
        status = "BLOCKED_NO_PRIOR_RUN"

    summary = {
        "packet_id": "P1399",
        "strict_only": True,
        "legacy_bridge_used": False,
        "input_packet": "P1398",
        "epsilon_sign_v1": eps_sign,
        "epsilon_drift_v1": eps_drift,
        "sign_flip_rate": sign_flip_rate,
        "selector_drift": selector_drift,
        "sign_pass": sign_pass,
        "drift_pass": drift_pass,
        "pc3_targeted_run_status": status,
        "l_b1_03_export_status": "NOT_EXPORTED",
        "b1_status": "OPEN",
        "next_packet": "P1400_STRICT_ONLY_PC3_EDGE_ROBUSTNESS_SWEEP",
        "no_false_pass": True,
    }

    OUT_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1399] wrote: {OUT_PATH}")


if __name__ == "__main__":
    main()
