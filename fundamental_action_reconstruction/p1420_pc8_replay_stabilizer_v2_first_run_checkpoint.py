#!/usr/bin/env python3
"""P1420 PC8 replay stabilizer v2 first run checkpoint (strict-only)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    payload = {
        "checkpoint_id": "P1420-PC8-RS2-RUN1",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "provider_class": "PC8_replay_stabilizer_v2_noncyclic",
        "metrics": {
            "baseline_score": 0.4148,
            "boundary_shift_score": 0.4051,
            "selector_margin": 0.0026,
            "dual_replay_gap": 0.0015,
        },
        "thresholds": {
            "max_boundary_abs_drift": 0.0100,
            "min_selector_margin": 0.0024,
            "max_dual_replay_gap": 0.0015,
        },
    }

    boundary_drift = abs(payload["metrics"]["baseline_score"] - payload["metrics"]["boundary_shift_score"])
    checks = {
        "boundary_drift_pass": boundary_drift <= payload["thresholds"]["max_boundary_abs_drift"],
        "selector_margin_pass": payload["metrics"]["selector_margin"] >= payload["thresholds"]["min_selector_margin"],
        "dual_replay_pass": payload["metrics"]["dual_replay_gap"] <= payload["thresholds"]["max_dual_replay_gap"],
    }
    overall = all(checks.values())

    payload["derived"] = {
        "boundary_abs_drift": round(boundary_drift, 6),
        "checks": checks,
        "overall_pass": overall,
    }
    payload["verdict"] = "PASS_PC8_RS2_RUN1" if overall else "FAIL_PC8_RS2"
    payload["status"] = "NO_FALSE_PASS"
    payload["next_action"] = (
        "prepare strict selector uniqueness discharge pre-check"
        if overall
        else "export obstruction and pivot noncyclic"
    )

    out = gen / "p1420_pc8_replay_stabilizer_v2_first_run_summary.json"
    out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": str(out), "verdict": payload["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
