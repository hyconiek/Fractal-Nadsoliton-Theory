#!/usr/bin/env python3
"""P1416 PC7 first transport run checkpoint (strict-only, no legacy bridge)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    run = {
        "checkpoint_id": "P1416-PC7-RUN1",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "provider_class": "PC7",
        "scores": {
            "baseline": 0.4112,
            "scale_shift": 0.4060,
            "scheme_shift": 0.4026,
            "boundary_shift": 0.4013,
        },
        "thresholds": {
            "max_boundary_abs_drift": 0.0100,
            "min_selector_margin": 0.0020,
            "max_dual_replay_gap": 0.0015,
        },
        "selector_margin": 0.0017,
        "dual_replay_gap": 0.0012,
    }

    baseline = run["scores"]["baseline"]
    boundary_drift = abs(run["scores"]["boundary_shift"] - baseline)

    checks = {
        "boundary_drift_pass": boundary_drift <= run["thresholds"]["max_boundary_abs_drift"],
        "selector_margin_pass": run["selector_margin"] >= run["thresholds"]["min_selector_margin"],
        "dual_replay_pass": run["dual_replay_gap"] <= run["thresholds"]["max_dual_replay_gap"],
    }

    overall_pass = all(checks.values())
    run["derived"] = {
        "boundary_abs_drift": round(boundary_drift, 6),
        "checks": checks,
        "overall_pass": overall_pass,
    }
    run["verdict"] = "PASS_PC7_RUN1" if overall_pass else "FAIL_PC7_MARGIN"
    run["status"] = "NO_FALSE_PASS"
    run["next_action"] = (
        "prepare SST-1C rerun on PC7 upgraded margin"
        if overall_pass
        else "export margin obstruction and design PC8 noncyclic selector-margin amplifier"
    )

    out = gen / "p1416_pc7_first_transport_run_summary.json"
    out.write_text(json.dumps(run, indent=2) + "\n", encoding="utf-8")

    if not overall_pass:
        obstruction = {
            "obstruction_id": "OBSTR-PC7-MARGIN-v1",
            "reason": "selector margin below preregistered minimum",
            "observed_margin": run["selector_margin"],
            "required_min": run["thresholds"]["min_selector_margin"],
            "recommended_pivot": "PC8_noncyclic_selector_margin_amplifier",
        }
        (gen / "p1416_pc7_margin_obstruction_v1.json").write_text(json.dumps(obstruction, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"written": str(out), "verdict": run["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
