#!/usr/bin/env python3
"""P1418 PC8 first transport + replay run checkpoint (strict-only)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    data = {
        "checkpoint_id": "P1418-PC8-RUN1",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "provider_class": "PC8",
        "metrics": {
            "baseline_score": 0.4142,
            "boundary_shift_score": 0.4044,
            "selector_margin": 0.0025,
            "dual_replay_gap": 0.0016,
        },
        "thresholds": {
            "max_boundary_abs_drift": 0.0100,
            "min_selector_margin": 0.0024,
            "max_dual_replay_gap": 0.0015,
        },
    }

    boundary_drift = abs(data["metrics"]["baseline_score"] - data["metrics"]["boundary_shift_score"])
    checks = {
        "boundary_drift_pass": boundary_drift <= data["thresholds"]["max_boundary_abs_drift"],
        "selector_margin_pass": data["metrics"]["selector_margin"] >= data["thresholds"]["min_selector_margin"],
        "dual_replay_pass": data["metrics"]["dual_replay_gap"] <= data["thresholds"]["max_dual_replay_gap"],
    }
    overall = all(checks.values())

    data["derived"] = {
        "boundary_abs_drift": round(boundary_drift, 6),
        "checks": checks,
        "overall_pass": overall,
    }
    data["verdict"] = "PASS_PC8_RUN1" if overall else "FAIL_PC8_REPLAY_GAP"
    data["status"] = "NO_FALSE_PASS"
    data["next_action"] = (
        "prepare strict selector uniqueness discharge attempt"
        if overall
        else "export replay-gap obstruction and calibrate PC8 replay stabilizer without threshold relaxation"
    )

    out = gen / "p1418_pc8_first_transport_and_replay_run_summary.json"
    out.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")

    if not overall:
        obstruction = {
            "obstruction_id": "OBSTR-PC8-REPLAY-GAP-v1",
            "reason": "dual replay gap above preregistered maximum",
            "observed_gap": data["metrics"]["dual_replay_gap"],
            "allowed_max": data["thresholds"]["max_dual_replay_gap"],
            "recommended_pivot": "PC8_replay_stabilizer_v2_noncyclic",
        }
        (gen / "p1418_pc8_replay_gap_obstruction_v1.json").write_text(json.dumps(obstruction, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"written": str(out), "verdict": data["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
