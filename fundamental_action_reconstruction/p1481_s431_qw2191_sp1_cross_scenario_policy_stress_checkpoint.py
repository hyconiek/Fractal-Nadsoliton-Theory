#!/usr/bin/env python3
"""P1481 S4.31: cross-scenario policy stress for SP1 in-policy sequences."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
POLICY = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"
SUMMARY = GEN / "p1481_s431_qw2191_sp1_cross_scenario_policy_stress_summary.json"
OBSTRUCTION = GEN / "p1481_s431_qw2191_sp1_cross_scenario_policy_stress_obstruction.json"

SCENARIOS = {
    "S1_descending": [-0.020, -0.021, -0.022],
    "S2_ascending": [-0.022, -0.021, -0.020],
    "S3_edge_then_recover": [-0.024, -0.023, -0.021],
}


def in_policy(sh: float, smin: float, smax: float) -> bool:
    return smin <= sh <= smax


def main() -> None:
    base = json.loads(P1468.read_text(encoding="utf-8"))
    policy = json.loads(POLICY.read_text(encoding="utf-8"))

    metric_a = float(base["arm_A_no_selector_premise_metric"])
    metric_b = float(base["arm_B_with_SP1_metric"])
    smin = float(policy["safe_shift_min"])
    smax = float(policy["safe_shift_max"])

    scenario_rows = []
    first_fail = None

    for name, shifts in SCENARIOS.items():
        rows = []
        for idx, sh in enumerate(shifts, start=1):
            if not in_policy(sh, smin, smax):
                row = {"iter": idx, "shift": sh, "status": "BLOCK"}
            else:
                delta = (metric_b + sh) - metric_a
                row = {"iter": idx, "shift": sh, "delta_metric_B_minus_A": delta, "status": "PASS" if delta > 0 else "FAIL"}
            rows.append(row)
            if row["status"] in {"FAIL", "BLOCK"} and first_fail is None:
                first_fail = {"scenario": name, **row}

        scenario_ok = all(r["status"] == "PASS" for r in rows)
        scenario_rows.append({"scenario": name, "rows": rows, "status": "PASS" if scenario_ok else "FAIL"})

    fail_count = sum(1 for s in scenario_rows if s["status"] == "FAIL")
    status = "PASS_SP1_CROSS_SCENARIO_LOCAL_ONLY" if fail_count == 0 else "FAIL_SP1_CROSS_SCENARIO"

    summary = {
        "packet": "P1481",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "policy_id": policy["policy_id"],
        "scenarios": scenario_rows,
        "fail_count": fail_count,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1481",
            "status": "FAIL_SP1_CROSS_SCENARIO",
            "first_fail": first_fail,
            "rule": "cross-scenario policy stress failed",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1481] status={status} fail_count={fail_count}")


if __name__ == "__main__":
    main()
