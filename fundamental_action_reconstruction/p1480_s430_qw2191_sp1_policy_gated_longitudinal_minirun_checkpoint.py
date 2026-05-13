#!/usr/bin/env python3
"""P1480 S4.30: policy-gated longitudinal mini-run for SP1 drift monitoring."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
POLICY = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"
SUMMARY = GEN / "p1480_s430_qw2191_sp1_policy_gated_longitudinal_minirun_summary.json"
OBSTRUCTION = GEN / "p1480_s430_qw2191_sp1_policy_gated_longitudinal_minirun_obstruction.json"

# all shifts are in-policy by construction
RUN_SHIFTS = [-0.022, -0.021, -0.020, -0.019]


def in_policy(shift: float, smin: float, smax: float) -> bool:
    return smin <= shift <= smax


def main() -> None:
    base = json.loads(P1468.read_text(encoding="utf-8"))
    policy = json.loads(POLICY.read_text(encoding="utf-8"))

    metric_a = float(base["arm_A_no_selector_premise_metric"])
    metric_b = float(base["arm_B_with_SP1_metric"])
    smin = float(policy["safe_shift_min"])
    smax = float(policy["safe_shift_max"])

    rows = []
    first_fail = None
    prev_delta = None
    for i, sh in enumerate(RUN_SHIFTS, start=1):
        if not in_policy(sh, smin, smax):
            row = {
                "iteration": i,
                "shift": sh,
                "status": "BLOCK",
                "reason": "out_of_policy",
            }
            rows.append(row)
            if first_fail is None:
                first_fail = row
            continue

        delta = (metric_b + sh) - metric_a
        drift = None if prev_delta is None else delta - prev_delta
        prev_delta = delta
        ok = delta > 0.0
        row = {
            "iteration": i,
            "shift": sh,
            "status": "PASS" if ok else "FAIL",
            "delta_metric_B_minus_A": delta,
            "delta_drift_vs_prev": drift,
        }
        rows.append(row)
        if (not ok) and first_fail is None:
            first_fail = row

    fail_count = sum(1 for r in rows if r["status"] in {"FAIL", "BLOCK"})
    status = "PASS_SP1_LONGITUDINAL_LOCAL_ONLY" if fail_count == 0 else "FAIL_SP1_LONGITUDINAL_LOCAL"

    summary = {
        "packet": "P1480",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "policy_id": policy["policy_id"],
        "rows": rows,
        "fail_count": fail_count,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1480",
            "status": "FAIL_SP1_LONGITUDINAL_LOCAL",
            "first_fail": first_fail,
            "rule": "policy-gated mini-run failed at local iteration",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1480] status={status} fail_count={fail_count}")


if __name__ == "__main__":
    main()
